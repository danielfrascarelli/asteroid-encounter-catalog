"""Detection pipeline — orchestrates prefilter → KD-tree scan → refinement.

Public entry point: :func:`detect_encounters`.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import polars as pl

from src.detect.kdtree_scan import scan_time_grid
from src.detect.parallel import scan_parallel
from src.detect.prefilter import compatible_pairs
from src.detect.refine import refine_candidates

logger = logging.getLogger(__name__)

# Above this N, np.triu_indices materialises O(N²) pairs (>35 GB at N=94k).
# Skip prefilter and rely on the KD-tree spatial query alone.
_PREFILTER_MAX_N = 5_000
PREFILTER_MAX_N = _PREFILTER_MAX_N  # public alias for provenance/auditing


def effective_prefilter_mode(n: int, *, enabled: bool, max_n: int = _PREFILTER_MAX_N) -> str:
    """Return how the orbital pair-list prefilter actually behaved for ``n`` bodies.

    The ``|Δa|``/``|Δi|`` pair-list prefilter (:func:`compatible_pairs`) is only
    built when ``n <= max_n``; above that it is silently skipped because
    materialising O(N²) pairs is infeasible (the KD-tree spatial query alone is
    used instead).  The declared config flag ``prefilter.enabled`` therefore does
    **not** tell you whether the ``|Δa|`` cut was applied to a given catalog.
    This helper makes the *effective* behaviour explicit for the provenance
    sidecar.

    Returns
    -------
    str
        ``"disabled"`` (config off), ``"applied"`` (pair-list built), or
        ``"skipped_large_n"`` (config on but ``n > max_n`` → KD-tree only, so the
        ``|Δa|`` cut did NOT shape the catalog).
    """
    if not enabled:
        return "disabled"
    return "applied" if n <= max_n else "skipped_large_n"


_SCHEMA = {
    "number_1": pl.Int32,
    "number_2": pl.Int32,
    "designation_1": pl.Utf8,
    "designation_2": pl.Utf8,
    "jd_tdb": pl.Float64,
    "dist_au": pl.Float64,
    "rel_vel_au_day": pl.Float64,
}

# Coarse-scan checkpoint: candidate tuples (idx_i, idx_j, t_coarse_jd, d_coarse_au)
# where idx_* are row indices into the *elements* DataFrame used for the scan.
# Persisted between the scan and the (long, memory-heavy) refinement so a crash in
# refinement does not lose the scan; a resume run reads it and refines directly.
# The indices are only valid for the SAME elements/subset that produced them.
_SCAN_CKPT_SCHEMA = {
    "idx_i": pl.Int64,
    "idx_j": pl.Int64,
    "t_coarse_jd": pl.Float64,
    "d_coarse_au": pl.Float64,
}


def write_scan_checkpoint(
    candidates: list[tuple[int, int, float, float]], path: str | Path
) -> None:
    """Persist coarse-scan candidates so refinement can be resumed after a crash."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if candidates:
        idx_i, idx_j, t_c, d_c = zip(*candidates)
    else:
        idx_i, idx_j, t_c, d_c = (), (), (), ()
    pl.DataFrame(
        {
            "idx_i": pl.Series(idx_i, dtype=pl.Int64),
            "idx_j": pl.Series(idx_j, dtype=pl.Int64),
            "t_coarse_jd": pl.Series(t_c, dtype=pl.Float64),
            "d_coarse_au": pl.Series(d_c, dtype=pl.Float64),
        },
        schema=_SCAN_CKPT_SCHEMA,
    ).write_parquet(path, compression="zstd")


def load_scan_checkpoint(path: str | Path) -> list[tuple[int, int, float, float]]:
    """Load coarse-scan candidates written by :func:`write_scan_checkpoint`."""
    df = pl.read_parquet(path)
    return list(
        zip(
            df["idx_i"].to_list(),
            df["idx_j"].to_list(),
            df["t_coarse_jd"].to_list(),
            df["d_coarse_au"].to_list(),
        )
    )


def detect_encounters(
    elements: pl.DataFrame,
    time_grid: np.ndarray,
    *,
    threshold_au: float,
    semimajor_diff_max_au: float,
    inclination_diff_max_deg: float,
    leaf_size: int,
    fine_step_seconds: float,
    window_hours: float,
    prefilter_enabled: bool,
    refinement_enabled: bool,
    n_workers: int | str,
    chunk_size_days: float,
    positions: np.ndarray | None = None,
    query_radius_au: float | None = None,
    force_kepler_refine: bool = False,
    scan_checkpoint_path: str | Path | None = None,
    resume_from_scan: str | Path | None = None,
) -> pl.DataFrame:
    """Detect close asteroid encounters over a time grid.

    Parameters
    ----------
    elements:
        Orbital elements DataFrame.  Required columns: ``number`` (Int32),
        ``designation`` (Utf8), ``a_au``, ``e``, ``i_deg``, ``Omega_deg``,
        ``omega_deg``, ``M_deg``, ``epoch_jd`` (all Float64).
    time_grid:
        JD TDB values to scan (e.g. from
        :func:`src.propagate.grid.make_time_grid`).
    threshold_au:
        Maximum closest-approach distance to record (AU).
    semimajor_diff_max_au:
        Prefilter: maximum |a₁ - a₂| (AU).
    inclination_diff_max_deg:
        Prefilter: maximum |i₁ - i₂| (degrees).
    leaf_size:
        ``cKDTree`` leaf_size parameter.
    fine_step_seconds:
        Fine grid time step for the refinement pass (seconds).
    window_hours:
        Half-width of the fine search window around each coarse epoch (hours).
    prefilter_enabled:
        Set to ``False`` to scan all N*(N-1)/2 pairs (for small test sets).
    refinement_enabled:
        Set to ``False`` to skip quadratic refinement (uses coarse results).
    positions:
        Optional ``(T, N, 3)`` pre-computed positions (e.g. from the N-body
        propagator or cache).  When supplied the coarse scan reads positions
        directly instead of re-propagating from Kepler elements.
    scan_checkpoint_path:
        If given, the coarse-scan candidates are written to this Parquet path
        *before* refinement (insurance: refinement is long and memory-heavy, so a
        crash there would otherwise lose the whole scan).  Ignored when resuming.
    resume_from_scan:
        If given, the prefilter and coarse scan are **skipped** and candidates are
        loaded from this checkpoint (written by a previous run's
        ``scan_checkpoint_path``); refinement then runs on them.  ``elements`` and
        ``time_grid`` must match the run that produced the checkpoint (the stored
        indices are row indices into ``elements``).  Lets a killed refinement be
        resumed later with a different ``n_workers``.

    Returns
    -------
    pl.DataFrame
        Columns: ``number_1``, ``number_2``, ``designation_1``,
        ``designation_2``, ``jd_tdb``, ``dist_au``, ``rel_vel_au_day``.
        Sorted by ``dist_au`` ascending.  Each pair appears at most once
        (the closest approach only).
    """
    n = len(elements)
    logger.info(
        "detect_encounters: %d asteroids | %d steps | threshold=%.5f AU",
        n,
        len(time_grid),
        threshold_au,
    )

    # --- Config validation (tribunal B1) ---
    # The Kepler refiner samples ±window_hours around each coarse sample. The
    # true minimum can be up to half a coarse step away from the nearest grid
    # point; a narrower window clips the epoch at the window edge and biases
    # the distance high (~60 % of the pre-2026-07 frozen catalog).
    if refinement_enabled and len(time_grid) > 1:
        grid_step_hours = float(time_grid[1] - time_grid[0]) * 24.0
        uses_kepler_refiner = force_kepler_refine or positions is None
        if uses_kepler_refiner and window_hours < grid_step_hours / 2.0:
            raise ValueError(
                f"detection.refinement.window_hours={window_hours:g} h is smaller than "
                f"half the coarse grid step ({grid_step_hours:g} h / 2 = "
                f"{grid_step_hours / 2.0:g} h): true minima between coarse samples "
                "would be clipped at the window edge (B1). Increase window_hours to "
                f"≥ {grid_step_hours / 2.0:g}."
            )

    # --- Resume path: skip prefilter + scan, load candidates from checkpoint ---
    if resume_from_scan is not None:
        candidates = load_scan_checkpoint(resume_from_scan)
        logger.info(
            "Resumed %d coarse candidates from scan checkpoint %s (prefilter + scan skipped)",
            len(candidates),
            resume_from_scan,
        )
        from src.detect.parallel import resolve_n_workers

        nw = resolve_n_workers(n_workers)
        return _refine_and_finalize(
            elements,
            candidates,
            threshold_au=threshold_au,
            fine_step_seconds=fine_step_seconds,
            window_hours=window_hours,
            refinement_enabled=refinement_enabled,
            nw=nw,
            positions=positions,
            time_grid=time_grid,
            force_kepler_refine=force_kepler_refine,
        )

    # --- Step 1: prefilter ---
    pairs: np.ndarray | None

    if prefilter_enabled:
        if n <= _PREFILTER_MAX_N:
            pairs = compatible_pairs(elements, semimajor_diff_max_au, inclination_diff_max_deg)
        else:
            logger.info(
                "N=%d > %d: skipping pair precomputation, KD-tree spatial filter only",
                n,
                _PREFILTER_MAX_N,
            )
            pairs = None
    else:
        if n <= _PREFILTER_MAX_N:
            ii, jj = np.triu_indices(n, k=1)
            pairs = np.stack([ii, jj], axis=1).astype(np.int32)
            logger.info("Prefilter disabled: %d pairs", len(pairs))
        else:
            logger.info("Prefilter disabled; N=%d — using KD-tree spatial filter only", n)
            pairs = None

    if pairs is not None and len(pairs) == 0:
        logger.info("No compatible pairs — catalog is empty")
        return pl.DataFrame(schema=_SCHEMA)

    # --- Step 2: KD-tree coarse scan ---
    from src.detect.parallel import resolve_n_workers

    nw = resolve_n_workers(n_workers)
    if query_radius_au is not None and query_radius_au > threshold_au:
        logger.info(
            "KD-tree query radius widened to %.5f AU (vs threshold %.5f AU) to "
            "compensate for coarse temporal sampling",
            query_radius_au,
            threshold_au,
        )
    if nw > 1:
        candidates = scan_parallel(
            elements,
            time_grid,
            pairs,
            threshold_au,
            leaf_size,
            n_workers,
            chunk_size_days,
            positions=positions,
            query_radius_au=query_radius_au,
        )
    else:
        candidates = scan_time_grid(
            elements,
            time_grid,
            pairs,
            threshold_au,
            leaf_size,
            positions=positions,
            query_radius_au=query_radius_au,
        )
    logger.info("%d coarse candidates after KD-tree scan", len(candidates))

    # --- Checkpoint: persist the scan before the long/memory-heavy refinement ---
    if scan_checkpoint_path is not None:
        write_scan_checkpoint(candidates, scan_checkpoint_path)
        logger.info(
            "Scan checkpoint written: %d candidates → %s (refinement resumable via "
            "resume_from_scan)",
            len(candidates),
            scan_checkpoint_path,
        )

    if not candidates:
        return pl.DataFrame(schema=_SCHEMA)

    # --- Step 3: refinement ---
    return _refine_and_finalize(
        elements,
        candidates,
        threshold_au=threshold_au,
        fine_step_seconds=fine_step_seconds,
        window_hours=window_hours,
        refinement_enabled=refinement_enabled,
        nw=nw,
        positions=positions,
        time_grid=time_grid,
        force_kepler_refine=force_kepler_refine,
    )


def _refine_and_finalize(
    elements: pl.DataFrame,
    candidates: list[tuple[int, int, float, float]],
    *,
    threshold_au: float,
    fine_step_seconds: float,
    window_hours: float,
    refinement_enabled: bool,
    nw: int,
    positions: np.ndarray | None,
    time_grid: np.ndarray | None,
    force_kepler_refine: bool,
) -> pl.DataFrame:
    """Refine coarse candidates (or pass them through) and return the sorted catalog.

    Shared by the normal path and the ``resume_from_scan`` path so both produce
    identical output.
    """
    if not candidates:
        return pl.DataFrame(schema=_SCHEMA)

    if refinement_enabled:
        # When the bulk cache is coarse (Strategy A), the quadratic-over-cache
        # refinement loses precision (3-point parabola over 12 h grid is way
        # less accurate than a 60 s Kepler scan). Pass force_kepler_refine=True
        # to skip the cache path inside refine_candidates and re-evaluate every
        # candidate with the 2-body propagator on a ±window_hours fine grid.
        refine_positions = None if force_kepler_refine else positions
        refine_time_grid = None if force_kepler_refine else time_grid
        result = refine_candidates(
            elements,
            candidates,
            threshold_au,
            fine_step_seconds,
            window_hours,
            positions=refine_positions,
            time_grid=refine_time_grid,
            n_workers=nw,
        )
    else:
        rows = [
            {
                "number_1": elements["number"][idx_i],
                "number_2": elements["number"][idx_j],
                "designation_1": elements["designation"][idx_i],
                "designation_2": elements["designation"][idx_j],
                "jd_tdb": t_jd,
                "dist_au": dist_au,
                "rel_vel_au_day": float("nan"),
            }
            for idx_i, idx_j, t_jd, dist_au in candidates
        ]
        result = pl.DataFrame(rows, schema=_SCHEMA)

    logger.info("Detection complete: %d encounters ≤ %.5f AU", len(result), threshold_au)
    return result.sort("dist_au")
