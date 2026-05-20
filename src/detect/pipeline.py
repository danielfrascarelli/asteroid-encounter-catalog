"""Detection pipeline — orchestrates prefilter → KD-tree scan → refinement.

Public entry point: :func:`detect_encounters`.
"""

from __future__ import annotations

import logging

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

_SCHEMA = {
    "number_1": pl.Int32,
    "number_2": pl.Int32,
    "designation_1": pl.Utf8,
    "designation_2": pl.Utf8,
    "jd_tdb": pl.Float64,
    "dist_au": pl.Float64,
    "rel_vel_au_day": pl.Float64,
}


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

    if not candidates:
        return pl.DataFrame(schema=_SCHEMA)

    # --- Step 3: refinement ---
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
