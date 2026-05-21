"""Close-encounter refinement — sub-grid search and quadratic interpolation.

For each coarse candidate ``(idx_i, idx_j, t_coarse, d_coarse)`` from the
KD-tree scan, this module:

1. Samples a fine-grained time window of ±``window_hours`` around ``t_coarse``
   at ``fine_step_seconds`` resolution.
2. Finds the true minimum-distance epoch with quadratic interpolation over the
   three grid points surrounding the argmin.
3. Computes the relative velocity at the minimum via centred finite differences
   (±1 fine step).
4. Drops encounters whose refined distance exceeds ``threshold_au``.
"""

from __future__ import annotations

import logging
import multiprocessing

import numpy as np
import polars as pl

from src.propagate.kepler import kepler_to_cartesian

logger = logging.getLogger(__name__)

_DEG = np.pi / 180.0
_SECONDS_PER_DAY = 86400.0


# ---------------------------------------------------------------------------
# Worker state for parallel Kepler refinement
# ---------------------------------------------------------------------------

_WORKER_ELEM_ROWS: list[dict] = []
_WORKER_THRESHOLD: float = 0.0
_WORKER_FINE_STEP: float = 0.0
_WORKER_HALF_WIN: float = 0.0


def _init_worker(
    elem_rows: list[dict],
    threshold: float,
    fine_step: float,
    half_win: float,
) -> None:
    global _WORKER_ELEM_ROWS, _WORKER_THRESHOLD, _WORKER_FINE_STEP, _WORKER_HALF_WIN
    _WORKER_ELEM_ROWS = elem_rows
    _WORKER_THRESHOLD = threshold
    _WORKER_FINE_STEP = fine_step
    _WORKER_HALF_WIN = half_win


def _refine_chunk(chunk: list[tuple[int, int, float, float]]) -> list[dict]:
    """Refine one chunk of Kepler-path candidates in a worker process."""
    rows: list[dict] = []
    for idx_i, idx_j, t_coarse, d_coarse in chunk:
        result = _refine_one_kepler(
            _WORKER_ELEM_ROWS[idx_i],
            _WORKER_ELEM_ROWS[idx_j],
            t_coarse,
            d_coarse,
            _WORKER_THRESHOLD,
            _WORKER_FINE_STEP,
            _WORKER_HALF_WIN,
        )
        if result is not None:
            rows.append(result)
    return rows


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _quadratic_min(
    t0: float,
    t1: float,
    t2: float,
    d0: float,
    d1: float,
    d2: float,
) -> tuple[float, float]:
    """Return (t_min, d_min) of the parabola through three equally-spaced points.

    Parameters
    ----------
    t0, t1, t2:
        Times (assumed t2 - t1 = t1 - t0 = h).
    d0, d1, d2:
        Distance values at those times.

    Returns
    -------
    (t_min, d_min)
        Vertex of the upward-opening parabola, clamped to [t0, t2].
        Falls back to the discrete argmin when the parabola opens downward.
    """
    h = t1 - t0
    denom = d0 - 2.0 * d1 + d2  # 2*A*h²; > 0 means parabola opens upward
    if denom <= 0.0:
        idx = int(np.argmin([d0, d1, d2]))
        return [t0, t1, t2][idx], [d0, d1, d2][idx]

    # Vertex: t_min = t1 + h*(d0 - d2) / (2*denom)
    dt = h * (d0 - d2) / (2.0 * denom)
    t_min = float(np.clip(t1 + dt, t0, t2))
    dt_c = t_min - t1  # clamped delta
    a_coef = denom / (2.0 * h * h)
    b_coef = (d2 - d0) / (2.0 * h)
    d_min = a_coef * dt_c**2 + b_coef * dt_c + d1
    # Numerical safety: never return a distance larger than the grid minimum
    d_min = min(float(d_min), d0, d1, d2)
    return t_min, d_min


def _propagate_pair(
    row_i: dict,
    row_j: dict,
    t_array: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return Cartesian positions for two asteroids at all *t_array* epochs.

    Parameters
    ----------
    row_i, row_j:
        Dicts with keys ``a_au, e, i_deg, Omega_deg, omega_deg, M_deg,
        epoch_jd`` (Python scalars).
    t_array:
        JD TDB values, shape ``(T,)``.

    Returns
    -------
    (pos_i, pos_j), each shape (T, 3)
    """

    def _pos(row: dict) -> np.ndarray:
        return kepler_to_cartesian(
            a_au=row["a_au"],
            e=row["e"],
            i_rad=float(row["i_deg"]) * _DEG,
            Omega_rad=float(row["Omega_deg"]) * _DEG,
            omega_rad=float(row["omega_deg"]) * _DEG,
            M0_rad=float(row["M_deg"]) * _DEG,
            epoch_jd=row["epoch_jd"],
            t_jd=t_array,
        )

    return _pos(row_i), _pos(row_j)


def _refine_one_kepler(
    row_i: dict,
    row_j: dict,
    t_coarse: float,
    d_coarse: float,
    threshold_au: float,
    fine_step_days: float,
    half_window_days: float,
) -> dict | None:
    """Refine a single Kepler-path candidate; returns None if distance exceeds threshold."""
    t_fine = np.arange(
        t_coarse - half_window_days,
        t_coarse + half_window_days + fine_step_days * 0.5,
        fine_step_days,
    )
    if len(t_fine) < 3:
        t_min = t_coarse
        d_min = d_coarse
    else:
        pos_i, pos_j = _propagate_pair(row_i, row_j, t_fine)
        dists = np.linalg.norm(pos_i - pos_j, axis=1)
        k = int(np.argmin(dists))
        if 0 < k < len(t_fine) - 1:
            t_min, d_min = _quadratic_min(
                float(t_fine[k - 1]),
                float(t_fine[k]),
                float(t_fine[k + 1]),
                float(dists[k - 1]),
                float(dists[k]),
                float(dists[k + 1]),
            )
        else:
            t_min = float(t_fine[k])
            d_min = float(dists[k])

    if d_min > threshold_au:
        return None

    dt = fine_step_days
    t_vel = np.array([t_min - dt, t_min + dt])
    pos_i_vel, pos_j_vel = _propagate_pair(row_i, row_j, t_vel)
    rel_vel_vec = (pos_i_vel[1] - pos_i_vel[0] - (pos_j_vel[1] - pos_j_vel[0])) / (2.0 * dt)
    return {
        "number_1": row_i["number"],
        "number_2": row_j["number"],
        "designation_1": row_i["designation"],
        "designation_2": row_j["designation"],
        "jd_tdb": t_min,
        "dist_au": d_min,
        "rel_vel_au_day": float(np.linalg.norm(rel_vel_vec)),
    }


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def refine_candidates(
    elements: pl.DataFrame,
    candidates: list[tuple[int, int, float, float]],
    threshold_au: float,
    fine_step_seconds: float = 60.0,
    window_hours: float = 2.0,
    positions: np.ndarray | None = None,
    time_grid: np.ndarray | None = None,
    n_workers: int = 1,
) -> pl.DataFrame:
    """Refine coarse KD-tree candidates to find true minimum-distance epochs.

    Parameters
    ----------
    elements:
        Orbital elements DataFrame.  Must include ``number`` (Int32),
        ``designation`` (Utf8), and the orbital columns.
    candidates:
        ``(idx_i, idx_j, t_coarse_jd, d_coarse_au)`` tuples from the
        KD-tree scan.
    threshold_au:
        Only encounters with refined distance ≤ this value are returned.
    fine_step_seconds:
        Time resolution of the fine search grid in seconds.
    window_hours:
        Half-width of the fine search window around each coarse epoch (hours).
    n_workers:
        Number of worker processes for the Kepler refinement path.  Has no
        effect when *positions* / *time_grid* are supplied (cache path is
        already fast and avoids memmap thrashing across processes).

    Returns
    -------
    pl.DataFrame
        Columns: ``number_1`` (Int32), ``number_2`` (Int32),
        ``designation_1`` (Utf8), ``designation_2`` (Utf8),
        ``jd_tdb`` (Float64), ``dist_au`` (Float64),
        ``rel_vel_au_day`` (Float64).
    """
    schema = {
        "number_1": pl.Int32,
        "number_2": pl.Int32,
        "designation_1": pl.Utf8,
        "designation_2": pl.Utf8,
        "jd_tdb": pl.Float64,
        "dist_au": pl.Float64,
        "rel_vel_au_day": pl.Float64,
    }

    if not candidates:
        return pl.DataFrame(schema=schema)

    fine_step_days = fine_step_seconds / _SECONDS_PER_DAY
    half_window_days = window_hours / 24.0

    # When a pre-computed N-body trajectory is supplied (rebound + cache), use
    # it for the per-candidate minimum instead of re-propagating with Kepler.
    # Falling back to Kepler here would overwrite N-body distances with 2-body
    # values, undoing the entire point of running rebound.
    use_cache = positions is not None and time_grid is not None
    cache_step_days = 0.0
    if use_cache:
        assert time_grid is not None  # mypy guard
        cache_step_days = float(time_grid[1] - time_grid[0]) if len(time_grid) > 1 else 0.0

    # Pre-materialise element rows for fast per-row access
    elem_rows = [{col: elements[col][k] for col in elements.columns} for k in range(len(elements))]

    rows: list[dict] = []

    # Parallel Kepler path: each candidate is independent, no shared mutable state.
    # Not used for the cache path — the cache is a large memmap that shouldn't be
    # re-opened across many processes (page-fault thrashing).
    if n_workers > 1 and not use_cache:
        chunk_size = max(500, len(candidates) // (n_workers * 4))
        chunks = [candidates[i : i + chunk_size] for i in range(0, len(candidates), chunk_size)]
        logger.info(
            "Parallel refinement: %d candidates → %d chunks × %d workers",
            len(candidates),
            len(chunks),
            n_workers,
        )
        with multiprocessing.Pool(
            n_workers,
            initializer=_init_worker,
            initargs=(elem_rows, threshold_au, fine_step_days, half_window_days),
        ) as pool:
            for chunk_rows in pool.imap_unordered(_refine_chunk, chunks):
                rows.extend(chunk_rows)
    elif use_cache:
        # Batch all k0 lookups upfront and sort candidates by k0 so that sequential
        # candidates touch the same memmap page region, trading random I/O for
        # sequential I/O on the positions array.
        assert positions is not None and time_grid is not None
        t_coarses = np.array([c[2] for c in candidates])
        k0s = np.clip(np.searchsorted(time_grid, t_coarses), 1, len(time_grid) - 2).astype(int)
        order = np.argsort(k0s, kind="stable")

        _prev_k0 = -1
        _slab: np.ndarray | None = None

        for oi in order:
            idx_i, idx_j, t_coarse, _ = candidates[oi]
            k0 = k0s[oi]
            row_i = elem_rows[idx_i]
            row_j = elem_rows[idx_j]

            if k0 != _prev_k0:
                _slab = positions[k0 - 1 : k0 + 2]  # memmap view; sort ensures locality
                _prev_k0 = k0
            assert _slab is not None  # _prev_k0=-1 forces the first branch

            p_i = _slab[:, idx_i]
            p_j = _slab[:, idx_j]
            d3 = np.linalg.norm(p_i - p_j, axis=1)
            t_min, d_min = _quadratic_min(
                float(time_grid[k0 - 1]),
                float(time_grid[k0]),
                float(time_grid[k0 + 1]),
                float(d3[0]),
                float(d3[1]),
                float(d3[2]),
            )

            if d_min > threshold_au:
                continue

            rel_vel_vec = (p_i[2] - p_i[0] - (p_j[2] - p_j[0])) / (2.0 * cache_step_days)
            rows.append(
                {
                    "number_1": row_i["number"],
                    "number_2": row_j["number"],
                    "designation_1": row_i["designation"],
                    "designation_2": row_j["designation"],
                    "jd_tdb": t_min,
                    "dist_au": d_min,
                    "rel_vel_au_day": float(np.linalg.norm(rel_vel_vec)),
                }
            )
    else:
        for idx_i, idx_j, t_coarse, d_coarse in candidates:
            row_i = elem_rows[idx_i]
            row_j = elem_rows[idx_j]
            result = _refine_one_kepler(
                row_i, row_j, t_coarse, d_coarse, threshold_au, fine_step_days, half_window_days
            )
            if result is None:
                continue
            rows.append(result)

    if not rows:
        return pl.DataFrame(schema=schema)

    return pl.DataFrame(rows, schema=schema)
