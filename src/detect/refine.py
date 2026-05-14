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

import numpy as np
import polars as pl

from src.propagate.kepler import kepler_to_cartesian

logger = logging.getLogger(__name__)

_DEG = np.pi / 180.0
_SECONDS_PER_DAY = 86400.0


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


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def refine_candidates(
    elements: pl.DataFrame,
    candidates: list[tuple[int, int, float, float]],
    threshold_au: float,
    fine_step_seconds: float = 60.0,
    window_hours: float = 2.0,
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

    # Pre-materialise element rows for fast per-row access
    elem_rows = [{col: elements[col][k] for col in elements.columns} for k in range(len(elements))]

    rows: list[dict] = []

    for idx_i, idx_j, t_coarse, _d_coarse in candidates:
        row_i = elem_rows[idx_i]
        row_j = elem_rows[idx_j]

        t_start = t_coarse - half_window_days
        t_end = t_coarse + half_window_days
        t_fine = np.arange(t_start, t_end + fine_step_days * 0.5, fine_step_days)

        if len(t_fine) < 3:
            # Window too narrow for interpolation — use coarse values
            t_min = t_coarse
            d_min = _d_coarse
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
            continue

        # Relative velocity: centred finite differences at t_min
        dt = fine_step_days
        t_vel = np.array([t_min - dt, t_min + dt])
        pos_i_vel, pos_j_vel = _propagate_pair(row_i, row_j, t_vel)
        rel_vel_vec = (pos_i_vel[1] - pos_i_vel[0] - (pos_j_vel[1] - pos_j_vel[0])) / (2.0 * dt)
        rel_vel_au_day = float(np.linalg.norm(rel_vel_vec))

        rows.append(
            {
                "number_1": row_i["number"],
                "number_2": row_j["number"],
                "designation_1": row_i["designation"],
                "designation_2": row_j["designation"],
                "jd_tdb": t_min,
                "dist_au": d_min,
                "rel_vel_au_day": rel_vel_au_day,
            }
        )

    if not rows:
        return pl.DataFrame(schema=schema)

    return pl.DataFrame(rows, schema=schema)
