"""Time-grid utilities for the propagation step."""

from __future__ import annotations

from collections.abc import Iterator

import numpy as np
import polars as pl

from src.propagate.kepler import propagate_df


def make_time_grid(start_jd: float, end_jd: float, step_hours: float) -> np.ndarray:
    """Return evenly-spaced JD TDB values from *start_jd* to *end_jd* (inclusive).

    Parameters
    ----------
    start_jd:
        Start of the grid in JD TDB.
    end_jd:
        End of the grid in JD TDB (included if reachable within half a step).
    step_hours:
        Time step in hours.

    Returns
    -------
    np.ndarray of float64, shape (N,)
    """
    step_days = step_hours / 24.0
    return np.arange(start_jd, end_jd + step_days * 0.5, step_days)


def propagate_grid(
    elements: pl.DataFrame,
    time_grid: np.ndarray,
) -> Iterator[tuple[float, np.ndarray]]:
    """Yield (t_jd, positions) for each step in *time_grid*.

    Positions are computed step-by-step to keep memory bounded at O(N_asteroids)
    rather than materialising the full (N_times × N_asteroids × 3) array.

    Parameters
    ----------
    elements:
        Orbital elements with columns: a_au, e, i_deg, Omega_deg, omega_deg,
        M_deg, epoch_jd.
    time_grid:
        JD TDB values to evaluate (e.g. from :func:`make_time_grid`).

    Yields
    ------
    tuple[float, np.ndarray]
        ``(t_jd, positions)`` where *positions* has shape ``(N_asteroids, 3)`` in AU.

    Notes
    -----
    Full-grid caching (for repeated detection runs on the same elements and
    time window) is deferred to Phase 4 (parallelisation).
    """
    for t in time_grid:
        yield float(t), propagate_df(elements, float(t))
