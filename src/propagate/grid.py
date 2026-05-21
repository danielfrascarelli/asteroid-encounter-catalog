"""Time-grid utilities for the propagation step.

The :func:`propagate_grid` iterator is the canonical step-by-step view of the
propagated catalogue and is consumed by the KD-tree scan.  When ``method`` is
``"rebound"`` the underlying N-body trajectory is computed once (eagerly or
from a memory-mapped cache) and replayed step-by-step; ``"kepler"`` falls back
to the analytic 2-body propagator that re-evaluates positions on demand.
"""

from __future__ import annotations

import logging
from collections.abc import Iterator

import numpy as np
import polars as pl

from src.propagate.kepler import propagate_df

logger = logging.getLogger(__name__)


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
    positions: np.ndarray | None = None,
) -> Iterator[tuple[float, np.ndarray]]:
    """Yield (t_jd, positions) for each step in *time_grid*.

    Two modes
    ---------
    - If *positions* is None, the iterator computes positions step-by-step
      via the analytic Kepler propagator.  Memory: O(N_asteroids).
    - If *positions* is given (shape ``(T, N, 3)``), the iterator simply
      replays its slices.  Used for the N-body branch and the cached
      trajectory in Phase 2.

    Parameters
    ----------
    elements:
        Orbital elements with columns: a_au, e, i_deg, Omega_deg, omega_deg,
        M_deg, epoch_jd.
    time_grid:
        JD TDB values to evaluate (e.g. from :func:`make_time_grid`).
    positions:
        Pre-computed positions of shape ``(T, N, 3)``.  When supplied,
        *elements* is unused (and may be empty).

    Yields
    ------
    tuple[float, np.ndarray]
        ``(t_jd, positions)`` where *positions* has shape ``(N_asteroids, 3)`` in AU.
    """
    if positions is not None:
        if positions.shape[0] != len(time_grid):
            raise ValueError(
                f"positions has T={positions.shape[0]} but time_grid has {len(time_grid)} steps"
            )
        for k, t in enumerate(time_grid):
            yield float(t), positions[k]
        return

    for t in time_grid:
        yield float(t), propagate_df(elements, float(t))


def propagate_full_grid(
    elements: pl.DataFrame,
    time_grid: np.ndarray,
    *,
    method: str = "kepler",
    rebound_kwargs: dict | None = None,
    cache_dir: str | None = None,
    cache_key: str | None = None,
    cache_format: str = "zarr",
) -> np.ndarray | None:
    """Pre-compute the full ``(T, N, 3)`` trajectory or return None for streaming Kepler.

    Used by the detection pipeline to materialise the N-body trajectory once
    (so parallel workers can share a memory-mapped or in-memory copy instead of
    each re-integrating from the snapshot epoch).

    Parameters
    ----------
    elements:
        Orbital elements DataFrame.
    time_grid:
        JD TDB epochs.
    method:
        ``"kepler"`` — return ``None`` (workers will use the streaming
        analytic propagator, which is cheap and stateless).
        ``"rebound"`` — integrate the full grid with planetary perturbers
        and return the trajectory array.
    rebound_kwargs:
        Extra keyword arguments forwarded to
        :func:`src.propagate.nbody.propagate_grid_nbody`.
    cache_dir:
        Optional directory under which to store / load the trajectory
        (Phase 2).  When supplied together with *cache_key*, the trajectory
        is memmapped from disk after first computation.
    cache_key:
        Unique key identifying this trajectory in the cache (see
        :mod:`src.propagate.cache`).

    Returns
    -------
    np.ndarray or None
        Trajectory of shape ``(T, N, 3)`` for the N-body method, otherwise
        ``None`` (streaming Kepler).
    """
    method = method.lower()
    if method == "kepler":
        return None
    if method == "rebound":
        rebound_kwargs = dict(rebound_kwargs or {})
        if cache_dir is not None and cache_key is not None:
            from src.propagate.cache import load_or_compute_trajectory

            return load_or_compute_trajectory(
                elements=elements,
                time_grid=time_grid,
                cache_dir=cache_dir,
                cache_key=cache_key,
                rebound_kwargs=rebound_kwargs,
                cache_format=cache_format,
            )
        from src.propagate.nbody import propagate_grid_nbody

        return propagate_grid_nbody(elements, time_grid, **rebound_kwargs)
    raise ValueError(f"Unknown propagation method: {method!r}")
