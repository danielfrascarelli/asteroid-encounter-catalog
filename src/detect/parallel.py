"""Parallel temporal scan — splits the time grid into chunks for multiprocessing.Pool.

Each chunk is an independent slice of the time grid; workers find the best
(minimum-distance) epoch per pair within their slice.  Results are merged
across chunks before the refinement step.

Public entry point: :func:`scan_parallel`.
"""

from __future__ import annotations

import logging
import multiprocessing as mp
import os

import numpy as np
import polars as pl
from tqdm import tqdm

from src.detect.kdtree_scan import scan_time_grid

logger = logging.getLogger(__name__)

# Per-worker globals — set once by the pool initializer, never mutated.
_G_ELEMENTS: pl.DataFrame | None = None
_G_PAIRS: np.ndarray | None = None


def _init_worker(elements: pl.DataFrame, pairs: np.ndarray | None) -> None:
    global _G_ELEMENTS, _G_PAIRS
    _G_ELEMENTS = elements
    _G_PAIRS = pairs


def _scan_chunk(
    args: tuple[np.ndarray, float, int],
) -> list[tuple[int, int, float, float]]:
    chunk, threshold_au, leaf_size = args
    assert _G_ELEMENTS is not None
    return scan_time_grid(_G_ELEMENTS, chunk, _G_PAIRS, threshold_au, leaf_size)


def _make_chunks(time_grid: np.ndarray, chunk_size_days: float) -> list[np.ndarray]:
    """Split *time_grid* into contiguous sub-arrays of ~*chunk_size_days*."""
    step_days = float(time_grid[1] - time_grid[0]) if len(time_grid) > 1 else 1.0
    steps_per_chunk = max(1, int(round(chunk_size_days / step_days)))
    return [time_grid[i : i + steps_per_chunk] for i in range(0, len(time_grid), steps_per_chunk)]


def _merge_candidates(
    results: list[list[tuple[int, int, float, float]]],
) -> list[tuple[int, int, float, float]]:
    """Merge per-chunk lists, keeping the minimum-distance epoch per pair."""
    best: dict[tuple[int, int], tuple[float, float]] = {}
    for chunk_result in results:
        for idx_i, idx_j, t_jd, dist in chunk_result:
            key = (idx_i, idx_j)
            if key not in best or dist < best[key][1]:
                best[key] = (t_jd, dist)
    return [(k[0], k[1], v[0], v[1]) for k, v in best.items()]


def resolve_n_workers(n_workers: int | str) -> int:
    """Resolve ``"auto"`` to ``os.cpu_count()``, otherwise cast to int."""
    if n_workers == "auto":
        return os.cpu_count() or 1
    return int(n_workers)


def scan_parallel(
    elements: pl.DataFrame,
    time_grid: np.ndarray,
    pairs: np.ndarray | None,
    threshold_au: float,
    leaf_size: int = 30,
    n_workers: int | str = "auto",
    chunk_size_days: float = 30.0,
) -> list[tuple[int, int, float, float]]:
    """Parallel drop-in replacement for :func:`~src.detect.kdtree_scan.scan_time_grid`.

    Splits *time_grid* into chunks of ~*chunk_size_days* days and distributes
    them across *n_workers* processes.  The per-chunk results are merged by
    keeping the minimum-distance epoch per asteroid pair.

    Parameters
    ----------
    elements:
        Orbital elements DataFrame (same as for :func:`scan_time_grid`).
    time_grid:
        Full JD TDB array to scan.
    pairs:
        Prefilter-compatible pair indices, or ``None`` to accept all spatial
        neighbours from the KD-tree query.
    threshold_au:
        Maximum inter-asteroid distance to record (AU).
    leaf_size:
        ``cKDTree`` leaf_size parameter.
    n_workers:
        Number of worker processes.  ``"auto"`` uses ``os.cpu_count()``.
    chunk_size_days:
        Approximate length of each time chunk in days.

    Returns
    -------
    list of (idx_i, idx_j, best_t_jd, best_dist_au)
        One entry per pair whose minimum observed distance is ≤ *threshold_au*.
    """
    nw = resolve_n_workers(n_workers)
    chunks = _make_chunks(time_grid, chunk_size_days)
    tasks = [(chunk, threshold_au, leaf_size) for chunk in chunks]

    logger.info(
        "Parallel scan: %d workers | %d chunks (~%.0f days each) | %d total steps",
        nw,
        len(chunks),
        chunk_size_days,
        len(time_grid),
    )

    # Limit numpy/BLAS/OpenMP to 1 thread per worker to prevent oversubscription.
    # Without this, each of the N workers spawns ~N numpy threads → N² threads
    # competing for N CPUs, causing most cores to idle on context-switch overhead.
    # Must be set before Pool creation: forkserver spawns a fresh interpreter that
    # inherits the parent's env, so workers see these values before importing numpy.
    for _var in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_MAX_THREADS",
    ):
        os.environ.setdefault(_var, "1")

    # "forkserver" avoids the deadlock risk from forking a multi-threaded
    # parent (polars uses Arrow thread pools) while still being faster than
    # "spawn" on repeated runs inside Docker.
    ctx = mp.get_context("forkserver")
    with ctx.Pool(
        processes=nw,
        initializer=_init_worker,
        initargs=(elements, pairs),
    ) as pool:
        all_results: list[list[tuple[int, int, float, float]]] = []
        with tqdm(total=len(chunks), desc="Scanning chunks", unit="chunk") as pbar:
            for result in pool.imap_unordered(_scan_chunk, tasks):
                all_results.append(result)
                pbar.update(1)

    candidates = _merge_candidates(all_results)
    logger.info("Parallel scan complete: %d candidate pairs", len(candidates))
    return candidates
