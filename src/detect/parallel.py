"""Parallel temporal scan — splits the time grid into chunks for multiprocessing.Pool.

Each chunk is an independent slice of the time grid; workers find the best
(minimum-distance) epoch per pair within their slice.  Results are merged
across chunks before the refinement step.

When pre-computed positions are supplied (N-body branch, or a cached
trajectory), workers receive only chunk-relative slices instead of
re-propagating — this is the path used by Phases 1–3.

Public entry point: :func:`scan_parallel`.
"""

from __future__ import annotations

import logging
import multiprocessing as mp
import os
import tempfile
from pathlib import Path
from typing import Any

import numpy as np
import polars as pl
from tqdm import tqdm

from src.detect.kdtree_scan import scan_time_grid

logger = logging.getLogger(__name__)

# Above this pickled size, the elements DataFrame is spilled to a Parquet
# tempfile instead of passed through the forkserver initargs pipe (which
# deadlocks at production scale). Module-level so tests can force the spill.
_ELEMENTS_SPILL_BYTES = 1_000_000

# Per-worker globals — set once by the pool initializer, never mutated.
_G_ELEMENTS: pl.DataFrame | None = None
_G_PAIRS: np.ndarray | None = None
_G_POSITIONS: Any | None = None


class _PositionsWindow:
    """Time-window adapter over a full trajectory array.

    ``np.memmap[start:stop]`` is cheap (returns a view), so memmap-backed
    trajectories are kept lazy and indexed one slab at a time.  Zarr-backed
    trajectories are *not* cheap to index one slab at a time — every read
    decompresses the underlying chunk regardless of any encoded-chunk cache.
    For those, eagerly slicing the full window once decompresses each
    overlapping chunk only once; subsequent slab reads are numpy views into
    the materialised block.
    """

    def __init__(self, positions: Any, start: int, stop: int) -> None:
        self._start = int(start)
        self.shape = (int(stop - start), positions.shape[1], positions.shape[2])
        # Eagerly materialise zarr-backed windows (TrajectoryView) into a
        # contiguous numpy block — one decompression per chunk vs one per slab.
        from src.propagate.cache import TrajectoryView

        if isinstance(positions, TrajectoryView):
            self._positions = np.asarray(positions[self._start : int(stop)])
            self._eager = True
        else:
            self._positions = positions
            self._eager = False

    def __getitem__(self, idx: int) -> np.ndarray:
        if self._eager:
            return np.asarray(self._positions[int(idx)])
        return np.asarray(self._positions[self._start + int(idx)])


def _init_worker(
    elements: pl.DataFrame | None,
    elements_path: str | None,
    pairs: np.ndarray | None,
    pairs_path: str | None,
    positions: np.ndarray | None,
    positions_memmap: tuple[str, tuple[int, int, int], str] | None,
    positions_zarr_path: str | None,
) -> None:
    """Initialise per-worker globals.

    Large objects (``elements`` DataFrame, ``pairs`` array, ``positions``
    trajectory) are passed by path when their pickle size would otherwise
    overwhelm the multiprocessing queue.  ``elements_path`` points to a Parquet
    file with the orbital elements (the full numbered population is ~32 MB —
    passing it through forkserver ``initargs`` to every worker deadlocks the
    forkserver pipe; observed at production scale hanging the scan at ~5/35
    chunks). ``pairs_path`` points to an ``.npy`` file with the prefilter pair
    indices; ``positions_memmap`` is the ``(filename, shape, dtype)`` triple
    used to re-open a disk-backed memmap; ``positions_zarr_path`` is a Zarr
    DirectoryStore path opened through :func:`src.propagate.cache.
    open_trajectory_for_worker` so the worker keeps recently-decompressed
    chunks resident in an LRU cache. Workers re-load these locally — the OS
    shares the underlying chunk files across processes via the page cache.
    """
    global _G_ELEMENTS, _G_PAIRS, _G_POSITIONS
    _G_ELEMENTS = pl.read_parquet(elements_path) if elements_path is not None else elements
    if pairs_path is not None:
        _G_PAIRS = np.load(pairs_path, mmap_mode="r")
    else:
        _G_PAIRS = pairs
    if positions_zarr_path is not None:
        from src.propagate.cache import open_trajectory_for_worker

        _G_POSITIONS = open_trajectory_for_worker(positions_zarr_path)
    elif positions_memmap is not None:
        filename, shape, dtype_str = positions_memmap
        _G_POSITIONS = np.memmap(filename, dtype=np.dtype(dtype_str), mode="r", shape=shape)
    else:
        _G_POSITIONS = positions


def _scan_chunk(
    args: tuple[np.ndarray, np.ndarray, float, int, float],
) -> list[tuple[int, int, float, float]]:
    chunk_times, chunk_indices, threshold_au, leaf_size, query_radius_au = args
    assert _G_ELEMENTS is not None
    positions_chunk: Any | None = None
    if _G_POSITIONS is not None:
        positions_chunk = _PositionsWindow(
            _G_POSITIONS,
            int(chunk_indices[0]),
            int(chunk_indices[-1]) + 1,
        )
    return scan_time_grid(
        _G_ELEMENTS,
        chunk_times,
        _G_PAIRS,
        threshold_au,
        leaf_size,
        positions=positions_chunk,
        query_radius_au=query_radius_au,
    )


def _make_chunks(
    time_grid: np.ndarray, chunk_size_days: float
) -> list[tuple[np.ndarray, np.ndarray]]:
    """Split *time_grid* into contiguous (chunk_times, chunk_indices) pairs."""
    step_days = float(time_grid[1] - time_grid[0]) if len(time_grid) > 1 else 1.0
    steps_per_chunk = max(1, int(round(chunk_size_days / step_days)))
    chunks: list[tuple[np.ndarray, np.ndarray]] = []
    for i in range(0, len(time_grid), steps_per_chunk):
        end = min(i + steps_per_chunk, len(time_grid))
        chunks.append((time_grid[i:end], np.arange(i, end, dtype=np.int64)))
    return chunks


def _merge_candidates(
    results: list[list[tuple[int, int, float, float]]],
) -> list[tuple[int, int, float, float]]:
    """Merge per-chunk lists, keeping the minimum-distance epoch per pair.

    Ties on distance are broken by the earlier epoch (lexicographic (dist, t)),
    so the merge is deterministic regardless of chunk arrival order
    (``imap_unordered``): float-equal distances at two epochs previously kept
    whichever chunk happened to arrive first.
    """
    best: dict[tuple[int, int], tuple[float, float]] = {}
    for chunk_result in results:
        for idx_i, idx_j, t_jd, dist in chunk_result:
            key = (idx_i, idx_j)
            if key not in best or (dist, t_jd) < (best[key][1], best[key][0]):
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
    positions: np.ndarray | None = None,
    query_radius_au: float | None = None,
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
    positions:
        Optional pre-computed ``(T, N, 3)`` positions array (e.g. N-body output
        or a memmapped cache).  When supplied, workers consume the
        corresponding chunk slice instead of re-propagating.

    Returns
    -------
    list of (idx_i, idx_j, best_t_jd, best_dist_au)
        One entry per pair whose minimum observed distance is ≤ *threshold_au*.
    """
    nw = resolve_n_workers(n_workers)
    chunks = _make_chunks(time_grid, chunk_size_days)
    qr = query_radius_au if query_radius_au is not None else threshold_au
    tasks = [
        (chunk_times, chunk_indices, threshold_au, leaf_size, qr)
        for chunk_times, chunk_indices in chunks
    ]

    # Memmap-backed trajectories (cache hits) are passed as (filename, shape,
    # dtype) so each worker re-opens its own read-only map.  Zarr-backed caches
    # come wrapped in a TrajectoryView; we pass the underlying DirectoryStore
    # path so each worker opens its own LRU-cached view.  Pickling either
    # object directly via initargs would force a full copy per worker — fatal
    # at the production scale of ~30 GB raw.
    positions_memmap: tuple[str, tuple[int, int, int], str] | None = None
    positions_zarr_path: str | None = None
    positions_inmem: np.ndarray | None = None
    if positions is not None:
        # Avoid importing TrajectoryView at module import time so this file
        # has no hard zarr dependency for the streaming-Kepler path.
        from src.propagate.cache import TrajectoryView

        if isinstance(positions, TrajectoryView):
            store_path = positions.zarr_path
            if store_path is not None:
                positions_zarr_path = str(store_path)
            else:
                positions_inmem = positions  # in-memory zarr / unsupported store
        elif isinstance(positions, np.memmap) and positions.filename:
            positions_memmap = (
                str(positions.filename),
                tuple(positions.shape),  # type: ignore[assignment]
                positions.dtype.str,
            )
        else:
            positions_inmem = positions

    if positions is None:
        mode = "streaming"
    elif positions_zarr_path is not None:
        mode = f"zarr:{Path(positions_zarr_path).name}"
    elif positions_memmap is not None:
        mode = f"memmap:{Path(positions_memmap[0]).name}"
    else:
        mode = "in-memory"

    # The prefilter pair list can grow to tens of MB at moderate N (e.g. 2 000
    # asteroids → ~1.3 M pairs). Passing such an array through Pool initargs
    # pickles it once per worker and has been observed to deadlock the
    # forkserver. Spill to a tempfile and let each worker mmap it instead.
    pairs_path: str | None = None
    pairs_inmem: np.ndarray | None = pairs
    pairs_tmp_dir: tempfile.TemporaryDirectory[str] | None = None
    if pairs is not None and pairs.nbytes > 1_000_000:  # >1 MB → spill
        pairs_tmp_dir = tempfile.TemporaryDirectory(prefix="gaia_pairs_")
        pairs_path = str(Path(pairs_tmp_dir.name) / "pairs.npy")
        np.save(pairs_path, pairs)
        pairs_inmem = None
        logger.info(
            "Pairs array (%.1f MB) spilled to tempfile %s for worker sharing",
            pairs.nbytes / 1e6,
            pairs_path,
        )

    # The elements DataFrame is passed to every worker. At production scale the
    # full numbered population is ~32 MB, and pickling it once per worker through
    # the forkserver initargs pipe deadlocks the forkserver (observed hanging the
    # scan at ~5/35 chunks with all workers idle in futex_wait). Spill it to a
    # Parquet tempfile and let each worker read it locally, mirroring the pairs
    # spill above. Workers only need `elements` for propagation in the streaming
    # path; the positions (cache) path ignores it, but spilling is correct for both.
    elements_path: str | None = None
    elements_inmem: pl.DataFrame | None = elements
    elements_tmp_dir: tempfile.TemporaryDirectory[str] | None = None
    if elements.estimated_size() > _ELEMENTS_SPILL_BYTES:  # spill above threshold
        elements_tmp_dir = tempfile.TemporaryDirectory(prefix="gaia_elements_")
        elements_path = str(Path(elements_tmp_dir.name) / "elements.parquet")
        elements.write_parquet(elements_path)
        elements_inmem = None
        logger.info(
            "Elements DataFrame (%.1f MB) spilled to tempfile %s for worker sharing",
            elements.estimated_size() / 1e6,
            elements_path,
        )

    logger.info(
        "Parallel scan: %d workers | %d chunks (~%.0f days each) | %d total steps | positions=%s",
        nw,
        len(chunks),
        chunk_size_days,
        len(time_grid),
        mode,
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
        # Polars/Arrow keeps its own Rust (rayon) thread pool that OMP/BLAS vars
        # do NOT control. Pin it to 1 per worker: 20 workers × a multi-thread
        # polars pool oversubscribes, and the pool is not fork-safe.
        "POLARS_MAX_THREADS",
    ):
        os.environ.setdefault(_var, "1")

    # Start method = "spawn". forkserver preloads ``__main__`` (run_pipeline),
    # which imports numpy/scipy/polars and initialises their native thread pools
    # in the forkserver process; forked workers then inherit those pools with
    # their locks held by threads that no longer exist → every worker deadlocks
    # in futex_wait the moment it touches polars/BLAS (observed at production
    # scale: the scan froze at ~5-6/35 chunks, all workers idle in futex_wait,
    # CPU ~0 %). "spawn" starts each worker as a fresh interpreter that imports
    # the libraries itself, so no locked thread pool is inherited. initargs are
    # small now (elements/pairs/positions are passed by path, not value), so
    # spawn's per-worker pickling cost is negligible against a multi-hour scan.
    ctx = mp.get_context("spawn")
    with ctx.Pool(
        processes=nw,
        initializer=_init_worker,
        initargs=(
            elements_inmem,
            elements_path,
            pairs_inmem,
            pairs_path,
            positions_inmem,
            positions_memmap,
            positions_zarr_path,
        ),
    ) as pool:
        all_results: list[list[tuple[int, int, float, float]]] = []
        with tqdm(total=len(chunks), desc="Scanning chunks", unit="chunk") as pbar:
            for result in pool.imap_unordered(_scan_chunk, tasks):
                all_results.append(result)
                pbar.update(1)

    candidates = _merge_candidates(all_results)
    logger.info("Parallel scan complete: %d candidate pairs", len(candidates))
    if pairs_tmp_dir is not None:
        pairs_tmp_dir.cleanup()
    if elements_tmp_dir is not None:
        elements_tmp_dir.cleanup()
    return candidates
