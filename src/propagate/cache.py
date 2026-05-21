"""Trajectory cache — persist the ``(T, N, 3)`` propagated grid to disk.

The N-body propagation is the dominant cost of repeated pipeline runs.  This
module stores the full trajectory either as a raw ``float32`` memory-mapped
file (legacy) or as a transposed, chunked, Delta + BitRound + Blosc-zstd-
bit-shuffled Zarr v2 directory (default).  Subsequent runs whose manifest key
matches stream the trajectory back at near-zero cost: a memory-mapped view
for ``memmap``, a chunked decompressed view for ``zarr``.  Both backends
share their underlying files across worker processes via the OS page cache.

Why transpose?
--------------
The KD-tree scan consumes one ``(N, 3)`` time slab at a time, which makes
``chunks=(1, N, 3)`` the natural read-pattern choice.  Empirically, however,
that layout puts uncorrelated asteroid positions adjacent in the chunk's
flat buffer, so Blosc bit-shuffle finds no repeated patterns and compression
caps at ~1.06× — essentially nothing.  Adding ``BitRound`` lifts the ceiling
to ~2× by zero-padding low mantissa bits but no more.

Storing the array transposed as ``(3, N, T)`` puts the per-asteroid, per-coord
time series contiguous in memory.  Trajectories are smooth in *time*, so
``numcodecs.Delta`` followed by Blosc-zstd-bitshuffle yields ~6.7× compression
losslessly with ``BitRound(keepbits=16)`` — max position error 30 µAU ≈ 4.5
km, six orders of magnitude tighter than the 0.05-AU encounter threshold.

Consumers still see the logical ``(T, N, 3)`` interface via
:class:`TrajectoryView`, an adapter that maps ``positions[t]`` and
``positions[k0:k1]`` to the appropriate transposed slice + transpose-back.
Workers wrap the zarr DirectoryStore in an ``LRUStoreCache`` so the
decompressed chunk for a given T-window is reused across consecutive slab
reads — turning the layout's "many chunks per slab" cost from a per-read
hit into a one-time cost per T-window.

Manifest key
------------
The cache key is derived from inputs that fully determine the trajectory:

* MPCORB snapshot SHA-256 (positions and elements depend on the orbital data)
* time grid bounds + step + length
* propagator method and per-method options (integrator, planets, major
  asteroids flag, dt)
* asteroid count (guard against changes in subset / ordering)

A change in any of the above invalidates the cache.  The on-disk *format*
(memmap vs zarr) is intentionally NOT part of the key — caches in different
formats live in side-by-side files (``.npy`` vs ``.zarr/`` directory) and
do not collide.

Files written
-------------
- ``<cache_dir>/trajectory_<key>.npy``      — raw float32, shape (T, N, 3).
- ``<cache_dir>/trajectory_<key>.zarr/``    — Zarr v2 directory store with
                                              shape (3, N, T), chunks
                                              (3, N, T_chunk), filters
                                              [BitRound, Delta], compressor
                                              Blosc(zstd, BITSHUFFLE, 5).
- ``<cache_dir>/trajectory_<key>.json``     — manifest sidecar (format-agnostic).
"""

from __future__ import annotations

import hashlib
import json
import logging
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np
import polars as pl

logger = logging.getLogger(__name__)

# Bumped to 2 when the zarr backend switched from (T, N, 3) chunks=(1, N, 3)
# to transposed (3, N, T) chunks=(3, N, T_chunk) with Delta+BitRound filters.
# v1 zarr caches are still readable as raw arrays but will be re-computed
# because the manifest version check no longer matches.
_CACHE_VERSION = 2

# Time-axis chunk size for the transposed zarr layout.  At T=2193 (3 yr at 12h)
# this gives 9 chunks total; large enough that the per-chunk metadata cost is
# negligible, small enough that a worker's 30-day scan window touches only
# 1–2 chunks (each ~50 MiB raw / ~7 MiB on disk at typical compression).
_DEFAULT_T_CHUNK = 256

# Mantissa bits kept by the BitRound filter.  For values up to ~5 AU (exponent
# bit = 2), keepbits=16 gives max error 2^(2-16) = 6e-5 AU ≈ 9 km — five
# orders of magnitude tighter than the 0.05-AU encounter detection threshold.
_DEFAULT_BITROUND_KEEPBITS = 16

# Default per-worker LRU cache size.  At chunks=(3, N=150k, 256) ≈ 440 MiB raw
# (~60 MiB on disk after compression), 512 MiB fits several decompressed
# chunks resident — enough for a worker to keep its full time window plus
# overlap with neighbouring chunks.
_DEFAULT_WORKER_LRU_BYTES = 512 * 1024 * 1024


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def build_cache_key(
    *,
    snapshot_sha: Path | str,
    time_grid: np.ndarray,
    method: str,
    rebound_kwargs: dict[str, Any],
    n_asteroids: int,
) -> str:
    """Compute a short content-addressed key for a trajectory.

    Parameters
    ----------
    snapshot_sha:
        Either the SHA-256 string or a Path to the MPCORB.DAT file.
    time_grid:
        The JD TDB grid being propagated.
    method:
        ``"kepler"`` or ``"rebound"``.
    rebound_kwargs:
        Options forwarded to the N-body propagator.
    n_asteroids:
        Number of asteroids in the subset (also folded into the key).

    Returns
    -------
    str
        16-character hex digest unique to this combination.
    """
    if isinstance(snapshot_sha, Path):
        sha = _sha256_file(snapshot_sha)
    else:
        sha = str(snapshot_sha)

    payload = {
        "snapshot_sha": sha,
        "t_start": float(time_grid[0]),
        "t_end": float(time_grid[-1]),
        "t_n": int(len(time_grid)),
        "t_step": float(time_grid[1] - time_grid[0]) if len(time_grid) > 1 else 0.0,
        "method": method,
        "rebound": {
            k: rebound_kwargs.get(k)
            for k in (
                "integrator",
                "dt_days",
                "include_planets",
                "include_major_asteroids",
            )
        },
        "n_asteroids": int(n_asteroids),
        "version": _CACHE_VERSION,
    }
    blob = json.dumps(payload, sort_keys=True, default=str).encode()
    return hashlib.sha256(blob).hexdigest()[:16]


def _paths(cache_dir: str | Path, key: str) -> tuple[Path, Path]:
    """Return (npy_path, manifest_path) for the legacy memmap backend."""
    base = Path(cache_dir)
    return base / f"trajectory_{key}.npy", base / f"trajectory_{key}.json"


def _zarr_paths(cache_dir: str | Path, key: str) -> tuple[Path, Path]:
    """Return (zarr_dir, manifest_path) for the Zarr backend."""
    base = Path(cache_dir)
    return base / f"trajectory_{key}.zarr", base / f"trajectory_{key}.json"


def _zarr_dir_size_bytes(zarr_dir: Path) -> int:
    """Sum on-disk size of all files under a Zarr DirectoryStore."""
    return sum(f.stat().st_size for f in zarr_dir.rglob("*") if f.is_file())


def _validate_manifest(manifest_path: Path, expected_shape: tuple[int, int, int]) -> bool:
    """Return True if the manifest exists and matches expectations."""
    if not manifest_path.is_file():
        return False
    try:
        meta = json.loads(manifest_path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        logger.warning("Cache manifest %s unreadable (%s) — recomputing", manifest_path, exc)
        return False
    if meta.get("version") != _CACHE_VERSION:
        logger.info(
            "Cache version mismatch (%s vs %s) — recomputing",
            meta.get("version"),
            _CACHE_VERSION,
        )
        return False
    if tuple(meta.get("shape", ())) != expected_shape:
        logger.info(
            "Cache shape mismatch %s vs %s — recomputing",
            tuple(meta.get("shape", ())),
            expected_shape,
        )
        return False
    return True


def _write_manifest(
    manifest_path: Path,
    *,
    expected_shape: tuple[int, int, int],
    cache_key: str,
    rebound_kwargs: dict[str, Any],
    time_grid: np.ndarray,
    fmt: str,
) -> None:
    manifest_path.write_text(
        json.dumps(
            {
                "version": _CACHE_VERSION,
                "shape": list(expected_shape),
                "dtype": "float32",
                "format": fmt,
                "cache_key": cache_key,
                "rebound_kwargs": {
                    k: rebound_kwargs.get(k)
                    for k in (
                        "integrator",
                        "dt_days",
                        "include_planets",
                        "include_major_asteroids",
                    )
                },
                "t_start_jd": float(time_grid[0]),
                "t_end_jd": float(time_grid[-1]),
                "n_asteroids": int(expected_shape[1]),
                "n_steps": int(len(time_grid)),
                "computed_at": datetime.now(UTC).isoformat(),
            },
            indent=2,
        )
    )


# ---------------------------------------------------------------------------
# TrajectoryView — (T, N, 3) facade over a (3, N, T) zarr array
# ---------------------------------------------------------------------------


class TrajectoryView:
    """Adapter that presents the logical ``(T, N, 3)`` interface over the
    transposed ``(3, N, T)`` on-disk zarr layout.

    The transposed storage layout is what unlocks the Delta filter's ability
    to exploit time-axis smoothness for compression — but every consumer of
    the cache (KD-tree scan, refinement) was written against the natural
    ``(T, N, 3)`` shape.  This class makes the layout choice transparent.

    Decompressed-chunk cache
    ------------------------
    The KD-tree scan reads slabs sequentially within a worker's time window.
    Without help, ``view[t]`` followed by ``view[t+1]`` would decompress the
    same underlying ``(3, N, T_chunk)`` zarr chunk **twice** — zarr's
    ``LRUStoreCache`` caches the *compressed* bytes but still calls Blosc on
    every access.  We keep the last ``decompressed_chunk_lru`` chunks
    resident as decompressed ndarrays so consecutive slab reads inside a
    T-chunk window decompress only once.  At production scale (N=150k,
    T_chunk=256) each cached chunk is ~440 MiB raw; the default cache of 2
    chunks (~880 MiB per worker) fits a worker's 30-day scan window plus
    overlap with neighbouring chunks.

    Supported indexing:

    * ``view[t]`` for an integer ``t`` → returns the ``(N, 3)`` slab at
      time index ``t``.
    * ``view[start:stop]`` → returns the ``(k, N, 3)`` slab covering the
      slice on the time axis.
    """

    def __init__(self, z: Any, decompressed_chunk_lru: int = 2) -> None:
        if z.ndim != 3 or z.shape[0] != 3:
            raise ValueError(f"TrajectoryView expects on-disk shape (3, N, T); got {z.shape}")
        self._z = z
        # Logical shape exposed to consumers
        self.shape: tuple[int, int, int] = (int(z.shape[2]), int(z.shape[1]), 3)
        self.dtype = np.dtype(z.dtype)
        self.ndim = 3
        # Decompressed-chunk LRU; key = chunk-along-T index, value = (3, N, T_chunk) array.
        self._t_chunk_size: int = int(z.chunks[2])
        self._chunk_lru_max: int = max(1, int(decompressed_chunk_lru))
        self._chunk_cache: dict[int, np.ndarray] = {}

    def _get_t_chunk(self, ci: int) -> np.ndarray:
        """Return the decompressed (3, N, k) array for T-chunk index *ci*.

        Inserts into the LRU on miss; evicts the oldest entry past
        ``_chunk_lru_max``.  Python ``dict`` preserves insertion order, which
        makes the "oldest" entry the first key.
        """
        cached = self._chunk_cache.get(ci)
        if cached is not None:
            # Refresh LRU position: re-insert at the end.
            del self._chunk_cache[ci]
            self._chunk_cache[ci] = cached
            return cached
        t_start = ci * self._t_chunk_size
        t_end = min(t_start + self._t_chunk_size, self.shape[0])
        chunk = np.asarray(self._z[:, :, t_start:t_end])
        self._chunk_cache[ci] = chunk
        # Evict the oldest entries until we're back at the cap.
        while len(self._chunk_cache) > self._chunk_lru_max:
            oldest_key = next(iter(self._chunk_cache))
            if oldest_key == ci:
                break  # never evict what we just put in
            del self._chunk_cache[oldest_key]
        return chunk

    def __getitem__(self, idx: Any) -> np.ndarray:
        if isinstance(idx, (int, np.integer)):
            t = int(idx)
            if t < 0:
                t += self.shape[0]
            ci = t // self._t_chunk_size
            local_t = t - ci * self._t_chunk_size
            chunk = self._get_t_chunk(ci)  # (3, N, k)
            # chunk[:, :, local_t] is (3, N); transpose to (N, 3).
            return np.ascontiguousarray(chunk[:, :, local_t].T)
        if isinstance(idx, slice):
            start, stop, step = idx.indices(self.shape[0])
            if step != 1:
                # Uncommon — fall back to direct zarr slicing.
                sl = np.asarray(self._z[:, :, idx])  # (3, N, k)
                return np.ascontiguousarray(sl.transpose(2, 1, 0))
            if stop <= start:
                return np.empty((0, self.shape[1], 3), dtype=self.dtype)
            ci_lo = start // self._t_chunk_size
            ci_hi = (stop - 1) // self._t_chunk_size
            if ci_lo == ci_hi:
                chunk = self._get_t_chunk(ci_lo)
                lo = start - ci_lo * self._t_chunk_size
                hi = stop - ci_lo * self._t_chunk_size
                # chunk[:, :, lo:hi] → (3, N, k); transpose → (k, N, 3).
                return np.ascontiguousarray(chunk[:, :, lo:hi].transpose(2, 1, 0))
            # Slice spans multiple T-chunks: assemble piece-wise so we still
            # benefit from the decompressed-chunk LRU.
            parts: list[np.ndarray] = []
            for ci in range(ci_lo, ci_hi + 1):
                chunk = self._get_t_chunk(ci)
                chunk_t_start = ci * self._t_chunk_size
                chunk_t_end = chunk_t_start + chunk.shape[2]
                lo = max(start, chunk_t_start) - chunk_t_start
                hi = min(stop, chunk_t_end) - chunk_t_start
                parts.append(chunk[:, :, lo:hi].transpose(2, 1, 0))
            return np.ascontiguousarray(np.concatenate(parts, axis=0))
        raise TypeError(
            f"TrajectoryView only supports int or slice indexing on the time axis; "
            f"got {type(idx).__name__}"
        )

    @property
    def zarr_path(self) -> str | None:
        """Filesystem path of the underlying DirectoryStore, if any."""
        store = self._z.store
        # If wrapped in LRUStoreCache, traverse to the inner store
        inner = getattr(store, "_store", store)
        return getattr(inner, "path", None) or getattr(inner, "dir_path", None)

    def __array__(self, dtype: Any = None) -> np.ndarray:
        full = np.asarray(self._z).transpose(2, 1, 0)  # (3,N,T) → (T,N,3)
        if dtype is not None:
            full = full.astype(dtype, copy=False)
        return np.ascontiguousarray(full)


def open_trajectory_for_worker(
    zarr_path: str | Path,
    max_cache_bytes: int = _DEFAULT_WORKER_LRU_BYTES,
) -> TrajectoryView:
    """Open a cached trajectory for read-only worker access with chunk caching.

    Wraps the DirectoryStore in :class:`zarr.storage.LRUStoreCache` so
    consecutive slab reads within the same time-chunk window pay the Blosc
    decompression cost only once.  Returns a :class:`TrajectoryView` adapter
    exposing the logical ``(T, N, 3)`` interface.
    """
    import zarr
    from zarr.storage import DirectoryStore, LRUStoreCache

    store = DirectoryStore(str(zarr_path))
    cached = LRUStoreCache(store, max_size=int(max_cache_bytes))
    z = zarr.open(cached, mode="r")
    return TrajectoryView(z)


# ---------------------------------------------------------------------------
# memmap backend (legacy)
# ---------------------------------------------------------------------------


def _load_or_compute_memmap(
    *,
    elements: pl.DataFrame,
    time_grid: np.ndarray,
    cache_dir: Path,
    cache_key: str,
    rebound_kwargs: dict[str, Any],
) -> np.ndarray:
    npy_path, manifest_path = _paths(cache_dir, cache_key)
    expected_shape = (len(time_grid), len(elements), 3)

    if npy_path.is_file() and _validate_manifest(manifest_path, expected_shape):
        size_gb = npy_path.stat().st_size / 1e9
        logger.info("Cache HIT (memmap): %s (%.2f GB) — memory-mapping", npy_path, size_gb)
        return np.memmap(npy_path, dtype=np.float32, mode="r", shape=expected_shape)

    logger.info("Cache MISS (memmap): computing trajectory → %s", npy_path)
    from src.propagate.nbody import propagate_grid_nbody

    tmp = npy_path.with_suffix(".npy.tmp")
    out_mm = np.memmap(tmp, dtype=np.float32, mode="w+", shape=expected_shape)
    propagate_grid_nbody(elements, time_grid, out=out_mm, **rebound_kwargs)
    out_mm.flush()
    del out_mm
    tmp.replace(npy_path)

    _write_manifest(
        manifest_path,
        expected_shape=expected_shape,
        cache_key=cache_key,
        rebound_kwargs=rebound_kwargs,
        time_grid=time_grid,
        fmt="memmap",
    )
    logger.info(
        "Cache WRITE (memmap): %s (%.2f GB) + manifest",
        npy_path,
        npy_path.stat().st_size / 1e9,
    )
    return np.memmap(npy_path, dtype=np.float32, mode="r", shape=expected_shape)


# ---------------------------------------------------------------------------
# zarr backend (default)
# ---------------------------------------------------------------------------


def _transcode_memmap_to_zarr(
    src_mm: np.memmap,
    z: Any,
    t_chunk: int,
) -> None:
    """Stream a ``(T, N, 3)`` memmap into a ``(3, N, T)`` zarr in T_chunk slabs.

    Each write spans exactly one zarr chunk along the time axis, so the
    Blosc compressor sees a well-formed input shape and the Delta filter
    operates on the per-asteroid, per-coord time series — the contiguous
    layout that makes Delta effective.
    """
    T = src_mm.shape[0]
    for t_start in range(0, T, t_chunk):
        t_end = min(t_start + t_chunk, T)
        slab = np.array(src_mm[t_start:t_end], copy=True)  # (k, N, 3)
        transposed = np.ascontiguousarray(np.transpose(slab, (2, 1, 0)))  # (3, N, k)
        z[:, :, t_start:t_end] = transposed


def _load_or_compute_zarr(
    *,
    elements: pl.DataFrame,
    time_grid: np.ndarray,
    cache_dir: Path,
    cache_key: str,
    rebound_kwargs: dict[str, Any],
    t_chunk: int = _DEFAULT_T_CHUNK,
    bitround_keepbits: int = _DEFAULT_BITROUND_KEEPBITS,
) -> TrajectoryView:
    import shutil

    import zarr
    from numcodecs import BitRound, Blosc, Delta

    zarr_dir, manifest_path = _zarr_paths(cache_dir, cache_key)
    expected_shape = (len(time_grid), len(elements), 3)  # logical
    storage_shape = (3, len(elements), len(time_grid))  # on-disk transposed

    if zarr_dir.is_dir() and _validate_manifest(manifest_path, expected_shape):
        size_gb = _zarr_dir_size_bytes(zarr_dir) / 1e9
        raw_gb = float(np.prod(expected_shape)) * 4 / 1e9
        logger.info(
            "Cache HIT (zarr): %s (%.2f GB on disk, %.2f GB uncompressed, ratio %.2f×)",
            zarr_dir,
            size_gb,
            raw_gb,
            raw_gb / max(size_gb, 1e-9),
        )
        return TrajectoryView(zarr.open(str(zarr_dir), mode="r"))

    logger.info("Cache MISS (zarr): computing trajectory → %s", zarr_dir)
    from src.propagate.nbody import propagate_grid_nbody

    # Integrate into a temp memmap with the natural (T, N, 3) layout the
    # nbody propagator expects, then transcode to the transposed zarr
    # layout.  Peak disk = raw memmap (≈4 × T × N × 3 bytes) + final
    # compressed zarr (≈raw / 6).
    tmp_npy = zarr_dir.with_suffix(".transcode-tmp.npy")
    if tmp_npy.exists():
        tmp_npy.unlink()
    tmp_zarr = zarr_dir.with_suffix(".zarr.tmp")
    if tmp_zarr.exists():
        shutil.rmtree(tmp_zarr)

    out_mm = np.memmap(tmp_npy, dtype=np.float32, mode="w+", shape=expected_shape)
    try:
        propagate_grid_nbody(elements, time_grid, out=out_mm, **rebound_kwargs)
        out_mm.flush()

        # Filter order matters: BitRound first zeros out the lowest mantissa
        # bits on the raw positions; Delta second turns rounded values into
        # tiny diffs that share their high bits; Blosc(zstd, BITSHUFFLE)
        # finishes the job.  Reversing the order (Delta-first, BitRound on
        # the diffs) compounds rounding error through the delta chain and
        # blows up the position error.
        chunks_used = (3, expected_shape[1], min(t_chunk, max(1, expected_shape[0])))
        z = zarr.open(
            str(tmp_zarr),
            mode="w",
            shape=storage_shape,
            chunks=chunks_used,
            dtype="float32",
            compressor=Blosc(cname="zstd", clevel=5, shuffle=Blosc.BITSHUFFLE),
            filters=[BitRound(keepbits=bitround_keepbits), Delta(dtype="float32")],
        )
        _transcode_memmap_to_zarr(out_mm, z, t_chunk=chunks_used[2])
    except Exception:
        shutil.rmtree(tmp_zarr, ignore_errors=True)
        raise
    finally:
        del out_mm
        if tmp_npy.exists():
            tmp_npy.unlink()

    if zarr_dir.exists():
        shutil.rmtree(zarr_dir)
    tmp_zarr.replace(zarr_dir)

    _write_manifest(
        manifest_path,
        expected_shape=expected_shape,
        cache_key=cache_key,
        rebound_kwargs=rebound_kwargs,
        time_grid=time_grid,
        fmt="zarr",
    )
    size_gb = _zarr_dir_size_bytes(zarr_dir) / 1e9
    raw_gb = float(np.prod(expected_shape)) * 4 / 1e9
    logger.info(
        "Cache WRITE (zarr): %s (%.2f GB on disk, %.2f GB uncompressed, ratio %.2f×) + manifest",
        zarr_dir,
        size_gb,
        raw_gb,
        raw_gb / max(size_gb, 1e-9),
    )
    return TrajectoryView(zarr.open(str(zarr_dir), mode="r"))


def load_or_compute_trajectory(
    *,
    elements: pl.DataFrame,
    time_grid: np.ndarray,
    cache_dir: str | Path,
    cache_key: str,
    rebound_kwargs: dict[str, Any],
    cache_format: str = "zarr",
) -> np.ndarray | TrajectoryView:
    """Return a logical ``(T, N, 3)`` float32 trajectory, materialising it if necessary.

    Look-up order:

    1. If a cached trajectory in *cache_format* exists and its manifest
       validates against ``(len(time_grid), len(elements), 3)``, open it
       read-only (``np.memmap`` for memmap, :class:`TrajectoryView` over a
       transposed Zarr ``(3, N, T)`` array for zarr) and return it.
    2. Otherwise compute via :func:`src.propagate.nbody.propagate_grid_nbody`
       into a temp memmap, transcode to the chosen on-disk format, and
       return a read-only view.

    Parameters
    ----------
    elements:
        Orbital elements DataFrame.
    time_grid:
        JD TDB epochs.
    cache_dir:
        Directory under which to store / read the trajectory.
    cache_key:
        Identifier built by :func:`build_cache_key`.
    rebound_kwargs:
        Keyword arguments forwarded to the propagator on a cache miss.
    cache_format:
        ``"zarr"`` (default) for the transposed compressed layout, or
        ``"memmap"`` for the legacy raw float32 layout.

    Returns
    -------
    np.ndarray or TrajectoryView
        Read-only view of the trajectory.  Both backends duck-type as
        ``(T, N, 3)`` arrays for slab access (``positions[step_idx]``).
    """
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    fmt = cache_format.lower()
    if fmt == "memmap":
        return _load_or_compute_memmap(
            elements=elements,
            time_grid=time_grid,
            cache_dir=cache_dir,
            cache_key=cache_key,
            rebound_kwargs=rebound_kwargs,
        )
    if fmt == "zarr":
        return _load_or_compute_zarr(
            elements=elements,
            time_grid=time_grid,
            cache_dir=cache_dir,
            cache_key=cache_key,
            rebound_kwargs=rebound_kwargs,
        )
    raise ValueError(f"Unknown cache_format: {cache_format!r} (expected 'zarr' or 'memmap')")
