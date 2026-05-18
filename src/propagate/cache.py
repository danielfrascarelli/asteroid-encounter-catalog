"""Trajectory cache — persist the ``(T, N, 3)`` propagated grid to disk.

The N-body propagation is the dominant cost of repeated pipeline runs.  This
module stores the full trajectory as a raw ``float32`` memory-mapped file
alongside a JSON sidecar manifest.  Subsequent runs whose manifest key matches
load the trajectory in milliseconds (via ``np.memmap``) and share it across
worker processes through OS page-cache mappings — no copies required.

Manifest key
------------
The cache key is derived from inputs that fully determine the trajectory:

* MPCORB snapshot SHA-256 (positions and elements depend on the orbital data)
* time grid bounds + step + length
* propagator method and per-method options (integrator, planets, major
  asteroids flag, dt)
* asteroid count (guard against changes in subset / ordering)

A change in any of the above invalidates the cache.

Files written
-------------
- ``<cache_dir>/trajectory_<key>.npy``       — raw float32, shape (T, N, 3).
- ``<cache_dir>/trajectory_<key>.json``      — manifest sidecar.
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

_CACHE_VERSION = 1  # bump when on-disk format changes


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
    base = Path(cache_dir)
    return base / f"trajectory_{key}.npy", base / f"trajectory_{key}.json"


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


def load_or_compute_trajectory(
    *,
    elements: pl.DataFrame,
    time_grid: np.ndarray,
    cache_dir: str | Path,
    cache_key: str,
    rebound_kwargs: dict[str, Any],
) -> np.ndarray:
    """Return a ``(T, N, 3)`` float32 trajectory, materialising it if necessary.

    Look-up order:

    1. If ``<cache_dir>/trajectory_<key>.{npy,json}`` exist and the manifest
       validates, memory-map the array and return it.
    2. Otherwise compute via :func:`src.propagate.nbody.propagate_grid_nbody`,
       write to disk, then return a memory-mapped view.

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

    Returns
    -------
    np.ndarray
        Either a fresh in-memory array or a ``np.memmap`` (read-only).
    """
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    npy_path, manifest_path = _paths(cache_dir, cache_key)
    expected_shape = (len(time_grid), len(elements), 3)

    if npy_path.is_file() and _validate_manifest(manifest_path, expected_shape):
        size_gb = npy_path.stat().st_size / 1e9
        logger.info("Cache HIT: %s (%.2f GB) — memory-mapping", npy_path, size_gb)
        return np.memmap(npy_path, dtype=np.float32, mode="r", shape=expected_shape)

    logger.info("Cache MISS: computing trajectory → %s", npy_path)
    from src.propagate.nbody import propagate_grid_nbody

    trajectory = propagate_grid_nbody(elements, time_grid, **rebound_kwargs)
    if trajectory.dtype != np.float32:
        trajectory = trajectory.astype(np.float32)

    # Write to a temp file then atomic-rename so a crash mid-write can't leave
    # a half-baked cache that future runs would accept on shape alone.
    tmp = npy_path.with_suffix(".npy.tmp")
    memmap = np.memmap(tmp, dtype=np.float32, mode="w+", shape=expected_shape)
    memmap[:] = trajectory
    memmap.flush()
    del memmap
    tmp.replace(npy_path)

    manifest_path.write_text(
        json.dumps(
            {
                "version": _CACHE_VERSION,
                "shape": list(expected_shape),
                "dtype": "float32",
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
                "n_asteroids": int(len(elements)),
                "n_steps": int(len(time_grid)),
                "computed_at": datetime.now(UTC).isoformat(),
            },
            indent=2,
        )
    )

    logger.info(
        "Cache WRITE: %s (%.2f GB) + manifest",
        npy_path,
        npy_path.stat().st_size / 1e9,
    )
    return np.memmap(npy_path, dtype=np.float32, mode="r", shape=expected_shape)
