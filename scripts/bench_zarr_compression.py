"""One-off bench: compare on-disk size and IO speed between the memmap and
zarr backends of the trajectory cache at production-scale payload sizes.

Run inside Docker:

    docker compose run --rm pipeline python -m scripts.bench_zarr_compression \
        --n-asteroids 20000 --years 3 --step-hours 12 --n-workers 4
"""

from __future__ import annotations

import argparse
import logging
import shutil
import time
from pathlib import Path

import numpy as np

from src.ingest.mpcorb import parse_mpcorb
from src.ingest.mpcorb_archive import discover_snapshots, select_for_window
from src.propagate.cache import (
    _zarr_dir_size_bytes,
    build_cache_key,
    load_or_compute_trajectory,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
)
log = logging.getLogger("bench_zarr")


def _make_grid(t_start_jd: float, years: float, step_hours: float) -> np.ndarray:
    step_days = step_hours / 24.0
    n = int(round(years * 365.25 / step_days)) + 1
    return np.array([t_start_jd + k * step_days for k in range(n)], dtype=np.float64)


def _stream_read_time(positions, label: str) -> tuple[float, float]:
    """Step-by-step read: emulates the KD-tree scan access pattern."""
    t0 = time.monotonic()
    checksum = 0.0
    for k in range(positions.shape[0]):
        slab = np.asarray(positions[k])
        # touch the data to defeat lazy paging
        checksum += float(slab[0, 0])
    dt = time.monotonic() - t0
    log.info(
        "[%s] step-by-step read %d slabs in %.2fs (%.1f slabs/s)",
        label,
        positions.shape[0],
        dt,
        positions.shape[0] / dt,
    )
    return dt, checksum


def _stream_read_lru(zarr_path: str, n_steps: int, label: str) -> tuple[float, float]:
    """Step-by-step read via the production worker path (LRUStoreCache)."""
    from src.propagate.cache import open_trajectory_for_worker

    view = open_trajectory_for_worker(zarr_path)
    t0 = time.monotonic()
    checksum = 0.0
    for k in range(n_steps):
        slab = np.asarray(view[k])
        checksum += float(slab[0, 0])
    dt = time.monotonic() - t0
    log.info(
        "[%s] step-by-step LRU read %d slabs in %.2fs (%.1f slabs/s)",
        label,
        n_steps,
        dt,
        n_steps / dt,
    )
    return dt, checksum


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--n-asteroids", type=int, default=20_000)
    p.add_argument("--years", type=float, default=3.0)
    p.add_argument("--step-hours", type=float, default=12.0)
    p.add_argument(
        "--dt-hours",
        type=float,
        default=1.0,
        help="REBOUND WHFast integrator step (independent of grid)",
    )
    p.add_argument("--semimajor-min-au", type=float, default=2.0)
    p.add_argument("--semimajor-max-au", type=float, default=3.5)
    p.add_argument("--n-workers", type=int, default=4)
    p.add_argument("--cache-dir", type=str, default="data/cache/bench_zarr")
    p.add_argument(
        "--keep-cache",
        action="store_true",
        help="Don't wipe the cache dir before running (skip computation if cached)",
    )
    args = p.parse_args()

    cache_dir = Path(args.cache_dir)
    if not args.keep_cache and cache_dir.exists():
        log.info("Wiping existing %s", cache_dir)
        shutil.rmtree(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    # --- Load elements ---
    snapshots = discover_snapshots(Path("data/raw"))
    if not snapshots:
        log.error("No MPCORB snapshots found under data/raw")
        return 1
    # Pick snapshot whose epoch is closest to the window centre.
    # Window: 2015-06-01 ± years/2.
    centre_jd = 2_457_174.5  # JD TDB for 2015-06-01
    t_start_jd = centre_jd - args.years * 365.25 / 2
    snap = select_for_window(snapshots, t_start_jd, t_start_jd + args.years * 365.25)
    log.info("Using snapshot %s", snap.path.name)

    elements = parse_mpcorb(
        snap.path,
        only_numbered=True,
        semimajor_min_au=args.semimajor_min_au,
        semimajor_max_au=args.semimajor_max_au,
    )
    log.info(
        "Loaded %d asteroids from MPCORB after a∈[%.2f,%.2f]",
        len(elements),
        args.semimajor_min_au,
        args.semimajor_max_au,
    )
    if len(elements) > args.n_asteroids:
        elements = elements.head(args.n_asteroids)
    log.info("Using subset of %d asteroids", len(elements))

    grid = _make_grid(t_start_jd, args.years, args.step_hours)
    log.info("Grid: %d steps Δt=%.2fh from JD %.3f", len(grid), args.step_hours, t_start_jd)

    raw_bytes = len(grid) * len(elements) * 3 * 4
    log.info(
        "Uncompressed payload: %.2f GB (T=%d N=%d float32)",
        raw_bytes / 1e9,
        len(grid),
        len(elements),
    )

    rebound_kwargs = dict(
        include_planets=["sun", "jupiter", "saturn"],
        include_major_asteroids=False,
        integrator="whfast",
        dt_days=args.dt_hours / 24.0,
        n_workers=args.n_workers,
    )

    cache_key = build_cache_key(
        snapshot_sha=snap.path,
        time_grid=grid,
        method="rebound",
        rebound_kwargs=rebound_kwargs,
        n_asteroids=len(elements),
    )
    log.info("Cache key: %s", cache_key)

    # ============================================================
    # Format 1: memmap
    # ============================================================
    log.info("=" * 60)
    log.info("RUN 1: cache_format=memmap")
    log.info("=" * 60)
    t0 = time.monotonic()
    pos_mm = load_or_compute_trajectory(
        elements=elements,
        time_grid=grid,
        cache_dir=cache_dir,
        cache_key=cache_key,
        rebound_kwargs=rebound_kwargs,
        cache_format="memmap",
    )
    mm_write_time = time.monotonic() - t0
    mm_path = cache_dir / f"trajectory_{cache_key}.npy"
    mm_bytes = mm_path.stat().st_size
    log.info("[memmap] write+integrate: %.2fs, on-disk %.2f GB", mm_write_time, mm_bytes / 1e9)

    # ============================================================
    # Format 2: zarr
    # ============================================================
    log.info("=" * 60)
    log.info("RUN 2: cache_format=zarr")
    log.info("=" * 60)
    t0 = time.monotonic()
    pos_zarr = load_or_compute_trajectory(
        elements=elements,
        time_grid=grid,
        cache_dir=cache_dir,
        cache_key=cache_key,
        rebound_kwargs=rebound_kwargs,
        cache_format="zarr",
    )
    zarr_write_time = time.monotonic() - t0
    zarr_dir = cache_dir / f"trajectory_{cache_key}.zarr"
    zarr_bytes = _zarr_dir_size_bytes(zarr_dir)
    log.info("[zarr] write+integrate: %.2fs, on-disk %.2f GB", zarr_write_time, zarr_bytes / 1e9)

    # ============================================================
    # Read benches
    # ============================================================
    log.info("=" * 60)
    log.info("READ BENCHMARK (sequential slab access)")
    log.info("=" * 60)
    mm_read_time, mm_checksum = _stream_read_time(pos_mm, "memmap")
    zarr_read_time, zarr_checksum = _stream_read_time(pos_zarr, "zarr (raw)")
    # Production read path: workers use open_trajectory_for_worker which wraps
    # the store in LRUStoreCache. Chunks decompressed once per worker time window.
    zarr_dir = cache_dir / f"trajectory_{cache_key}.zarr"
    zarr_lru_time, _ = _stream_read_lru(str(zarr_dir), len(grid), "zarr+LRU")

    # ============================================================
    # Bit-exact?
    # ============================================================
    log.info("=" * 60)
    log.info("CORRECTNESS")
    log.info("=" * 60)
    # Compare a few slabs spread through the time axis to keep memory bounded.
    sample_idxs = np.linspace(0, len(grid) - 1, num=20).astype(int)
    max_diff = 0.0
    bit_exact = True
    for k in sample_idxs:
        a = np.asarray(pos_mm[k])
        b = np.asarray(pos_zarr[k])
        if not np.array_equal(a, b):
            bit_exact = False
            max_diff = max(max_diff, float(np.max(np.abs(a - b))))
    if bit_exact:
        log.info("memmap[k] == zarr[k] bit-exact at %d sampled timesteps ✓", len(sample_idxs))
    else:
        log.warning("Mismatch detected: max abs diff over samples = %.3e", max_diff)

    # ============================================================
    # Summary
    # ============================================================
    ratio = raw_bytes / max(zarr_bytes, 1)
    log.info("=" * 60)
    log.info("SUMMARY  (N=%d × T=%d, raw=%.2f GB)", len(elements), len(grid), raw_bytes / 1e9)
    log.info("=" * 60)
    log.info(
        "  memmap on-disk:    %8.2f GB   write %7.1fs   read         %6.2fs",
        mm_bytes / 1e9,
        mm_write_time,
        mm_read_time,
    )
    log.info(
        "  zarr on-disk:      %8.2f GB   write %7.1fs   read (raw)   %6.2fs",
        zarr_bytes / 1e9,
        zarr_write_time,
        zarr_read_time,
    )
    log.info("                                                 read (LRU)   %6.2fs", zarr_lru_time)
    log.info("  compression ratio (raw/zarr):         %.2fx", ratio)
    log.info(
        "  zarr/memmap write overhead:           %+.1f%%",
        100.0 * (zarr_write_time - mm_write_time) / max(mm_write_time, 1e-9),
    )
    log.info(
        "  zarr+LRU/memmap read overhead:        %+.1f%%",
        100.0 * (zarr_lru_time - mm_read_time) / max(mm_read_time, 1e-9),
    )
    log.info(
        "  bit-exact:                            %s",
        "yes" if bit_exact else "NO (lossy via BitRound)",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
