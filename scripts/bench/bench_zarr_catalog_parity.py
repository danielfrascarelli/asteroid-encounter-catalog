"""End-to-end parity check: detect_encounters with memmap vs zarr cache.

Runs the detection pipeline twice on the same subset, once per cache backend,
and compares the resulting encounter catalogs.  Acceptable difference is
bounded by the BitRound precision of the zarr backend (~4,488 km half-quantum
at keepbits=16) — about three orders of magnitude tighter than the 0.05-AU
(~7.48 × 10⁶ km) threshold, which makes it safe for the coarse scan but NOT
for sub-micro-AU validation work.
"""

from __future__ import annotations

import argparse
import logging
import shutil
import time
from pathlib import Path

import numpy as np
import polars as pl

from src.detect.pipeline import detect_encounters
from src.ingest.mpcorb import parse_mpcorb
from src.ingest.mpcorb_archive import discover_snapshots, select_for_window
from src.propagate.cache import build_cache_key
from src.propagate.grid import make_time_grid, propagate_full_grid

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
log = logging.getLogger("parity")


def _run(elements, grid, rebound_kwargs, snap_path, cache_dir, cache_format):
    key = build_cache_key(
        snapshot_sha=snap_path,
        time_grid=grid,
        method="rebound",
        rebound_kwargs=rebound_kwargs,
        n_asteroids=len(elements),
    )
    t0 = time.monotonic()
    positions = propagate_full_grid(
        elements,
        grid,
        method="rebound",
        rebound_kwargs=rebound_kwargs,
        cache_dir=str(cache_dir),
        cache_key=key,
        cache_format=cache_format,
    )
    log.info(
        "[%s] propagation %.1fs  shape=%s",
        cache_format,
        time.monotonic() - t0,
        positions.shape if positions is not None else None,
    )

    # Wider query radius for Strategy A coarse grid (12h × 25 km/s margin).
    v_max_au_per_day = 25.0 * 86_400.0 / 1.495_978_707e8
    widen_au = v_max_au_per_day * (12.0 / 24.0)
    query_radius_au = 0.05 + widen_au

    t0 = time.monotonic()
    cat = detect_encounters(
        elements,
        grid,
        threshold_au=0.05,
        semimajor_diff_max_au=0.5,
        inclination_diff_max_deg=30.0,
        leaf_size=30,
        fine_step_seconds=60.0,
        window_hours=2.0,
        prefilter_enabled=True,
        refinement_enabled=True,
        n_workers=2,
        chunk_size_days=30.0,
        positions=positions,
        query_radius_au=query_radius_au,
        force_kepler_refine=True,
    )
    log.info("[%s] detection %.1fs  %d rows", cache_format, time.monotonic() - t0, len(cat))
    return cat


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--n-asteroids", type=int, default=500)
    p.add_argument("--years", type=float, default=0.5)
    p.add_argument("--cache-dir", default="data/cache/bench_zarr_parity")
    args = p.parse_args()

    cache_dir = Path(args.cache_dir)
    if cache_dir.exists():
        shutil.rmtree(cache_dir)
    cache_dir.mkdir(parents=True)

    centre_jd = 2_457_174.5
    t_start_jd = centre_jd - args.years * 365.25 / 2
    snapshots = discover_snapshots(Path("data/raw"))
    snap = select_for_window(snapshots, t_start_jd, t_start_jd + args.years * 365.25)
    log.info("Snapshot %s", snap.path.name)

    elements = parse_mpcorb(
        snap.path, only_numbered=True, semimajor_min_au=2.0, semimajor_max_au=3.5
    )
    elements = elements.head(args.n_asteroids)
    log.info("Using %d asteroids", len(elements))

    grid = make_time_grid(t_start_jd, t_start_jd + args.years * 365.25, step_hours=12.0)
    log.info("Grid: %d steps Δt=12h", len(grid))

    rebound_kwargs = dict(
        include_planets=["sun", "jupiter", "saturn"],
        include_major_asteroids=False,
        integrator="whfast",
        dt_days=1.0 / 24.0,
        n_workers=2,
    )

    cat_mm = _run(elements, grid, rebound_kwargs, snap.path, cache_dir, "memmap")
    cat_zr = _run(elements, grid, rebound_kwargs, snap.path, cache_dir, "zarr")

    # Compare: same pairs detected?
    def pair_key(row):
        return (row["number_1"], row["number_2"])

    set_mm = set(map(pair_key, cat_mm.iter_rows(named=True)))
    set_zr = set(map(pair_key, cat_zr.iter_rows(named=True)))
    only_mm = set_mm - set_zr
    only_zr = set_zr - set_mm
    common = set_mm & set_zr

    log.info("=" * 60)
    log.info("CATALOG PARITY  (zarr ratio expected ≥ 4×; precision loss << 0.05 AU)")
    log.info("=" * 60)
    log.info("memmap pairs: %d", len(set_mm))
    log.info("zarr   pairs: %d", len(set_zr))
    log.info(
        "common:       %d  (memmap-only %d, zarr-only %d)", len(common), len(only_mm), len(only_zr)
    )

    # Per-pair distance difference
    if common:
        d_mm = {pair_key(r): r["dist_au"] for r in cat_mm.iter_rows(named=True)}
        d_zr = {pair_key(r): r["dist_au"] for r in cat_zr.iter_rows(named=True)}
        diffs = np.array([d_mm[p] - d_zr[p] for p in common])
        log.info(
            "distance Δ:   mean %.2e AU  max |Δ| %.2e AU  std %.2e AU",
            float(diffs.mean()),
            float(np.max(np.abs(diffs))),
            float(diffs.std()),
        )

    # Major bodies sanity (will only show if subset is large enough)
    for n in (1, 2, 4, 10):
        m = cat_mm.filter((pl.col("number_1") == n) | (pl.col("number_2") == n))
        z = cat_zr.filter((pl.col("number_1") == n) | (pl.col("number_2") == n))
        if len(m) or len(z):
            log.info("body (%d):   memmap %d, zarr %d", n, len(m), len(z))

    return 0 if set_mm == set_zr else 2


if __name__ == "__main__":
    raise SystemExit(main())
