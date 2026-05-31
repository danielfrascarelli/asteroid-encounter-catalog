"""Measure the Kepler-vs-N-body false-negative rate near the 0.05 AU threshold.

Track C Stage 2. The frozen catalog (and the Kepler-bias note,
``docs/kepler_threshold_bias_paper.md``) can only observe pairs with
Kepler ``dist < 0.05`` AU, so it measures *upward* crossings (Kepler `<0.05`,
N-body `≥0.05`) but **cannot** measure *downward* crossings (Kepler `≥0.05`,
N-body `<0.05`) — those are censored, hence the false-negative rate of the
0.05 AU catalog is unmeasured.

This experiment removes the censoring by detecting at a **wider** Kepler
threshold (0.06 AU) on a body sample, isolating the pairs whose Kepler minimum
falls in the ``[0.05, 0.06)`` band, and re-refining each under full N-body. The
fraction whose true N-body minimum drops below 0.05 AU is the per-band
false-negative rate; combined with the band population it bounds how many real
< 0.05 AU encounters the frozen catalog misses for refinement reasons (distinct
from the prefilter-recall deficit in ``docs/prefilter_recall.md``).

Method
------
1. Sample ``--n-bodies`` numbered MBAs (a∈[1.5,4.0]) from the frozen MPCORB
   snapshot (seeded). The false-negative *rate* is a per-pair property, so a
   random body sample is unbiased for it.
2. Kepler-detect at threshold 0.06 AU (same prefilter/grid as the frozen run).
3. Keep pairs with ``dist_au ∈ [band_lo, band_hi)`` (default [0.05, 0.06)).
4. N-body refine each (``refine_pair_nbody``, ±``window_hours``, IAS15-grade)
   and count how many have N-body minimum < 0.05 AU (downward crossings).

Usage
-----
    docker compose run --rm pipeline python -m \
        scripts.validate.measure_threshold_false_negatives --n-bodies 6000
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import multiprocessing as mp
import time
from pathlib import Path

import polars as pl
from astropy.time import Time

from src.detect.pipeline import detect_encounters
from src.ingest.mpcorb import parse_mpcorb
from src.ingest.mpcorb_archive import discover_snapshots, select_for_window
from src.propagate.grid import make_time_grid
from src.utils.config import load_config

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
logger = logging.getLogger(__name__)

_AU_KM = 1.495_978_707e8
_ELEM_KEYS = ("a_au", "e", "i_deg", "Omega_deg", "omega_deg", "M_deg", "epoch_jd")


def _wilson_ci(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    denom = 1.0 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = (z / denom) * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (max(0.0, centre - half), min(1.0, centre + half))


# Worker globals (set per process) so the elements lookup is not re-pickled per task.
_ELEM_BY_NUMBER: dict[int, dict] = {}
_WINDOW_HOURS = 12.0


def _init_worker(elem_by_number: dict[int, dict], window_hours: float) -> None:
    global _ELEM_BY_NUMBER, _WINDOW_HOURS
    _ELEM_BY_NUMBER = elem_by_number
    _WINDOW_HOURS = window_hours


def _refine_one(task: tuple[int, int, float]) -> dict:
    """Refine one band pair under N-body. Returns dict with nbody distance."""
    from scripts.validate.refine_pair_nbody import refine_pair_nbody

    n1, n2, t_center = task
    e1 = dict(_ELEM_BY_NUMBER[n1], number=n1)
    e2 = dict(_ELEM_BY_NUMBER[n2], number=n2)
    try:
        res = refine_pair_nbody(
            elements_1=e1,
            elements_2=e2,
            t_center_jd=t_center,
            window_hours=_WINDOW_HOURS,
            sample_dt_seconds=60.0,
            include_planets=("sun", "jupiter", "saturn"),
            include_major_asteroids=True,
        )
        # near-boundary: the true minimum may lie outside the window if the
        # refined epoch sits within ~1 sample step of either edge.
        near_boundary = abs(res.t_min_nbody_jd - t_center) >= (_WINDOW_HOURS / 24.0) - 1e-3
        return {
            "number_1": n1,
            "number_2": n2,
            "dist_au_nbody": float(res.dist_au_nbody),
            "t_min_nbody_jd": float(res.t_min_nbody_jd),
            "near_boundary": bool(near_boundary),
            "converged": bool(res.converged),
            "error": None,
        }
    except Exception as exc:  # noqa: BLE001 — record and continue
        return {
            "number_1": n1,
            "number_2": n2,
            "dist_au_nbody": float("nan"),
            "t_min_nbody_jd": float("nan"),
            "near_boundary": False,
            "converged": False,
            "error": str(exc),
        }


def main() -> int:
    parser = argparse.ArgumentParser(description="Measure threshold false-negative rate")
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--n-bodies", type=int, default=6000)
    parser.add_argument("--band-lo", type=float, default=0.05)
    parser.add_argument("--band-hi", type=float, default=0.06)
    parser.add_argument(
        "--max-pairs",
        type=int,
        default=8000,
        help="cap on band pairs to N-body refine (random sub-sample if exceeded)",
    )
    parser.add_argument("--window-hours", type=float, default=12.0)
    parser.add_argument("--n-workers", type=int, default=0, help="0 = cpu_count-2")
    parser.add_argument("--out-dir", default="data/output/kepler_false_negatives")
    args = parser.parse_args()

    cfg = load_config(args.config)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    seed = cfg.run.seed

    tw = cfg.time_window
    t_start = Time(tw.start, scale=tw.scale).tdb.jd
    t_end = Time(tw.end, scale=tw.scale).tdb.jd
    snap = select_for_window(discover_snapshots(Path(cfg.paths.raw)), t_start, t_end)
    logger.info("Snapshot %s", snap.path.name)

    sub = cfg.subset
    elements = parse_mpcorb(
        snap.path,
        only_numbered=True,
        semimajor_min_au=sub.semimajor_axis_au.min,
        semimajor_max_au=sub.semimajor_axis_au.max,
    )
    if args.n_bodies < len(elements):
        elements = elements.sample(n=args.n_bodies, seed=seed).sort("number")
    n_bodies = len(elements)
    logger.info("Sampled %d bodies (seed=%d)", n_bodies, seed)

    # --- Kepler detection at the widened threshold (band_hi) ---
    det = cfg.detection
    fine_step_hours = cfg.propagation.time_step_hours
    coarse_step_hours = cfg.propagation.coarse_step_hours or fine_step_hours
    grid = make_time_grid(t_start, t_end, step_hours=coarse_step_hours)
    v_max_au_day = det.max_relative_velocity_km_s * 86_400.0 / _AU_KM
    query_radius = args.band_hi + v_max_au_day * (coarse_step_hours / 24.0)

    logger.info("Kepler detection at threshold %.3f AU over %d bodies…", args.band_hi, n_bodies)
    t0 = time.monotonic()
    cat = detect_encounters(
        elements,
        grid,
        threshold_au=args.band_hi,
        semimajor_diff_max_au=det.prefilter.semimajor_diff_max_au,
        inclination_diff_max_deg=det.prefilter.inclination_diff_max_deg,
        leaf_size=det.kdtree.leaf_size,
        fine_step_seconds=det.refinement.fine_time_step_seconds,
        window_hours=det.refinement.window_hours,
        prefilter_enabled=True,
        refinement_enabled=True,
        n_workers=cfg.parallel.n_workers if cfg.parallel.enabled else 1,
        chunk_size_days=cfg.parallel.chunk_size_days,
        positions=None,
        query_radius_au=query_radius,
    )
    logger.info(
        "Detection: %d pairs <%.3f AU in %.1fs", len(cat), args.band_hi, time.monotonic() - t0
    )

    band = cat.filter((pl.col("dist_au") >= args.band_lo) & (pl.col("dist_au") < args.band_hi))
    n_band_total = len(band)
    n_below = len(cat.filter(pl.col("dist_au") < args.band_lo))
    logger.info(
        "Band [%.3f,%.3f): %d pairs (catalog <%.3f: %d)",
        args.band_lo,
        args.band_hi,
        n_band_total,
        args.band_lo,
        n_below,
    )
    if n_band_total == 0:
        logger.warning("No band pairs — increase --n-bodies.")
        return 1

    if n_band_total > args.max_pairs:
        band = band.sample(n=args.max_pairs, seed=seed)
        logger.info("Sub-sampled band to %d pairs for N-body refinement", args.max_pairs)

    # --- N-body refine each band pair ---
    elem_by_number = {
        int(r["number"]): {k: float(r[k]) for k in _ELEM_KEYS}
        for r in elements.select(["number", *_ELEM_KEYS]).iter_rows(named=True)
    }
    tasks = [
        (int(r["number_1"]), int(r["number_2"]), float(r["jd_tdb"]))
        for r in band.iter_rows(named=True)
    ]
    n_workers = args.n_workers or max(1, (mp.cpu_count() or 2) - 2)
    logger.info(
        "N-body refining %d band pairs on %d workers (±%.0f h)…",
        len(tasks),
        n_workers,
        args.window_hours,
    )
    t0 = time.monotonic()
    results = []  # type: list[dict]
    # 'spawn' (not the default 'fork'): the detection step above used threaded
    # BLAS, and forking after that deadlocks rebound/numpy workers (0% CPU hang).
    ctx = mp.get_context("spawn")
    with ctx.Pool(
        n_workers, initializer=_init_worker, initargs=(elem_by_number, args.window_hours)
    ) as pool:
        for i, r in enumerate(pool.imap_unordered(_refine_one, tasks, chunksize=8), 1):
            results.append(r)
            if i % 500 == 0:
                logger.info("  refined %d/%d", i, len(tasks))
    logger.info("N-body refinement done in %.1fs", time.monotonic() - t0)

    res_df = pl.DataFrame(results)
    band_keyed = band.rename({"dist_au": "dist_au_kepler", "jd_tdb": "t_min_kepler_jd"})
    merged = band_keyed.join(res_df, on=["number_1", "number_2"], how="inner")
    merged = merged.with_columns(
        (pl.col("dist_au_nbody") - pl.col("dist_au_kepler")).alias("delta_dist_au")
    )

    ok = merged.filter(pl.col("error").is_null() & pl.col("dist_au_nbody").is_finite())
    n_ref = len(ok)
    n_failed = len(merged) - n_ref
    n_cross_down = len(ok.filter(pl.col("dist_au_nbody") < args.band_lo))
    rate = n_cross_down / n_ref if n_ref else float("nan")
    ci_lo, ci_hi = _wilson_ci(n_cross_down, n_ref)

    merged.write_parquet(out_dir / "band_refined.parquet")
    summary = {
        "snapshot": snap.path.name,
        "n_bodies": n_bodies,
        "band": [args.band_lo, args.band_hi],
        "window_hours": args.window_hours,
        "n_below_threshold_kepler": n_below,
        "n_band_pairs_total": n_band_total,
        "n_band_pairs_refined": n_ref,
        "n_failed": n_failed,
        "n_crossing_down": n_cross_down,
        "false_negative_rate_in_band": rate,
        "rate_ci95": [ci_lo, ci_hi],
        "interpretation": (
            "Of Kepler pairs in [band_lo, band_hi), this fraction truly fall "
            "below band_lo under N-body — i.e. real <band_lo encounters the "
            "Kepler-threshold catalog censors. Combine with the upward-crossing "
            "rate (docs/kepler_threshold_bias_paper.md) for the full near-"
            "threshold confusion matrix."
        ),
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2))

    print("\n" + "=" * 70)
    print(
        f"KEPLER THRESHOLD FALSE-NEGATIVES — {n_bodies} bodies, band [{args.band_lo},{args.band_hi})"
    )
    print("=" * 70)
    print(f"Band pairs (Kepler in band):     {n_band_total}")
    print(f"N-body refined:                  {n_ref}  (failed {n_failed})")
    print(f"Crossed DOWN (N-body < {args.band_lo}):  {n_cross_down}")
    print(
        f"False-negative rate in band:     {100 * rate:.2f}%  [95% CI {100*ci_lo:.2f}, {100*ci_hi:.2f}]"
    )
    print("=" * 70)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
