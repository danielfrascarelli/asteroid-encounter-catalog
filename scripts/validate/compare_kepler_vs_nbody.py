"""Parallel comparator: Kepler vs N-body refinement on a sample of catalog pairs.

Reads ``data/cache/nbody_validation/sample_1000.parquet`` (produced by
``sample_for_nbody_check.py``), runs :func:`refine_pair_nbody` over each pair
in a multiprocessing pool, and writes
``data/output/kepler_vs_nbody_comparison.parquet`` with the Kepler/N-body
values side-by-side plus orbital metadata for downstream analysis (A.4).

The comparator is a *pure measurement* — it does not modify the frozen catalog;
its sole purpose is to characterise the error of the Kepler refinement step
documented in [FROZEN_RUN.md].

Output schema
-------------
    number_1, number_2, bin_id
    a_1, e_1, i_1, q_1, a_2, e_2, i_2, q_2, delta_a_au
    dist_au_kepler, dist_au_nbody, delta_dist_au
    t_min_kepler_jd, t_min_nbody_jd, delta_t_min_hours
    rel_vel_kepler, rel_vel_nbody, delta_rel_vel_au_day
    nbody_converged, nbody_energy_drift, n_samples
    error_message  (NULL on success, str on failure)

Usage
-----
    docker compose run --rm pipeline python -m scripts.validate.compare_kepler_vs_nbody \\
        --workers 24 --window-hours 6.0 --sample-dt-seconds 60.0
"""

from __future__ import annotations

import argparse
import logging
import multiprocessing as mp
import os
import time
import traceback
from dataclasses import dataclass
from pathlib import Path

import polars as pl
from tqdm import tqdm

from scripts.validate.refine_pair_nbody import refine_pair_nbody

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
logger = logging.getLogger(__name__)

SAMPLE = Path("data/cache/nbody_validation/sample_1000.parquet")
OUT = Path("data/output/kepler_vs_nbody_comparison.parquet")


@dataclass(frozen=True)
class _PairTask:
    """Inputs to refine a single pair (carries everything; no parquet IO in workers)."""

    number_1: int
    number_2: int
    bin_id: str
    jd_tdb_kepler: float
    dist_au_kepler: float
    rel_vel_kepler: float
    a_1: float
    e_1: float
    i_1: float
    Omega_1: float
    omega_1: float
    M_1: float
    epoch_1: float
    a_2: float
    e_2: float
    i_2: float
    Omega_2: float
    omega_2: float
    M_2: float
    epoch_2: float
    q_1: float
    q_2: float
    delta_a_au: float


_WORKER_CONFIG: dict = {}


def _init_worker(window_hours: float, sample_dt_seconds: float,
                 warmup_dt_seconds: float, include_major_asteroids: bool) -> None:
    """Pool initializer: pin one BLAS thread per worker and cache refiner kwargs."""
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    _WORKER_CONFIG["window_hours"] = window_hours
    _WORKER_CONFIG["sample_dt_seconds"] = sample_dt_seconds
    _WORKER_CONFIG["warmup_dt_seconds"] = warmup_dt_seconds
    _WORKER_CONFIG["include_major_asteroids"] = include_major_asteroids


def _refine_one(task: _PairTask) -> dict:
    """Run one N-body refinement and return a result row (or an error row)."""
    elements_1 = dict(
        a_au=task.a_1, e=task.e_1, i_deg=task.i_1,
        Omega_deg=task.Omega_1, omega_deg=task.omega_1, M_deg=task.M_1,
        epoch_jd=task.epoch_1,
    )
    elements_2 = dict(
        a_au=task.a_2, e=task.e_2, i_deg=task.i_2,
        Omega_deg=task.Omega_2, omega_deg=task.omega_2, M_deg=task.M_2,
        epoch_jd=task.epoch_2,
    )
    row = {
        "number_1": task.number_1,
        "number_2": task.number_2,
        "bin_id": task.bin_id,
        "a_1": task.a_1, "e_1": task.e_1, "i_1": task.i_1, "q_1": task.q_1,
        "a_2": task.a_2, "e_2": task.e_2, "i_2": task.i_2, "q_2": task.q_2,
        "delta_a_au": task.delta_a_au,
        "dist_au_kepler": task.dist_au_kepler,
        "t_min_kepler_jd": task.jd_tdb_kepler,
        "rel_vel_kepler": task.rel_vel_kepler,
        "dist_au_nbody": None,
        "t_min_nbody_jd": None,
        "rel_vel_nbody": None,
        "delta_dist_au": None,
        "delta_t_min_hours": None,
        "delta_rel_vel_au_day": None,
        "nbody_converged": False,
        "nbody_energy_drift": None,
        "n_samples": None,
        "near_boundary": False,
        "error_message": None,
    }
    try:
        result = refine_pair_nbody(
            elements_1=elements_1,
            elements_2=elements_2,
            t_center_jd=task.jd_tdb_kepler,
            window_hours=_WORKER_CONFIG["window_hours"],
            sample_dt_seconds=_WORKER_CONFIG["sample_dt_seconds"],
            warmup_dt_seconds=_WORKER_CONFIG["warmup_dt_seconds"],
            include_major_asteroids=_WORKER_CONFIG["include_major_asteroids"],
        )
        row["dist_au_nbody"] = result.dist_au_nbody
        row["t_min_nbody_jd"] = result.t_min_nbody_jd
        row["rel_vel_nbody"] = result.rel_vel_au_day
        row["delta_dist_au"] = result.dist_au_nbody - task.dist_au_kepler
        delta_t_h = (result.t_min_nbody_jd - task.jd_tdb_kepler) * 24.0
        row["delta_t_min_hours"] = delta_t_h
        row["delta_rel_vel_au_day"] = result.rel_vel_au_day - task.rel_vel_kepler
        row["nbody_converged"] = result.converged
        row["nbody_energy_drift"] = result.energy_drift
        row["n_samples"] = result.n_samples
        window_h = _WORKER_CONFIG["window_hours"]
        # Flag minima within 5% of either boundary — those may be clipped
        # (true minimum could lie outside the integration window).
        row["near_boundary"] = bool(abs(delta_t_h) > 0.95 * window_h)
    except Exception as exc:  # capture per-pair failures, don't kill the pool
        row["error_message"] = f"{type(exc).__name__}: {exc}\n{traceback.format_exc(limit=2)}"
    return row


def _build_tasks(sample: pl.DataFrame) -> list[_PairTask]:
    tasks: list[_PairTask] = []
    for r in sample.iter_rows(named=True):
        tasks.append(
            _PairTask(
                number_1=int(r["number_1"]),
                number_2=int(r["number_2"]),
                bin_id=str(r["bin_id"]),
                jd_tdb_kepler=float(r["jd_tdb"]),
                dist_au_kepler=float(r["dist_au"]),
                rel_vel_kepler=float(r["rel_vel_au_day"]),
                a_1=float(r["a_1"]), e_1=float(r["e_1"]), i_1=float(r["i_1"]),
                Omega_1=float(r["Omega_1"]), omega_1=float(r["omega_1"]),
                M_1=float(r["M_1"]), epoch_1=float(r["epoch_1"]),
                a_2=float(r["a_2"]), e_2=float(r["e_2"]), i_2=float(r["i_2"]),
                Omega_2=float(r["Omega_2"]), omega_2=float(r["omega_2"]),
                M_2=float(r["M_2"]), epoch_2=float(r["epoch_2"]),
                q_1=float(r["q_1"]), q_2=float(r["q_2"]),
                delta_a_au=float(r["delta_a_au"]),
            )
        )
    return tasks


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", type=Path, default=SAMPLE)
    parser.add_argument("--out", type=Path, default=OUT)
    parser.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 4) - 4))
    parser.add_argument("--window-hours", type=float, default=12.0,
                        help="Half-width of the N-body refinement window (hours). "
                             "Default 12h: a smoke run showed Kepler t_min offsets up "
                             "to 4h, so a 6h window risks clipping the true minimum.")
    parser.add_argument("--sample-dt-seconds", type=float, default=60.0)
    parser.add_argument("--warmup-dt-seconds", type=float, default=600.0)
    parser.add_argument("--no-major-asteroids", action="store_true",
                        help="Disable Ceres/Pallas/Vesta/Hygiea as massive bodies.")
    parser.add_argument("--limit", type=int, default=None,
                        help="Optional: process only the first N pairs (smoke testing).")
    args = parser.parse_args()

    logger.info("Reading sample: %s", args.sample)
    sample = pl.read_parquet(args.sample)
    logger.info("Sample rows: %d", len(sample))
    if args.limit:
        sample = sample.head(args.limit)
        logger.info("Limit applied: %d rows", len(sample))

    tasks = _build_tasks(sample)
    logger.info("Built %d tasks; spawning %d workers", len(tasks), args.workers)

    t0 = time.monotonic()
    rows: list[dict] = []

    include_major = not args.no_major_asteroids
    ctx = mp.get_context("spawn")  # avoid fork+REBOUND fingerprint issues
    with ctx.Pool(
        processes=args.workers,
        initializer=_init_worker,
        initargs=(args.window_hours, args.sample_dt_seconds, args.warmup_dt_seconds, include_major),
    ) as pool:
        for row in tqdm(
            pool.imap_unordered(_refine_one, tasks, chunksize=4),
            total=len(tasks),
            desc="N-body refine",
        ):
            rows.append(row)

    elapsed = time.monotonic() - t0
    logger.info("Pool finished in %.1f s (%.2f s/pair avg)", elapsed, elapsed / max(1, len(rows)))

    df = pl.from_dicts(rows)
    n_ok = int(df.filter(pl.col("error_message").is_null()).height)
    n_failed = int(df.filter(pl.col("error_message").is_not_null()).height)
    logger.info("Results: ok=%d  failed=%d", n_ok, n_failed)

    if n_ok:
        ok = df.filter(pl.col("error_message").is_null())
        dd = ok["delta_dist_au"].abs()
        dt = ok["delta_t_min_hours"].abs()
        logger.info(
            "|Δdist_au|   median=%.6f  p95=%.6f  p99=%.6f  max=%.6f",
            float(dd.median()), float(dd.quantile(0.95)),
            float(dd.quantile(0.99)), float(dd.max()),
        )
        logger.info(
            "|Δt_min_h|   median=%.4f  p95=%.4f  p99=%.4f  max=%.4f",
            float(dt.median()), float(dt.quantile(0.95)),
            float(dt.quantile(0.99)), float(dt.max()),
        )
        drift = ok["nbody_energy_drift"]
        logger.info("energy_drift max=%.2e", float(drift.max()))
        n_boundary = int(ok.filter(pl.col("near_boundary")).height)
        logger.info(
            "near_boundary (|Δt| > 0.95 × window): %d / %d (%.1f%%) — "
            "true minimum may lie outside the integration window for these pairs",
            n_boundary, len(ok), 100.0 * n_boundary / max(1, len(ok)),
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.write_parquet(args.out)
    logger.info("Wrote %s (%d rows)", args.out, len(df))
    return 0 if n_failed == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
