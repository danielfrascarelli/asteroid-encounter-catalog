"""Run Track 1 Stage B N-body refinement for a shard of the selected subset.

This script consumes the subset produced by ``select_stageb_nbody_subset.py``
and writes one checkpoint parquet per processed shard. It intentionally refuses
to process the full 8M+ row subset unless the caller explicitly asks for
``--all``; production runs should use ``--shard-index`` and ``--shard-size``.

Example smoke run
-----------------
    docker compose run --rm pipeline python -m scripts.validate.refine_stageb_nbody \\
        --input data/cache/nbody_validation/stageb_subset_smoke.parquet \\
        --limit 25 --workers 4 --output-dir data/output/stageb_nbody_smoke

Example production shard
------------------------
    docker compose run --rm pipeline python -m scripts.validate.refine_stageb_nbody \\
        --shard-index 0 --shard-size 10000 --workers 24

Example resumable batch
-----------------------
    docker compose run --rm pipeline python -m scripts.validate.refine_stageb_nbody \\
        --shard-index 0 --shard-size 10000 --num-shards 10 --workers 24
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

SUBSET = Path("data/cache/nbody_validation/stageb_selective_subset.parquet")
OUTPUT_DIR = Path("data/output/stageb_nbody_shards")


@dataclass(frozen=True)
class _Task:
    row_index: int
    number_1: int
    number_2: int
    designation_1: str
    designation_2: str
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
    q_1: float
    a_2: float
    e_2: float
    i_2: float
    Omega_2: float
    omega_2: float
    M_2: float
    epoch_2: float
    q_2: float
    delta_a_au: float
    e_max: float
    i_max: float
    q_min: float
    stageb_reason: str


_WORKER_CONFIG: dict = {}


def _init_worker(
    window_hours: float,
    sample_dt_seconds: float,
    warmup_dt_seconds: float,
    include_major_asteroids: bool,
) -> None:
    for var in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_MAX_THREADS",
    ):
        os.environ.setdefault(var, "1")
    _WORKER_CONFIG["window_hours"] = window_hours
    _WORKER_CONFIG["sample_dt_seconds"] = sample_dt_seconds
    _WORKER_CONFIG["warmup_dt_seconds"] = warmup_dt_seconds
    _WORKER_CONFIG["include_major_asteroids"] = include_major_asteroids


def _task_from_row(row_index: int, row: dict) -> _Task:
    return _Task(
        row_index=row_index,
        number_1=int(row["number_1"]),
        number_2=int(row["number_2"]),
        designation_1=str(row["designation_1"]),
        designation_2=str(row["designation_2"]),
        jd_tdb_kepler=float(row["jd_tdb"]),
        dist_au_kepler=float(row["dist_au"]),
        rel_vel_kepler=float(row["rel_vel_au_day"]),
        a_1=float(row["a_1"]),
        e_1=float(row["e_1"]),
        i_1=float(row["i_1"]),
        Omega_1=float(row["node_1"]),
        omega_1=float(row["argperi_1"]),
        M_1=float(row["mean_anomaly_1"]),
        epoch_1=float(row["epoch_1"]),
        q_1=float(row["q_1"]),
        a_2=float(row["a_2"]),
        e_2=float(row["e_2"]),
        i_2=float(row["i_2"]),
        Omega_2=float(row["node_2"]),
        omega_2=float(row["argperi_2"]),
        M_2=float(row["mean_anomaly_2"]),
        epoch_2=float(row["epoch_2"]),
        q_2=float(row["q_2"]),
        delta_a_au=float(row["delta_a_au"]),
        e_max=float(row["e_max"]),
        i_max=float(row["i_max"]),
        q_min=float(row["q_min"]),
        stageb_reason=str(row["stageb_reason"]),
    )


def _base_row(task: _Task) -> dict:
    return {
        "row_index": task.row_index,
        "number_1": task.number_1,
        "number_2": task.number_2,
        "designation_1": task.designation_1,
        "designation_2": task.designation_2,
        "stageb_reason": task.stageb_reason,
        "a_1": task.a_1,
        "e_1": task.e_1,
        "i_1": task.i_1,
        "q_1": task.q_1,
        "a_2": task.a_2,
        "e_2": task.e_2,
        "i_2": task.i_2,
        "q_2": task.q_2,
        "delta_a_au": task.delta_a_au,
        "e_max": task.e_max,
        "i_max": task.i_max,
        "q_min": task.q_min,
        "dist_au_kepler": task.dist_au_kepler,
        "t_min_kepler_jd": task.jd_tdb_kepler,
        "rel_vel_kepler": task.rel_vel_kepler,
        "dist_au_nbody": None,
        "t_min_nbody_jd": None,
        "rel_vel_nbody": None,
        "delta_dist_au": None,
        "delta_t_min_hours": None,
        "delta_rel_vel_au_day": None,
        "refinement_method": "nbody",
        "nbody_converged": False,
        "nbody_energy_drift": None,
        "n_samples": None,
        "near_boundary": False,
        "error_message": None,
    }


def _refine_one(task: _Task) -> dict:
    row = _base_row(task)
    elements_1 = {
        "number": task.number_1,
        "a_au": task.a_1,
        "e": task.e_1,
        "i_deg": task.i_1,
        "Omega_deg": task.Omega_1,
        "omega_deg": task.omega_1,
        "M_deg": task.M_1,
        "epoch_jd": task.epoch_1,
    }
    elements_2 = {
        "number": task.number_2,
        "a_au": task.a_2,
        "e": task.e_2,
        "i_deg": task.i_2,
        "Omega_deg": task.Omega_2,
        "omega_deg": task.omega_2,
        "M_deg": task.M_2,
        "epoch_jd": task.epoch_2,
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
        row["near_boundary"] = bool(abs(delta_t_h) > 0.95 * _WORKER_CONFIG["window_hours"])
    except Exception as exc:
        row["error_message"] = f"{type(exc).__name__}: {exc}\n{traceback.format_exc(limit=2)}"
    return row


def _load_tasks(path: Path, *, offset: int, limit: int | None) -> list[_Task]:
    logger.info("Reading input shard: %s offset=%d limit=%s", path, offset, limit)
    lazy = pl.scan_parquet(path).with_row_index("row_index")
    if limit is None:
        frame = lazy.slice(offset).collect()
    else:
        frame = lazy.slice(offset, limit).collect()
    return [_task_from_row(int(row["row_index"]), row) for row in frame.iter_rows(named=True)]


def _resolve_windows(args: argparse.Namespace) -> list[tuple[int, int | None]]:
    if args.shard_index is not None:
        return [
            ((args.shard_index + shard_offset) * args.shard_size, args.shard_size)
            for shard_offset in range(args.num_shards)
        ]
    if args.limit is not None:
        return [(args.offset, args.limit)]
    if args.all:
        return [(args.offset, None)]
    raise SystemExit("Pass --shard-index, --limit, or --all; refusing accidental full run.")


def _out_path(output_dir: Path, offset: int, limit: int | None) -> Path:
    stop_label = "end" if limit is None else f"{offset + limit:09d}"
    return output_dir / f"stageb_nbody_{offset:09d}_{stop_label}.parquet"


def _run_window(
    *,
    pool: mp.pool.Pool,
    input_path: Path,
    output_dir: Path,
    offset: int,
    limit: int | None,
    force: bool,
    show_progress: bool,
) -> tuple[int, int, int]:
    out_path = _out_path(output_dir, offset, limit)
    if out_path.exists() and not force:
        logger.info("Shard already exists, skipping: %s", out_path)
        return 0, 0, 0

    tasks = _load_tasks(input_path, offset=offset, limit=limit)
    if not tasks:
        logger.warning("No rows to process for offset=%d limit=%s", offset, limit)
        return 0, 0, 0

    logger.info("Processing %d tasks for offset=%d limit=%s", len(tasks), offset, limit)
    rows: list[dict] = []
    t0 = time.monotonic()
    rows_iter = pool.imap_unordered(_refine_one, tasks, chunksize=4)
    if show_progress:
        rows_iter = tqdm(rows_iter, total=len(tasks), desc=f"Stage B {offset}")
    for row in rows_iter:
        rows.append(row)

    elapsed = time.monotonic() - t0
    df = pl.from_dicts(rows).sort("row_index")
    n_failed = df.filter(pl.col("error_message").is_not_null()).height
    n_unconverged = df.filter(~pl.col("nbody_converged")).height
    logger.info(
        "Finished offset=%d rows=%d in %.1f s (%.3f s/row avg); failed=%d unconverged=%d",
        offset,
        len(df),
        elapsed,
        elapsed / max(1, len(df)),
        n_failed,
        n_unconverged,
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    df.write_parquet(out_path)
    logger.info("Wrote %s", out_path)
    return len(df), n_failed, n_unconverged


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=SUBSET)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 4) - 4))
    parser.add_argument("--offset", type=int, default=0)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--all", action="store_true", help="Process all rows from --offset onward.")
    parser.add_argument("--shard-index", type=int, default=None)
    parser.add_argument("--shard-size", type=int, default=10_000)
    parser.add_argument(
        "--num-shards",
        type=int,
        default=1,
        help="With --shard-index, process this many consecutive shards in one pool.",
    )
    parser.add_argument("--force", action="store_true", help="Overwrite existing shard output.")
    parser.add_argument("--window-hours", type=float, default=12.0)
    parser.add_argument("--sample-dt-seconds", type=float, default=60.0)
    parser.add_argument("--warmup-dt-seconds", type=float, default=600.0)
    parser.add_argument("--no-major-asteroids", action="store_true")
    parser.add_argument("--no-progress", action="store_true", help="Disable tqdm progress bars.")
    args = parser.parse_args()

    windows = _resolve_windows(args)
    logger.info("Resolved %d window(s); spawning %d workers", len(windows), args.workers)

    t0_all = time.monotonic()
    total_rows = 0
    total_failed = 0
    total_unconverged = 0
    ctx = mp.get_context("spawn")
    with ctx.Pool(
        processes=args.workers,
        initializer=_init_worker,
        initargs=(
            args.window_hours,
            args.sample_dt_seconds,
            args.warmup_dt_seconds,
            not args.no_major_asteroids,
        ),
    ) as pool:
        for offset, limit in windows:
            n_rows, n_failed, n_unconverged = _run_window(
                pool=pool,
                input_path=args.input,
                output_dir=args.output_dir,
                offset=offset,
                limit=limit,
                force=args.force,
                show_progress=not args.no_progress,
            )
            total_rows += n_rows
            total_failed += n_failed
            total_unconverged += n_unconverged

    elapsed_all = time.monotonic() - t0_all
    logger.info(
        "Batch finished: rows=%d elapsed=%.1f s failed=%d unconverged=%d",
        total_rows,
        elapsed_all,
        total_failed,
        total_unconverged,
    )
    return 0 if total_failed == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
