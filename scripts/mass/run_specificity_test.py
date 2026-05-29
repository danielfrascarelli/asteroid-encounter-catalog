"""Stage 3 specificity test: real perturber vs N null perturbers.

For each candidate (target, real_perturber, encounter_date) listed in
``data/output/mass_followup_candidates.csv`` we

1. Sample ``--n-nulls`` MPCORB asteroids that share the orbital band of
   the real perturber but never come within 0.1 AU of the target inside
   the Gaia DR3 window (using the precomputed encounter catalog).
2. Re-run the joint orbit + mass fit on every null with the Stage 2
   Mahalanobis 2D likelihood, keeping the target and the encounter
   epoch fixed.
3. Compare the real ``M_fit``, ``log10_M_fit`` and ``chi2_red_joint``
   against the null distribution; report a one-sided p-value.

Usage
-----
    docker compose run --rm pipeline python -m scripts.mass.run_specificity_test \
        --targets 18105,44887,3294 --n-nulls 30 --workers 24
"""

from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import polars as pl

from scripts.mass.fit_mass_gaia_loo import _MPCORB_ARCHIVE_DIR, _best_mpcorb_snapshot
from src.ingest.mpcorb import parse_mpcorb
from src.mass.null_perturbers import sample_null_perturbers

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

_DEFAULT_CANDIDATES = Path("data/output/mass_followup_candidates.csv")
_DEFAULT_CATALOG = Path("data/output/encounters_catalog_hybrid_stageb.parquet")
_DEFAULT_OUTPUT_DIR = Path("data/output/specificity_v2")


def _run_one_fit(
    perturber: int,
    target: int,
    date_utc: str,
    output: Path,
    mpcorb: Path,
    config: str = "config.yaml",
) -> dict:
    cmd = [
        sys.executable,
        "-m",
        "scripts.mass.fit_mass_gaia_joint",
        "--config",
        config,
        "--perturber",
        str(perturber),
        "--target",
        str(target),
        "--date",
        date_utc,
        "--likelihood",
        "mahalanobis2d",
        "--mpcorb",
        str(mpcorb),
        "--output",
        str(output),
    ]
    t0 = time.monotonic()
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    elapsed = time.monotonic() - t0
    return {
        "perturber": perturber,
        "target": target,
        "returncode": result.returncode,
        "elapsed_s": elapsed,
        "stderr_tail": result.stderr[-400:] if result.stderr else "",
        "output": str(output),
    }


def _load_fit(path: Path) -> dict | None:
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return None


def _real_perturber_a(mpcorb_df: pl.DataFrame, number: int) -> float | None:
    row = mpcorb_df.filter(pl.col("number") == number)
    if row.height == 0:
        return None
    return float(row["a_au"][0])


def _percentile_rank(value: float, population: np.ndarray) -> float:
    """Fraction of population strictly less than ``value`` (0..1)."""
    if population.size == 0:
        return float("nan")
    return float(np.mean(population < value))


def _process_candidate(
    *,
    target: int,
    real_perturber: int,
    encounter_date: str,
    mpcorb_path: Path,
    mpcorb_df: pl.DataFrame,
    encounters_df: pl.DataFrame,
    output_dir: Path,
    n_nulls: int,
    a_window: float,
    h_window: float | None,
    workers: int,
    force: bool,
    seed: int,
) -> dict:
    a_real = _real_perturber_a(mpcorb_df, real_perturber)
    if a_real is None:
        logger.warning("Real perturber %d not in MPCORB %s", real_perturber, mpcorb_path.name)
        return {
            "target": target,
            "real_perturber": real_perturber,
            "status": "skip_missing_real_perturber",
        }
    nulls = sample_null_perturbers(
        target_number=target,
        real_perturber_number=real_perturber,
        real_perturber_a_au=a_real,
        mpcorb=mpcorb_df,
        encounters=encounters_df,
        n_nulls=n_nulls,
        a_window_au=a_window,
        min_separation_au=0.1,
        h_window_mag=h_window,
        seed=seed + target,
    )
    if not nulls:
        return {
            "target": target,
            "real_perturber": real_perturber,
            "status": "no_nulls_eligible",
        }
    logger.info(
        "Target %d real %d: %d nulls sampled, dispatching joint fits...",
        target,
        real_perturber,
        len(nulls),
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    real_out = output_dir / f"real_{real_perturber:06d}_{target:06d}.json"
    null_results: list[dict] = []
    tasks: list[tuple[int, Path]] = [(real_perturber, real_out)] + [
        (null, output_dir / f"null_{null:06d}_{target:06d}.json") for null in nulls
    ]

    with ProcessPoolExecutor(max_workers=workers) as pool:
        future_to_task = {}
        for perturber, out_path in tasks:
            if out_path.exists() and not force:
                logger.info("  skip existing %s", out_path.name)
                null_results.append(
                    {
                        "perturber": perturber,
                        "target": target,
                        "returncode": 0,
                        "elapsed_s": 0.0,
                        "stderr_tail": "",
                        "output": str(out_path),
                    }
                )
                continue
            fut = pool.submit(
                _run_one_fit,
                perturber,
                target,
                encounter_date,
                out_path,
                mpcorb_path,
            )
            future_to_task[fut] = (perturber, out_path)
        for fut in as_completed(future_to_task):
            res = fut.result()
            if res["returncode"] != 0:
                logger.warning(
                    "  fit failed perturber=%d rc=%d tail=%s",
                    res["perturber"],
                    res["returncode"],
                    res["stderr_tail"][-200:],
                )
            null_results.append(res)

    chi2_nulls: list[float] = []
    mass_nulls: list[float] = []
    log10_nulls: list[float] = []
    null_meta: list[dict] = []
    real_fit_path = real_out
    real_fit = _load_fit(real_fit_path)
    for null in nulls:
        path = output_dir / f"null_{null:06d}_{target:06d}.json"
        fit = _load_fit(path)
        if fit is None:
            null_meta.append({"perturber": null, "status": "no_output"})
            continue
        if not fit.get("joint_success"):
            null_meta.append({"perturber": null, "status": "fit_failed"})
            continue
        chi2_nulls.append(float(fit["chi2_red_joint"]))
        mass_nulls.append(float(fit["mass_kg"]))
        log10_nulls.append(float(fit["log10_mass"]))
        null_meta.append(
            {
                "perturber": null,
                "status": "ok",
                "chi2_red_joint": fit["chi2_red_joint"],
                "mass_kg": fit["mass_kg"],
                "log10_mass": fit["log10_mass"],
            }
        )

    chi2_real = float(real_fit["chi2_red_joint"]) if real_fit else float("nan")
    mass_real = float(real_fit["mass_kg"]) if real_fit else float("nan")
    log10_real = float(real_fit["log10_mass"]) if real_fit else float("nan")

    chi2_nulls_arr = np.asarray(chi2_nulls, dtype=float)
    mass_nulls_arr = np.asarray(mass_nulls, dtype=float)
    log10_nulls_arr = np.asarray(log10_nulls, dtype=float)

    summary = {
        "target": target,
        "real_perturber": real_perturber,
        "encounter_date": encounter_date,
        "n_nulls_requested": n_nulls,
        "n_nulls_eligible": len(nulls),
        "n_nulls_success": int(mass_nulls_arr.size),
        "chi2_red_real": chi2_real,
        "mass_real_kg": mass_real,
        "log10_mass_real": log10_real,
        "chi2_red_null_median": (
            float(np.median(chi2_nulls_arr)) if chi2_nulls_arr.size else float("nan")
        ),
        "chi2_red_null_p10": (
            float(np.quantile(chi2_nulls_arr, 0.1)) if chi2_nulls_arr.size else float("nan")
        ),
        "chi2_red_null_p90": (
            float(np.quantile(chi2_nulls_arr, 0.9)) if chi2_nulls_arr.size else float("nan")
        ),
        "mass_null_median_kg": (
            float(np.median(mass_nulls_arr)) if mass_nulls_arr.size else float("nan")
        ),
        "log10_mass_null_median": (
            float(np.median(log10_nulls_arr)) if log10_nulls_arr.size else float("nan")
        ),
        "p_value_chi2_below": (
            _percentile_rank(chi2_real, chi2_nulls_arr) if chi2_nulls_arr.size else float("nan")
        ),
        "p_value_mass_above": (
            1.0 - _percentile_rank(mass_real, mass_nulls_arr)
            if mass_nulls_arr.size
            else float("nan")
        ),
        "null_details": null_meta,
    }
    summary_path = output_dir / f"specificity_{real_perturber:06d}_{target:06d}.json"
    summary_path.write_text(json.dumps(summary, indent=2))
    logger.info(
        "Target %d done: chi2_real=%.3f vs chi2_null_median=%.3f; mass_real=%.3e vs mass_null_median=%.3e (p_chi2_below=%.3f, p_mass_above=%.3f)",
        target,
        chi2_real,
        summary["chi2_red_null_median"],
        mass_real,
        summary["mass_null_median_kg"],
        summary["p_value_chi2_below"],
        summary["p_value_mass_above"],
    )
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--candidates", type=Path, default=_DEFAULT_CANDIDATES)
    parser.add_argument("--catalog", type=Path, default=_DEFAULT_CATALOG)
    parser.add_argument(
        "--targets",
        type=str,
        default=None,
        help="Comma-separated target numbers; if omitted, use all viable candidates.",
    )
    parser.add_argument("--n-nulls", type=int, default=30)
    parser.add_argument("--a-window", type=float, default=0.5)
    parser.add_argument("--h-window", type=float, default=1.5)
    parser.add_argument("--workers", type=int, default=24)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output-dir", type=Path, default=_DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--summary",
        type=Path,
        default=Path("data/output/specificity_test_v2.csv"),
    )
    parser.add_argument("--mpcorb", type=Path, default=None)
    args = parser.parse_args()

    candidates = pl.read_csv(args.candidates)
    if "viable_obs" in candidates.columns:
        candidates = candidates.filter(pl.col("viable_obs"))
    if args.targets:
        wanted = {int(x) for x in args.targets.split(",") if x.strip()}
        candidates = candidates.filter(pl.col("target_number").is_in(list(wanted)))
    if candidates.height == 0:
        logger.error("No candidates matched the target list.")
        return 1
    logger.info("Loaded %d candidate(s) for specificity test", candidates.height)

    logger.info("Loading encounter catalog from %s", args.catalog)
    encounters = pl.read_parquet(args.catalog, columns=["number_1", "number_2", "dist_au"])
    logger.info("Catalog rows: %d", encounters.height)

    if args.mpcorb is None:
        from astropy.time import Time

        dates_sorted = sorted(candidates["date_utc"].to_list())
        median_date = dates_sorted[len(dates_sorted) // 2]
        median_jd = float(Time(median_date, scale="utc").tdb.jd)
        args.mpcorb = _best_mpcorb_snapshot(_MPCORB_ARCHIVE_DIR, median_jd)
        logger.info("Auto-selected MPCORB snapshot: %s", args.mpcorb.name)

    mpcorb_df = parse_mpcorb(str(args.mpcorb), semimajor_min_au=0.0, semimajor_max_au=50.0)

    rows = []
    for row in candidates.iter_rows(named=True):
        target = int(row["target_number"])
        real_perturber = int(row["perturber_number"])
        encounter_date = str(row["date_utc"])
        result = _process_candidate(
            target=target,
            real_perturber=real_perturber,
            encounter_date=encounter_date,
            mpcorb_path=args.mpcorb,
            mpcorb_df=mpcorb_df,
            encounters_df=encounters,
            output_dir=args.output_dir,
            n_nulls=args.n_nulls,
            a_window=args.a_window,
            h_window=args.h_window,
            workers=args.workers,
            force=args.force,
            seed=args.seed,
        )
        row_summary = {k: v for k, v in result.items() if k != "null_details"}
        rows.append(row_summary)

    out_df = pl.DataFrame(rows)
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    out_df.write_csv(args.summary)
    logger.info("Wrote specificity summary to %s", args.summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
