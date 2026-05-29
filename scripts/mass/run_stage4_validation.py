"""Stage 4: validation against literature masses (Ceres / Vesta / Pallas / Hygiea).

For each calibrator perturber, pick its closest approaches in the hybrid
encounter catalog, refit the joint mass model with the Stage 2 Mahalanobis
2D likelihood, and report ``z = (M_fit - M_lit) / sqrt(sigma_fit^2 +
sigma_lit^2)``. A pass means ``|z| < 3``.

Usage
-----
    docker compose run --rm pipeline python -m scripts.mass.run_stage4_validation \\
        --top-per-perturber 5 --workers 12
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

import polars as pl
from astropy.time import Time

from scripts.mass.fit_mass_gaia_loo import _MPCORB_ARCHIVE_DIR, _best_mpcorb_snapshot

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

_DEFAULT_CATALOG = Path("data/output/encounters_catalog_hybrid_stageb.parquet")
_DEFAULT_OUTPUT_DIR = Path("data/output/stage4_validation")


@dataclass(frozen=True)
class Calibrator:
    number: int
    name: str
    mass_kg: float
    mass_sigma_kg: float
    source: str


_CALIBRATORS: list[Calibrator] = [
    Calibrator(1, "Ceres", 4.71e20, 0.04e20, "DAWN (Park+ 2016)"),
    Calibrator(4, "Vesta", 2.59e20, 0.01e20, "DAWN (Russell+ 2012)"),
    Calibrator(2, "Pallas", 2.05e20, 0.05e20, "Goffin (2014)"),
    Calibrator(10, "Hygiea", 8.3e19, 0.4e19, "Vernazza+ (2020)"),
]


def _select_top_targets(
    catalog: pl.DataFrame,
    perturber: int,
    top_n: int,
    max_target_number: int = 100_000,
    max_dist_au: float = 0.05,
) -> pl.DataFrame:
    """Pick the closest distinct targets where ``perturber`` is the heavy body.

    Restricts to numbered targets below ``max_target_number`` because higher
    numbers in MPCORB are generally fainter / younger and have few or no
    Gaia DR3 observations.
    """
    sub = catalog.filter((pl.col("number_1") == perturber) | (pl.col("number_2") == perturber))
    sub = sub.with_columns(
        pl.when(pl.col("number_1") == perturber)
        .then(pl.col("number_2"))
        .otherwise(pl.col("number_1"))
        .alias("target")
    )
    sub = sub.filter((pl.col("target") < max_target_number) & (pl.col("dist_au") < max_dist_au))
    sub = sub.sort("dist_au").unique(subset=["target"], keep="first")
    return sub.head(top_n)


def _run_one_fit(
    perturber: int,
    target: int,
    date_utc: str,
    output: Path,
    mpcorb: Path,
) -> dict:
    cmd = [
        sys.executable,
        "-m",
        "scripts.mass.fit_mass_gaia_joint",
        "--config",
        "config.yaml",
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
    return {
        "perturber": perturber,
        "target": target,
        "returncode": result.returncode,
        "elapsed_s": time.monotonic() - t0,
        "stderr_tail": result.stderr[-300:] if result.stderr else "",
        "output": str(output),
    }


def _load_fit(path: Path) -> dict | None:
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return None


def _z_score(
    m_fit: float | None,
    sigma_fit: float | None,
    m_lit: float,
    sigma_lit: float,
) -> float | None:
    if m_fit is None or sigma_fit is None or not math.isfinite(m_fit):
        return None
    denom = math.sqrt(max(sigma_fit, 0.0) ** 2 + sigma_lit**2)
    if denom <= 0:
        return None
    return (m_fit - m_lit) / denom


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--catalog", type=Path, default=_DEFAULT_CATALOG)
    parser.add_argument("--top-per-perturber", type=int, default=5)
    parser.add_argument("--workers", type=int, default=12)
    parser.add_argument("--output-dir", type=Path, default=_DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--summary",
        type=Path,
        default=Path("data/output/stage4_validation_summary.csv"),
    )
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Loading catalog from %s", args.catalog)
    catalog = pl.read_parquet(args.catalog, columns=["number_1", "number_2", "dist_au", "jd_tdb"])
    logger.info("Catalog rows: %d", catalog.height)

    tasks: list[tuple[Calibrator, int, str, Path, Path]] = []
    for cal in _CALIBRATORS:
        top = _select_top_targets(catalog, cal.number, args.top_per_perturber)
        if top.height == 0:
            logger.warning("No encounters under threshold for %s (%d)", cal.name, cal.number)
            continue
        for row in top.iter_rows(named=True):
            target = int(row["target"])
            date_utc = Time(float(row["jd_tdb"]), format="jd", scale="tdb").utc.iso[:10]
            out_path = args.output_dir / f"fit_{cal.number:06d}_{target:06d}_stage4.json"
            mpcorb = _best_mpcorb_snapshot(_MPCORB_ARCHIVE_DIR, float(row["jd_tdb"]))
            tasks.append((cal, target, date_utc, out_path, mpcorb))

    logger.info("Dispatching %d Stage 4 fits", len(tasks))
    futures = {}
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for cal, target, date_utc, out_path, mpcorb in tasks:
            if out_path.exists() and not args.force:
                logger.info("  skip existing %s", out_path.name)
                continue
            fut = pool.submit(_run_one_fit, cal.number, target, date_utc, out_path, mpcorb)
            futures[fut] = (cal, target, date_utc, out_path)
        for fut in as_completed(futures):
            res = fut.result()
            if res["returncode"] != 0:
                logger.warning(
                    "  fit failed perturber=%d target=%d tail=%s",
                    res["perturber"],
                    res["target"],
                    res["stderr_tail"][-200:],
                )

    rows = []
    for cal, target, date_utc, out_path, _mpcorb in tasks:
        fit = _load_fit(out_path)
        if fit is None:
            rows.append(
                {
                    "perturber": cal.number,
                    "perturber_name": cal.name,
                    "target": target,
                    "encounter_date": date_utc,
                    "status": "no_output",
                    "mass_lit_kg": cal.mass_kg,
                    "mass_lit_sigma_kg": cal.mass_sigma_kg,
                    "literature_source": cal.source,
                }
            )
            continue
        if not fit.get("joint_success"):
            rows.append(
                {
                    "perturber": cal.number,
                    "perturber_name": cal.name,
                    "target": target,
                    "encounter_date": date_utc,
                    "status": "fit_failed",
                    "mass_lit_kg": cal.mass_kg,
                    "mass_lit_sigma_kg": cal.mass_sigma_kg,
                    "literature_source": cal.source,
                }
            )
            continue
        m_fit = float(fit["mass_kg"])
        sigma_fit = float(fit.get("mass_sigma_kg") or float("nan"))
        z = _z_score(m_fit, sigma_fit, cal.mass_kg, cal.mass_sigma_kg)
        rows.append(
            {
                "perturber": cal.number,
                "perturber_name": cal.name,
                "target": target,
                "encounter_date": date_utc,
                "status": "ok",
                "n_joint": fit.get("n_joint"),
                "chi2_red_joint": fit.get("chi2_red_joint"),
                "mass_fit_kg": m_fit,
                "mass_fit_sigma_kg": sigma_fit,
                "log10_mass_fit": fit.get("log10_mass"),
                "mass_lit_kg": cal.mass_kg,
                "mass_lit_sigma_kg": cal.mass_sigma_kg,
                "literature_source": cal.source,
                "z_score": z,
                "ratio_fit_over_lit": m_fit / cal.mass_kg if cal.mass_kg > 0 else None,
            }
        )

    df = pl.DataFrame(rows).sort(["perturber", "target"])
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    df.write_csv(args.summary)
    logger.info("Wrote %d rows to %s", df.height, args.summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
