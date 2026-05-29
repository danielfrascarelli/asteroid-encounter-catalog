"""Run joint orbit+mass fits over mass follow-up candidates.

The runner is resumable: by default it skips candidates whose
``fit_<perturber>_<target>_joint.json`` output already exists.

Usage
-----
    docker compose run --rm pipeline python -m scripts.mass.run_joint_batch --top 5
"""

from __future__ import annotations

import argparse
import logging
import subprocess
import sys
import time
from pathlib import Path

import polars as pl

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

_DEFAULT_CANDIDATES = Path("data/output/mass_followup_candidates.csv")
_DEFAULT_OUTPUT_DIR = Path("data/output")


def _fit_path(output_dir: Path, perturber: int, target: int) -> Path:
    return output_dir / f"fit_{perturber:06d}_{target:06d}_joint.json"


def _command(args: argparse.Namespace, row: dict, out_path: Path) -> list[str]:
    cmd = [
        sys.executable,
        "-m",
        "scripts.mass.fit_mass_gaia_joint",
        "--config",
        args.config,
        "--perturber",
        str(int(row["perturber_number"])),
        "--target",
        str(int(row["target_number"])),
        "--date",
        str(row["date_utc"]),
        "--loo-window-days",
        str(args.loo_window_days),
        "--blackout-days",
        str(args.blackout_days),
        "--dt-days",
        str(args.dt_days),
        "--integrator",
        args.integrator,
        "--background-n",
        str(args.background_n),
        "--loo-max-nfev",
        str(args.loo_max_nfev),
        "--max-nfev",
        str(args.max_nfev),
        "--output",
        str(out_path),
    ]
    if args.mpcorb is not None:
        cmd.extend(["--mpcorb", str(args.mpcorb)])
    return cmd


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--candidates", type=Path, default=_DEFAULT_CANDIDATES)
    parser.add_argument("--output-dir", type=Path, default=_DEFAULT_OUTPUT_DIR)
    parser.add_argument("--top", type=int, default=None)
    parser.add_argument("--include-nonviable", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--mpcorb", type=Path, default=None)
    parser.add_argument("--loo-window-days", type=float, default=180.0)
    parser.add_argument("--blackout-days", type=float, default=7.0)
    parser.add_argument("--dt-days", type=float, default=1.0)
    parser.add_argument("--integrator", default="whfast", choices=["whfast", "ias15"])
    parser.add_argument("--background-n", type=int, default=20)
    parser.add_argument("--loo-max-nfev", type=int, default=800)
    parser.add_argument("--max-nfev", type=int, default=800)
    parser.add_argument("--sleep-seconds", type=float, default=0.0)
    parser.add_argument(
        "--report",
        type=Path,
        default=Path("data/output/joint_batch_run_report.csv"),
    )
    args = parser.parse_args()

    candidates = pl.read_csv(args.candidates)
    if not args.include_nonviable and "viable_obs" in candidates.columns:
        candidates = candidates.filter(pl.col("viable_obs"))
    candidates = candidates.sort("rank") if "rank" in candidates.columns else candidates
    if args.top is not None:
        candidates = candidates.head(args.top)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    for i, row in enumerate(candidates.iter_rows(named=True), start=1):
        perturber = int(row["perturber_number"])
        target = int(row["target_number"])
        out_path = _fit_path(args.output_dir, perturber, target)
        if out_path.exists() and not args.force:
            logger.info("[%d/%d] Skipping existing %s", i, candidates.height, out_path.name)
            rows.append(
                {
                    "perturber": perturber,
                    "target": target,
                    "date_utc": row["date_utc"],
                    "status": "skipped_existing",
                    "returncode": 0,
                    "output": str(out_path),
                }
            )
            continue

        logger.info(
            "[%d/%d] Fitting (%d, %d) %s", i, candidates.height, perturber, target, row["date_utc"]
        )
        t0 = time.monotonic()
        result = subprocess.run(_command(args, row, out_path), check=False)
        elapsed = time.monotonic() - t0
        status = "ok" if result.returncode == 0 else "failed"
        rows.append(
            {
                "perturber": perturber,
                "target": target,
                "date_utc": row["date_utc"],
                "status": status,
                "returncode": result.returncode,
                "elapsed_s": elapsed,
                "output": str(out_path),
            }
        )
        if args.sleep_seconds > 0:
            time.sleep(args.sleep_seconds)

    report = pl.DataFrame(rows)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    report.write_csv(args.report)
    n_failed = report.filter(pl.col("returncode") != 0).height
    logger.info("Wrote batch report to %s; failed=%d", args.report, n_failed)
    return 1 if n_failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
