"""Null test for the step-model detection: shift the step epoch off the real date.

If the step detected by ``step_model_test.py`` is genuinely caused by the
asteroid encounter, then placing the step at the WRONG epoch (e.g. 60 days
earlier or later) should significantly reduce the evidence (ΔBIC drops
toward zero or negative).

This script runs the same step-vs-drift comparison but with t_step displaced
from the encounter epoch by ``--offset-days`` days.

A signal that is truly impulsive at the encounter epoch will fail this test
(low ΔBIC). A "detection" that survives the test was probably reflecting
slow orbital drift, not the encounter.

Output: ``data/output/step_model_null.csv``
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import polars as pl

from scripts.dev.step_model_test import _fit_linear, _fit_linear_step, bic

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)


_BLACKOUT_DAYS = 7.0


def analyze_with_offset(residual_csv: Path, offset_days: float, axis: str = "dra_mas") -> dict:
    import numpy as np

    df = pl.read_csv(residual_csv)
    days = df["days_from_encounter"].to_numpy()
    y = df[axis].to_numpy()
    mask = np.abs(days) >= _BLACKOUT_DAYS
    t = days[mask]
    y = y[mask]
    if len(t) < 10:
        return {"delta_bic": float("nan"), "step_amp": float("nan")}

    _, _, chi2_lin = _fit_linear(t, y)
    _, _, c_step, chi2_step = _fit_linear_step(t, y, t_step=offset_days)

    bic_lin = bic(chi2_lin, len(t), 2)
    bic_step = bic(chi2_step, len(t), 3)
    return {"delta_bic": bic_lin - bic_step, "step_amp": c_step, "n": len(t)}


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument(
        "--detections",
        type=Path,
        default=Path("data/output/deflection_detections.csv"),
    )
    p.add_argument(
        "--residuals-dir",
        type=Path,
        default=Path("data/output/deflection_residuals"),
    )
    p.add_argument(
        "--offset-days",
        type=float,
        default=60.0,
        help="Days to shift the step epoch away from the real encounter.",
    )
    p.add_argument(
        "--output",
        type=Path,
        default=Path("data/output/step_model_null.csv"),
    )
    args = p.parse_args()

    det = pl.read_csv(args.detections)
    logger.info(
        "Running null test with step offset = %+.0f days on %d candidates",
        args.offset_days,
        det.height,
    )

    rows: list[dict] = []
    for r in det.iter_rows(named=True):
        target_no = r["target_number"]
        perturber = int(r["perturber_number"])
        if target_no is None:
            continue
        target_no = int(target_no)
        resid_path = args.residuals_dir / f"{perturber:06d}_{target_no:06d}.csv"
        if not resid_path.exists():
            continue

        real = analyze_with_offset(resid_path, offset_days=0.0, axis="dra_mas")
        fake = analyze_with_offset(resid_path, offset_days=args.offset_days, axis="dra_mas")
        real_dec = analyze_with_offset(resid_path, offset_days=0.0, axis="ddec_mas")
        fake_dec = analyze_with_offset(resid_path, offset_days=args.offset_days, axis="ddec_mas")

        rows.append(
            {
                "perturber_number": perturber,
                "perturber_name": r["perturber_name"],
                "target_designation": r["target_designation"],
                "dbic_real_ra": real["delta_bic"],
                "dbic_fake_ra": fake["delta_bic"],
                "step_real_ra": real["step_amp"],
                "step_fake_ra": fake["step_amp"],
                "dbic_real_dec": real_dec["delta_bic"],
                "dbic_fake_dec": fake_dec["delta_bic"],
            }
        )

    out = pl.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.write_csv(args.output)
    logger.info("Wrote %d rows to %s", out.height, args.output)

    # Summary
    finite = out.filter(pl.col("dbic_real_ra").is_finite())
    n_real_strong = int(
        finite.filter((pl.col("dbic_real_ra") > 6) | (pl.col("dbic_real_dec") > 6)).height
    )
    n_fake_strong = int(
        finite.filter((pl.col("dbic_fake_ra") > 6) | (pl.col("dbic_fake_dec") > 6)).height
    )
    logger.info("")
    logger.info("NULL TEST RESULT (step offset = %+.0f days):", args.offset_days)
    logger.info("  Real encounter epoch:  %d / %d strong (ΔBIC>6)", n_real_strong, finite.height)
    logger.info("  Offset epoch:          %d / %d strong (ΔBIC>6)", n_fake_strong, finite.height)
    if n_fake_strong < n_real_strong:
        excess = n_real_strong - n_fake_strong
        logger.info(
            "  → REAL excess of %d detections (%.0f%%) — encounter is the right epoch ✓",
            excess,
            100 * excess / max(1, finite.height),
        )
    else:
        logger.warning("  → No real-vs-fake excess — the step is not encounter-specific")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
