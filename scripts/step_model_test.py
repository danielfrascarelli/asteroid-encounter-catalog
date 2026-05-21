"""Linear-drift + impulsive-step model fit to deflection residuals.

For each candidate, fits two competing models to the per-transit residuals
(observed_RA − Horizons_RA in mas, vs. time):

  Model A:  residual(t) = a + b·(t − t_enc)            (linear drift only)
  Model B:  residual(t) = a + b·(t − t_enc) + c·H(t)   (drift + step at encounter)

where H(t) = 0 for t < t_enc and 1 for t > t_enc. We compute ΔBIC = BIC_A − BIC_B
to decide whether the step is justified by the data.

ΔBIC > 6 means strong evidence for the step (i.e., a real impulsive
perturbation at the encounter epoch); ΔBIC < 2 means no evidence; in between
is borderline.

Output: ``data/output/step_model_results.csv``

Usage
-----
    docker compose run --rm pipeline python -m scripts.step_model_test
"""

from __future__ import annotations

import argparse
import logging
import math
from pathlib import Path

import numpy as np
import polars as pl

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)


_BLACKOUT_DAYS = 7.0


def _fit_linear(t: np.ndarray, y: np.ndarray) -> tuple[float, float, float]:
    """Least-squares fit y = a + b·t. Returns (a, b, chi2)."""
    A = np.column_stack([np.ones_like(t), t])
    coef, _, _, _ = np.linalg.lstsq(A, y, rcond=None)
    a, b = float(coef[0]), float(coef[1])
    resid = y - (a + b * t)
    chi2 = float(np.sum(resid * resid))
    return a, b, chi2


def _fit_linear_step(
    t: np.ndarray, y: np.ndarray, t_step: float
) -> tuple[float, float, float, float]:
    """Fit y = a + b·t + c·H(t > t_step). Returns (a, b, c, chi2)."""
    H = (t > t_step).astype(float)
    A = np.column_stack([np.ones_like(t), t, H])
    coef, _, _, _ = np.linalg.lstsq(A, y, rcond=None)
    a, b, c = float(coef[0]), float(coef[1]), float(coef[2])
    resid = y - (a + b * t + c * H)
    chi2 = float(np.sum(resid * resid))
    return a, b, c, chi2


def bic(chi2: float, n: int, k: int) -> float:
    """Bayesian Information Criterion for a least-squares fit.

    Assumes Gaussian residuals with known dispersion estimated from chi2/(n-k).
    """
    if n <= k:
        return float("inf")
    sigma2 = chi2 / (n - k)
    if sigma2 <= 0.0:
        return float("inf")
    return n * math.log(sigma2) + k * math.log(n)


def analyze_residuals(
    residual_csv: Path,
    encounter_date: str,
    axis: str = "dra_mas",
) -> dict:
    """Run the step-model test on one residual file.

    Parameters
    ----------
    residual_csv:
        Path to a per-candidate residuals CSV produced by detect_deflections.
    encounter_date:
        ISO UTC date of the encounter (used to centre the step).
    axis:
        Which residual axis to analyse (``"dra_mas"`` or ``"ddec_mas"``).

    Returns
    -------
    dict
        Contains ``delta_bic``, ``step_amplitude_mas``, ``step_significance``.
    """
    df = pl.read_csv(residual_csv)
    days = df["days_from_encounter"].to_numpy()
    y = df[axis].to_numpy()

    # Apply the same ±blackout exclusion to avoid the encounter epoch itself
    mask = np.abs(days) >= _BLACKOUT_DAYS
    t = days[mask]
    y = y[mask]

    if len(t) < 10:
        return {
            "delta_bic": float("nan"),
            "step_amplitude_mas": float("nan"),
            "step_chi2_red": float("nan"),
            "drift_rate_mas_per_day": float("nan"),
            "n_used": len(t),
        }

    a_lin, b_lin, chi2_lin = _fit_linear(t, y)
    a_step, b_step, c_step, chi2_step = _fit_linear_step(t, y, t_step=0.0)

    bic_lin = bic(chi2_lin, len(t), 2)
    bic_step = bic(chi2_step, len(t), 3)
    dbic = bic_lin - bic_step

    n_dof_step = max(1, len(t) - 3)
    chi2_red_step = chi2_step / n_dof_step

    # Step significance ≈ c / σ_c
    # σ_c ≈ σ_resid × sqrt(n / (n_before · n_after))
    n_before = int(np.sum(t < 0))
    n_after = int(np.sum(t > 0))
    if n_before > 0 and n_after > 0:
        sigma_resid = math.sqrt(chi2_step / n_dof_step)
        sigma_c = sigma_resid * math.sqrt(len(t) / (n_before * n_after))
        signif = c_step / sigma_c if sigma_c > 0 else float("nan")
    else:
        signif = float("nan")

    return {
        "delta_bic": dbic,
        "step_amplitude_mas": c_step,
        "step_significance": signif,
        "step_chi2_red": chi2_red_step,
        "linear_chi2_red": chi2_lin / max(1, len(t) - 2),
        "drift_rate_mas_per_day": b_step,
        "n_used": len(t),
        "n_before": n_before,
        "n_after": n_after,
    }


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
        "--output",
        type=Path,
        default=Path("data/output/step_model_results.csv"),
    )
    args = p.parse_args()

    det = pl.read_csv(args.detections)
    logger.info("Analysing step model on %d candidates", det.height)

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

        stats_ra = analyze_residuals(resid_path, r["date_utc"], axis="dra_mas")
        stats_dec = analyze_residuals(resid_path, r["date_utc"], axis="ddec_mas")

        rows.append(
            {
                "perturber_number": perturber,
                "perturber_name": r["perturber_name"],
                "target_number": target_no,
                "target_designation": r["target_designation"],
                "date_utc": r["date_utc"],
                "expected_muas": r.get("expected_muas"),
                "dbic_ra": stats_ra["delta_bic"],
                "step_ra_mas": stats_ra["step_amplitude_mas"],
                "step_signif_ra": stats_ra["step_significance"],
                "drift_ra_per_day": stats_ra["drift_rate_mas_per_day"],
                "dbic_dec": stats_dec["delta_bic"],
                "step_dec_mas": stats_dec["step_amplitude_mas"],
                "step_signif_dec": stats_dec["step_significance"],
                "drift_dec_per_day": stats_dec["drift_rate_mas_per_day"],
                "n_used": stats_ra["n_used"],
                "encounter_detected": (
                    max(abs(stats_ra["delta_bic"] or 0), abs(stats_dec["delta_bic"] or 0)) > 6.0
                ),
            }
        )

    out = pl.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.write_csv(args.output)
    logger.info("Wrote %d rows to %s", out.height, args.output)

    # Summary
    finite = out.filter(pl.col("dbic_ra").is_finite())
    n_strong = int(finite.filter((pl.col("dbic_ra") > 6.0) | (pl.col("dbic_dec") > 6.0)).height)
    n_marginal = int(
        finite.filter(
            ((pl.col("dbic_ra") > 2.0) & (pl.col("dbic_ra") <= 6.0))
            | ((pl.col("dbic_dec") > 2.0) & (pl.col("dbic_dec") <= 6.0))
        ).height
    )
    logger.info("")
    logger.info("STEP-MODEL TEST SUMMARY:")
    logger.info(
        "  Strong evidence for step (ΔBIC > 6 in RA or Dec):  %d / %d",
        n_strong,
        finite.height,
    )
    logger.info(
        "  Marginal evidence (2 < ΔBIC ≤ 6):                  %d / %d",
        n_marginal,
        finite.height,
    )

    # Top candidates by step significance
    sig_col = pl.max_horizontal(pl.col("dbic_ra").abs(), pl.col("dbic_dec").abs())
    top = finite.with_columns(sig_col.alias("max_dbic")).sort("max_dbic", descending=True).head(15)
    logger.info("")
    logger.info("Top 15 by ΔBIC (best step evidence):")
    header = f"  {'Perturber':<22} {'Target':<18} {'ΔBIC_RA':>10} {'step_RA':>10} {'ΔBIC_Dec':>10} {'step_Dec':>10}"
    logger.info(header)
    logger.info("-" * len(header))
    for r in top.iter_rows(named=True):
        logger.info(
            "  (%-3d) %-15s %-18s %10.1f %10.1f %10.1f %10.1f",
            r["perturber_number"],
            (r["perturber_name"] or "")[:15],
            (r["target_designation"] or "")[:18],
            r["dbic_ra"],
            r["step_ra_mas"],
            r["dbic_dec"],
            r["step_dec_mas"],
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
