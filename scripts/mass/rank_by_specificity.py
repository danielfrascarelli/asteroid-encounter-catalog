"""Rank candidates by encounter-specificity of the detected signal.

For each candidate, computes ΔBIC at the real encounter epoch AND at a grid
of offset epochs (every ±15 days from −150 to +150). A candidate where the
real epoch is the GLOBAL maximum of ΔBIC has a perturbation-like signature
specific to the encounter. A candidate whose ΔBIC is similar (or higher) at
a wrong epoch is dominated by drift / structure, not the encounter.

Specificity metric:
    specificity = ΔBIC_real − max(ΔBIC at offsets ≠ 0)

A large positive value means the real epoch is uniquely best.
A negative value means a different epoch fits even better — likely systematic.

Output: ``data/output/specificity_ranking.csv``
"""

from __future__ import annotations

import argparse
import logging
import math
from pathlib import Path

import numpy as np
import polars as pl

from scripts.dev.step_model_test import _fit_linear, _fit_linear_step, bic

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)


_BLACKOUT_DAYS = 7.0
_OFFSETS = list(range(-150, 151, 15))  # 21 points incl. 0


def specificity_per_axis(t: np.ndarray, y: np.ndarray) -> dict:
    """Compute ΔBIC at each offset; return real value, best offset, specificity."""
    if len(t) < 10:
        return {
            "dbic_real": float("nan"),
            "dbic_best_offset": float("nan"),
            "best_offset_days": float("nan"),
            "specificity": float("nan"),
        }
    _, _, chi2_lin = _fit_linear(t, y)
    bic_lin = bic(chi2_lin, len(t), 2)

    dbics: list[tuple[int, float]] = []
    for off in _OFFSETS:
        # require at least 3 points on each side of the step
        n_left = int(np.sum(t < off))
        n_right = int(np.sum(t > off))
        if n_left < 3 or n_right < 3:
            continue
        _, _, _, chi2_step = _fit_linear_step(t, y, t_step=float(off))
        dbic = bic_lin - bic(chi2_step, len(t), 3)
        dbics.append((off, dbic))

    if not dbics:
        return {
            "dbic_real": float("nan"),
            "dbic_best_offset": float("nan"),
            "best_offset_days": float("nan"),
            "specificity": float("nan"),
        }

    dbic_real_lookup = {off: d for off, d in dbics if off == 0}
    dbic_real = dbic_real_lookup.get(0, float("nan"))
    others = [(off, d) for off, d in dbics if off != 0]
    if not others:
        return {
            "dbic_real": dbic_real,
            "dbic_best_offset": float("nan"),
            "best_offset_days": float("nan"),
            "specificity": float("nan"),
        }
    best_off, best_d = max(others, key=lambda p: p[1])
    return {
        "dbic_real": dbic_real,
        "dbic_best_offset": best_d,
        "best_offset_days": float(best_off),
        "specificity": dbic_real - best_d,
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
        default=Path("data/output/specificity_ranking.csv"),
    )
    args = p.parse_args()

    det = pl.read_csv(args.detections)
    logger.info(
        "Computing specificity for %d candidates over %d offsets", det.height, len(_OFFSETS)
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

        df = pl.read_csv(resid_path)
        days = df["days_from_encounter"].to_numpy()
        mask = np.abs(days) >= _BLACKOUT_DAYS

        ra = specificity_per_axis(days[mask], df["dra_mas"].to_numpy()[mask])
        dec = specificity_per_axis(days[mask], df["ddec_mas"].to_numpy()[mask])

        spec_ra = ra["specificity"]
        spec_dec = dec["specificity"]
        max_spec = max(
            spec_ra if math.isfinite(spec_ra) else -math.inf,
            spec_dec if math.isfinite(spec_dec) else -math.inf,
        )

        rows.append(
            {
                "perturber_number": perturber,
                "perturber_name": r["perturber_name"],
                "target_number": target_no,
                "target_designation": r["target_designation"],
                "date_utc": r["date_utc"],
                "expected_muas": r.get("expected_muas"),
                "dbic_real_ra": ra["dbic_real"],
                "best_offset_ra": ra["best_offset_days"],
                "specificity_ra": spec_ra,
                "dbic_real_dec": dec["dbic_real"],
                "best_offset_dec": dec["best_offset_days"],
                "specificity_dec": spec_dec,
                "max_specificity": max_spec if math.isfinite(max_spec) else float("nan"),
                "encounter_specific": max_spec > 2.0,
            }
        )

    out = pl.DataFrame(rows).sort("max_specificity", descending=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.write_csv(args.output)
    logger.info("Wrote %d rows to %s", out.height, args.output)

    # Summary
    n_specific = int(out["encounter_specific"].sum())
    logger.info("")
    logger.info(
        "Encounter-specific detections (specificity > 2 BIC units): %d / %d", n_specific, out.height
    )

    logger.info("")
    logger.info("Top 20 by max_specificity:")
    header = (
        f"  {'Perturber':<20} {'Target':<18} "
        f"{'ΔBIC_RA':>9} {'spec_RA':>8} {'ΔBIC_Dec':>9} {'spec_Dec':>8}"
    )
    logger.info(header)
    logger.info("-" * len(header))
    for r in out.head(20).iter_rows(named=True):
        tag = "★" if r["encounter_specific"] else " "
        logger.info(
            "%s (%-4d) %-13s %-18s %9.1f %8.1f %9.1f %8.1f",
            tag,
            r["perturber_number"],
            (r["perturber_name"] or "")[:13],
            (r["target_designation"] or "")[:18],
            r["dbic_real_ra"],
            r["specificity_ra"],
            r["dbic_real_dec"],
            r["specificity_dec"],
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
