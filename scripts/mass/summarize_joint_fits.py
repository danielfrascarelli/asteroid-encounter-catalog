"""Consolidate joint orbit+mass fit JSON files into a summary CSV.

Usage:
    docker compose run --rm pipeline python -m scripts.mass.summarize_joint_fits
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import polars as pl

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

_OUTPUT_DIR = Path("data/output")


def main() -> int:
    fit_files = sorted(_OUTPUT_DIR.glob("fit_*_joint.json"))
    logger.info("Found %d joint fit JSON files", len(fit_files))

    rows = []
    for path in fit_files:
        try:
            data = json.loads(path.read_text())
        except (json.JSONDecodeError, OSError) as exc:
            logger.warning("Skipping %s: %s", path.name, exc)
            continue
        rows.append(
            {
                "fit_file": path.stem,
                "perturber": data.get("perturber"),
                "target": data.get("target"),
                "encounter_date": data.get("encounter_date", ""),
                "joint_success": data.get("joint_success"),
                "joint_nfev": data.get("joint_nfev"),
                "n_joint": data.get("n_joint"),
                "n_pre": data.get("n_pre"),
                "n_post": data.get("n_post"),
                "loo_orbit_chi2_red": data.get("loo_orbit_chi2_red"),
                "chi2_red_joint": data.get("chi2_red_joint"),
                "jtj_condition": data.get("jtj_condition"),
                "active_bounds": ",".join(data.get("active_bounds", [])),
                "mass_kg": data.get("mass_kg"),
                "mass_sigma_kg": data.get("mass_sigma_kg"),
                "log10_mass": data.get("log10_mass"),
                "log10_mass_sigma": data.get("log10_mass_sigma"),
                "da_rel": data.get("da_rel"),
                "de": data.get("de"),
                "di_deg": data.get("di_deg"),
                "dOmega_deg": data.get("dOmega_deg"),
                "domega_deg": data.get("domega_deg"),
                "dM_deg": data.get("dM_deg"),
            }
        )

    if not rows:
        logger.error("No valid joint fit files found.")
        return 1

    df = pl.DataFrame(rows).sort("chi2_red_joint", nulls_last=True)
    out_path = _OUTPUT_DIR / "loo_batch_results_joint.csv"
    df.write_csv(out_path)
    logger.info("Wrote %d rows to %s", df.height, out_path)

    print("\n=== JOINT ORBIT+MASS FIT SUMMARY ===\n")
    print(
        df.select(
            [
                "fit_file",
                "joint_success",
                "n_joint",
                "chi2_red_joint",
                "mass_kg",
                "mass_sigma_kg",
                "active_bounds",
            ]
        )
        .to_pandas()
        .to_string(index=False)
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
