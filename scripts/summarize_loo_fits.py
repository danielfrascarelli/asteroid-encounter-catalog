"""Consolidate all data/output/fit_*.json files into a ranked summary CSV.

Usage:
    docker compose run --rm pipeline python -m scripts.summarize_loo_fits
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
    fit_files = sorted(_OUTPUT_DIR.glob("fit_*.json"))
    logger.info("Found %d fit JSON files", len(fit_files))

    rows = []
    for path in fit_files:
        try:
            d = json.loads(path.read_text())
        except (json.JSONDecodeError, OSError) as exc:
            logger.warning("Skipping %s: %s", path.name, exc)
            continue

        mass = d.get("mass_kg")
        if mass is None:
            continue

        rows.append(
            {
                "fit_file": path.stem,
                "perturber": d.get("perturber"),
                "target": d.get("target"),
                "encounter_date": d.get("encounter_date", ""),
                "n_loo_orbit": d.get("n_loo_orbit", 0),
                "n_pre": d.get("n_pre", 0),
                "n_post": d.get("n_post", 0),
                "loo_orbit_al_rms_mas": d.get("loo_orbit_al_rms_mas"),
                "loo_orbit_chi2_red": d.get("loo_orbit_chi2_red"),
                "shift_al_mas": d.get("shift_al_mas"),
                "rms_pre_al_mas": d.get("rms_pre_al_mas"),
                "rms_post_al_mas": d.get("rms_post_al_mas"),
                "mass_kg": mass,
                "mass_sigma_inflated_kg": d.get("mass_sigma_inflated_kg"),
                "log10_mass": d.get("log10_mass"),
                "log10_mass_sigma_inflated": d.get("log10_mass_sigma_inflated"),
                "chi2_red_window": d.get("chi2_red"),
                "background_n": d.get("background_n", 0),
                "loo_window_days": d.get("loo_window_days", 180.0),
            }
        )

    if not rows:
        logger.error("No valid fit files found.")
        return 1

    df = (
        pl.DataFrame(rows)
        .sort("chi2_red_window", nulls_last=True)
    )

    out_path = _OUTPUT_DIR / "loo_batch_results.csv"
    df.write_csv(out_path)
    logger.info("Wrote %d rows → %s", len(df), out_path)

    # Print summary to stdout
    print("\n=== LOO BATCH FIT SUMMARY (sorted by window chi²) ===\n")
    print(
        df.select(
            [
                "fit_file",
                "n_loo_orbit",
                "loo_orbit_al_rms_mas",
                "shift_al_mas",
                "mass_kg",
                "mass_sigma_inflated_kg",
                "chi2_red_window",
            ]
        )
        .to_pandas()
        .to_string(index=False)
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
