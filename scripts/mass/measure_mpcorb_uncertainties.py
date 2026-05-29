"""Measure formal orbital-element uncertainties from JPL SBDB.

For each calibrator (target + perturber) used in Stage 4 validation,
query the JPL Small-Body Database (SBDB) API and extract the formal
sigma on the six classical elements (a, e, i, Omega, omega, M).

This is the empirical input for the ``TightPriors`` configuration of
the joint mass fit (Track A Stage 1 of the follow-up plan).

Output: ``data/output/mpcorb_uncertainties_per_element.csv``
"""

from __future__ import annotations

import argparse
import json
import logging
import time
from pathlib import Path

import polars as pl
import requests

logger = logging.getLogger(__name__)

SBDB_URL = "https://ssd-api.jpl.nasa.gov/sbdb.api"

# Stage 4 calibrators (extracted from data/output/stage4_validation_summary.csv).
PERTURBERS: tuple[int, ...] = (1, 2, 4, 10)
TARGETS: tuple[int, ...] = (
    # Ceres
    18937,
    # Pallas
    28036,
    47563,
    59882,
    60093,
    73243,
    # Hygiea
    4803,
    16772,
    45989,
    47605,
    58775,
)

# Map SBDB element names → our internal column names.
ELEMENT_MAP: dict[str, str] = {
    "a": "a_au",
    "e": "e",
    "i": "i_deg",
    "om": "Omega_deg",
    "w": "omega_deg",
    "ma": "M_deg",
}


def query_sbdb(designation: str | int, *, timeout: float = 30.0) -> dict:
    """Query JPL SBDB for an object; return full JSON response."""
    params = {"sstr": str(designation), "full-prec": 1}
    response = requests.get(SBDB_URL, params=params, timeout=timeout)
    response.raise_for_status()
    return response.json()


def extract_sigmas(sbdb_json: dict) -> dict[str, dict[str, float | str]]:
    """Extract value + sigma for each of the six classical elements."""
    out: dict[str, dict[str, float | str]] = {}
    for element in sbdb_json["orbit"]["elements"]:
        name = element.get("name")
        if name not in ELEMENT_MAP:
            continue
        internal_name = ELEMENT_MAP[name]
        sigma_str = element.get("sigma")
        if sigma_str in (None, "n.a.", ""):
            sigma_val: float | None = None
        else:
            sigma_val = float(sigma_str)
        out[internal_name] = {
            "value": float(element["value"]),
            "sigma": sigma_val,
            "units": element.get("units") or "",
        }
    return out


def collect_uncertainties(numbers: list[int]) -> pl.DataFrame:
    """Query SBDB for each number and tabulate sigmas per element."""
    rows: list[dict] = []
    for n in numbers:
        logger.info("Querying SBDB for (%d)", n)
        try:
            payload = query_sbdb(n)
        except requests.RequestException as exc:
            logger.warning("SBDB query failed for %d: %s", n, exc)
            continue

        obj = payload.get("object", {})
        orbit = payload.get("orbit", {})
        sigmas = extract_sigmas(payload)
        row = {
            "number": n,
            "fullname": obj.get("fullname"),
            "orbit_class": obj.get("orbit_class", {}).get("name"),
            "condition_code": orbit.get("condition_code"),
            "n_obs_used": orbit.get("n_obs_used"),
            "data_arc_days": orbit.get("data_arc"),
            "rms_arcsec": orbit.get("rms"),
            "epoch_jd": orbit.get("epoch"),
        }
        for internal_name in ELEMENT_MAP.values():
            block = sigmas.get(internal_name, {})
            row[f"{internal_name}_value"] = block.get("value")
            row[f"{internal_name}_sigma"] = block.get("sigma")
        rows.append(row)
        time.sleep(0.3)  # be polite to the API

    return pl.DataFrame(rows)


def summarize(df: pl.DataFrame, *, label: str) -> None:
    """Log percentile summary of sigmas across the population."""
    logger.info("--- summary (%s, n=%d) ---", label, df.height)
    sigma_cols = [
        "a_au_sigma",
        "e_sigma",
        "i_deg_sigma",
        "Omega_deg_sigma",
        "omega_deg_sigma",
        "M_deg_sigma",
    ]
    for col in sigma_cols:
        s = df[col].drop_nulls()
        if s.is_empty():
            continue
        logger.info(
            "  %s: median=%.3e  p90=%.3e  max=%.3e",
            col,
            s.median(),
            s.quantile(0.90),
            s.max(),
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/output/mpcorb_uncertainties_per_element.csv"),
    )
    parser.add_argument("--log-level", default="INFO")
    args = parser.parse_args()
    logging.basicConfig(
        level=args.log_level,
        format="%(asctime)s | %(levelname)s | %(message)s",
    )

    all_numbers = list(PERTURBERS) + list(TARGETS)
    df = collect_uncertainties(all_numbers)

    perturbers_df = df.filter(pl.col("number").is_in(list(PERTURBERS)))
    targets_df = df.filter(pl.col("number").is_in(list(TARGETS)))
    summarize(perturbers_df, label="perturbers")
    summarize(targets_df, label="targets")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    df.write_csv(args.output)
    logger.info("Wrote %s (%d rows)", args.output, df.height)

    summary = {
        "n_queried": len(all_numbers),
        "n_returned": df.height,
        "targets_median_sigma": {
            col: float(targets_df[col].drop_nulls().median())
            for col in (
                "a_au_sigma",
                "e_sigma",
                "i_deg_sigma",
                "Omega_deg_sigma",
                "omega_deg_sigma",
                "M_deg_sigma",
            )
            if not targets_df[col].drop_nulls().is_empty()
        },
    }
    summary_path = args.output.with_suffix(".summary.json")
    summary_path.write_text(json.dumps(summary, indent=2))
    logger.info("Wrote %s", summary_path)


if __name__ == "__main__":
    main()
