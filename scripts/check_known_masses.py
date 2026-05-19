"""Cross-check detected-perturbation candidates against JPL Small-Body Database.

For each perturber in ``data/output/deflection_detections.csv``, queries the
JPL Small-Body Database (SBDB) for a published GM value.  Identifies which
of our detected perturbers ALREADY have a published mass — those are NOT
new mass determinations.

The candidates without a published GM are the ones that could yield a
genuinely new asteroid mass with proper orbit-fit analysis.

Output: ``data/output/perturbers_known_masses.csv``

Usage
-----
    docker compose run --rm pipeline python -m scripts.check_known_masses
"""

from __future__ import annotations

import argparse
import logging
import time
from pathlib import Path

import polars as pl
from astroquery.jplsbdb import SBDB

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)


_G = 6.674e-11
_KM3_S2_TO_KG = 1e9 / _G       # GM in km³/s² → M in kg


def query_sbdb(number: int) -> dict:
    """Query JPL SBDB for asteroid *number*. Returns dict with name + gm if available."""
    try:
        data = SBDB.query(str(number), phys=True)
    except Exception as exc:  # noqa: BLE001
        return {"error": str(exc)}

    out: dict = {"name": None, "gm": None, "gm_sig": None, "diameter": None}
    obj = data.get("object", {})
    out["name"] = obj.get("fullname")

    phys = data.get("phys_par", {})
    if isinstance(phys, dict):
        for key in ("GM", "Mass"):
            if key in phys:
                try:
                    val = phys[key]
                    if hasattr(val, "value"):
                        out["gm"] = float(val.value)
                    else:
                        out["gm"] = float(val)
                except (TypeError, ValueError):
                    pass
        for key in ("GM_sig", "Mass_sig"):
            if key in phys:
                try:
                    val = phys[key]
                    if hasattr(val, "value"):
                        out["gm_sig"] = float(val.value)
                    else:
                        out["gm_sig"] = float(val)
                except (TypeError, ValueError):
                    pass
        if "diameter" in phys:
            try:
                val = phys["diameter"]
                out["diameter"] = float(val.value) if hasattr(val, "value") else float(val)
            except (TypeError, ValueError):
                pass
    return out


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument(
        "--detections",
        type=Path,
        default=Path("data/output/deflection_detections.csv"),
    )
    p.add_argument(
        "--output",
        type=Path,
        default=Path("data/output/perturbers_known_masses.csv"),
    )
    p.add_argument(
        "--rate-limit",
        type=float,
        default=0.5,
        help="Seconds between SBDB queries.",
    )
    args = p.parse_args()

    if not args.detections.exists():
        logger.error("File not found: %s", args.detections)
        return 1

    det = pl.read_csv(args.detections)
    perturbers = det.select(["perturber_number", "perturber_name", "detection"]).unique(
        subset=["perturber_number"]
    )
    logger.info("Querying SBDB for %d unique perturbers…", perturbers.height)

    rows: list[dict] = []
    for i, r in enumerate(perturbers.iter_rows(named=True), start=1):
        num = int(r["perturber_number"])
        logger.info("[%d/%d] (%d) %s…", i, perturbers.height, num, r["perturber_name"])
        info = query_sbdb(num)
        time.sleep(args.rate_limit)

        gm = info.get("gm")
        mass_kg = (gm * _KM3_S2_TO_KG) if (gm is not None) else None

        rows.append(
            {
                "perturber_number": num,
                "perturber_name": r["perturber_name"],
                "sbdb_name": info.get("name"),
                "gm_km3_s2": gm,
                "mass_kg": mass_kg,
                "mass_kg_sig": (
                    info["gm_sig"] * _KM3_S2_TO_KG if info.get("gm_sig") is not None else None
                ),
                "diameter_km": info.get("diameter"),
                "has_published_mass": gm is not None,
                "detected_in_our_pipeline": r["detection"] == "yes",
                "novel_candidate": (gm is None) and (r["detection"] == "yes"),
                "sbdb_error": info.get("error"),
            }
        )

    out = pl.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.write_csv(args.output)
    logger.info("Wrote %d rows to %s", out.height, args.output)

    # Summary
    n_total = out.height
    n_with_mass = int(out["has_published_mass"].sum())
    n_detected = int(out["detected_in_our_pipeline"].sum())
    n_novel = int(out["novel_candidate"].sum())
    logger.info("")
    logger.info("Summary of %d unique perturbers in detection catalog:", n_total)
    logger.info("  Already have published mass in SBDB: %d", n_with_mass)
    logger.info("  Detected by our pipeline (≥3σ shift): %d", n_detected)
    logger.info("  >>> NOVEL MASS CANDIDATES (detected AND no published mass): %d <<<",
                n_novel)

    novel = out.filter(pl.col("novel_candidate")).sort("perturber_number")
    if novel.height > 0:
        logger.info("")
        logger.info("Novel candidates (the asteroids we could newly weigh):")
        for r in novel.iter_rows(named=True):
            label = f"({r['perturber_number']}) {r['perturber_name']}"
            d_str = f"D~{r['diameter_km']:.0f} km" if r['diameter_km'] is not None else ""
            logger.info("  • %-22s  %s", label, d_str)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
