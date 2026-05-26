"""Regenerate the JPL Horizons fixture used by `test_refine_pair_matches_jpl_horizons`.

Picks a small set of representative pairs from the catalog (diverse in e, q, i),
runs the N-body refiner to obtain ``t_min_nbody_jd``, then queries Horizons for
the heliocentric ecliptic state vector of each asteroid at that epoch.  The
test compares the refiner's `dist_au_nbody` to the Horizons-derived value
(``|p1 − p2|``) within 1×10⁻⁴ AU.

Requires network to reach https://ssd.jpl.nasa.gov/api/horizons.api ; run
manually whenever the catalog, MPCORB snapshot, or refiner changes:

    docker compose run --rm pipeline python -m scripts.validate.refresh_horizons_fixture

The fixture is committed at ``tests/fixtures/jpl_horizons_pairs.json``.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np
import polars as pl
from astropy.time import Time
from astroquery.jplhorizons import Horizons

from scripts.validate.refine_pair_nbody import refine_pair_nbody
from src.ingest.mpcorb import parse_mpcorb

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
logger = logging.getLogger(__name__)

MPCORB = Path("data/raw/mpcorb_archive/MPCORB_20160217.DAT")
COMPARISON = Path("data/output/kepler_vs_nbody_comparison.parquet")
FIXTURE = Path("tests/fixtures/jpl_horizons_pairs.json")

# Hand-picked pairs covering different corners of orbital space, all sourced
# from ``data/output/kepler_vs_nbody_comparison.parquet`` (964-pair sample).
# The ``t_center_kepler_jd`` used as the refiner anchor is **read at runtime**
# from the parquet's ``t_min_kepler_jd`` for each pair, so the fixture stays
# consistent with whatever A.3 sample is currently committed — no hardcoded
# epochs that can drift when the sample is regenerated.
#
# Format: (number_1, number_2, label).
_PAIRS: list[tuple[int, int, str]] = [
    (96599, 188686, "low-e main-belt, inclined"),
    (34717, 153349, "very high-e (NEA-like q < 1 AU)"),
    (47784, 310584, "low-e high-i (i > 15°)"),
    (128402, 229088, "moderate-e Mars-grazing"),
    (26, 305582, "canonical low-e low-i main-belt"),
]


def _horizons_helio_xyz(num: int, jd_tdb: float) -> np.ndarray:
    """Heliocentric ecliptic position (AU) of asteroid `num` at `jd_tdb`."""
    t = Time(jd_tdb, format="jd", scale="tdb")
    # The Horizons API requires UTC time strings.  We give a 2-min window
    # straddling the epoch with a 1-min step and pick the closer sample.
    delta_min = 1.0
    h = Horizons(
        id=str(num),
        location="@10",  # @10 = Sun (heliocentric)
        epochs={
            "start": (t - delta_min / 1440.0).utc.iso[:19],
            "stop": (t + delta_min / 1440.0).utc.iso[:19],
            "step": "1m",
        },
        id_type="smallbody",
    )
    tbl = h.vectors(refplane="ecliptic")
    # Find the row whose datetime_jd is closest to jd_tdb.
    times = np.array([float(x) for x in tbl["datetime_jd"]])
    k = int(np.argmin(np.abs(times - jd_tdb)))
    return np.array([float(tbl["x"][k]), float(tbl["y"][k]), float(tbl["z"][k])])


def main() -> int:
    logger.info("Loading MPCORB: %s", MPCORB)
    elements = parse_mpcorb(MPCORB, only_numbered=True)
    rows = {int(r["number"]): r for r in elements.to_dicts()}

    logger.info("Loading comparison parquet: %s", COMPARISON)
    comparison = pl.read_parquet(COMPARISON)

    fixture: dict = {
        "mpcorb_snapshot": MPCORB.name,
        "comparison_snapshot": COMPARISON.name,
        "refiner_settings": {
            "window_hours": 12.0,
            "sample_dt_seconds": 60.0,
            "warmup_dt_seconds": 600.0,
            "include_major_asteroids": True,
        },
        "pairs": [],
    }

    for n1, n2, label in _PAIRS:
        if n1 not in rows or n2 not in rows:
            logger.warning("Skipping (%d, %d): missing from MPCORB", n1, n2)
            continue
        match = comparison.filter((pl.col("number_1") == n1) & (pl.col("number_2") == n2))
        if match.height == 0:
            logger.warning("Skipping (%d, %d): not present in %s", n1, n2, COMPARISON.name)
            continue
        t_kepler = float(match["t_min_kepler_jd"][0])
        logger.info("Refining (%d, %d) — %s @ t_kepler=%.6f", n1, n2, label, t_kepler)
        result = refine_pair_nbody(
            elements_1=rows[n1],
            elements_2=rows[n2],
            t_center_jd=t_kepler,
            window_hours=fixture["refiner_settings"]["window_hours"],
            sample_dt_seconds=fixture["refiner_settings"]["sample_dt_seconds"],
            warmup_dt_seconds=fixture["refiner_settings"]["warmup_dt_seconds"],
            include_major_asteroids=fixture["refiner_settings"]["include_major_asteroids"],
        )

        logger.info(
            "  refiner: dist=%.6f AU  t_min=%.6f JD  drift=%.2e",
            result.dist_au_nbody,
            result.t_min_nbody_jd,
            result.energy_drift,
        )

        logger.info("  querying Horizons @ jd=%.6f for #%d and #%d", result.t_min_nbody_jd, n1, n2)
        p1 = _horizons_helio_xyz(n1, result.t_min_nbody_jd)
        p2 = _horizons_helio_xyz(n2, result.t_min_nbody_jd)
        horizons_dist = float(np.linalg.norm(p1 - p2))
        delta = abs(result.dist_au_nbody - horizons_dist)
        logger.info(
            "  horizons: dist=%.6f AU  Δ(refiner−horizons)=%.6e AU",
            horizons_dist,
            delta,
        )

        fixture["pairs"].append(
            {
                "number_1": n1,
                "number_2": n2,
                "label": label,
                "t_center_kepler_jd": t_kepler,
                "t_min_nbody_jd": float(result.t_min_nbody_jd),
                "dist_au_nbody": float(result.dist_au_nbody),
                "horizons_p1_au": [float(p1[0]), float(p1[1]), float(p1[2])],
                "horizons_p2_au": [float(p2[0]), float(p2[1]), float(p2[2])],
                "horizons_dist_au": horizons_dist,
                "delta_refiner_minus_horizons_au": float(result.dist_au_nbody - horizons_dist),
            }
        )

    FIXTURE.parent.mkdir(parents=True, exist_ok=True)
    FIXTURE.write_text(json.dumps(fixture, indent=2))
    logger.info("Wrote %s with %d pairs", FIXTURE, len(fixture["pairs"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
