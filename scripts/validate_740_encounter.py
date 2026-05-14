"""Validate close-encounter detection for (740) Cantabia and its orbital neighbor (7999) Nesvorny.

Both asteroids have nearly identical semi-major axis (~3.0515 AU) and inclination (~11 deg),
making them a good test pair. The Gaia DR3 orbital elements have osc_epoch within the
observation window (2014–2017), so propagation is at most ~2 years, not 10.

Usage:
    docker compose run --rm pipeline python -m scripts.validate_740_encounter
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import polars as pl
from astropy.time import Time

from src.detect.pipeline import detect_encounters
from src.ingest.gaia_orbits import load_gaia_orbits
from src.propagate.grid import make_time_grid

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)

_GAIA_START = "2014-07-25T00:00:00"
_GAIA_END = "2017-05-28T00:00:00"
_THRESHOLD_AU = 0.05  # wider than pipeline default to catch near-misses
_PAIR = [740, 7999]


def main() -> None:
    orbits_path = Path("data/raw/gaia_orbits.parquet")
    all_orbits = load_gaia_orbits(orbits_path)

    pair = all_orbits.filter(pl.col("number").is_in(_PAIR))
    if len(pair) < 2:
        found = pair["number"].to_list()
        missing = [n for n in _PAIR if n not in found]
        logger.error("Missing from gaia_orbits: %s", missing)
        sys.exit(1)

    logger.info("Pair orbital elements:")
    for row in pair.iter_rows(named=True):
        logger.info(
            "  (%d) %-15s  a=%.4f AU  e=%.4f  i=%.2f°  epoch_jd=%.1f",
            row["number"],
            row["designation"],
            row["a_au"],
            row["e"],
            row["i_deg"],
            row["epoch_jd"],
        )

    t_start = Time(_GAIA_START, scale="utc").tdb.jd
    t_end = Time(_GAIA_END, scale="utc").tdb.jd
    grid = make_time_grid(t_start, t_end, step_hours=1.0)
    logger.info(
        "Time grid: %s → %s  (%d steps)",
        _GAIA_START,
        _GAIA_END,
        len(grid),
    )

    results = detect_encounters(
        pair,
        grid,
        threshold_au=_THRESHOLD_AU,
        prefilter_enabled=False,  # only 2 asteroids, skip prefilter
        refinement_enabled=True,
    )

    if len(results) == 0:
        logger.info(
            "No encounter within %.3f AU between (%d) and (%d) during Gaia DR3 window.",
            _THRESHOLD_AU,
            *_PAIR,
        )
    else:
        logger.info("Encounters found:")
        for row in results.iter_rows(named=True):
            t = Time(row["jd_tdb"], format="jd", scale="tdb")
            logger.info(
                "  (%d)–(%d)  dist=%.6f AU  v_rel=%.4f AU/day  epoch=%s",
                row["number_1"],
                row["number_2"],
                row["dist_au"],
                row["rel_vel_au_day"],
                t.utc.iso,
            )


if __name__ == "__main__":
    sys.exit(main())
