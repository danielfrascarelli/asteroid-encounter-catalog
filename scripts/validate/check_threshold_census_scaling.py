"""Verify the N^2 scaling used by the threshold-censoring extrapolation (M4).

The Kepler-threshold false-negative census
(:mod:`scripts.validate.measure_threshold_false_negatives`) draws a *uniform
random* sample of ``n_bodies`` numbered asteroids and counts the pairs it finds
below the catalogue threshold. Because pair counts in a fixed distance band
scale as ``N^2`` for a random draw, the sample's sub-threshold pair count,
scaled by ``(N_universe / n_bodies)^2``, must reproduce the full catalogue's
row count. This script performs that check and is the published, reproducible
gate for the census-design finding (tribunal R2, M4): if the sample were
stratified rather than random the scaling would fail.

Usage
-----
docker compose run --rm pipeline python -m scripts.validate.check_threshold_census_scaling \
    --summary data/output/kepler_false_negatives/summary.json \
    --catalog data/output/encounters_catalog_rebound_005au_b1fix.parquet
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
logger = logging.getLogger(__name__)

_DEFAULT_SUMMARY = Path("data/output/kepler_false_negatives/summary.json")
_DEFAULT_CATALOG = Path("data/output/encounters_catalog_rebound_005au_b1fix.parquet")
_DEFAULT_N_UNIVERSE = 449_454


def _catalog_count(catalog: Path) -> int:
    """Return the number of rows in the encounter catalogue via DuckDB."""
    import duckdb

    con = duckdb.connect()
    try:
        (n,) = con.execute("SELECT count(*) FROM read_parquet(?)", [str(catalog)]).fetchone()
    finally:
        con.close()
    return int(n)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, default=_DEFAULT_SUMMARY)
    parser.add_argument(
        "--catalog",
        type=Path,
        default=_DEFAULT_CATALOG,
        help="Encounter catalogue whose row count is the N^2 target.",
    )
    parser.add_argument(
        "--catalog-count",
        type=int,
        default=None,
        help="Skip the DuckDB count and use this value for the catalogue size.",
    )
    parser.add_argument("--n-universe", type=int, default=_DEFAULT_N_UNIVERSE)
    parser.add_argument("--tol", type=float, default=0.02, help="Max allowed relative error.")
    args = parser.parse_args()

    summary = json.loads(args.summary.read_text())
    n_bodies = int(summary["n_bodies"])
    n_below = int(summary["n_below_threshold_kepler"])

    catalog_count = (
        args.catalog_count if args.catalog_count is not None else _catalog_count(args.catalog)
    )

    scale = (args.n_universe / n_bodies) ** 2
    predicted = n_below * scale
    rel_err = abs(predicted - catalog_count) / catalog_count

    logger.info("Census sample: n_bodies=%d, n_below_threshold=%d", n_bodies, n_below)
    logger.info("N^2 factor (%d/%d)^2 = %.1f", args.n_universe, n_bodies, scale)
    logger.info("Predicted catalogue size: %.3e", predicted)
    logger.info("Actual catalogue size:    %.3e (%d)", catalog_count, catalog_count)
    logger.info(
        "Relative error: %.4f (%.2f %%); tolerance %.2f %%", rel_err, 100 * rel_err, 100 * args.tol
    )

    ok = rel_err <= args.tol
    logger.info("N^2 scaling check: %s", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
