"""Select the Track 1 Stage B subset that needs N-body refinement.

Stage A found that Kepler-refined distances are least defensible for pairs
with low perihelion or high eccentricity. This script applies that criterion
to the frozen 72M-row catalog by joining each pair to the MPCORB snapshot and
writing a reproducible subset parquet with all elements needed by the N-body
runner.

Default criterion
-----------------
    q_min < 1.8 AU OR e_max > 0.3

where:
    q_min = min(a_1 * (1 - e_1), a_2 * (1 - e_2))
    e_max = max(e_1, e_2)

Usage
-----
Summary only, safe first pass:

    docker compose run --rm pipeline python -m scripts.validate.select_stageb_nbody_subset \\
        --summary-only

Materialize the full subset:

    docker compose run --rm pipeline python -m scripts.validate.select_stageb_nbody_subset

Smoke output:

    docker compose run --rm pipeline python -m scripts.validate.select_stageb_nbody_subset \\
        --limit 1000 --out data/cache/nbody_validation/stageb_subset_smoke.parquet
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import duckdb
import polars as pl

from src.ingest.mpcorb import parse_mpcorb

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
logger = logging.getLogger(__name__)

CATALOG = Path("data/output/encounters_catalog_rebound_005au.parquet")
MPCORB = Path("data/raw/mpcorb_archive/MPCORB_20160217.DAT")
ELEMENTS_CACHE = Path("data/cache/nbody_validation/mpcorb_stageb_elements.parquet")
OUT = Path("data/cache/nbody_validation/stageb_selective_subset.parquet")


def _sql_path(path: Path) -> str:
    """Return a DuckDB-safe single-quoted path literal."""
    return "'" + str(path).replace("'", "''") + "'"


def _ensure_elements_cache(mpcorb: Path, elements_cache: Path, *, force: bool) -> None:
    """Write the small MPCORB element table used by DuckDB joins."""
    if elements_cache.exists() and not force:
        logger.info("Reusing MPCORB elements cache: %s", elements_cache)
        return

    logger.info("Parsing MPCORB: %s", mpcorb)
    elements = parse_mpcorb(mpcorb, only_numbered=True).select(
        [
            "number",
            "a_au",
            "e",
            "i_deg",
            "Omega_deg",
            "omega_deg",
            "M_deg",
            "epoch_jd",
            "H",
        ]
    )
    # DuckDB identifiers are case-insensitive, so MPCORB's Omega/omega pair
    # must be physically renamed before writing the cache. Otherwise DuckDB can
    # read both as the same column and silently poison the state vectors.
    elements = elements.rename(
        {
            "Omega_deg": "node_deg",
            "omega_deg": "argperi_deg",
            "M_deg": "mean_anomaly_deg",
            "H": "h_mag",
        }
    )
    elements = elements.with_columns((pl.col("a_au") * (1.0 - pl.col("e"))).alias("q_au"))

    elements_cache.parent.mkdir(parents=True, exist_ok=True)
    elements.write_parquet(elements_cache)
    logger.info("Wrote %s (%d rows)", elements_cache, len(elements))


def _joined_cte(catalog: Path, elements_cache: Path) -> str:
    """Common DuckDB CTE used by the summary and materialization queries."""
    cat = _sql_path(catalog)
    elem = _sql_path(elements_cache)
    return f"""
WITH
cat AS (
    SELECT number_1, number_2, designation_1, designation_2,
           jd_tdb, dist_au, rel_vel_au_day
    FROM read_parquet({cat})
),
e1 AS (
    SELECT number AS number_1,
           a_au AS a_1, e AS e_1, i_deg AS i_1,
           node_deg AS node_1, argperi_deg AS argperi_1, mean_anomaly_deg AS mean_anomaly_1,
           epoch_jd AS epoch_1, h_mag AS h_1, q_au AS q_1
    FROM read_parquet({elem})
),
e2 AS (
    SELECT number AS number_2,
           a_au AS a_2, e AS e_2, i_deg AS i_2,
           node_deg AS node_2, argperi_deg AS argperi_2, mean_anomaly_deg AS mean_anomaly_2,
           epoch_jd AS epoch_2, h_mag AS h_2, q_au AS q_2
    FROM read_parquet({elem})
),
joined AS (
    SELECT
        cat.*,
        e1.a_1, e1.e_1, e1.i_1, e1.node_1, e1.argperi_1, e1.mean_anomaly_1,
        e1.epoch_1, e1.h_1, e1.q_1,
        e2.a_2, e2.e_2, e2.i_2, e2.node_2, e2.argperi_2, e2.mean_anomaly_2,
        e2.epoch_2, e2.h_2, e2.q_2,
        abs(e1.a_1 - e2.a_2) AS delta_a_au,
        greatest(e1.e_1, e2.e_2) AS e_max,
        greatest(e1.i_1, e2.i_2) AS i_max,
        least(e1.q_1, e2.q_2) AS q_min
    FROM cat
    JOIN e1 USING (number_1)
    JOIN e2 USING (number_2)
),
classified AS (
    SELECT
        *,
        CASE
            WHEN q_min < 1.8 AND e_max > 0.3 THEN 'q_min,e_max'
            WHEN q_min < 1.8 THEN 'q_min'
            WHEN e_max > 0.3 THEN 'e_max'
            ELSE 'none'
        END AS stageb_reason
    FROM joined
)
"""


def _summary_query(catalog: Path, elements_cache: Path) -> str:
    return _joined_cte(catalog, elements_cache) + """
SELECT
    count(*) AS n_total,
    sum(CASE WHEN stageb_reason != 'none' THEN 1 ELSE 0 END) AS n_stageb,
    sum(CASE WHEN q_min < 1.8 THEN 1 ELSE 0 END) AS n_q_min,
    sum(CASE WHEN e_max > 0.3 THEN 1 ELSE 0 END) AS n_e_max,
    sum(CASE WHEN q_min < 1.8 AND e_max > 0.3 THEN 1 ELSE 0 END) AS n_both,
    min(q_min) AS q_min_min,
    max(e_max) AS e_max_max,
    avg(CASE WHEN stageb_reason != 'none' THEN 1.0 ELSE 0.0 END) AS frac_stageb
FROM classified;
"""


def _select_query(catalog: Path, elements_cache: Path, *, limit: int | None) -> str:
    limit_clause = f"\nLIMIT {int(limit)}" if limit is not None else ""
    return _joined_cte(catalog, elements_cache) + f"""
SELECT
    number_1, number_2, designation_1, designation_2,
    jd_tdb, dist_au, rel_vel_au_day,
    a_1, e_1, i_1, node_1, argperi_1, mean_anomaly_1, epoch_1, h_1, q_1,
    a_2, e_2, i_2, node_2, argperi_2, mean_anomaly_2, epoch_2, h_2, q_2,
    delta_a_au, e_max, i_max, q_min, stageb_reason
FROM classified
WHERE stageb_reason != 'none'
ORDER BY number_1, number_2
{limit_clause}
"""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", type=Path, default=CATALOG)
    parser.add_argument("--mpcorb", type=Path, default=MPCORB)
    parser.add_argument("--elements-cache", type=Path, default=ELEMENTS_CACHE)
    parser.add_argument("--out", type=Path, default=OUT)
    parser.add_argument("--summary-only", action="store_true")
    parser.add_argument("--limit", type=int, default=None, help="Write only first N selected rows.")
    parser.add_argument("--force-elements-cache", action="store_true")
    args = parser.parse_args()

    _ensure_elements_cache(args.mpcorb, args.elements_cache, force=args.force_elements_cache)

    con = duckdb.connect()
    logger.info("Computing Stage B subset summary")
    summary = con.execute(_summary_query(args.catalog, args.elements_cache)).fetchdf()
    logger.info("\n%s", summary.to_string(index=False))

    if args.summary_only:
        return 0

    args.out.parent.mkdir(parents=True, exist_ok=True)
    logger.info("Writing Stage B subset to %s", args.out)
    select_sql = _select_query(args.catalog, args.elements_cache, limit=args.limit)
    con.execute(f"COPY ({select_sql}) TO {_sql_path(args.out)} (FORMAT PARQUET)")
    n_out = con.execute(f"SELECT count(*) FROM read_parquet({_sql_path(args.out)})").fetchone()[0]
    logger.info("Wrote %s (%d rows)", args.out, n_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
