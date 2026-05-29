"""Assemble the Track 1 Stage B hybrid Kepler/N-body catalog.

Inputs
------
- Frozen detection catalog: Kepler-refined geometric values.
- Stage B N-body shards produced by ``refine_stageb_nbody.py``.

Output
------
A hybrid parquet where:

- rows present in the Stage B shard set use N-body values in the canonical
  ``jd_tdb``, ``dist_au`` and ``rel_vel_au_day`` columns;
- rows not present in the shard set keep the original Kepler values;
- audit columns preserve both values and identify the method:
  ``refinement_method``, ``*_kepler``, ``*_nbody``, ``delta_*``,
  ``nbody_converged``, ``nbody_energy_drift``, ``near_boundary``.

By default the script is strict: it verifies that every row in the selected
Stage B subset has a successful N-body shard result before writing the full
hybrid catalog. Use ``--allow-incomplete`` only for smoke tests.
"""

from __future__ import annotations

import argparse
import json
import logging
from datetime import UTC, datetime
from pathlib import Path

import duckdb

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
logger = logging.getLogger(__name__)

CATALOG = Path("data/output/encounters_catalog_rebound_005au.parquet")
SUBSET = Path("data/cache/nbody_validation/stageb_selective_subset.parquet")
SHARDS_GLOB = "data/output/stageb_nbody_shards/*.parquet"
OUT = Path("data/output/encounters_catalog_hybrid_stageb.parquet")


def _sql_path(path: str | Path) -> str:
    return "'" + str(path).replace("'", "''") + "'"


def _read_glob(glob: str) -> str:
    return f"read_parquet({_sql_path(glob)}, union_by_name=true)"


def _deduped_nbody_cte(shards_glob: str) -> str:
    shards = _read_glob(shards_glob)
    return f"""
raw_nbody AS (
    SELECT * FROM {shards}
),
nbody AS (
    SELECT * EXCLUDE (dedupe_rank)
    FROM (
        SELECT
            raw_nbody.*,
            row_number() OVER (
                PARTITION BY number_1, number_2, t_min_kepler_jd
                ORDER BY
                    CASE WHEN error_message IS NULL THEN 0 ELSE 1 END,
                    CASE WHEN nbody_converged THEN 0 ELSE 1 END,
                    row_index DESC
            ) AS dedupe_rank
        FROM raw_nbody
    )
    WHERE dedupe_rank = 1
)"""


def _shard_summary_query(shards_glob: str) -> str:
    nbody_cte = _deduped_nbody_cte(shards_glob)
    return f"""
WITH
{nbody_cte}
SELECT
    (SELECT count(*) FROM raw_nbody) AS n_rows_raw,
    count(*) AS n_rows,
    (SELECT count(*) FROM raw_nbody) - count(*) AS n_duplicate_rows,
    sum(CASE WHEN error_message IS NOT NULL THEN 1 ELSE 0 END) AS n_failed,
    sum(CASE WHEN NOT nbody_converged THEN 1 ELSE 0 END) AS n_unconverged,
    max(nbody_energy_drift) AS max_energy_drift,
    max(abs(delta_dist_au)) AS max_abs_delta_dist_au,
    sum(CASE WHEN near_boundary THEN 1 ELSE 0 END) AS n_near_boundary
FROM nbody;
"""


def _coverage_query(subset: Path, shards_glob: str) -> str:
    subset_sql = _sql_path(subset)
    nbody_cte = _deduped_nbody_cte(shards_glob)
    return f"""
WITH
subset AS (
    SELECT number_1, number_2, jd_tdb
    FROM read_parquet({subset_sql})
),
{nbody_cte}
SELECT
    (SELECT count(*) FROM subset) AS n_subset,
    count(nbody.number_1) AS n_matched,
    sum(CASE WHEN nbody.error_message IS NOT NULL THEN 1 ELSE 0 END) AS n_failed,
    sum(CASE WHEN nbody.number_1 IS NOT NULL AND NOT nbody.nbody_converged THEN 1 ELSE 0 END)
        AS n_unconverged,
    (SELECT count(*) FROM subset) - count(nbody.number_1) AS n_missing
FROM subset
LEFT JOIN nbody
    ON subset.number_1 = nbody.number_1
   AND subset.number_2 = nbody.number_2
   AND subset.jd_tdb = nbody.t_min_kepler_jd;
"""


def _assemble_query(catalog: Path, shards_glob: str) -> str:
    catalog_sql = _sql_path(catalog)
    nbody_cte = _deduped_nbody_cte(shards_glob)
    return f"""
WITH
cat AS (
    SELECT * FROM read_parquet({catalog_sql})
),
{nbody_cte}
SELECT
    cat.number_1,
    cat.number_2,
    cat.designation_1,
    cat.designation_2,
    COALESCE(nbody.t_min_nbody_jd, cat.jd_tdb) AS jd_tdb,
    COALESCE(nbody.dist_au_nbody, cat.dist_au) AS dist_au,
    COALESCE(nbody.rel_vel_nbody, cat.rel_vel_au_day) AS rel_vel_au_day,
    CASE WHEN nbody.number_1 IS NULL THEN 'kepler' ELSE 'nbody' END AS refinement_method,
    cat.jd_tdb AS t_min_kepler_jd,
    cat.dist_au AS dist_au_kepler,
    cat.rel_vel_au_day AS rel_vel_kepler,
    nbody.t_min_nbody_jd,
    nbody.dist_au_nbody,
    nbody.rel_vel_nbody,
    nbody.delta_dist_au,
    nbody.delta_t_min_hours,
    nbody.delta_rel_vel_au_day,
    nbody.nbody_converged,
    nbody.nbody_energy_drift,
    nbody.n_samples,
    COALESCE(nbody.near_boundary, false) AS near_boundary,
    nbody.error_message
FROM cat
LEFT JOIN nbody
    ON cat.number_1 = nbody.number_1
   AND cat.number_2 = nbody.number_2
   AND cat.jd_tdb = nbody.t_min_kepler_jd
"""


def _write_sidecar(out: Path, summary: dict, coverage: dict, args: argparse.Namespace) -> None:
    sidecar = out.parent / (out.stem + "_provenance.json")
    meta = {
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "catalog": str(args.catalog),
        "subset": str(args.subset),
        "shards_glob": args.shards_glob,
        "output": str(out),
        "allow_incomplete": args.allow_incomplete,
        "shard_summary": summary,
        "stageb_coverage": coverage,
    }
    sidecar.write_text(json.dumps(meta, indent=2, default=str))
    logger.info("Wrote provenance sidecar: %s", sidecar)


def _fetch_one_dict(con: duckdb.DuckDBPyConnection, query: str) -> dict:
    """Execute *query* and return the first row as a plain dict."""
    cursor = con.execute(query)
    columns = [desc[0] for desc in cursor.description]
    row = cursor.fetchone()
    if row is None:
        raise RuntimeError("Query returned no rows")
    return dict(zip(columns, row, strict=True))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", type=Path, default=CATALOG)
    parser.add_argument("--subset", type=Path, default=SUBSET)
    parser.add_argument("--shards-glob", default=SHARDS_GLOB)
    parser.add_argument("--out", type=Path, default=OUT)
    parser.add_argument("--allow-incomplete", action="store_true")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    if args.out.exists() and not args.force:
        raise SystemExit(f"Output exists; pass --force to overwrite: {args.out}")

    con = duckdb.connect()
    logger.info("Reading shard summary from %s", args.shards_glob)
    summary = _fetch_one_dict(con, _shard_summary_query(args.shards_glob))
    logger.info("Shard summary: %s", summary)

    coverage = _fetch_one_dict(con, _coverage_query(args.subset, args.shards_glob))
    logger.info("Stage B coverage: %s", coverage)

    incomplete = (
        coverage["n_missing"] != 0 or coverage["n_failed"] != 0 or coverage["n_unconverged"] != 0
    )
    if incomplete and not args.allow_incomplete:
        raise SystemExit(
            "Stage B shards are incomplete or contain failed/unconverged rows: "
            f"{coverage}. Pass --allow-incomplete only for smoke/testing."
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    logger.info("Writing hybrid catalog: %s", args.out)
    con.execute(
        f"COPY ({_assemble_query(args.catalog, args.shards_glob)}) TO {_sql_path(args.out)} (FORMAT PARQUET)"
    )
    n_out = con.execute(f"SELECT count(*) FROM read_parquet({_sql_path(args.out)})").fetchone()[0]
    logger.info("Wrote %s (%d rows)", args.out, n_out)
    _write_sidecar(args.out, summary, coverage, args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
