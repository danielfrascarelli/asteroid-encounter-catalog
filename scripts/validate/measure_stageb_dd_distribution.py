"""Characterise the Kepler <-> N-body refined-distance discrepancy (tribunal M6).

Catalogue membership requires the Kepler-refined distance ``dist_au_kepler``
to fall below the catalogue threshold (0.05 AU), but a pair only becomes a
*candidate* in the first place if the coarse N-body scan placed it within the
query radius ``r_q = 0.0572 AU``. The margin between the two,

    margin = r_q - threshold = 0.0572 - 0.05 = 7.2 mAU,

is the safety buffer against candidates being missed due to Kepler/N-body
disagreement. If the discrepancy

    delta_dist_au = dist_au_nbody - dist_au_kepler

exceeds that margin for a pair whose true (N-body) distance is just inside
the threshold, the pair could have been dropped before ever being scored by
Kepler refinement.

This script streams the Stage-B N-body re-refinement shards (one row per
re-refined pair, produced by ``scripts/validate/refine_stageb_nbody.py`` /
the Stage-B pipeline) through DuckDB and reports the distribution of
``delta_dist_au`` for converged rows, both overall and restricted to the
near-threshold band ``dist_au_kepler in [0.045, 0.05)`` AU that matters most
for M6. It is read-only: shard files are only ever read via
``read_parquet(glob)`` and are never opened for writing, so it is safe to run
concurrently against shards still being produced by an in-progress Stage-B
run.

Usage
-----
    docker compose run --rm pipeline python -m scripts.validate.measure_stageb_dd_distribution \\
        --shards-glob "data/output/stageb_nbody_shards_b1fix/*.parquet"
"""

from __future__ import annotations

import argparse
import glob
import logging

import duckdb

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
logger = logging.getLogger(__name__)

DEFAULT_SHARDS_GLOB = "data/output/stageb_nbody_shards_b1fix/*.parquet"
CATALOG_THRESHOLD_AU = 0.05
DEFAULT_RQ_AU = 0.0572
NEAR_THRESHOLD_LO_AU = 0.045
NEAR_THRESHOLD_HI_AU = 0.05


def _stats_query(con: duckdb.DuckDBPyConnection, shards_glob: str, extra_where: str) -> dict:
    """Run one streaming aggregation over the shard glob and return summary stats.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open DuckDB connection (in-memory; no data is persisted).
    shards_glob : str
        Glob pattern passed straight to DuckDB's ``read_parquet``.
    extra_where : str
        Additional SQL boolean expression ANDed onto the base
        ``nbody_converged AND delta_dist_au IS NOT NULL`` filter (e.g. to
        restrict to the near-threshold Kepler-distance band). Must reference
        only columns present in the shard schema.

    Returns
    -------
    dict
        Keys: ``n``, ``mean``, ``std``, ``median``, ``p50_abs``, ``p90_abs``,
        ``p99_abs``, ``p999_abs``, ``max_abs``, ``n_exceed_margin``. All
        aggregates are computed by DuckDB in a single streaming pass over the
        parquet row groups (no full materialisation into Python).
    """
    query = f"""
        SELECT
            COUNT(*) AS n,
            AVG(delta_dist_au) AS mean,
            STDDEV_SAMP(delta_dist_au) AS std,
            MEDIAN(delta_dist_au) AS median,
            QUANTILE_CONT(ABS(delta_dist_au), 0.50) AS p50_abs,
            QUANTILE_CONT(ABS(delta_dist_au), 0.90) AS p90_abs,
            QUANTILE_CONT(ABS(delta_dist_au), 0.99) AS p99_abs,
            QUANTILE_CONT(ABS(delta_dist_au), 0.999) AS p999_abs,
            MAX(ABS(delta_dist_au)) AS max_abs,
            SUM(CASE WHEN ABS(delta_dist_au) > $margin THEN 1 ELSE 0 END) AS n_exceed_margin
        FROM read_parquet($shards_glob)
        WHERE nbody_converged
          AND delta_dist_au IS NOT NULL
          AND ({extra_where})
    """
    row = con.execute(query, {"margin": MARGIN_AU, "shards_glob": shards_glob}).fetchone()
    cols = [d[0] for d in con.description]
    return dict(zip(cols, row))


def _row_counts(con: duckdb.DuckDBPyConnection, shards_glob: str) -> tuple[int, int]:
    """Return (n_total_rows, n_converged_rows) across all matched shards."""
    query = """
        SELECT
            COUNT(*) AS n_total,
            SUM(CASE WHEN nbody_converged THEN 1 ELSE 0 END) AS n_converged
        FROM read_parquet($shards_glob)
    """
    n_total, n_converged = con.execute(query, {"shards_glob": shards_glob}).fetchone()
    return int(n_total), int(n_converged or 0)


def _log_stats(label: str, stats: dict, margin_au: float) -> None:
    """Log one stats block in a consistent, human-readable format."""
    n = int(stats["n"])
    if n == 0:
        logger.info("[%s] n=0 rows (no converged rows in this subset)", label)
        return
    frac_exceed = stats["n_exceed_margin"] / n
    logger.info(
        "[%s] n=%d  mean(delta)=%+.6f AU  std(delta)=%.6f AU  median(delta)=%+.6f AU",
        label,
        n,
        stats["mean"],
        stats["std"],
        stats["median"],
    )
    logger.info(
        "[%s] |delta| percentiles: p50=%.6f  p90=%.6f  p99=%.6f  p99.9=%.6f  max=%.6f  (AU)",
        label,
        stats["p50_abs"],
        stats["p90_abs"],
        stats["p99_abs"],
        stats["p999_abs"],
        stats["max_abs"],
    )
    logger.info(
        "[%s] fraction |delta| > margin (%.4f AU): %d / %d = %.6f (%.4f%%)",
        label,
        margin_au,
        int(stats["n_exceed_margin"]),
        n,
        frac_exceed,
        100.0 * frac_exceed,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--shards-glob",
        type=str,
        default=DEFAULT_SHARDS_GLOB,
        help="Glob pattern for Stage-B N-body shard parquet files (read-only).",
    )
    parser.add_argument(
        "--threshold-au",
        type=float,
        default=CATALOG_THRESHOLD_AU,
        help="Catalogue membership threshold on Kepler-refined distance (AU).",
    )
    parser.add_argument(
        "--rq-au",
        type=float,
        default=DEFAULT_RQ_AU,
        help="Coarse N-body scan query radius r_q (AU); margin = rq - threshold.",
    )
    parser.add_argument(
        "--near-lo-au",
        type=float,
        default=NEAR_THRESHOLD_LO_AU,
        help="Lower bound (inclusive) of the near-threshold Kepler-distance band (AU).",
    )
    parser.add_argument(
        "--near-hi-au",
        type=float,
        default=NEAR_THRESHOLD_HI_AU,
        help="Upper bound (exclusive) of the near-threshold Kepler-distance band (AU).",
    )
    args = parser.parse_args()

    global MARGIN_AU
    MARGIN_AU = args.rq_au - args.threshold_au

    n_files = len(glob.glob(args.shards_glob))
    logger.info(
        "shards_glob=%r matched %d files -- PRELIMINARY (partial shard set; Stage-B run may still be in progress)",
        args.shards_glob,
        n_files,
    )
    if n_files == 0:
        logger.error("No shard files matched %r; aborting.", args.shards_glob)
        return 1

    logger.info(
        "threshold=%.4f AU  r_q=%.4f AU  margin=r_q-threshold=%.6f AU (%.1f mAU)",
        args.threshold_au,
        args.rq_au,
        MARGIN_AU,
        1000.0 * MARGIN_AU,
    )
    logger.info(
        "near-threshold band: dist_au_kepler in [%.4f, %.4f) AU", args.near_lo_au, args.near_hi_au
    )

    con = duckdb.connect()

    n_total, n_converged = _row_counts(con, args.shards_glob)
    logger.info(
        "row counts across matched shards: total=%d  converged=%d (%.2f%%)",
        n_total,
        n_converged,
        100.0 * n_converged / max(1, n_total),
    )

    overall = _stats_query(con, args.shards_glob, extra_where="TRUE")
    _log_stats("ALL converged", overall, MARGIN_AU)

    near_where = (
        f"dist_au_kepler >= {args.near_lo_au!r} AND dist_au_kepler < {args.near_hi_au!r}"
    )
    near = _stats_query(con, args.shards_glob, extra_where=near_where)
    _log_stats(
        f"near-threshold [{args.near_lo_au},{args.near_hi_au})",
        near,
        MARGIN_AU,
    )

    logger.info("PRELIMINARY result -- partial shard set (%d files); re-run once Stage-B completes.", n_files)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
