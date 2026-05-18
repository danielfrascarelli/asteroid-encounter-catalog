"""Cross-match the generated encounter catalog against Galád & Gray (2002).

Loads:
  - data/raw/galad_2002_encounters.parquet  (predicted encounters, parsed from
    the published HTML of A&A 391, 1115)
  - data/output/encounters_catalog.parquet  (our pipeline output)

Procedure
---------
1. Restrict Galád events to those with a parseable date falling inside the
   Gaia DR3 observation window.
2. Partition into two groups:

   * ``r <= threshold`` (default 0.01 AU)  — events that *should* be in our
     catalog given the configured detection threshold.
   * ``r >  threshold``                    — wider flybys reported for
     transparency but not expected in our catalog.

3. For each "expected" event, look for a matching (perturber, target) pair
   in our catalog whose JD TDB lies within a tolerance window of the Galád
   epoch.  Tolerance defaults to ±31 days for day-precision dates and
   ±365 days for year-only dates (Galád lists some entries with only a
   year — the parser sets July-1st as a midpoint and tags the precision).

4. Print per-event hit/miss outcomes and a summary detection rate.

Outputs
-------
    data/output/galad_2002_matches.csv
    data/output/galad_2002_misses.csv

Usage
-----
    docker compose run --rm pipeline python -m scripts.validate_galad_2002
    docker compose run --rm pipeline python -m scripts.validate_galad_2002 \\
        --impact-threshold-au 0.01 --tolerance-days 31 --year-tolerance-days 365
"""

from __future__ import annotations

import argparse
import datetime as dt
import logging
import sys
from pathlib import Path

import polars as pl
from astropy.time import Time

from src.utils.config import load_config

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)


def _date_to_jd_tdb(date_obj: dt.date) -> float:
    """Convert a calendar date (UTC midnight) to JD in TDB scale."""
    return float(Time(date_obj.isoformat(), scale="utc").tdb.jd)


def _load_galad(path: Path, win_start: str, win_end: str) -> pl.DataFrame:
    """Return Galád rows with a parseable date inside the Gaia window.

    Both perturber and target numbers must be present — pure "name-only"
    rows cannot be cross-matched against the integer-keyed catalog.
    """
    df = pl.read_parquet(path).filter(
        pl.col("date_parsed").is_not_null()
        & pl.col("perturber_number").is_not_null()
        & pl.col("target_number").is_not_null()
        & pl.col("r_au").is_not_null()
    )
    start = dt.date.fromisoformat(win_start)
    end = dt.date.fromisoformat(win_end)
    return df.filter((pl.col("date_parsed") >= start) & (pl.col("date_parsed") <= end))


def _find_match(
    catalog: pl.DataFrame, pert: int, targ: int, epoch_jd: float, tol_days: float
) -> pl.DataFrame:
    """Return catalog rows for (perturber, target) within ``tol_days`` of ``epoch_jd``."""
    return catalog.filter(
        (
            ((pl.col("number_1") == pert) & (pl.col("number_2") == targ))
            | ((pl.col("number_1") == targ) & (pl.col("number_2") == pert))
        )
        & ((pl.col("jd_tdb") - epoch_jd).abs() <= tol_days)
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument(
        "--impact-threshold-au",
        type=float,
        default=None,
        help="r cutoff above which Galád events are reported but not "
        "expected to appear in our catalog. Defaults to config detection.threshold_au.",
    )
    parser.add_argument(
        "--tolerance-days",
        type=float,
        default=31.0,
        help="Date-match window around the Galád epoch for day-precision entries. "
        "Default: 31.",
    )
    parser.add_argument(
        "--year-tolerance-days",
        type=float,
        default=365.0,
        help="Date-match window for year-only Galád entries. Default: 365.",
    )
    args = parser.parse_args()

    cfg = load_config(args.config)
    impact_thresh = (
        args.impact_threshold_au
        if args.impact_threshold_au is not None
        else cfg.detection.threshold_au
    )

    galad_path = Path(cfg.paths.raw) / cfg.sources.galad_2002.output_filename
    catalog_path = Path(cfg.paths.output) / f"{cfg.output.filename}.{cfg.output.format}"
    if not galad_path.exists():
        logger.error("Galád catalog not found: %s — run download_galad_2002.", galad_path)
        return 1
    if not catalog_path.exists():
        logger.error("Pipeline catalog not found: %s — run run_pipeline.", catalog_path)
        return 1

    catalog = pl.read_parquet(catalog_path)
    logger.info("Pipeline catalog: %d encounters", len(catalog))

    galad = _load_galad(galad_path, cfg.time_window.start[:10], cfg.time_window.end[:10])
    logger.info(
        "Galád 2002 events in Gaia window (%s → %s): %d",
        cfg.time_window.start[:10],
        cfg.time_window.end[:10],
        len(galad),
    )

    if len(galad) == 0:
        logger.warning("No Galád events fall inside the Gaia window — nothing to match.")
        return 0

    expected = galad.filter(pl.col("r_au") <= impact_thresh).sort("r_au")
    wider = galad.filter(pl.col("r_au") > impact_thresh)
    logger.info(
        "Partition by r ≤ %.4f AU: %d expected in our catalog, %d wider flybys (not expected)",
        impact_thresh,
        len(expected),
        len(wider),
    )

    matched: list[dict] = []
    missed: list[dict] = []

    logger.info("")
    logger.info("=== Cross-match: expected events (r ≤ %.4f AU) ===", impact_thresh)
    for row in expected.iter_rows(named=True):
        pert = int(row["perturber_number"])
        targ = int(row["target_number"])
        date_obj = row["date_parsed"]
        precision = row["date_precision"]
        date_str = date_obj.isoformat()
        epoch_jd = _date_to_jd_tdb(date_obj)
        tol = args.year_tolerance_days if precision == "year" else args.tolerance_days
        precision_tag = "" if precision == "day" else "  [low-precision date]"

        hits = _find_match(catalog, pert, targ, epoch_jd, tol)

        if len(hits) == 0:
            logger.warning(
                "  ✗ MISS  (%d, %d)  Galád=%s  r=%.5f AU  v=%.2f km/s  P=%.1f km/s%s",
                pert,
                targ,
                date_str,
                row["r_au"],
                row["v_km_s"] if row["v_km_s"] is not None else float("nan"),
                row["p_km_s"] if row["p_km_s"] is not None else float("nan"),
                precision_tag,
            )
            missed.append(
                {
                    "perturber": pert,
                    "perturber_name": row["perturber_name"],
                    "target": targ,
                    "target_name": row["target_name"],
                    "galad_date": date_str,
                    "date_precision": precision,
                    "galad_r_au": row["r_au"],
                    "galad_v_km_s": row["v_km_s"],
                    "galad_p_km_s": row["p_km_s"],
                    "source_table": row["source_table"],
                }
            )
        else:
            best = hits.sort("dist_au").head(1).row(0, named=True)
            our_date = Time(best["jd_tdb"], format="jd", scale="tdb").utc.iso[:10]
            delta_days = best["jd_tdb"] - epoch_jd
            logger.info(
                "  ✓ HIT   (%d, %d)  Galád=%s  Ours=%s (Δ=%+.1fd)  "
                "r_G=%.5f AU  dist_ours=%.5f AU%s",
                pert,
                targ,
                date_str,
                our_date,
                delta_days,
                row["r_au"],
                best["dist_au"],
                precision_tag,
            )
            matched.append(
                {
                    "perturber": pert,
                    "perturber_name": row["perturber_name"],
                    "target": targ,
                    "target_name": row["target_name"],
                    "galad_date": date_str,
                    "date_precision": precision,
                    "galad_r_au": row["r_au"],
                    "galad_v_km_s": row["v_km_s"],
                    "galad_p_km_s": row["p_km_s"],
                    "source_table": row["source_table"],
                    "our_date": our_date,
                    "our_dist_au": best["dist_au"],
                    "delta_days": delta_days,
                }
            )

    n_exp = len(expected)
    n_hit = len(matched)
    rate = (100.0 * n_hit / n_exp) if n_exp > 0 else float("nan")
    logger.info("")
    logger.info("=== Summary ===")
    logger.info("Galád events in Gaia window:          %d", len(galad))
    logger.info("  with day-precision dates:           %d", galad.filter(pl.col("date_precision") == "day").height)
    logger.info("  with year-only dates:               %d", galad.filter(pl.col("date_precision") == "year").height)
    logger.info("Expected events (r ≤ %.4f AU):     %d", impact_thresh, n_exp)
    logger.info("Matched in our catalog:               %d (%.1f%%)", n_hit, rate)
    logger.info("Missed:                               %d", len(missed))
    logger.info("Wider flybys reported only:           %d", len(wider))

    if matched:
        offsets = sorted(abs(m["delta_days"]) for m in matched)
        logger.info(
            "Date offset |Δt|: median=%.1f d, max=%.1f d",
            offsets[len(offsets) // 2],
            max(offsets),
        )

    # Save full match report for downstream analysis.  CSVs are always
    # emitted (possibly with only a header row) so consumers can rely on
    # their existence after a successful run.
    report_dir = Path(cfg.paths.output)
    matches_schema = {
        "perturber": pl.Int64,
        "perturber_name": pl.Utf8,
        "target": pl.Int64,
        "target_name": pl.Utf8,
        "galad_date": pl.Utf8,
        "date_precision": pl.Utf8,
        "galad_r_au": pl.Float64,
        "galad_v_km_s": pl.Float64,
        "galad_p_km_s": pl.Float64,
        "source_table": pl.Utf8,
        "our_date": pl.Utf8,
        "our_dist_au": pl.Float64,
        "delta_days": pl.Float64,
    }
    misses_schema = {
        "perturber": pl.Int64,
        "perturber_name": pl.Utf8,
        "target": pl.Int64,
        "target_name": pl.Utf8,
        "galad_date": pl.Utf8,
        "date_precision": pl.Utf8,
        "galad_r_au": pl.Float64,
        "galad_v_km_s": pl.Float64,
        "galad_p_km_s": pl.Float64,
        "source_table": pl.Utf8,
    }
    pl.DataFrame(matched, schema=matches_schema).write_csv(
        report_dir / "galad_2002_matches.csv"
    )
    pl.DataFrame(missed, schema=misses_schema).write_csv(
        report_dir / "galad_2002_misses.csv"
    )
    logger.info("Match/miss CSVs written to %s/", report_dir)

    return 0


if __name__ == "__main__":
    sys.exit(main())
