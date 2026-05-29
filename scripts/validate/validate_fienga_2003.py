"""Cross-match the generated encounter catalog against Fienga et al. (2003).

Loads:
  - data/raw/fienga_2003_encounters.parquet  (predicted encounters, VizieR J/A+A/406/751)
  - data/output/encounters_catalog.parquet   (our pipeline output)

Procedure
---------
1. Restrict Fienga events to those with explicit ``Epoch.MGE`` (tables A.1, A.3).
2. Keep only events whose epoch falls inside the Gaia DR3 observation window.
3. Partition into two groups:

   * ``Impact <= threshold`` (default 0.01 AU) — events that *should* be in our
     catalog given the configured detection threshold.
   * ``Impact >  threshold`` — flyby predictions whose geometric distance is
     wider than our threshold; not expected in our catalog, reported for
     transparency.

4. For each "expected" event, look for a matching (Perturber, Target) pair in
   our catalog whose JD TDB lies within ``MATCH_TOLERANCE_DAYS`` of the
   predicted epoch.  Fienga epochs are at the 1st of the month, so the actual
   encounter can fall anywhere inside that month → tolerance ≥ 31 days.

5. Print per-event match/miss outcomes and a summary detection rate.

Usage
-----
    docker compose run --rm pipeline python -m scripts.validate_fienga_2003
    docker compose run --rm pipeline python -m scripts.validate_fienga_2003 \\
        --impact-threshold-au 0.01 --tolerance-days 31
    docker compose run --rm pipeline python -m scripts.validate_fienga_2003 \\
        --catalog data/output/encounters_catalog_hybrid_stageb.parquet
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


def _date_to_jd_tdb(date_str: str) -> float:
    """Convert an ISO date string (UTC midnight) to JD in TDB scale."""
    return float(Time(date_str, scale="utc").tdb.jd)


def _load_fienga(path: Path, win_start: str, win_end: str) -> pl.DataFrame:
    """Return Fienga rows with non-null ``Epoch.MGE`` inside the Gaia window."""
    df = pl.read_parquet(path).filter(pl.col("Epoch.MGE").is_not_null())
    df = df.with_columns(pl.col("Epoch.MGE").str.to_date("%Y-%m-%d").alias("epoch_date"))
    start = dt.date.fromisoformat(win_start)
    end = dt.date.fromisoformat(win_end)
    return df.filter((pl.col("epoch_date") >= start) & (pl.col("epoch_date") <= end))


def _find_match(
    catalog: pl.DataFrame, pert: int, targ: int, epoch_jd: float, tol_days: float
) -> pl.DataFrame:
    """Return rows of *catalog* matching (Perturber, Target) within tol_days of epoch_jd."""
    return catalog.filter(
        (
            ((pl.col("number_1") == pert) & (pl.col("number_2") == targ))
            | ((pl.col("number_1") == targ) & (pl.col("number_2") == pert))
        )
        & ((pl.col("jd_tdb") - epoch_jd).abs() <= tol_days)
    )


def _load_relevant_catalog(
    catalog_path: Path,
    expected: pl.DataFrame,
    tolerance_days: float,
) -> pl.DataFrame:
    """Load only catalog rows that can match the expected Fienga events."""
    schema = {
        "number_1": pl.Int64,
        "number_2": pl.Int64,
        "jd_tdb": pl.Float64,
        "dist_au": pl.Float64,
    }
    if len(expected) == 0:
        return pl.DataFrame(schema=schema)

    numbers: set[int] = set()
    epoch_jds: list[float] = []
    for row in expected.iter_rows(named=True):
        numbers.add(int(row["Perturber"]))
        numbers.add(int(row["Target"]))
        epoch_jds.append(_date_to_jd_tdb(row["epoch_date"].isoformat()))

    min_jd = min(epoch_jds) - tolerance_days
    max_jd = max(epoch_jds) + tolerance_days
    number_list = sorted(numbers)
    catalog = (
        pl.scan_parquet(catalog_path)
        .select(["number_1", "number_2", "jd_tdb", "dist_au"])
        .filter(
            (pl.col("number_1").is_in(number_list) | pl.col("number_2").is_in(number_list))
            & (pl.col("jd_tdb") >= min_jd)
            & (pl.col("jd_tdb") <= max_jd)
        )
        .collect()
    )
    logger.info("Relevant catalog slice: %d encounters from %s", len(catalog), catalog_path)
    return catalog


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument(
        "--catalog",
        type=Path,
        default=None,
        help="Encounter catalog parquet to validate. Defaults to the catalog path in config.",
    )
    parser.add_argument(
        "--impact-threshold-au",
        type=float,
        default=None,
        help="Impact parameter cutoff above which Fienga events are reported but not "
        "expected to appear in our catalog. Defaults to config detection.threshold_au.",
    )
    parser.add_argument(
        "--tolerance-days",
        type=float,
        default=31.0,
        help="Date-match window around the Fienga predicted epoch (Fienga gives "
        "monthly resolution, so ≥31 days makes sense). Default: 31.",
    )
    parser.add_argument(
        "--distance-tolerance-au",
        type=float,
        default=1e-4,
        help="Tolerance used to report distance agreement against Fienga. Default: 1e-4 AU.",
    )
    parser.add_argument(
        "--require-distance-tolerance",
        action="store_true",
        help="Return non-zero if any expected event is missing or outside --distance-tolerance-au.",
    )
    args = parser.parse_args()

    cfg = load_config(args.config)
    impact_thresh = (
        args.impact_threshold_au
        if args.impact_threshold_au is not None
        else cfg.detection.threshold_au
    )

    fienga_path = Path(cfg.paths.raw) / cfg.sources.fienga_2003.output_filename
    catalog_path = (
        args.catalog
        if args.catalog is not None
        else Path(cfg.paths.output) / f"{cfg.output.filename}.{cfg.output.format}"
    )
    if not fienga_path.exists():
        logger.error("Fienga catalog not found: %s — run download_fienga_2003.", fienga_path)
        return 1
    if not catalog_path.exists():
        logger.error("Pipeline catalog not found: %s — run run_pipeline.", catalog_path)
        return 1

    fienga = _load_fienga(fienga_path, cfg.time_window.start[:10], cfg.time_window.end[:10])
    logger.info(
        "Fienga 2003 events in Gaia window (%s → %s): %d",
        cfg.time_window.start[:10],
        cfg.time_window.end[:10],
        len(fienga),
    )

    expected = fienga.filter(pl.col("Impact") <= impact_thresh).sort("Impact")
    wider = fienga.filter(pl.col("Impact") > impact_thresh)
    logger.info(
        "Partition by Impact ≤ %.4f AU: %d expected in our catalog, %d wider flybys "
        "(not expected)",
        impact_thresh,
        len(expected),
        len(wider),
    )
    catalog = _load_relevant_catalog(catalog_path, expected, args.tolerance_days)

    matched: list[dict] = []
    missed: list[dict] = []

    logger.info("")
    logger.info("=== Cross-match: expected events (Impact ≤ %.4f AU) ===", impact_thresh)
    for row in expected.iter_rows(named=True):
        pert, targ = row["Perturber"], row["Target"]
        date_str = row["epoch_date"].isoformat()
        epoch_jd = _date_to_jd_tdb(date_str)
        hits = _find_match(catalog, pert, targ, epoch_jd, args.tolerance_days)

        if len(hits) == 0:
            logger.warning(
                "  ✗ MISS  (%d, %d)  Fienga=%s  Impact=%.5f AU  Vel=%.2f km/s",
                pert,
                targ,
                date_str,
                row["Impact"],
                row["Vel"],
            )
            missed.append(
                {
                    "perturber": pert,
                    "target": targ,
                    "fienga_date": date_str,
                    "fienga_impact_au": row["Impact"],
                    "fienga_vel_km_s": row["Vel"],
                }
            )
        else:
            best = hits.sort("dist_au").head(1).row(0, named=True)
            our_date = Time(best["jd_tdb"], format="jd", scale="tdb").utc.iso[:10]
            delta_days = best["jd_tdb"] - epoch_jd
            delta_dist_au = float(best["dist_au"]) - float(row["Impact"])
            abs_delta_dist_au = abs(delta_dist_au)
            within_distance_tolerance = abs_delta_dist_au <= args.distance_tolerance_au
            logger.info(
                "  ✓ HIT   (%d, %d)  Fienga=%s  Ours=%s (Δ=%+.1fd)  "
                "Impact_F=%.5f AU  dist_ours=%.5f AU  Δdist=%+.6f AU",
                pert,
                targ,
                date_str,
                our_date,
                delta_days,
                row["Impact"],
                best["dist_au"],
                delta_dist_au,
            )
            matched.append(
                {
                    "perturber": pert,
                    "target": targ,
                    "fienga_date": date_str,
                    "fienga_impact_au": row["Impact"],
                    "our_date": our_date,
                    "our_dist_au": best["dist_au"],
                    "delta_days": delta_days,
                    "delta_dist_au": delta_dist_au,
                    "abs_delta_dist_au": abs_delta_dist_au,
                    "within_distance_tolerance": within_distance_tolerance,
                }
            )

    n_exp = len(expected)
    n_hit = len(matched)
    rate = (100.0 * n_hit / n_exp) if n_exp > 0 else float("nan")
    logger.info("")
    logger.info("=== Summary ===")
    logger.info("Expected events (Impact ≤ %.4f AU): %d", impact_thresh, n_exp)
    logger.info("Matched in our catalog:              %d (%.1f%%)", n_hit, rate)
    logger.info("Missed:                              %d", len(missed))
    logger.info("Wider flybys reported only:          %d", len(wider))

    if matched:
        offsets = [abs(m["delta_days"]) for m in matched]
        logger.info(
            "Date offset |Δt|: median=%.1f d, max=%.1f d",
            sorted(offsets)[len(offsets) // 2],
            max(offsets),
        )
        dist_offsets = sorted(abs(m["delta_dist_au"]) for m in matched)
        n_within = sum(bool(m["within_distance_tolerance"]) for m in matched)
        logger.info(
            "Distance offset |Δdist|: median=%.6f AU, max=%.6f AU, within %.1e AU: %d/%d",
            dist_offsets[len(dist_offsets) // 2],
            max(dist_offsets),
            args.distance_tolerance_au,
            n_within,
            len(matched),
        )
        distance_failures = [m for m in matched if not bool(m["within_distance_tolerance"])]
        if distance_failures:
            logger.warning(
                "Distance tolerance failures (> %.1e AU): %d",
                args.distance_tolerance_au,
                len(distance_failures),
            )
    else:
        distance_failures = []

    # Save full match report for downstream analysis
    report_dir = Path(cfg.paths.output)
    if matched:
        pl.DataFrame(matched).write_csv(report_dir / "fienga_2003_matches.csv")
    if missed:
        pl.DataFrame(missed).write_csv(report_dir / "fienga_2003_misses.csv")
    logger.info("Match/miss CSVs written to %s/", report_dir)

    if args.require_distance_tolerance and (missed or distance_failures):
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
