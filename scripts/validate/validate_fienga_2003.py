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


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--config", default="config.yaml")
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
    args = parser.parse_args()

    cfg = load_config(args.config)
    impact_thresh = (
        args.impact_threshold_au
        if args.impact_threshold_au is not None
        else cfg.detection.threshold_au
    )

    fienga_path = Path(cfg.paths.raw) / cfg.sources.fienga_2003.output_filename
    catalog_path = Path(cfg.paths.output) / f"{cfg.output.filename}.{cfg.output.format}"
    if not fienga_path.exists():
        logger.error("Fienga catalog not found: %s — run download_fienga_2003.", fienga_path)
        return 1
    if not catalog_path.exists():
        logger.error("Pipeline catalog not found: %s — run run_pipeline.", catalog_path)
        return 1

    catalog = pl.read_parquet(catalog_path)
    logger.info("Pipeline catalog: %d encounters", len(catalog))

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
            logger.info(
                "  ✓ HIT   (%d, %d)  Fienga=%s  Ours=%s (Δ=%+.1fd)  "
                "Impact_F=%.5f AU  dist_ours=%.5f AU",
                pert,
                targ,
                date_str,
                our_date,
                delta_days,
                row["Impact"],
                best["dist_au"],
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

    # Save full match report for downstream analysis
    report_dir = Path(cfg.paths.output)
    if matched:
        pl.DataFrame(matched).write_csv(report_dir / "fienga_2003_matches.csv")
    if missed:
        pl.DataFrame(missed).write_csv(report_dir / "fienga_2003_misses.csv")
    logger.info("Match/miss CSVs written to %s/", report_dir)

    return 0


if __name__ == "__main__":
    sys.exit(main())
