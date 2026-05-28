"""Cross-match the generated encounter catalog against Goffin (2014).

Loads:
  - data/raw/goffin_2014_encounters.parquet  (encounters used for mass determinations,
    VizieR J/A+A/565/A56)
  - data/output/encounters_catalog.parquet   (our pipeline output)

Procedure
---------
1. Inspect the downloaded table for epoch and distance columns (VizieR column
   names vary between catalog versions; the script adapts automatically).
2. Keep only events whose epoch falls inside the Gaia DR3 observation window.
3. Partition into two groups:

   * ``dist <= threshold`` (default 0.01 AU) — events that *should* be in our
     catalog given the configured detection threshold.
   * ``dist >  threshold`` — encounters whose geometric distance exceeds our
     threshold; not expected in our catalog, reported for transparency.

4. For each "expected" event, look for a matching (Perturber, Target) pair in
   our catalog whose JD TDB lies within ``MATCH_TOLERANCE_DAYS`` of the
   predicted epoch.

5. Print per-event match/miss outcomes and a summary detection rate.

Usage
-----
    docker compose run --rm pipeline python -m scripts.validate_goffin_2014
    docker compose run --rm pipeline python -m scripts.validate_goffin_2014 \\
        --dist-threshold-au 0.01 --tolerance-days 31
    docker compose run --rm pipeline python -m scripts.validate_goffin_2014 \\
        --catalog data/output/encounters_catalog_hybrid_stageb.parquet
"""

from __future__ import annotations

import argparse
import datetime
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

# Candidate column names for each semantic field, in preference order.
# VizieR column names vary slightly between catalog versions.
_PERTURBER_COLS = ["Pert", "Perturber", "pert", "perturber", "Num1", "num1"]
_TARGET_COLS = ["Targ", "Target", "targ", "target", "Num2", "num2"]
_EPOCH_COLS = ["Date", "Epoch", "date", "epoch", "Date1", "TDB"]
_DIST_COLS = ["Dist", "Delta", "dist", "delta", "dmin", "Dmin", "Impact", "b"]
_VEL_COLS = ["Vrel", "Vel", "vrel", "vel", "V", "v"]


def _pick_col(df: pl.DataFrame, candidates: list[str]) -> str | None:
    """Return the first candidate column name that exists in *df*, or None."""
    for c in candidates:
        if c in df.columns:
            return c
    return None


def _date_to_jd_tdb(date_str: str) -> float:
    """Convert an ISO date string (UTC midnight) to JD in TDB scale."""
    return float(Time(date_str, scale="utc").tdb.jd)


def _load_goffin(path: Path, win_start: str, win_end: str) -> tuple[pl.DataFrame, dict[str, str]]:
    """Return Goffin rows inside the Gaia window, plus a column-name map.

    Returns
    -------
    (filtered_df, col_map)
        col_map keys: 'perturber', 'target', 'epoch', 'dist'  (all present)
        and optionally 'vel'.
    """
    df = pl.read_parquet(path)
    logger.info("Goffin raw columns: %s", df.columns)

    perturber_col = _pick_col(df, _PERTURBER_COLS)
    target_col = _pick_col(df, _TARGET_COLS)
    epoch_col = _pick_col(df, _EPOCH_COLS)
    dist_col = _pick_col(df, _DIST_COLS)
    vel_col = _pick_col(df, _VEL_COLS)

    missing = [
        name
        for name, col in [
            ("perturber", perturber_col),
            ("target", target_col),
            ("epoch", epoch_col),
            ("dist", dist_col),
        ]
        if col is None
    ]
    if missing:
        raise RuntimeError(
            f"Cannot find required columns {missing} in Goffin parquet. " f"Available: {df.columns}"
        )

    assert perturber_col and target_col and epoch_col and dist_col  # mypy

    col_map: dict[str, str] = {
        "perturber": perturber_col,
        "target": target_col,
        "epoch": epoch_col,
        "dist": dist_col,
    }
    if vel_col:
        col_map["vel"] = vel_col

    # Drop rows with null epoch or perturber/target
    df = df.filter(
        pl.col(epoch_col).is_not_null()
        & pl.col(perturber_col).is_not_null()
        & pl.col(target_col).is_not_null()
    )

    # Parse epoch — may be a string date or a numeric year
    epoch_series = df[epoch_col]
    if epoch_series.dtype == pl.Utf8:
        df = df.with_columns(
            pl.col(epoch_col).str.to_date("%Y-%m-%d", strict=False).alias("_epoch_date")
        )
    elif epoch_series.dtype in (pl.Float32, pl.Float64, pl.Int32, pl.Int64):
        # Assume decimal year (e.g., 2015.32) → convert to approximate date
        def _decimal_year_to_date(y: float) -> datetime.date:
            year = int(y)
            frac = y - year
            start = datetime.date(year, 1, 1)
            end = datetime.date(year + 1, 1, 1)
            delta = (end - start).days * frac
            return start + datetime.timedelta(days=delta)

        dates = [
            _decimal_year_to_date(float(v)) if v is not None else None
            for v in epoch_series.to_list()
        ]
        df = df.with_columns(pl.Series("_epoch_date", dates, dtype=pl.Date))
    else:
        raise RuntimeError(f"Unrecognised epoch column dtype: {epoch_series.dtype}")

    start = datetime.date.fromisoformat(win_start)
    end = datetime.date.fromisoformat(win_end)
    df = df.filter(
        pl.col("_epoch_date").is_not_null()
        & (pl.col("_epoch_date") >= start)
        & (pl.col("_epoch_date") <= end)
    )
    return df, col_map


def _find_match(
    catalog: pl.DataFrame, pert: int, targ: int, epoch_jd: float, tol_days: float
) -> pl.DataFrame:
    """Return catalog rows matching (Perturber, Target) within tol_days of epoch_jd."""
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
    pert_col: str,
    targ_col: str,
    tolerance_days: float,
) -> pl.DataFrame:
    """Load only catalog rows that can match the expected Goffin events."""
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
        numbers.add(int(row[pert_col]))
        numbers.add(int(row[targ_col]))
        epoch_jds.append(_date_to_jd_tdb(row["_epoch_date"].isoformat()))

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
        "--dist-threshold-au",
        type=float,
        default=None,
        help="Distance cutoff: Goffin events closer than this are expected in our catalog. "
        "Defaults to config detection.threshold_au.",
    )
    parser.add_argument(
        "--tolerance-days",
        type=float,
        default=31.0,
        help="Date-match window around the Goffin predicted epoch (days). Default: 31.",
    )
    parser.add_argument(
        "--distance-tolerance-au",
        type=float,
        default=1e-4,
        help="Tolerance used to report distance agreement against Goffin. Default: 1e-4 AU.",
    )
    parser.add_argument(
        "--require-distance-tolerance",
        action="store_true",
        help="Return non-zero if any expected event is missing or outside --distance-tolerance-au.",
    )
    args = parser.parse_args()

    cfg = load_config(args.config)
    dist_thresh = (
        args.dist_threshold_au if args.dist_threshold_au is not None else cfg.detection.threshold_au
    )

    goffin_path = Path(cfg.paths.raw) / cfg.sources.goffin_2014.output_filename
    catalog_path = (
        args.catalog
        if args.catalog is not None
        else Path(cfg.paths.output) / f"{cfg.output.filename}.{cfg.output.format}"
    )

    if not goffin_path.exists():
        logger.error("Goffin catalog not found: %s — run download_goffin_2014.", goffin_path)
        return 1
    if not catalog_path.exists():
        logger.error("Pipeline catalog not found: %s — run run_pipeline.", catalog_path)
        return 1

    try:
        goffin, col_map = _load_goffin(
            goffin_path, cfg.time_window.start[:10], cfg.time_window.end[:10]
        )
    except RuntimeError as exc:
        logger.error("Failed to parse Goffin catalog: %s", exc)
        return 1

    logger.info(
        "Goffin 2014 events in Gaia window (%s → %s): %d",
        cfg.time_window.start[:10],
        cfg.time_window.end[:10],
        len(goffin),
    )

    pert_col = col_map["perturber"]
    targ_col = col_map["target"]
    dist_col = col_map["dist"]

    expected = goffin.filter(pl.col(dist_col) <= dist_thresh).sort(dist_col)
    wider = goffin.filter(pl.col(dist_col) > dist_thresh)
    logger.info(
        "Partition by dist ≤ %.4f AU: %d expected in our catalog, %d wider (not expected)",
        dist_thresh,
        len(expected),
        len(wider),
    )
    catalog = _load_relevant_catalog(
        catalog_path, expected, pert_col, targ_col, args.tolerance_days
    )

    matched: list[dict] = []
    missed: list[dict] = []

    logger.info("")
    logger.info("=== Cross-match: expected events (dist ≤ %.4f AU) ===", dist_thresh)
    for row in expected.iter_rows(named=True):
        pert = int(row[pert_col])
        targ = int(row[targ_col])
        epoch_date = row["_epoch_date"].isoformat()
        epoch_jd = _date_to_jd_tdb(epoch_date)
        goffin_dist = float(row[dist_col])
        vel_str = f"  Vel={row[col_map['vel']]:.2f} km/s" if "vel" in col_map else ""

        hits = _find_match(catalog, pert, targ, epoch_jd, args.tolerance_days)

        if len(hits) == 0:
            logger.warning(
                "  ✗ MISS  (%d, %d)  Goffin=%s  dist=%.5f AU%s",
                pert,
                targ,
                epoch_date,
                goffin_dist,
                vel_str,
            )
            missed.append(
                {
                    "perturber": pert,
                    "target": targ,
                    "goffin_date": epoch_date,
                    "goffin_dist_au": goffin_dist,
                }
            )
        else:
            best = hits.sort("dist_au").head(1).row(0, named=True)
            our_date = Time(best["jd_tdb"], format="jd", scale="tdb").utc.iso[:10]
            delta_days = best["jd_tdb"] - epoch_jd
            delta_dist_au = float(best["dist_au"]) - goffin_dist
            abs_delta_dist_au = abs(delta_dist_au)
            within_distance_tolerance = abs_delta_dist_au <= args.distance_tolerance_au
            logger.info(
                "  ✓ HIT   (%d, %d)  Goffin=%s  Ours=%s (Δ=%+.1fd)  "
                "dist_G=%.5f AU  dist_ours=%.5f AU  Δdist=%+.6f AU",
                pert,
                targ,
                epoch_date,
                our_date,
                delta_days,
                goffin_dist,
                best["dist_au"],
                delta_dist_au,
            )
            matched.append(
                {
                    "perturber": pert,
                    "target": targ,
                    "goffin_date": epoch_date,
                    "goffin_dist_au": goffin_dist,
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
    logger.info("Expected events (dist ≤ %.4f AU): %d", dist_thresh, n_exp)
    logger.info("Matched in our catalog:           %d (%.1f%%)", n_hit, rate)
    logger.info("Missed:                           %d", len(missed))
    logger.info("Wider encounters reported only:   %d", len(wider))

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

    report_dir = Path(cfg.paths.output)
    if matched:
        pl.DataFrame(matched).write_csv(report_dir / "goffin_2014_matches.csv")
    if missed:
        pl.DataFrame(missed).write_csv(report_dir / "goffin_2014_misses.csv")
    logger.info("Match/miss CSVs written to %s/", report_dir)

    if args.require_distance_tolerance and (missed or distance_failures):
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
