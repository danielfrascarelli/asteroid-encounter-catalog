"""Read-only reconnaissance of the Gaia FPR SSO data model (FPR_INGEST_PLAN Stage 0).

Confirms — against the live Gaia archive — everything the FPR ingest must not
assume:

  1. Which ``gaiafpr`` (Focused Product Release) tables exist and which one holds
     per-transit asteroid astrometry (the FPR analogue of
     ``gaiadr3.sso_observation``).
  2. Its column set, diffed against the 14 columns the DR3 ingest uses today.
  3. The empirical observation window (MIN/MAX ``epoch``), row count, and number
     of distinct ``number_mp`` — to fix the temporal window and ``mp_max``.
  4. A side-by-side of (1) Ceres in DR3 vs FPR: epoch count and temporal span.
     FPR must show more epochs and a longer baseline, or our premise is wrong.

Nothing is written to the catalog. Output is a JSON report (machine-readable)
plus a human summary on stdout; the JSON feeds ``docs/gaia_fpr_data_model.md``.

Usage:
    docker run --rm -v "$PWD:/app" gaia-asteroid-encounters:latest \\
        python -m scripts.dev.probe_gaia_fpr
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

from astroquery.utils.tap.core import TapPlus

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

_ARCHIVE_URL = "https://gea.esac.esa.int/tap-server/tap"
_DR3_TABLE = "gaiadr3.sso_observation"

# The 14 columns the DR3 ingest + mass layer consume today (src/ingest/gaia_sso.py,
# scripts/mass/fit_mass_gaia_loo.py). Any FPR rename/absence shows up in the diff.
_DR3_COLUMNS_USED = {
    "solution_id",
    "source_id",
    "denomination",
    "number_mp",
    "epoch",
    "ra",
    "dec",
    "g_mag",
    "x_gaia",
    "y_gaia",
    "z_gaia",
    "ra_error_systematic",
    "dec_error_systematic",
    "ra_dec_correlation_systematic",
    "ra_error_random",
    "dec_error_random",
    "ra_dec_correlation_random",
    "position_angle_scan",
}

# DR3 epoch reference: epoch = JD_TCB - 2455197.5 (days since J2010.0 TCB).
_J2010_TCB_JD = 2455197.5


def _tap(archive_url: str) -> TapPlus:
    return TapPlus(url=archive_url)


def list_fpr_tables(tap: TapPlus) -> list[str]:
    """All table names whose schema is ``gaiafpr`` (or that mention fpr+sso)."""
    tables = tap.load_tables(only_names=True)
    names = [t.name for t in tables]
    fpr = sorted(n for n in names if n.lower().startswith("gaiafpr"))
    return fpr


def guess_sso_table(fpr_tables: list[str]) -> str | None:
    """Best guess for the per-transit SSO observation table within FPR."""
    cands = [t for t in fpr_tables if "sso" in t.lower() and "observation" in t.lower()]
    if cands:
        return cands[0]
    cands = [t for t in fpr_tables if "sso" in t.lower()]
    return cands[0] if cands else None


def describe_columns(tap: TapPlus, table: str) -> list[str]:
    """Column names for *table* via TAP metadata."""
    meta = tap.load_table(table)
    return [c.name for c in meta.columns]


def window_stats(tap: TapPlus, table: str) -> dict:
    adql = (
        "SELECT MIN(epoch) AS emin, MAX(epoch) AS emax, "
        "COUNT(*) AS nrows, COUNT(DISTINCT number_mp) AS nmp, "
        "MAX(number_mp) AS mpmax "
        f"FROM {table}"
    )
    job = tap.launch_job_async(adql)
    row = job.get_results()[0]
    emin = float(row["emin"])
    emax = float(row["emax"])
    return {
        "epoch_min_days": emin,
        "epoch_max_days": emax,
        "jd_tcb_min": emin + _J2010_TCB_JD,
        "jd_tcb_max": emax + _J2010_TCB_JD,
        "baseline_days": emax - emin,
        "baseline_months": (emax - emin) / 30.44,
        "n_rows": int(row["nrows"]),
        "n_distinct_number_mp": int(row["nmp"]),
        "max_number_mp": int(row["mpmax"]),
    }


def ceres_compare(tap: TapPlus, fpr_table: str) -> dict:
    out = {}
    for label, table in (("dr3", _DR3_TABLE), ("fpr", fpr_table)):
        adql = (
            "SELECT COUNT(*) AS n, MIN(epoch) AS emin, MAX(epoch) AS emax "
            f"FROM {table} WHERE number_mp = 1"
        )
        row = tap.launch_job_async(adql).get_results()[0]
        emin, emax = float(row["emin"]), float(row["emax"])
        out[label] = {
            "n_epochs": int(row["n"]),
            "epoch_min_days": emin,
            "epoch_max_days": emax,
            "baseline_days": emax - emin,
        }
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--archive-url", default=_ARCHIVE_URL)
    ap.add_argument("--table", default=None, help="Force FPR SSO table name (skip guessing)")
    ap.add_argument("--out", default="data/output/gaia_fpr_probe.json")
    args = ap.parse_args()

    tap = _tap(args.archive_url)
    report: dict = {"archive_url": args.archive_url, "dr3_table": _DR3_TABLE}

    logger.info("Listing gaiafpr tables...")
    fpr_tables = list_fpr_tables(tap)
    report["fpr_tables"] = fpr_tables
    logger.info("Found %d gaiafpr tables: %s", len(fpr_tables), fpr_tables)

    sso_table = args.table or guess_sso_table(fpr_tables)
    if sso_table is None:
        logger.error("Could not identify an FPR SSO table. Inspect fpr_tables above.")
        report["error"] = "no_sso_table_found"
        _write(report, args.out)
        return
    report["fpr_sso_table"] = sso_table
    logger.info("Using FPR SSO table: %s", sso_table)

    logger.info("Describing columns...")
    fpr_cols = describe_columns(tap, sso_table)
    report["fpr_columns"] = sorted(fpr_cols)
    fpr_col_set = set(fpr_cols)
    report["columns_missing_vs_dr3"] = sorted(_DR3_COLUMNS_USED - fpr_col_set)
    report["columns_present_vs_dr3"] = sorted(_DR3_COLUMNS_USED & fpr_col_set)

    logger.info("Querying window stats (may take a minute)...")
    report["window"] = window_stats(tap, sso_table)

    logger.info("Comparing (1) Ceres DR3 vs FPR...")
    report["ceres"] = ceres_compare(tap, sso_table)

    _write(report, args.out)
    _print_summary(report)


def _write(report: dict, out: str) -> None:
    p = Path(out)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(report, indent=2))
    logger.info("Wrote probe report to %s", p)


def _print_summary(r: dict) -> None:
    print("\n" + "=" * 70)
    print("GAIA FPR SSO DATA MODEL — PROBE SUMMARY")
    print("=" * 70)
    print(f"FPR SSO table       : {r.get('fpr_sso_table')}")
    miss = r.get("columns_missing_vs_dr3", [])
    print(f"Columns missing/renamed vs DR3 : {miss if miss else 'NONE — drop-in'}")
    w = r.get("window", {})
    if w:
        print(
            f"Window              : JD_TCB {w['jd_tcb_min']:.1f} … {w['jd_tcb_max']:.1f} "
            f"(~{w['baseline_months']:.1f} months, {w['n_rows']:,} rows, "
            f"{w['n_distinct_number_mp']:,} asteroids, mp_max={w['max_number_mp']})"
        )
    c = r.get("ceres", {})
    if c:
        d, f = c.get("dr3", {}), c.get("fpr", {})
        print(
            f"(1) Ceres epochs    : DR3={d.get('n_epochs')} (baseline {d.get('baseline_days', 0):.0f} d) "
            f"-> FPR={f.get('n_epochs')} (baseline {f.get('baseline_days', 0):.0f} d)"
        )
    print("=" * 70 + "\n")


if __name__ == "__main__":
    main()
