"""Download Galád & Gray (2002) asteroid-asteroid close-encounter catalog.

Reference:
    Galád A., Gray B. (2002)
    "Asteroid encounters suitable for mass determinations"
    A&A 391, 1115 — DOI: 10.1051/0004-6361:20020880

The catalog is **not** available on VizieR — the encounter tables are embedded
as HTML inside the published article:

    https://www.aanda.org/articles/aa/full/2002/33/aah3376/aah3376.right.html

The paper contains seven encounter tables:

* Table 1 — Asteroids heavily perturbed by (1) Ceres
* Table 2 — Asteroids heavily perturbed by (4) Vesta
* Table 3 — Asteroids perturbed by (2) Pallas
* Table 4 — Asteroids perturbed by (10) Hygiea, before 2000-01-01
* Table 5 — Asteroids perturbed by (10) Hygiea, after  2000-01-01
* Table 6 — Other large perturbers, before 1980-01-01
* Table 7 — Other large perturbers, 1980-1997

Tables 1-5 fix the perturber per-table; tables 6-7 list both perturber and
target on each row.  The columns are otherwise consistent:

    number  name  date(year/m/d)  r[AU]  v[km/s]  P[km/s]  obs(year)

The script fetches the HTML, parses every table with ``pandas.read_html`` (the
only realistic option since the rows are in <TR>/<TD> markup with no class
hooks), normalises into a single typed polars DataFrame, and writes:

    data/raw/galad_2002_encounters.parquet
    data/raw/galad_2002_metadata.json    (source URL, row counts, fetch time)

Usage:
    docker compose run --rm pipeline python -m scripts.download_galad_2002
    docker compose run --rm pipeline python -m scripts.download_galad_2002 --config config.yaml
"""

from __future__ import annotations

import argparse
import datetime as dt
import io
import json
import logging
import re
import sys
import urllib.request
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd
import polars as pl

from src.utils.config import load_config

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)


# Fixed perturbers for the per-target tables.  Indices match the order in
# which pandas.read_html returns the tables on the source page.
@dataclass(frozen=True)
class _FixedPerturber:
    number: int
    name: str
    label: str  # short string used in the ``source_table`` column


_FIXED_PERTURBERS: dict[int, _FixedPerturber] = {
    0: _FixedPerturber(1, "Ceres", "table1_ceres"),
    1: _FixedPerturber(4, "Vesta", "table2_vesta"),
    2: _FixedPerturber(2, "Pallas", "table3_pallas"),
    3: _FixedPerturber(10, "Hygiea", "table4_hygiea_pre2000"),
    4: _FixedPerturber(10, "Hygiea", "table5_hygiea_post2000"),
}

_MIXED_TABLES: dict[int, str] = {
    5: "table6_other_pre1980",
    6: "table7_other_1980_1997",
}


_NUMBER_RE = re.compile(r"\(?\s*(\d[\d\s ]*)\s*\)?")


def _clean_number(raw: object) -> int | None:
    """Strip parentheses and embedded whitespace, return int or None.

    Galád formats asteroid numbers as ``(14 375)`` with a non-breaking
    space; this helper handles both regular and U+00A0 spaces.
    """
    if raw is None:
        return None
    s = str(raw).replace(" ", " ").strip()
    if not s or s.lower() == "nan":
        return None
    m = _NUMBER_RE.search(s)
    if not m:
        return None
    digits = re.sub(r"\s+", "", m.group(1))
    try:
        return int(digits)
    except ValueError:
        return None


def _clean_str(raw: object) -> str | None:
    """Trim a string field, returning None for empty / NaN values."""
    if raw is None:
        return None
    s = str(raw).replace(" ", " ").strip()
    if not s or s.lower() == "nan":
        return None
    return s


def _clean_float(raw: object) -> float | None:
    """Coerce a numeric cell to float, returning None on failure."""
    if raw is None:
        return None
    s = str(raw).strip()
    if not s or s.lower() == "nan":
        return None
    try:
        return float(s)
    except ValueError:
        return None


def _parse_date(raw: object) -> tuple[dt.date | None, str | None]:
    """Parse a Galád ``year/m/d`` string.

    Returns
    -------
    tuple of (date, precision)
        ``precision`` is ``"day"`` when full y/m/d is given, ``"year"`` when
        only the year could be parsed, ``None`` if the cell is empty.
        For ``"year"`` precision a date of July 1st of that year is returned
        as a reasonable midpoint.
    """
    if raw is None:
        return None, None
    s = str(raw).replace(" ", " ").strip()
    if not s or s.lower() == "nan":
        return None, None

    parts = [p.strip() for p in s.split("/") if p.strip()]
    if len(parts) == 3:
        try:
            y, m, d = (int(p) for p in parts)
            return dt.date(y, m, d), "day"
        except ValueError:
            pass
    if len(parts) >= 1:
        try:
            y = int(parts[0])
            return dt.date(y, 7, 1), "year"
        except ValueError:
            pass
    return None, None


def _normalise_fixed(df: pd.DataFrame, fixed: _FixedPerturber) -> list[dict]:
    """Convert a table 1-5 DataFrame into the unified row schema."""
    # Two header rows precede the data; the actual data starts at row 2.
    rows: list[dict] = []
    for _, r in df.iloc[2:].iterrows():
        target_number = _clean_number(r.iloc[0])
        target_name = _clean_str(r.iloc[1])
        date_obj, precision = _parse_date(r.iloc[2])
        rows.append(
            {
                "perturber_number": fixed.number,
                "perturber_name": fixed.name,
                "perturber_diameter_km": None,
                "target_number": target_number,
                "target_name": target_name,
                "date_raw": _clean_str(r.iloc[2]),
                "date_parsed": date_obj,
                "date_precision": precision,
                "r_au": _clean_float(r.iloc[3]),
                "v_km_s": _clean_float(r.iloc[4]),
                "p_km_s": _clean_float(r.iloc[5]),
                "obs_year": _clean_number(r.iloc[6]),
                "source_table": fixed.label,
            }
        )
    return rows


def _normalise_mixed(df: pd.DataFrame, label: str) -> list[dict]:
    """Convert a table 6/7 DataFrame into the unified row schema."""
    rows: list[dict] = []
    for _, r in df.iloc[2:].iterrows():
        perturber_number = _clean_number(r.iloc[0])
        perturber_name = _clean_str(r.iloc[1])
        perturber_d = _clean_float(r.iloc[2])
        target_number = _clean_number(r.iloc[3])
        target_name = _clean_str(r.iloc[4])
        date_obj, precision = _parse_date(r.iloc[5])
        rows.append(
            {
                "perturber_number": perturber_number,
                "perturber_name": perturber_name,
                "perturber_diameter_km": perturber_d,
                "target_number": target_number,
                "target_name": target_name,
                "date_raw": _clean_str(r.iloc[5]),
                "date_parsed": date_obj,
                "date_precision": precision,
                "r_au": _clean_float(r.iloc[6]),
                "v_km_s": _clean_float(r.iloc[7]),
                "p_km_s": _clean_float(r.iloc[8]),
                "obs_year": _clean_number(r.iloc[9]),
                "source_table": label,
            }
        )
    return rows


def _fetch_html(url: str, timeout: float = 60.0) -> str:
    """Fetch *url* with a browser-like User-Agent and return decoded HTML."""
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read().decode("utf-8", errors="replace")


def download_galad_2002(config_path: str = "config.yaml") -> Path:
    """Fetch and parse the Galád & Gray (2002) encounter tables.

    Parameters
    ----------
    config_path:
        Path to the YAML configuration file.

    Returns
    -------
    Path
        Path to the unified parquet output.
    """
    cfg = load_config(config_path)
    raw_dir = Path(cfg.paths.raw)
    raw_dir.mkdir(parents=True, exist_ok=True)

    url = cfg.sources.galad_2002.source_url
    out_filename = cfg.sources.galad_2002.output_filename
    dest = raw_dir / out_filename
    meta_path = raw_dir / "galad_2002_metadata.json"

    logger.info("Fetching Galád & Gray (2002) HTML from %s …", url)
    html = _fetch_html(url)
    logger.info("Received %d bytes", len(html))

    tables = pd.read_html(io.StringIO(html), flavor="bs4")
    logger.info("Parsed %d <TABLE> blocks from the page", len(tables))
    if len(tables) < 7:
        raise RuntimeError(
            f"Expected ≥7 tables in Galád 2002 page, got {len(tables)}; "
            "the page structure may have changed."
        )

    all_rows: list[dict] = []
    per_table_counts: dict[str, int] = {}

    for idx, fixed in _FIXED_PERTURBERS.items():
        rows = _normalise_fixed(tables[idx], fixed)
        all_rows.extend(rows)
        per_table_counts[fixed.label] = len(rows)
        logger.info("  %-32s %d rows", fixed.label, len(rows))

    for idx, label in _MIXED_TABLES.items():
        rows = _normalise_mixed(tables[idx], label)
        all_rows.extend(rows)
        per_table_counts[label] = len(rows)
        logger.info("  %-32s %d rows", label, len(rows))

    schema = {
        "perturber_number": pl.Int64,
        "perturber_name": pl.Utf8,
        "perturber_diameter_km": pl.Float64,
        "target_number": pl.Int64,
        "target_name": pl.Utf8,
        "date_raw": pl.Utf8,
        "date_parsed": pl.Date,
        "date_precision": pl.Utf8,
        "r_au": pl.Float64,
        "v_km_s": pl.Float64,
        "p_km_s": pl.Float64,
        "obs_year": pl.Int64,
        "source_table": pl.Utf8,
    }
    df = pl.DataFrame(all_rows, schema=schema)
    df.write_parquet(dest, compression="zstd")
    logger.info("Wrote %d rows to %s", len(df), dest)

    # Sanity log: distribution by precision and Gaia-window count.
    win_start = dt.date.fromisoformat(cfg.time_window.start[:10])
    win_end = dt.date.fromisoformat(cfg.time_window.end[:10])
    in_win = df.filter(
        pl.col("date_parsed").is_not_null()
        & (pl.col("date_parsed") >= win_start)
        & (pl.col("date_parsed") <= win_end)
    )
    logger.info(
        "Date precision: %d day-precision, %d year-only, %d unparsed",
        df.filter(pl.col("date_precision") == "day").height,
        df.filter(pl.col("date_precision") == "year").height,
        df.filter(pl.col("date_precision").is_null()).height,
    )
    logger.info(
        "Encounters falling inside Gaia DR3 window (%s → %s): %d",
        win_start,
        win_end,
        len(in_win),
    )

    metadata = {
        "source_url": url,
        "downloaded_at": datetime.now(UTC).isoformat(),
        "filename": out_filename,
        "n_rows_total": len(df),
        "n_rows_per_table": per_table_counts,
        "size_bytes": dest.stat().st_size,
        "reference": (
            "Galád A., Gray B. (2002) "
            "A&A 391, 1115 — DOI: 10.1051/0004-6361:20020880"
        ),
    }
    meta_path.write_text(json.dumps(metadata, indent=2))
    logger.info("Metadata → %s", meta_path)
    return dest


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Download Galád & Gray (2002) encounter catalog from the published HTML"
    )
    parser.add_argument("--config", default="config.yaml", help="Path to config.yaml")
    args = parser.parse_args()
    download_galad_2002(args.config)
    return 0


if __name__ == "__main__":
    sys.exit(main())
