"""Download Goffin (2014) asteroid-asteroid close-encounter catalog from VizieR.

Reference:
    Goffin E. (2014)
    "New determination of asteroid masses from close encounters"
    A&A 565, A56 — DOI: 10.1051/0004-6361/201322766
    VizieR catalog: J/A+A/565/A56

The catalog contains the close encounters used to determine masses of 230
main-belt asteroids, spanning roughly 1900–2012.  Each row gives the
perturber number, test-asteroid number, encounter epoch, minimum distance,
and relative velocity.

Used as ground-truth validation for the encounter detection pipeline:
encounters in the Gaia DR3 window (2014-07-25 → 2017-05-28) that also
fall inside our detection threshold should appear in our generated catalog.

Usage:
    docker compose run --rm pipeline python -m scripts.download_goffin_2014
    docker compose run --rm pipeline python -m scripts.download_goffin_2014 --config config.yaml
"""

from __future__ import annotations

import argparse
import json
import logging
from datetime import UTC, datetime
from pathlib import Path

import polars as pl

from src.utils.config import load_config

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)


def _astropy_to_polars(table) -> pl.DataFrame:
    """Convert an astropy.table.Table to a polars DataFrame."""
    import numpy as np

    cols: dict[str, list] = {}
    for name in table.colnames:
        col = table[name]
        values: list = []
        for v in col:
            if hasattr(v, "mask") and bool(getattr(v, "mask", False)):
                values.append(None)
            elif isinstance(v, bytes):
                values.append(v.decode("ascii", errors="replace").strip())
            elif isinstance(v, np.generic):
                values.append(v.item())
            else:
                values.append(v)
        cols[name] = values
    return pl.DataFrame(cols)


def download_goffin_2014(config_path: str = "config.yaml") -> Path:
    """Fetch J/A+A/565/A56 from VizieR and write a single parquet.

    Returns
    -------
    Path
        Path to the combined parquet file.
    """
    from astroquery.vizier import Vizier

    cfg = load_config(config_path)
    raw_dir = Path(cfg.paths.raw)
    raw_dir.mkdir(parents=True, exist_ok=True)

    catalog_id = cfg.sources.goffin_2014.vizier_catalog
    out_filename = cfg.sources.goffin_2014.output_filename
    dest = raw_dir / out_filename
    meta_path = raw_dir / "goffin_2014_metadata.json"

    logger.info("Querying VizieR for catalog %s …", catalog_id)
    v = Vizier(row_limit=-1)
    table_list = v.get_catalogs(catalog_id)

    if len(table_list) == 0:
        raise RuntimeError(f"VizieR returned no tables for {catalog_id}")

    logger.info("Retrieved %d tables", len(table_list))

    parts: list[pl.DataFrame] = []
    per_table_counts: dict[str, int] = {}
    for tab in table_list:
        name = getattr(tab, "meta", {}).get("name", "unknown")
        df = _astropy_to_polars(tab)
        df = df.with_columns(pl.lit(name).alias("source_table"))
        parts.append(df)
        per_table_counts[name] = len(df)
        logger.info("  %-32s  %d rows  cols=%s", name, len(df), df.columns)

    combined = pl.concat(parts, how="diagonal")
    combined.write_parquet(dest, compression="zstd")
    logger.info("Wrote %d rows to %s", len(combined), dest)

    metadata = {
        "vizier_catalog": catalog_id,
        "downloaded_at": datetime.now(UTC).isoformat(),
        "filename": out_filename,
        "n_rows_total": len(combined),
        "n_rows_per_table": per_table_counts,
        "size_bytes": dest.stat().st_size,
        "reference": ("Goffin E. (2014) " "A&A 565, A56 — DOI: 10.1051/0004-6361/201322766"),
    }
    meta_path.write_text(json.dumps(metadata, indent=2))
    logger.info("Metadata → %s", meta_path)
    return dest


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Download Goffin (2014) encounter catalog from VizieR"
    )
    parser.add_argument("--config", default="config.yaml", help="Path to config.yaml")
    args = parser.parse_args()
    download_goffin_2014(args.config)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
