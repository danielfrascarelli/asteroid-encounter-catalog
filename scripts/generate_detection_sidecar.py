"""Backfill a detection-provenance sidecar for an existing catalog parquet.

Use this when a catalog parquet was generated before
:func:`src.catalog.writer.write_detection_sidecar` was wired into the pipeline.
The sidecar is reconstructed from the *current* config + git state, so the
output reflects what would be produced by re-running with the same config —
not necessarily what was actually used for the historical file.  Verify the
config dates match the parquet timestamp before trusting the sidecar.

Usage:
    docker compose run --rm pipeline python -m scripts.generate_detection_sidecar \\
        --catalog data/output/encounters_catalog_rebound_005au.parquet
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import polars as pl
from astropy.time import Time

from src.catalog.writer import write_detection_sidecar
from src.ingest.mpcorb_archive import discover_snapshots, select_for_window
from src.utils.config import load_config

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--catalog", required=True, help="Parquet path to backfill provenance for")
    p.add_argument("--config", default="config.yaml")
    args = p.parse_args()

    catalog_path = Path(args.catalog)
    if not catalog_path.exists():
        logger.error("Catalog not found: %s", catalog_path)
        return 1

    cfg = load_config(args.config)

    df = pl.scan_parquet(catalog_path)
    n = df.select(pl.len()).collect().item()
    logger.info("Catalog %s has %d rows", catalog_path.name, n)

    # Resolve MPCORB snapshot from config (mirror logic in run_pipeline)
    tw = cfg.time_window
    t_start_jd = Time(tw.start, scale=tw.scale).tdb.jd
    t_end_jd = Time(tw.end, scale=tw.scale).tdb.jd
    snapshots = discover_snapshots(cfg.paths.raw)
    snap = select_for_window(snapshots, t_start_jd, t_end_jd)
    if snap is None:
        logger.warning("No MPCORB snapshot resolvable; sidecar will omit MPCORB hash")
        mpcorb_path = None
    else:
        mpcorb_path = Path(snap.path)
        logger.info("Resolved MPCORB snapshot: %s", mpcorb_path.name)

    fine_step_hours = cfg.propagation.time_step_hours
    coarse_step_hours = cfg.propagation.coarse_step_hours or fine_step_hours
    use_tiered = coarse_step_hours > fine_step_hours

    # run_id encodes catalog mtime so re-runs of this script keep one run_id
    # per catalog version.
    mtime = catalog_path.stat().st_mtime
    run_id = "backfill_" + Time(mtime, format="unix").utc.iso.replace(" ", "T")[:19] + "Z"

    sidecar = write_detection_sidecar(
        catalog_path,
        run_id=run_id,
        n_encounters=n,
        cfg=cfg,
        mpcorb_path=mpcorb_path,
        coarse_step_hours=coarse_step_hours,
        fine_step_hours=fine_step_hours,
        use_tiered=use_tiered,
        force_kepler_refine=use_tiered,
    )
    logger.info(
        "Sidecar written: %s  (backfilled from current config — verify timestamps match)",
        sidecar,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
