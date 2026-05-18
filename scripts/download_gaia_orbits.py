"""Download Gaia DR3 SSO orbital elements from the Gaia Archive via TAP.

Saves to data/raw/gaia_orbits.parquet in pipeline schema (a_au, e, i_deg, …).
The osculating epoch (epoch_jd) is in TDB and falls within the Gaia DR3
observation window, eliminating the ~10-year backward extrapolation that
MPCORB requires.

Usage:
    docker compose run --rm pipeline python -m scripts.download_gaia_orbits
    docker compose run --rm pipeline python -m scripts.download_gaia_orbits --config config.yaml
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from src.ingest.gaia_orbits import download_gaia_orbits
from src.utils.config import load_config

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)


def main() -> int:
    parser = argparse.ArgumentParser(description="Download Gaia DR3 SSO orbital elements")
    parser.add_argument("--config", default="config.yaml", help="Path to config.yaml")
    args = parser.parse_args()

    cfg = load_config(args.config)
    dest = Path(cfg.paths.raw) / "gaia_orbits.parquet"
    cache_dir = Path(cfg.paths.cache) / "gaia_orbits_chunks"

    download_gaia_orbits(
        archive_url=cfg.sources.gaia_orbits.archive_url,
        dest=dest,
        batch_size=cfg.sources.gaia_orbits.batch_size,
        n_workers=cfg.sources.gaia_orbits.n_workers,
        cache_dir=cache_dir,
        max_retries=cfg.sources.gaia_orbits.max_retries,
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
