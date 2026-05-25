"""Download Gaia DR3 SSO observations from the Gaia Archive via TAP.

Saves to data/raw/gaia_sso.parquet.
The ``epoch`` column is **days since J2010.0 TCB** (``JD_TCB − 2455197.5``),
NOT a raw Julian Date. See the module docstring of ``src/ingest/gaia_sso.py``
for the full convention and conversion recipes.

Usage:
    docker compose run --rm pipeline python -m scripts.download_gaia_sso
    docker compose run --rm pipeline python -m scripts.download_gaia_sso --config config.yaml
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from src.ingest.gaia_sso import download_gaia_sso
from src.utils.config import load_config

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)


def main() -> int:
    parser = argparse.ArgumentParser(description="Download Gaia DR3 SSO observations")
    parser.add_argument("--config", default="config.yaml", help="Path to config.yaml")
    args = parser.parse_args()

    cfg = load_config(args.config)
    dest = Path(cfg.paths.raw) / "gaia_sso.parquet"

    cache_dir = Path(cfg.paths.cache) / "gaia_sso_chunks"
    download_gaia_sso(
        archive_url=cfg.sources.gaia_sso.archive_url,
        columns=cfg.sources.gaia_sso.columns,
        dest=dest,
        batch_size=cfg.sources.gaia_sso.batch_size,
        n_workers=cfg.sources.gaia_sso.n_workers,
        cache_dir=cache_dir,
        max_retries=cfg.sources.gaia_sso.max_retries,
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
