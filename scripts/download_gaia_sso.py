"""Download Gaia DR3 SSO observations from the Gaia Archive via TAP.

Saves to data/raw/gaia_sso.parquet.
The ``epoch`` column is in TCB (documented in the file and in the module
docstring of src/ingest/gaia_sso.py).

Usage:
    docker compose run --rm pipeline python -m scripts.download_gaia_sso
    docker compose run --rm pipeline python -m scripts.download_gaia_sso --config config.yaml
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from src.ingest.gaia_sso import download_gaia_sso
from src.utils.config import load_config

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)


def main() -> None:
    parser = argparse.ArgumentParser(description="Download Gaia DR3 SSO observations")
    parser.add_argument("--config", default="config.yaml", help="Path to config.yaml")
    args = parser.parse_args()

    cfg = load_config(args.config)
    dest = Path(cfg.paths.raw) / "gaia_sso.parquet"

    download_gaia_sso(
        archive_url=cfg.sources.gaia_sso.archive_url,
        columns=cfg.sources.gaia_sso.columns,
        dest=dest,
    )


if __name__ == "__main__":
    sys.exit(main())
