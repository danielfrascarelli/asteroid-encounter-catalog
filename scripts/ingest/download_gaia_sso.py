"""Download Gaia SSO observations from the Gaia Archive via TAP.

The Gaia release (DR3 or FPR) is chosen by ``sources.gaia_sso.release`` in the
config, overridable with ``--release``. Output and cache are scoped per release
so DR3 and FPR artefacts never mix:

    data/raw/gaia_sso_{release}.parquet
    data/cache/gaia_sso_chunks/{release}/

The ``epoch`` column is **days since J2010.0 TCB** (``JD_TCB − 2455197.5``),
NOT a raw Julian Date (same reference in both releases — see
docs/gaia_fpr_data_model.md). See the module docstring of
``src/ingest/gaia_sso.py`` for the full convention and conversion recipes.

Usage:
    docker compose run --rm pipeline python -m scripts.ingest.download_gaia_sso
    docker compose run --rm pipeline python -m scripts.ingest.download_gaia_sso --release fpr
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
    parser = argparse.ArgumentParser(description="Download Gaia SSO observations (DR3 or FPR)")
    parser.add_argument("--config", default="config.yaml", help="Path to config.yaml")
    parser.add_argument(
        "--release",
        default=None,
        help="Gaia release to download ('dr3' | 'fpr'). Defaults to config's gaia_sso.release.",
    )
    args = parser.parse_args()

    cfg = load_config(args.config)
    gaia = cfg.sources.gaia_sso
    if args.release is not None:
        gaia.release = args.release  # override the active release for this run
    rel = gaia.release
    release_cfg = gaia.active()

    dest = Path(cfg.paths.raw) / f"gaia_sso_{rel}.parquet"
    cache_dir = Path(cfg.paths.cache) / "gaia_sso_chunks" / rel

    logger.info(
        "Downloading Gaia release '%s' from table %s (mp_max=%d, window %s … %s)",
        rel,
        release_cfg.table,
        release_cfg.mp_max,
        release_cfg.window_start,
        release_cfg.window_end,
    )
    download_gaia_sso(
        archive_url=gaia.archive_url,
        columns=gaia.active_columns(),
        dest=dest,
        table=release_cfg.table,
        mp_max=release_cfg.mp_max,
        batch_size=gaia.batch_size,
        n_workers=gaia.n_workers,
        cache_dir=cache_dir,
        max_retries=gaia.max_retries,
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
