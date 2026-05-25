"""Download MPCORB.DAT from the Minor Planet Center.

Usage:
    docker compose run --rm pipeline python -m scripts.download_mpcorb
    docker compose run --rm pipeline python -m scripts.download_mpcorb --config config.yaml
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import logging
import urllib.request
from datetime import UTC, datetime
from pathlib import Path

from src.utils.config import load_config

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)


def _md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def download_mpcorb(config_path: str = "config.yaml") -> Path:
    """Download and decompress MPCORB.DAT, writing metadata alongside.

    Returns
    -------
    Path
        Path to the decompressed MPCORB.DAT file.
    """
    cfg = load_config(config_path)
    raw_dir = Path(cfg.paths.raw)
    raw_dir.mkdir(parents=True, exist_ok=True)

    url = cfg.sources.mpcorb.url
    dest = raw_dir / cfg.sources.mpcorb.local_filename
    meta_path = raw_dir / "mpcorb_metadata.json"

    logger.info("Downloading MPCORB from %s", url)
    with urllib.request.urlopen(url) as resp:  # noqa: S310
        compressed = resp.read()

    logger.info("Decompressing (%d bytes compressed)…", len(compressed))
    data = gzip.decompress(compressed)
    dest.write_bytes(data)
    logger.info("Saved %d bytes to %s", len(data), dest)

    metadata = {
        "url": url,
        "downloaded_at": datetime.now(UTC).isoformat(),
        "filename": cfg.sources.mpcorb.local_filename,
        "md5": _md5(dest),
        "size_bytes": dest.stat().st_size,
    }
    meta_path.write_text(json.dumps(metadata, indent=2))
    logger.info("Metadata → %s", meta_path)
    return dest


def main() -> int:
    parser = argparse.ArgumentParser(description="Download MPCORB.DAT from MPC")
    parser.add_argument("--config", default="config.yaml", help="Path to config.yaml")
    args = parser.parse_args()
    download_mpcorb(args.config)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
