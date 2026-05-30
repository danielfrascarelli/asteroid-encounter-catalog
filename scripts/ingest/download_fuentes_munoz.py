"""Download the Fuentes-Muñoz et al. (2025) machine-readable Table 5.

Fuentes-Muñoz, Farnocchia, Giorgini & Park (2025), "Asteroid Mass Estimation by
Mutual Perturbations during Close Encounters after Gaia FPR", AJ 170, 353
(DOI 10.3847/1538-3881/ae0cc9). Table 5 lists each perturber asteroid with the
pipe-delimited list of test asteroids that showed a mass signal — the
perturber→target encounter pairs used by
``scripts.validate.validate_fuentes_munoz_2025``.

Usage:
    docker compose run --rm pipeline python -m scripts.ingest.download_fuentes_munoz
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import urllib.request
from datetime import UTC, datetime
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
logger = logging.getLogger(__name__)

_URL = "https://content.cld.iop.org/journals/1538-3881/170/6/353/revision2/ajae0cc9t5_mrt.txt"
_OUT_DIR = Path("data/raw/fuentes_munoz_2025")
_FILENAME = "ajae0cc9t5_mrt.txt"


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def download(out_dir: Path = _OUT_DIR) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    dest = out_dir / _FILENAME
    logger.info("Downloading Fuentes-Muñoz 2025 Table 5 from %s", _URL)
    req = urllib.request.Request(_URL, headers={"User-Agent": "gaia-asteroid-encounters/1.0"})
    with urllib.request.urlopen(req) as resp:  # noqa: S310 (trusted official IOP URL)
        dest.write_bytes(resp.read())
    sha = _sha256(dest)
    meta = {
        "source_url": _URL,
        "reference": "Fuentes-Muñoz et al. 2025, AJ 170, 353 (DOI 10.3847/1538-3881/ae0cc9)",
        "table": "Table 5 — Asteroid initial masses and uncertainties",
        "downloaded_at": datetime.now(UTC).isoformat(),
        "filename": _FILENAME,
        "size_bytes": dest.stat().st_size,
        "sha256": sha,
    }
    (out_dir / "metadata.json").write_text(json.dumps(meta, indent=2))
    logger.info("Wrote %s (%d bytes, sha256 %s…)", dest, meta["size_bytes"], sha[:16])
    return dest


def main() -> int:
    parser = argparse.ArgumentParser(description="Download Fuentes-Muñoz 2025 Table 5")
    parser.add_argument("--out-dir", default=str(_OUT_DIR))
    args = parser.parse_args()
    download(Path(args.out_dir))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
