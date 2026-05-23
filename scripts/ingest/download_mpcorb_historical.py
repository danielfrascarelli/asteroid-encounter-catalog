"""Download a historical MPCORB.DAT snapshot via the Internet Archive Wayback Machine.

The MPC only serves the current MPCORB.DAT — older versions are not archived on
minorplanetcenter.net. The Wayback Machine, however, has captured the URL
several times between 2014 and 2016 (inside the Gaia DR3 observation window).

For accurate Kepler 2-body propagation of asteroids inside a target time window,
the snapshot's osculating epoch should sit close to that window. Otherwise the
neglected planetary perturbations accumulate, and the propagated positions
drift away from the true positions (already documented for the Doris–Geraldina
pair in the project notes).

Usage
-----
    # Closest snapshot to mid-2015
    docker compose run --rm pipeline python -m scripts.download_mpcorb_historical \\
        --year 2015 --month 6

    # Pick the closest snapshot to a specific date
    docker compose run --rm pipeline python -m scripts.download_mpcorb_historical \\
        --target-date 2015-07-01

Output
------
Writes ``data/raw/mpcorb_archive/MPCORB_<YYYYMMDD>.DAT`` (where ``<YYYYMMDD>``
is the Wayback capture timestamp) plus a JSON sidecar with metadata
(snapshot date, source URL, epoch_jd from the first record, SHA-256).
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

from src.ingest.mpcorb_archive import ARCHIVE_SUBDIR, write_sidecar
from src.utils.config import load_config

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)

MPC_URL = "http://www.minorplanetcenter.net/iau/MPCORB/MPCORB.DAT"
WAYBACK_AVAILABILITY = "http://archive.org/wayback/available"


def _query_wayback_closest(target_yyyymmdd: str) -> tuple[str, str]:
    """Ask the Wayback API for the closest snapshot to *target_yyyymmdd*.

    Returns
    -------
    (timestamp, raw_url)
        ``timestamp`` is the 14-digit Wayback identifier (YYYYMMDDhhmmss).
        ``raw_url`` is the full URL with the ``id_`` flag that returns the
        original bytes (no Wayback HTML wrapper).

    Raises
    ------
    RuntimeError
        If no snapshot is available near the requested date.
    """
    api = f"{WAYBACK_AVAILABILITY}?url={MPC_URL}&timestamp={target_yyyymmdd}"
    logger.info("Querying Wayback availability: %s", api)
    req = urllib.request.Request(api, headers={"User-Agent": "gaia-asteroid-pipeline/1.0"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req) as resp:  # noqa: S310
                payload = json.loads(resp.read().decode("utf-8"))
            break
        except urllib.error.HTTPError as exc:
            if exc.code == 429 and attempt < 3:
                wait = 2 ** (attempt + 1)
                logger.warning("Wayback API rate-limited (429); retrying in %ds…", wait)
                time.sleep(wait)
            else:
                raise
    else:
        raise RuntimeError("Wayback API rate-limited after 4 attempts")

    closest = payload.get("archived_snapshots", {}).get("closest")
    if not closest or not closest.get("available"):
        raise RuntimeError(f"No Wayback snapshot near {target_yyyymmdd}")

    timestamp = closest["timestamp"]  # YYYYMMDDhhmmss
    # Build the raw-byte URL ourselves (more reliable than `closest['url']`).
    raw_url = f"https://web.archive.org/web/{timestamp}id_/{MPC_URL}"
    return timestamp, raw_url


def _download_streaming(url: str, dest: Path) -> int:
    """Download *url* to *dest* with chunked writes. Returns bytes written."""
    logger.info("Downloading %s …", url)
    bytes_written = 0
    with urllib.request.urlopen(url) as resp:  # noqa: S310
        with dest.open("wb") as fh:
            while True:
                chunk = resp.read(1 << 20)  # 1 MiB
                if not chunk:
                    break
                fh.write(chunk)
                bytes_written += len(chunk)
    return bytes_written


def download_historical(
    config_path: str = "config.yaml",
    target_date: str = "2015-06-15",
) -> Path:
    """Download the Wayback snapshot closest to *target_date* (ISO).

    Returns the path to the downloaded ``.DAT`` file.
    """
    cfg = load_config(config_path)
    archive_dir = Path(cfg.paths.raw) / ARCHIVE_SUBDIR
    archive_dir.mkdir(parents=True, exist_ok=True)

    target_yyyymmdd = target_date.replace("-", "")
    timestamp, raw_url = _query_wayback_closest(target_yyyymmdd)

    capture_iso = f"{timestamp[0:4]}-{timestamp[4:6]}-{timestamp[6:8]}"
    dest = archive_dir / f"MPCORB_{timestamp[0:8]}.DAT"

    if dest.exists():
        logger.info("Snapshot already present: %s — skipping download", dest)
    else:
        n = _download_streaming(raw_url, dest)
        logger.info("Wrote %d bytes to %s (capture %s)", n, dest, capture_iso)

    sidecar = write_sidecar(dest, snapshot_date=capture_iso, source_url=raw_url)
    meta = json.loads(sidecar.read_text())
    logger.info(
        "Snapshot %s — epoch JD %.1f (TDB), size %d bytes",
        capture_iso,
        meta["epoch_jd_tdb"],
        meta["size_bytes"],
    )
    return dest


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--config", default="config.yaml")
    p.add_argument(
        "--target-date",
        help="ISO date (YYYY-MM-DD) — download the Wayback capture closest to this date.",
    )
    p.add_argument("--year", type=int, help="Shortcut: pick year, use --month or default 6.")
    p.add_argument("--month", type=int, default=6, help="Month used with --year (default: 6).")
    args = p.parse_args()

    if args.target_date and args.year:
        p.error("Use either --target-date or --year, not both.")

    if args.target_date:
        target = args.target_date
    elif args.year:
        target = f"{args.year:04d}-{args.month:02d}-15"
    else:
        p.error("Provide --target-date or --year.")
        return 2  # unreachable, silence type checker

    download_historical(args.config, target_date=target)
    return 0


if __name__ == "__main__":
    sys.exit(main())
