"""Discover and select among multiple MPCORB.DAT snapshots indexed by epoch.

The MPC publishes MPCORB with a rolling osculating epoch (typically updated
every few months). For accurate Kepler 2-body propagation, the snapshot's
epoch should be close to the target propagation window — otherwise neglected
planetary perturbations accumulate over the propagation interval.

This module manages an archive of MPCORB snapshots stored under
``data/raw/mpcorb_archive/``, each accompanied by a JSON sidecar with metadata
(snapshot date, source URL, epoch_jd of the contained records, SHA-256).
The current MPCORB at ``data/raw/MPCORB.DAT`` is also recognised if present.

Public API
----------
- ``MpcorbSnapshot``           — dataclass describing one snapshot.
- ``read_first_epoch_jd``      — extract epoch_jd from the first data line.
- ``write_sidecar``            — write the JSON metadata next to a .DAT file.
- ``discover_snapshots``       — find all available snapshots in ``data/raw``.
- ``select_for_window``        — pick the snapshot best suited for a JD window.
"""

from __future__ import annotations

import hashlib
import json
import logging
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

from src.ingest.mpcorb import unpack_epoch

logger = logging.getLogger(__name__)

ARCHIVE_SUBDIR = "mpcorb_archive"


@dataclass(frozen=True)
class MpcorbSnapshot:
    """One MPCORB.DAT file with metadata about its osculating epoch."""

    path: Path
    snapshot_date: str  # ISO date (YYYY-MM-DD) when the snapshot was captured
    epoch_jd: float     # JD TDB of the osculating epoch of records inside


def sha256(path: Path) -> str:
    """Hex digest of a file's SHA-256 (streaming, low memory)."""
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def read_first_epoch_jd(path: Path) -> float:
    """Return the epoch_jd (TDB) of the first parseable data record in *path*.

    MPCORB has a global header followed by one record per asteroid; all
    asteroids in a given file share the same osculating epoch in practice.
    Reading the first valid record is enough to identify the snapshot epoch.
    """
    with path.open(encoding="ascii", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if len(line) < 103:
                continue
            try:
                float(line[92:103])  # semimajor axis must be parseable
                return unpack_epoch(line[20:25])
            except (KeyError, ValueError, IndexError):
                continue
    raise ValueError(f"No parseable data records found in {path}")


def write_sidecar(
    dat_path: Path, *, snapshot_date: str, source_url: str, epoch_jd: float | None = None
) -> Path:
    """Write a JSON sidecar describing *dat_path*. Returns the sidecar path."""
    if epoch_jd is None:
        epoch_jd = read_first_epoch_jd(dat_path)
    meta = {
        "filename": dat_path.name,
        "snapshot_date": snapshot_date,
        "epoch_jd_tdb": epoch_jd,
        "source_url": source_url,
        "downloaded_at": datetime.now(UTC).isoformat(),
        "size_bytes": dat_path.stat().st_size,
        "sha256": sha256(dat_path),
    }
    sidecar = dat_path.with_suffix(".json")
    sidecar.write_text(json.dumps(meta, indent=2))
    return sidecar


def _load_sidecar(dat_path: Path) -> MpcorbSnapshot | None:
    """Return an MpcorbSnapshot if a sidecar exists for *dat_path*, else None."""
    sidecar = dat_path.with_suffix(".json")
    if not sidecar.is_file():
        return None
    try:
        meta = json.loads(sidecar.read_text())
        return MpcorbSnapshot(
            path=dat_path,
            snapshot_date=meta["snapshot_date"],
            epoch_jd=float(meta["epoch_jd_tdb"]),
        )
    except (KeyError, ValueError, json.JSONDecodeError) as exc:
        logger.warning("Ignoring malformed sidecar %s: %s", sidecar, exc)
        return None


def discover_snapshots(raw_dir: Path) -> list[MpcorbSnapshot]:
    """Find all MPCORB snapshots under *raw_dir*.

    Looks for:
      * ``raw_dir/MPCORB.DAT``                   (the legacy/current file)
      * ``raw_dir/mpcorb_archive/MPCORB_*.DAT``  (historical archive)

    Each must have a sidecar JSON with ``epoch_jd_tdb``. For the legacy file
    the sidecar is generated on the fly if absent.

    Returns
    -------
    list of MpcorbSnapshot, sorted by ``epoch_jd`` ascending.
    """
    raw_dir = Path(raw_dir)
    snapshots: list[MpcorbSnapshot] = []

    legacy = raw_dir / "MPCORB.DAT"
    if legacy.is_file():
        snap = _load_sidecar(legacy)
        if snap is None:
            # Generate sidecar on the fly using the file's epoch
            epoch_jd = read_first_epoch_jd(legacy)
            write_sidecar(
                legacy,
                snapshot_date="unknown",
                source_url="https://www.minorplanetcenter.net/iau/MPCORB/MPCORB.DAT.gz",
                epoch_jd=epoch_jd,
            )
            snap = _load_sidecar(legacy)
        if snap is not None:
            snapshots.append(snap)

    archive = raw_dir / ARCHIVE_SUBDIR
    if archive.is_dir():
        for dat in sorted(archive.glob("MPCORB_*.DAT")):
            snap = _load_sidecar(dat)
            if snap is None:
                logger.warning("Snapshot %s has no sidecar — skipping", dat)
                continue
            snapshots.append(snap)

    snapshots.sort(key=lambda s: s.epoch_jd)
    return snapshots


def select_for_window(
    snapshots: list[MpcorbSnapshot], t_start_jd: float, t_end_jd: float
) -> MpcorbSnapshot:
    """Pick the snapshot whose ``epoch_jd`` is closest to the window midpoint.

    Raises
    ------
    ValueError
        If *snapshots* is empty.
    """
    if not snapshots:
        raise ValueError("No MPCORB snapshots available")
    target = 0.5 * (t_start_jd + t_end_jd)
    return min(snapshots, key=lambda s: abs(s.epoch_jd - target))
