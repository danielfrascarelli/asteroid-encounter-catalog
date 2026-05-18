"""Tests for src.ingest.mpcorb_archive."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from src.ingest.mpcorb_archive import (
    MpcorbSnapshot,
    discover_snapshots,
    select_for_window,
    write_sidecar,
)


def _make_minimal_mpcorb_line(epoch_packed: str) -> str:
    """Build a single fixed-width MPCORB.DAT record with a parseable epoch + a_au.

    Only the columns used by the archive helpers (epoch field at 20:25 and
    semi-major axis at 92:103) are filled meaningfully; everything else is
    padding.
    """
    line = " " * 220
    # Insert the packed epoch at 20:25 (5 chars)
    line = line[:20] + epoch_packed + line[25:]
    # Insert semimajor axis at 92:103 (11 chars, e.g. "  2.7691330")
    a_str = "  2.7691330"
    line = line[:92] + a_str + line[103:]
    return line


def test_write_sidecar_then_load(tmp_path: Path) -> None:
    """A sidecar written next to a .DAT file must be discoverable as a snapshot."""
    dat = tmp_path / "MPCORB_20150524.DAT"
    dat.write_text(_make_minimal_mpcorb_line("K156R") + "\n")  # 2015-06-27 TT

    sidecar = write_sidecar(dat, snapshot_date="2015-05-24", source_url="http://example/x")
    assert sidecar.is_file()

    meta = json.loads(sidecar.read_text())
    assert meta["filename"] == "MPCORB_20150524.DAT"
    assert meta["snapshot_date"] == "2015-05-24"
    # Epoch K156R = 2015-06-27 TT ≈ JD 2457200.5 TDB
    assert 2457199.0 < meta["epoch_jd_tdb"] < 2457202.0
    assert len(meta["sha256"]) == 64


def test_discover_returns_sorted_by_epoch(tmp_path: Path) -> None:
    """Snapshots must be sorted by ascending epoch_jd."""
    archive = tmp_path / "mpcorb_archive"
    archive.mkdir()
    for name, packed in [
        ("MPCORB_20150524.DAT", "K156R"),  # 2015-06-27
        ("MPCORB_20130301.DAT", "K134H"),  # 2013-04-17
        ("MPCORB_20251201.DAT", "K25BF"),  # 2025-11-15
    ]:
        p = archive / name
        p.write_text(_make_minimal_mpcorb_line(packed) + "\n")
        write_sidecar(p, snapshot_date=name.split("_")[1][:4], source_url="x")

    snaps = discover_snapshots(tmp_path)
    assert len(snaps) == 3
    epochs = [s.epoch_jd for s in snaps]
    assert epochs == sorted(epochs)


def test_select_for_window_picks_closest_midpoint(tmp_path: Path) -> None:
    """For a target window centred on 2015-12, the 2015 snapshot must win over 2013/2025."""
    archive = tmp_path / "mpcorb_archive"
    archive.mkdir()
    samples = [
        ("MPCORB_20130301.DAT", "K134H"),  # epoch ≈ JD 2456400 (2013-04-17)
        ("MPCORB_20150524.DAT", "K156R"),  # epoch ≈ JD 2457200 (2015-06-27)
        ("MPCORB_20251201.DAT", "K25BF"),  # epoch ≈ JD 2461000 (2025-11-15)
    ]
    for name, packed in samples:
        p = archive / name
        p.write_text(_make_minimal_mpcorb_line(packed) + "\n")
        write_sidecar(p, snapshot_date="x", source_url="x")

    snaps = discover_snapshots(tmp_path)
    # Gaia window 2014-07-25 → 2017-05-28 in JD TDB
    t_start = 2456864.0
    t_end = 2457902.0
    chosen = select_for_window(snaps, t_start, t_end)
    assert chosen.path.name == "MPCORB_20150524.DAT"


def test_select_for_window_empty_raises() -> None:
    with pytest.raises(ValueError):
        select_for_window([], 2457000.0, 2458000.0)


def test_select_with_single_snapshot(tmp_path: Path) -> None:
    """If only one snapshot is available, it must be selected regardless of window."""
    snap = MpcorbSnapshot(path=Path("/tmp/dummy.dat"), snapshot_date="x", epoch_jd=2400000.0)
    chosen = select_for_window([snap], 2500000.0, 2600000.0)
    assert chosen is snap
