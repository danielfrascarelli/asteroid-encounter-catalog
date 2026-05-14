"""Tests for src/ingest/gaia_orbits.py.

The TAP download itself is not tested here (requires network + Gaia archive).
Tests cover: column mapping, angle conversion, parallel mock download,
load/save round-trip, chunk caching/resume, and retry logic.
"""

from __future__ import annotations

import math
from pathlib import Path
from unittest.mock import MagicMock, patch

import polars as pl
import pytest

from src.ingest.gaia_orbits import (
    _chunk_path,
    _fetch_chunk,
    _to_pipeline_schema,
    download_gaia_orbits,
    load_gaia_orbits,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_DEG = math.pi / 180.0


def _make_raw_df() -> pl.DataFrame:
    """Minimal DataFrame mimicking a raw gaiadr3.sso_orbits TAP response."""
    return pl.DataFrame(
        {
            "number_mp": [1, 2, 4],
            "denomination": ["Ceres", "Pallas", "Vesta"],
            "semi_major_axis": [2.769, 2.774, 2.362],
            "eccentricity": [0.076, 0.231, 0.089],
            "inclination": [10.6 * _DEG, 34.8 * _DEG, 7.1 * _DEG],
            "arg_perihelion": [73.6 * _DEG, 310.0 * _DEG, 151.2 * _DEG],
            "long_asc_node": [80.3 * _DEG, 173.1 * _DEG, 103.9 * _DEG],
            "mean_anomaly": [77.4 * _DEG, 78.5 * _DEG, 20.9 * _DEG],
            "osc_epoch": [1666.0, 1666.0, 1666.0],  # days since J2010.0; 1666+2455197.5=2456863.5
        }
    )


def _tap_mock_for(df: pl.DataFrame) -> MagicMock:
    """Return a TapPlus mock whose sync job yields *df*."""
    from astropy.table import Table

    table = Table.from_pandas(df.to_pandas())
    job = MagicMock()
    job.get_results.return_value = table

    tap = MagicMock()
    tap.launch_job.return_value = job
    return tap


# ---------------------------------------------------------------------------
# Chunk path naming
# ---------------------------------------------------------------------------


def test_chunk_path_naming(tmp_path: Path) -> None:
    assert _chunk_path(tmp_path, 1, 5000).name == "orb_0000001_0005000.parquet"
    assert _chunk_path(tmp_path, 150001, 155000).name == "orb_0150001_0155000.parquet"


# ---------------------------------------------------------------------------
# _to_pipeline_schema: column renaming and unit conversion
# ---------------------------------------------------------------------------


def test_to_pipeline_schema_has_expected_columns() -> None:
    result = _to_pipeline_schema(_make_raw_df())
    expected = {
        "number",
        "designation",
        "a_au",
        "e",
        "i_deg",
        "omega_deg",
        "Omega_deg",
        "M_deg",
        "epoch_jd",
    }
    assert set(result.columns) == expected


def test_to_pipeline_schema_drops_radian_columns() -> None:
    result = _to_pipeline_schema(_make_raw_df())
    for col in ("inclination", "arg_perihelion", "long_asc_node", "mean_anomaly"):
        assert col not in result.columns


def test_to_pipeline_schema_converts_inclination() -> None:
    df = _make_raw_df()
    result = _to_pipeline_schema(df)
    expected = df["inclination"][0] * (180.0 / math.pi)
    assert abs(result["i_deg"][0] - expected) < 1e-10


def test_to_pipeline_schema_converts_arg_perihelion() -> None:
    df = _make_raw_df()
    result = _to_pipeline_schema(df)
    expected = df["arg_perihelion"][0] * (180.0 / math.pi)
    assert abs(result["omega_deg"][0] - expected) < 1e-10


def test_to_pipeline_schema_converts_long_asc_node() -> None:
    df = _make_raw_df()
    result = _to_pipeline_schema(df)
    expected = df["long_asc_node"][0] * (180.0 / math.pi)
    assert abs(result["Omega_deg"][0] - expected) < 1e-10


def test_to_pipeline_schema_converts_mean_anomaly() -> None:
    df = _make_raw_df()
    result = _to_pipeline_schema(df)
    expected = df["mean_anomaly"][0] * (180.0 / math.pi)
    assert abs(result["M_deg"][0] - expected) < 1e-10


def test_to_pipeline_schema_casts_number_to_int32() -> None:
    result = _to_pipeline_schema(_make_raw_df())
    assert result["number"].dtype == pl.Int32


def test_to_pipeline_schema_preserves_epoch_jd() -> None:
    df = _make_raw_df()
    result = _to_pipeline_schema(df)
    assert result["epoch_jd"][0] == pytest.approx(2456863.5)


# ---------------------------------------------------------------------------
# download_gaia_orbits — mocked parallel TAP
# ---------------------------------------------------------------------------


@patch("src.ingest.gaia_orbits.TapPlus")
def test_download_writes_parquet(mock_tap_cls: MagicMock, tmp_path: Path) -> None:
    mock_tap_cls.return_value = _tap_mock_for(_make_raw_df())
    dest = tmp_path / "gaia_orbits.parquet"
    download_gaia_orbits(
        archive_url="http://fake",
        dest=dest,
        mp_max=10,
        n_workers=2,
        cache_dir=tmp_path / "chunks",
    )
    assert dest.exists()


@patch("src.ingest.gaia_orbits.TapPlus")
def test_download_returns_dataframe(mock_tap_cls: MagicMock, tmp_path: Path) -> None:
    mock_tap_cls.return_value = _tap_mock_for(_make_raw_df())
    df = download_gaia_orbits(
        archive_url="http://fake",
        dest=tmp_path / "out.parquet",
        mp_max=10,
        n_workers=2,
        cache_dir=tmp_path / "chunks",
    )
    assert isinstance(df, pl.DataFrame)


@patch("src.ingest.gaia_orbits.TapPlus")
def test_download_returns_pipeline_columns(mock_tap_cls: MagicMock, tmp_path: Path) -> None:
    mock_tap_cls.return_value = _tap_mock_for(_make_raw_df())
    df = download_gaia_orbits(
        archive_url="http://fake",
        dest=tmp_path / "out.parquet",
        mp_max=10,
        n_workers=2,
        cache_dir=tmp_path / "chunks",
    )
    for col in (
        "number",
        "designation",
        "a_au",
        "e",
        "i_deg",
        "Omega_deg",
        "omega_deg",
        "M_deg",
        "epoch_jd",
    ):
        assert col in df.columns, f"Missing column: {col}"


@patch("src.ingest.gaia_orbits.TapPlus")
def test_download_angles_in_degrees(mock_tap_cls: MagicMock, tmp_path: Path) -> None:
    """All angular columns must be in [0, 360) — confirming rad→deg conversion happened."""
    mock_tap_cls.return_value = _tap_mock_for(_make_raw_df())
    df = download_gaia_orbits(
        archive_url="http://fake",
        dest=tmp_path / "out.parquet",
        mp_max=10,
        n_workers=2,
        cache_dir=tmp_path / "chunks",
    )
    for col in ("i_deg", "omega_deg", "Omega_deg", "M_deg"):
        vals = df[col].to_list()
        for v in vals:
            assert 0.0 <= v < 360.0, f"{col}={v} not in [0, 360) — still in radians?"


# ---------------------------------------------------------------------------
# Resume: skip cached ranges
# ---------------------------------------------------------------------------


@patch("src.ingest.gaia_orbits.TapPlus")
def test_resume_skips_cached_ranges(mock_tap_cls: MagicMock, tmp_path: Path) -> None:
    """TapPlus is NOT called for ranges that already have a chunk file."""
    mock_tap_cls.return_value = _tap_mock_for(_make_raw_df())
    cache_dir = tmp_path / "chunks"
    cache_dir.mkdir()
    # Pre-populate first batch (batch_size=5, mp_max=10 → ranges: (1,5), (6,10))
    _make_raw_df().write_parquet(_chunk_path(cache_dir, 1, 5))
    download_gaia_orbits(
        archive_url="http://fake",
        dest=tmp_path / "out.parquet",
        mp_max=10,
        n_workers=2,
        batch_size=5,
        cache_dir=cache_dir,
    )
    # (1,5) cached → only (6,10) is a TAP call
    assert mock_tap_cls.call_count == 1


# ---------------------------------------------------------------------------
# Retry on failure
# ---------------------------------------------------------------------------


@patch("src.ingest.gaia_orbits.time.sleep", return_value=None)
@patch("src.ingest.gaia_orbits.TapPlus")
def test_retry_on_failure(mock_tap_cls: MagicMock, mock_sleep: MagicMock, tmp_path: Path) -> None:
    """_fetch_chunk retries on failure and succeeds on the Nth attempt."""
    good_tap = _tap_mock_for(_make_raw_df())
    mock_tap_cls.side_effect = [RuntimeError("boom"), RuntimeError("boom"), good_tap]
    cache_dir = tmp_path / "chunks"
    cache_dir.mkdir()
    cp = _chunk_path(cache_dir, 1, 5)
    _fetch_chunk("http://fake", 1, 5, cache_path=cp, max_retries=3)
    assert cp.exists()
    assert mock_tap_cls.call_count == 3
    assert mock_sleep.call_count == 2  # slept before attempt 2 and 3


# ---------------------------------------------------------------------------
# load_gaia_orbits
# ---------------------------------------------------------------------------


@patch("src.ingest.gaia_orbits.TapPlus")
def test_load_roundtrip(mock_tap_cls: MagicMock, tmp_path: Path) -> None:
    mock_tap_cls.return_value = _tap_mock_for(_make_raw_df())
    dest = tmp_path / "gaia_orbits.parquet"
    download_gaia_orbits(
        archive_url="http://fake",
        dest=dest,
        mp_max=10,
        n_workers=1,
        cache_dir=tmp_path / "chunks",
    )
    loaded = load_gaia_orbits(dest)
    assert isinstance(loaded, pl.DataFrame)
    assert "epoch_jd" in loaded.columns
    assert "a_au" in loaded.columns


def test_load_missing_file_raises(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        load_gaia_orbits(tmp_path / "nonexistent.parquet")
