"""Tests for src/ingest/gaia_sso.py.

The TAP download itself is not tested here (requires network + Gaia archive).
Tests cover: column validation, parallel mock download, load/save round-trip,
and epoch scale documentation.
"""

from pathlib import Path
from unittest.mock import MagicMock, patch

import polars as pl
import pytest

from src.ingest.gaia_sso import _REQUIRED_COLUMNS, download_gaia_sso, load_gaia_sso

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_sso_df() -> pl.DataFrame:
    """Minimal DataFrame that mimics a real Gaia SSO download."""
    return pl.DataFrame(
        {
            "solution_id": [1, 1, 1],
            "source_id": [100, 200, 300],
            "denomination": ["Ceres", "Ceres", "Vesta"],
            "number_mp": [1, 1, 4],
            "transit_id": [10, 11, 20],
            "observation_id": [1001, 1002, 2001],
            "epoch": [2456863.5, 2456864.0, 2456865.0],
            "epoch_utc": [2456863.4, 2456863.9, 2456864.9],
            "ra": [180.1, 180.2, 90.3],
            "dec": [5.1, 5.2, -3.4],
            "g_mag": [7.5, 7.6, 6.8],
            "x_gaia": [0.1, 0.11, 0.2],
            "y_gaia": [0.9, 0.91, -0.5],
            "z_gaia": [0.05, 0.06, 0.3],
        }
    )


COLUMNS = [
    "solution_id", "source_id", "denomination", "number_mp",
    "transit_id", "observation_id", "epoch", "epoch_utc",
    "ra", "dec", "g_mag", "x_gaia", "y_gaia", "z_gaia",
]


def _tap_mock_for(df: pl.DataFrame) -> MagicMock:
    """Return a TapPlus mock whose async job yields *df*."""
    from astropy.table import Table

    table = Table.from_pandas(df.to_pandas())
    job = MagicMock()
    job.get_results.return_value = table

    tap = MagicMock()
    tap.launch_job_async.return_value = job
    return tap


# ---------------------------------------------------------------------------
# Column validation
# ---------------------------------------------------------------------------


def test_missing_required_column_raises(tmp_path: Path) -> None:
    bad_columns = [c for c in COLUMNS if c != "epoch"]
    with pytest.raises(ValueError, match="Missing required columns"):
        download_gaia_sso(
            archive_url="http://fake",
            columns=bad_columns,
            dest=tmp_path / "out.parquet",
            mp_max=10,
        )


def test_required_columns_set_is_nonempty() -> None:
    assert len(_REQUIRED_COLUMNS) >= 4


# ---------------------------------------------------------------------------
# download_gaia_sso — mocked parallel TAP
#
# Pass mp_max=10, n_workers=2 so the function creates 3 ranges:
#   (1,5), (6,10), (None,None) — all answered by the same mock.
# ---------------------------------------------------------------------------


@patch("src.ingest.gaia_sso.TapPlus")
def test_download_writes_parquet(mock_tap_cls: MagicMock, tmp_path: Path) -> None:
    mock_tap_cls.return_value = _tap_mock_for(_make_sso_df())
    dest = tmp_path / "gaia_sso.parquet"
    download_gaia_sso(
        archive_url="http://fake", columns=COLUMNS, dest=dest,
        mp_max=10, n_workers=2,
    )
    assert dest.exists()


@patch("src.ingest.gaia_sso.TapPlus")
def test_download_returns_dataframe(mock_tap_cls: MagicMock, tmp_path: Path) -> None:
    mock_tap_cls.return_value = _tap_mock_for(_make_sso_df())
    df = download_gaia_sso(
        archive_url="http://fake", columns=COLUMNS,
        dest=tmp_path / "out.parquet",
        mp_max=10, n_workers=2,
    )
    assert isinstance(df, pl.DataFrame)


@patch("src.ingest.gaia_sso.TapPlus")
def test_ceres_observations_present(mock_tap_cls: MagicMock, tmp_path: Path) -> None:
    mock_tap_cls.return_value = _tap_mock_for(_make_sso_df())
    df = download_gaia_sso(
        archive_url="http://fake", columns=COLUMNS,
        dest=tmp_path / "out.parquet",
        mp_max=10, n_workers=2,
    )
    # Each range returns the same mock df with Ceres rows
    assert len(df.filter(pl.col("number_mp") == 1)) > 0


@patch("src.ingest.gaia_sso.TapPlus")
def test_progress_logging(mock_tap_cls: MagicMock, tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    """Verify that progress lines with percentage and timing are logged."""
    mock_tap_cls.return_value = _tap_mock_for(_make_sso_df())
    import logging
    with caplog.at_level(logging.INFO, logger="src.ingest.gaia_sso"):
        download_gaia_sso(
            archive_url="http://fake", columns=COLUMNS,
            dest=tmp_path / "out.parquet",
            mp_max=10, n_workers=2,
        )
    progress_lines = [r for r in caplog.records if "%" in r.message and "s" in r.message]
    assert len(progress_lines) > 0, "Expected progress log lines with % and timing"


# ---------------------------------------------------------------------------
# load_gaia_sso
# ---------------------------------------------------------------------------


def test_load_roundtrip(tmp_path: Path) -> None:
    df = _make_sso_df()
    dest = tmp_path / "gaia_sso.parquet"
    df.write_parquet(dest)
    loaded = load_gaia_sso(dest)
    assert loaded.shape == df.shape
    assert set(loaded.columns) == set(df.columns)


def test_load_missing_file_raises(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        load_gaia_sso(tmp_path / "nonexistent.parquet")


def test_epoch_column_present(tmp_path: Path) -> None:
    df = _make_sso_df()
    dest = tmp_path / "gaia_sso.parquet"
    df.write_parquet(dest)
    assert "epoch" in load_gaia_sso(dest).columns
