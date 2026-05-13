"""Tests for src/ingest/gaia_sso.py.

The TAP download itself is not tested here (requires network + Gaia archive).
Tests cover: column validation, load/save round-trip, and epoch scale docs.
"""

from pathlib import Path
from unittest.mock import MagicMock, patch

import polars as pl
import pytest

from src.ingest.gaia_sso import _REQUIRED_COLUMNS, download_gaia_sso, load_gaia_sso

# ---------------------------------------------------------------------------
# Fixtures
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
            "epoch": [2456863.5, 2456864.0, 2456865.0],  # TCB (fake)
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
        )


def test_required_columns_set_is_nonempty() -> None:
    assert len(_REQUIRED_COLUMNS) >= 4


# ---------------------------------------------------------------------------
# download_gaia_sso — mocked TAP
# ---------------------------------------------------------------------------


def _mock_tap_job(df: pl.DataFrame) -> MagicMock:
    """Return a mock TapPlus that yields *df* on the first call, empty on second."""
    import pandas as pd
    from astropy.table import Table

    pandas_df = df.to_pandas()
    astropy_table = Table.from_pandas(pandas_df)
    empty_table = Table.from_pandas(pd.DataFrame(columns=pandas_df.columns))

    job_full = MagicMock()
    job_full.get_results.return_value = astropy_table

    job_empty = MagicMock()
    job_empty.get_results.return_value = empty_table

    tap = MagicMock()
    tap.launch_job.side_effect = [job_full, job_empty]
    return tap


@patch("src.ingest.gaia_sso.TapPlus")
def test_download_writes_parquet(mock_tap_cls: MagicMock, tmp_path: Path) -> None:
    mock_tap_cls.return_value = _mock_tap_job(_make_sso_df())
    dest = tmp_path / "gaia_sso.parquet"
    df = download_gaia_sso(archive_url="http://fake", columns=COLUMNS, dest=dest)
    assert dest.exists()
    assert len(df) == 3


@patch("src.ingest.gaia_sso.TapPlus")
def test_download_returns_dataframe(mock_tap_cls: MagicMock, tmp_path: Path) -> None:
    mock_tap_cls.return_value = _mock_tap_job(_make_sso_df())
    df = download_gaia_sso(
        archive_url="http://fake",
        columns=COLUMNS,
        dest=tmp_path / "out.parquet",
    )
    assert isinstance(df, pl.DataFrame)


@patch("src.ingest.gaia_sso.TapPlus")
def test_ceres_observations_present(mock_tap_cls: MagicMock, tmp_path: Path) -> None:
    mock_tap_cls.return_value = _mock_tap_job(_make_sso_df())
    df = download_gaia_sso(
        archive_url="http://fake",
        columns=COLUMNS,
        dest=tmp_path / "out.parquet",
    )
    ceres_rows = df.filter(pl.col("number_mp") == 1)
    assert len(ceres_rows) > 0, "Expected observations for Ceres (number_mp=1)"


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
    """epoch column must be present — it holds TCB Julian Dates."""
    df = _make_sso_df()
    dest = tmp_path / "gaia_sso.parquet"
    df.write_parquet(dest)
    loaded = load_gaia_sso(dest)
    assert "epoch" in loaded.columns
