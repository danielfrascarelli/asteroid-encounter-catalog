"""Tests for src/ingest/gaia_sso.py.

The TAP download itself is not tested here (requires network + Gaia archive).
Tests cover: column validation, parallel mock download, load/save round-trip,
chunk caching/resume, cross-batch reuse, and retry logic.
"""

from pathlib import Path
from unittest.mock import MagicMock, patch

import polars as pl
import pytest

from src.ingest.gaia_sso import (
    _REQUIRED_COLUMNS,
    _build_chunk_from_cache,
    _chunk_path,
    _fetch_range,
    download_gaia_sso,
    load_gaia_sso,
)

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
# Pass mp_max=10, n_workers=2 so the function creates 2 ranges:
#   (1,10), (None,None) — both answered by the same mock.
# ---------------------------------------------------------------------------


@patch("src.ingest.gaia_sso.TapPlus")
def test_download_writes_parquet(mock_tap_cls: MagicMock, tmp_path: Path) -> None:
    mock_tap_cls.return_value = _tap_mock_for(_make_sso_df())
    dest = tmp_path / "gaia_sso.parquet"
    download_gaia_sso(
        archive_url="http://fake", columns=COLUMNS, dest=dest,
        mp_max=10, n_workers=2, cache_dir=tmp_path / "chunks",
    )
    assert dest.exists()


@patch("src.ingest.gaia_sso.TapPlus")
def test_download_returns_dataframe(mock_tap_cls: MagicMock, tmp_path: Path) -> None:
    mock_tap_cls.return_value = _tap_mock_for(_make_sso_df())
    df = download_gaia_sso(
        archive_url="http://fake", columns=COLUMNS,
        dest=tmp_path / "out.parquet",
        mp_max=10, n_workers=2, cache_dir=tmp_path / "chunks",
    )
    assert isinstance(df, pl.DataFrame)


@patch("src.ingest.gaia_sso.TapPlus")
def test_ceres_observations_present(mock_tap_cls: MagicMock, tmp_path: Path) -> None:
    mock_tap_cls.return_value = _tap_mock_for(_make_sso_df())
    df = download_gaia_sso(
        archive_url="http://fake", columns=COLUMNS,
        dest=tmp_path / "out.parquet",
        mp_max=10, n_workers=2, cache_dir=tmp_path / "chunks",
    )
    assert len(df.filter(pl.col("number_mp") == 1)) > 0


@patch("src.ingest.gaia_sso.TapPlus")
def test_progress_logging(
    mock_tap_cls: MagicMock, tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """Verify that progress lines with percentage and timing are logged."""
    mock_tap_cls.return_value = _tap_mock_for(_make_sso_df())
    import logging
    with caplog.at_level(logging.INFO, logger="src.ingest.gaia_sso"):
        download_gaia_sso(
            archive_url="http://fake", columns=COLUMNS,
            dest=tmp_path / "out.parquet",
            mp_max=10, n_workers=2, cache_dir=tmp_path / "chunks",
        )
    progress_lines = [r for r in caplog.records if "%" in r.message and "s" in r.message]
    assert len(progress_lines) > 0, "Expected progress log lines with % and timing"


# ---------------------------------------------------------------------------
# Chunk path naming
# ---------------------------------------------------------------------------


def test_chunk_path_naming(tmp_path: Path) -> None:
    assert _chunk_path(tmp_path, 1, 5000).name == "mp_0000001_0005000.parquet"
    assert _chunk_path(tmp_path, None, None).name == "unnumbered.parquet"


# ---------------------------------------------------------------------------
# Resume: skip cached ranges
# ---------------------------------------------------------------------------


@patch("src.ingest.gaia_sso.TapPlus")
def test_resume_skips_cached_ranges(mock_tap_cls: MagicMock, tmp_path: Path) -> None:
    """TapPlus is NOT called for ranges that already have a chunk file."""
    mock_tap_cls.return_value = _tap_mock_for(_make_sso_df())
    cache_dir = tmp_path / "chunks"
    cache_dir.mkdir()
    # Pre-populate first batch (batch_size=5, mp_max=10 → ranges: (1,5), (6,10), unnumbered)
    _make_sso_df().write_parquet(_chunk_path(cache_dir, 1, 5))
    download_gaia_sso(
        archive_url="http://fake", columns=COLUMNS,
        dest=tmp_path / "out.parquet", mp_max=10, n_workers=2,
        batch_size=5, cache_dir=cache_dir,
    )
    # (1,5) is cached; (6,10) and unnumbered are pending → 2 TAP calls, not 3
    assert mock_tap_cls.call_count == 2


# ---------------------------------------------------------------------------
# Retry on failure
# ---------------------------------------------------------------------------


@patch("src.ingest.gaia_sso.time.sleep", return_value=None)
@patch("src.ingest.gaia_sso.TapPlus")
def test_retry_on_failure(
    mock_tap_cls: MagicMock, mock_sleep: MagicMock, tmp_path: Path
) -> None:
    """_fetch_range retries on failure and succeeds on the Nth attempt."""
    good_tap = _tap_mock_for(_make_sso_df())
    mock_tap_cls.side_effect = [RuntimeError("boom"), RuntimeError("boom"), good_tap]
    cache_dir = tmp_path / "chunks"
    cache_dir.mkdir()
    cp = _chunk_path(cache_dir, 1, 5)
    _fetch_range("http://fake", ", ".join(COLUMNS), 1, 5, cache_path=cp, max_retries=3)
    assert cp.exists()
    assert mock_tap_cls.call_count == 3
    assert mock_sleep.call_count == 2  # slept before attempt 2 and 3


# ---------------------------------------------------------------------------
# Cross-batch reuse
# ---------------------------------------------------------------------------


def test_cross_batch_reuse_builds_chunk(tmp_path: Path) -> None:
    """_build_chunk_from_cache merges two old chunks to cover a wider range."""
    cache_dir = tmp_path / "chunks"
    cache_dir.mkdir()

    # Two prior-run chunks: mp 1–5 and mp 6–10
    df_low = pl.DataFrame({"number_mp": [1, 2, 3], "source_id": [10, 20, 30]})
    df_high = pl.DataFrame({"number_mp": [6, 7], "source_id": [60, 70]})
    df_low.write_parquet(_chunk_path(cache_dir, 1, 5))
    df_high.write_parquet(_chunk_path(cache_dir, 6, 10))

    dest = cache_dir / "mp_0000001_0010000.parquet"
    result = _build_chunk_from_cache(cache_dir, 1, 10, dest)

    assert result is True
    assert dest.exists()
    loaded = pl.read_parquet(dest)
    assert set(loaded["number_mp"].to_list()) == {1, 2, 3, 6, 7}


def test_cross_batch_reuse_gap_returns_false(tmp_path: Path) -> None:
    """_build_chunk_from_cache returns False when existing chunks leave a gap."""
    cache_dir = tmp_path / "chunks"
    cache_dir.mkdir()

    # Only mp 1–5 exists; we need 1–10
    pl.DataFrame({"number_mp": [1, 2], "source_id": [10, 20]}).write_parquet(
        _chunk_path(cache_dir, 1, 5)
    )
    dest = cache_dir / "mp_0000001_0010000.parquet"
    assert _build_chunk_from_cache(cache_dir, 1, 10, dest) is False
    assert not dest.exists()


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
