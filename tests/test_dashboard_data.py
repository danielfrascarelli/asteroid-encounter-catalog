"""Tests for the dashboard data layer (memory-safe catalog access)."""

from __future__ import annotations

from pathlib import Path

import polars as pl

from src.dashboard import data as dd


def _make_catalog(tmp_path: Path, n: int) -> Path:
    df = pl.DataFrame(
        {
            "number_1": list(range(n)),
            "number_2": list(range(100, 100 + n)),
            "dist_au": [0.05 - i * (0.04 / max(n, 1)) for i in range(n)],  # descending
            "gaia_observable": [i % 2 == 0 for i in range(n)],
            "class_1": ["MBA"] * n,
            "class_2": ["MBA"] * n,
        }
    )
    path = tmp_path / "cat.parquet"
    df.write_parquet(path)
    return path


def test_stats_are_global(tmp_path: Path) -> None:
    path = _make_catalog(tmp_path, 100)
    stats = dd.catalog_stats(path)
    assert stats["n_total"] == 100
    assert stats["n_gaia_observable"] == 50
    # smallest dist is the last row: 0.05 - 99*(0.04/100)
    assert stats["dist_min_au"] < 0.05


def test_load_uncapped_below_cap(tmp_path: Path) -> None:
    path = _make_catalog(tmp_path, 100)
    df, capped = dd.load_catalog_display(path, cap=300)
    assert not capped
    assert len(df) == 100


def test_load_capped_returns_closest(tmp_path: Path) -> None:
    path = _make_catalog(tmp_path, 1000)
    df, capped = dd.load_catalog_display(path, cap=200)
    assert capped
    assert len(df) == 200
    # capped slice must be the 200 closest encounters
    assert df["dist_au"].max() <= 0.05
    full = pl.read_parquet(path)
    assert df["dist_au"].max() <= full["dist_au"].sort()[199] + 1e-12


def test_resolve_prefers_full(tmp_path: Path, monkeypatch) -> None:
    full = tmp_path / "full.parquet"
    small = tmp_path / "small.parquet"
    _make_catalog(tmp_path, 5).rename(full)
    pl.read_parquet(full).write_parquet(small)
    monkeypatch.setattr(dd, "_FULL_CATALOG", full)
    monkeypatch.setattr(dd, "_SMALL_CATALOG", small)
    assert dd.resolve_catalog_path() == full
    # falls back to small when full is absent
    monkeypatch.setattr(dd, "_FULL_CATALOG", tmp_path / "missing.parquet")
    assert dd.resolve_catalog_path() == small
