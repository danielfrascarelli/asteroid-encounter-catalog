"""Tests for src.propagate.cache (trajectory cache key generation + manifest)."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from src.propagate.cache import build_cache_key


def _kwargs(integrator: str = "whfast", planets: list[str] | None = None) -> dict:
    return {
        "integrator": integrator,
        "dt_days": 1.0,
        "include_planets": planets or ["sun", "jupiter", "saturn"],
        "include_major_asteroids": False,
    }


def test_cache_key_is_deterministic() -> None:
    """Calling build_cache_key twice with the same inputs must yield the same key."""
    grid = np.arange(2457000.0, 2458000.0, 1.0)
    k1 = build_cache_key(
        snapshot_sha="abc123",
        time_grid=grid,
        method="rebound",
        rebound_kwargs=_kwargs(),
        n_asteroids=10,
    )
    k2 = build_cache_key(
        snapshot_sha="abc123",
        time_grid=grid,
        method="rebound",
        rebound_kwargs=_kwargs(),
        n_asteroids=10,
    )
    assert k1 == k2
    assert len(k1) == 16
    assert all(c in "0123456789abcdef" for c in k1)


def test_cache_key_changes_with_snapshot_sha() -> None:
    grid = np.arange(2457000.0, 2458000.0, 1.0)
    k1 = build_cache_key(
        snapshot_sha="abc123",
        time_grid=grid,
        method="rebound",
        rebound_kwargs=_kwargs(),
        n_asteroids=10,
    )
    k2 = build_cache_key(
        snapshot_sha="def456",
        time_grid=grid,
        method="rebound",
        rebound_kwargs=_kwargs(),
        n_asteroids=10,
    )
    assert k1 != k2


def test_cache_key_changes_with_planets() -> None:
    grid = np.arange(2457000.0, 2458000.0, 1.0)
    k1 = build_cache_key(
        snapshot_sha="abc123",
        time_grid=grid,
        method="rebound",
        rebound_kwargs=_kwargs(planets=["sun", "jupiter"]),
        n_asteroids=10,
    )
    k2 = build_cache_key(
        snapshot_sha="abc123",
        time_grid=grid,
        method="rebound",
        rebound_kwargs=_kwargs(planets=["sun", "jupiter", "saturn"]),
        n_asteroids=10,
    )
    assert k1 != k2


def test_cache_key_changes_with_n_asteroids() -> None:
    grid = np.arange(2457000.0, 2458000.0, 1.0)
    k1 = build_cache_key(
        snapshot_sha="abc123",
        time_grid=grid,
        method="rebound",
        rebound_kwargs=_kwargs(),
        n_asteroids=10,
    )
    k2 = build_cache_key(
        snapshot_sha="abc123",
        time_grid=grid,
        method="rebound",
        rebound_kwargs=_kwargs(),
        n_asteroids=100,
    )
    assert k1 != k2


def test_cache_key_changes_with_time_grid() -> None:
    grid_a = np.arange(2457000.0, 2458000.0, 1.0)
    grid_b = np.arange(2457000.0, 2459000.0, 1.0)  # different end
    k1 = build_cache_key(
        snapshot_sha="abc123",
        time_grid=grid_a,
        method="rebound",
        rebound_kwargs=_kwargs(),
        n_asteroids=10,
    )
    k2 = build_cache_key(
        snapshot_sha="abc123",
        time_grid=grid_b,
        method="rebound",
        rebound_kwargs=_kwargs(),
        n_asteroids=10,
    )
    assert k1 != k2


def test_cache_key_accepts_path_for_sha(tmp_path: Path) -> None:
    """When snapshot_sha is a Path, the file is hashed."""
    f = tmp_path / "tiny.bin"
    f.write_bytes(b"hello world")
    grid = np.arange(2457000.0, 2458000.0, 1.0)
    k = build_cache_key(
        snapshot_sha=f,
        time_grid=grid,
        method="rebound",
        rebound_kwargs=_kwargs(),
        n_asteroids=10,
    )
    assert len(k) == 16
