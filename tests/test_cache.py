"""Tests for src.propagate.cache (trajectory cache key generation + manifest)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import polars as pl
import pytest

from src.propagate.cache import build_cache_key, load_or_compute_trajectory


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


def test_zarr_cache_writes_transposed_compressed_chunked_array(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Zarr cache stores (3, N, T) transposed with BitRound+Delta+Blosc-zstd."""
    pytest.importorskip("zarr")
    from numcodecs import BitRound, Blosc, Delta

    from src.propagate.cache import TrajectoryView
    import src.propagate.nbody as nbody

    elements = pl.DataFrame({"number": [1, 2]})
    time_grid = np.array([2457000.5, 2457001.5, 2457002.5], dtype=np.float64)
    # Use a smooth signal across time so the Delta filter has something to
    # exploit; small mantissas survive BitRound(16) exactly.
    base = np.arange(18, dtype=np.float32).reshape(3, 2, 3)
    expected = base
    calls = {"n": 0}

    def fake_propagate_grid_nbody(
        elements: pl.DataFrame,
        time_grid: np.ndarray,
        *,
        out,
        **kwargs,
    ):
        calls["n"] += 1
        out[:] = expected
        return out

    monkeypatch.setattr(nbody, "propagate_grid_nbody", fake_propagate_grid_nbody)

    view = load_or_compute_trajectory(
        elements=elements,
        time_grid=time_grid,
        cache_dir=tmp_path,
        cache_key="abc123",
        rebound_kwargs=_kwargs(),
        cache_format="zarr",
    )

    assert calls["n"] == 1
    assert isinstance(view, TrajectoryView)
    # Consumers see the logical (T, N, 3) shape regardless of on-disk layout.
    assert view.shape == expected.shape
    assert view.dtype == np.dtype("float32")

    # Inspect the underlying zarr array directly for filter/compressor checks.
    inner = view._z
    T, N = expected.shape[0], expected.shape[1]
    assert inner.shape == (3, N, T)  # transposed on disk
    # T_chunk gets clamped to T when shorter than the default.
    assert inner.chunks == (3, N, T)
    assert inner.dtype == np.dtype("float32")
    assert inner.compressor.cname == "zstd"
    assert inner.compressor.clevel == 5
    assert inner.compressor.shuffle == Blosc.BITSHUFFLE
    assert inner.filters is not None and len(inner.filters) == 2
    assert isinstance(inner.filters[0], BitRound)
    assert inner.filters[0].keepbits == 16
    assert isinstance(inner.filters[1], Delta)
    np.testing.assert_array_equal(np.asarray(view[:]), expected)

    def fail_if_called(*args, **kwargs):
        raise AssertionError("cache hit should not recompute")

    monkeypatch.setattr(nbody, "propagate_grid_nbody", fail_if_called)
    view_hit = load_or_compute_trajectory(
        elements=elements,
        time_grid=time_grid,
        cache_dir=tmp_path,
        cache_key="abc123",
        rebound_kwargs=_kwargs(),
        cache_format="zarr",
    )
    assert isinstance(view_hit, TrajectoryView)
    np.testing.assert_array_equal(np.asarray(view_hit[:]), expected)


def test_trajectory_view_indexing(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """TrajectoryView maps ``[t]`` and ``[k0:k1]`` against the transposed zarr."""
    pytest.importorskip("zarr")
    from src.propagate.cache import TrajectoryView
    import src.propagate.nbody as nbody

    elements = pl.DataFrame({"number": [1, 2, 3]})
    time_grid = np.array(
        [2457000.5, 2457001.5, 2457002.5, 2457003.5, 2457004.5], dtype=np.float64
    )
    expected = np.arange(45, dtype=np.float32).reshape(5, 3, 3)

    def fake(elements, time_grid, *, out, **kwargs):
        out[:] = expected
        return out

    monkeypatch.setattr(nbody, "propagate_grid_nbody", fake)

    view = load_or_compute_trajectory(
        elements=elements,
        time_grid=time_grid,
        cache_dir=tmp_path,
        cache_key="view-idx",
        rebound_kwargs=_kwargs(),
        cache_format="zarr",
    )
    assert isinstance(view, TrajectoryView)

    # Single-step indexing
    for t in range(5):
        np.testing.assert_array_equal(view[t], expected[t])

    # Slice indexing
    np.testing.assert_array_equal(view[1:4], expected[1:4])

    # Round-trip through __array__
    np.testing.assert_array_equal(np.asarray(view), expected)


def test_open_trajectory_for_worker_round_trip(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A worker re-opening the cached zarr sees an equivalent TrajectoryView."""
    pytest.importorskip("zarr")
    from src.propagate.cache import TrajectoryView, open_trajectory_for_worker
    import src.propagate.nbody as nbody

    elements = pl.DataFrame({"number": [1, 2]})
    time_grid = np.array([2457000.5, 2457001.5, 2457002.5], dtype=np.float64)
    expected = np.arange(18, dtype=np.float32).reshape(3, 2, 3)

    def fake(elements, time_grid, *, out, **kwargs):
        out[:] = expected
        return out

    monkeypatch.setattr(nbody, "propagate_grid_nbody", fake)

    view = load_or_compute_trajectory(
        elements=elements,
        time_grid=time_grid,
        cache_dir=tmp_path,
        cache_key="worker-roundtrip",
        rebound_kwargs=_kwargs(),
        cache_format="zarr",
    )
    assert isinstance(view, TrajectoryView)
    zarr_path = view.zarr_path
    assert zarr_path is not None

    worker_view = open_trajectory_for_worker(zarr_path)
    assert isinstance(worker_view, TrajectoryView)
    assert worker_view.shape == expected.shape
    np.testing.assert_array_equal(np.asarray(worker_view[:]), expected)
    # zarr_path must still be discoverable even when the store is wrapped in
    # LRUStoreCache — otherwise scan_parallel can't pass the path on to its
    # own workers.
    assert worker_view.zarr_path is not None
