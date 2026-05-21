"""Integration test for the memmap-backed parallel scan path.

The fix in src/detect/parallel.py detects ``np.memmap`` ``positions`` and passes
the underlying file path to each worker (instead of pickling the array via
initargs). This test verifies that the memmap branch produces identical
results to the in-memory branch on a tiny example, and that the workers
correctly re-open the memmap from disk.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import polars as pl
import pytest

from src.detect.parallel import scan_parallel


def _make_elements(n: int) -> pl.DataFrame:
    """Tiny synthetic catalog: nearly-circular co-planar orbits ~2.7 AU."""
    rng = np.random.default_rng(seed=42)
    a = 2.7 + rng.normal(0.0, 0.001, size=n)
    e = np.zeros(n)
    i_deg = np.zeros(n)
    omega_asc_deg = np.zeros(n)
    omega_deg = np.zeros(n)
    # Spread the mean anomaly so some pairs come close at different times.
    mean_anom_deg = rng.uniform(0.0, 360.0, size=n)
    return pl.DataFrame(
        {
            "number": np.arange(1, n + 1, dtype=np.int32),
            "designation": [f"({i+1})" for i in range(n)],
            "a_au": a,
            "e": e,
            "i_deg": i_deg,
            "Omega_deg": omega_asc_deg,
            "omega_deg": omega_deg,
            "M_deg": mean_anom_deg,
            "epoch_jd": np.full(n, 2457200.5),
        },
        schema_overrides={"number": pl.Int32, "designation": pl.Utf8},
    )


def _propagate_grid_inmemory(elements: pl.DataFrame, time_grid: np.ndarray) -> np.ndarray:
    """Naive in-memory propagation: yields (T, N, 3) float32 trajectory."""
    from src.propagate.kepler import propagate_df

    out = np.empty((len(time_grid), len(elements), 3), dtype=np.float32)
    for k, t in enumerate(time_grid):
        out[k] = propagate_df(elements, float(t))
    return out


def test_scan_parallel_memmap_matches_streaming(tmp_path: Path) -> None:
    """The memmap-backed scan must return the same pair set as the streaming scan."""
    n_ast = 25
    elements = _make_elements(n_ast)
    time_grid = np.linspace(2457200.5, 2457230.5, 31)  # 30 days, daily

    # 1) Streaming Kepler (positions=None).
    result_streaming = scan_parallel(
        elements,
        time_grid,
        pairs=None,
        threshold_au=0.5,  # Loose threshold so we always get some hits.
        leaf_size=10,
        n_workers=2,
        chunk_size_days=10.0,
        positions=None,
    )

    # 2) Pre-computed positions handed in directly (in-memory).
    positions_inmem = _propagate_grid_inmemory(elements, time_grid)
    result_inmem = scan_parallel(
        elements,
        time_grid,
        pairs=None,
        threshold_au=0.5,
        leaf_size=10,
        n_workers=2,
        chunk_size_days=10.0,
        positions=positions_inmem,
    )

    # 3) Same positions, but persisted to disk + memmapped back.
    npy_path = tmp_path / "traj.npy"
    mm = np.memmap(npy_path, dtype=np.float32, mode="w+", shape=positions_inmem.shape)
    mm[:] = positions_inmem
    mm.flush()
    del mm
    positions_mm = np.memmap(npy_path, dtype=np.float32, mode="r", shape=positions_inmem.shape)
    result_memmap = scan_parallel(
        elements,
        time_grid,
        pairs=None,
        threshold_au=0.5,
        leaf_size=10,
        n_workers=2,
        chunk_size_days=10.0,
        positions=positions_mm,
    )

    # The pair set should be identical across all three branches.
    pairs_streaming = {(r[0], r[1]) for r in result_streaming}
    pairs_inmem = {(r[0], r[1]) for r in result_inmem}
    pairs_memmap = {(r[0], r[1]) for r in result_memmap}

    assert pairs_streaming == pairs_inmem, "in-memory should match streaming"
    assert pairs_inmem == pairs_memmap, "memmap should match in-memory"

    # And the minimum distances should agree to float32 precision.
    by_pair_inmem = {(r[0], r[1]): r[3] for r in result_inmem}
    by_pair_memmap = {(r[0], r[1]): r[3] for r in result_memmap}
    for pair, d in by_pair_inmem.items():
        assert abs(d - by_pair_memmap[pair]) < 1e-6


def test_scan_parallel_zarr_matches_inmemory(tmp_path: Path) -> None:
    """The transposed (3,N,T) zarr-backed scan must match the in-memory scan."""
    zarr = pytest.importorskip("zarr")

    from src.propagate.cache import TrajectoryView

    n_ast = 25
    elements = _make_elements(n_ast)
    time_grid = np.linspace(2457200.5, 2457230.5, 31)
    positions_inmem = _propagate_grid_inmemory(elements, time_grid)

    result_inmem = scan_parallel(
        elements,
        time_grid,
        pairs=None,
        threshold_au=0.5,
        leaf_size=10,
        n_workers=2,
        chunk_size_days=10.0,
        positions=positions_inmem,
    )

    # Build a transposed (3, N, T) zarr directly — mirrors what cache.py
    # produces — and wrap it in TrajectoryView for the scan.
    T, N, _ = positions_inmem.shape
    zarr_path = tmp_path / "traj.zarr"
    z = zarr.open(
        str(zarr_path),
        mode="w",
        shape=(3, N, T),
        chunks=(3, N, min(16, T)),  # forces multiple chunks across T
        dtype="float32",
    )
    z[:] = np.ascontiguousarray(np.transpose(positions_inmem, (2, 1, 0)))
    positions_view = TrajectoryView(zarr.open(str(zarr_path), mode="r"))

    result_zarr = scan_parallel(
        elements,
        time_grid,
        pairs=None,
        threshold_au=0.5,
        leaf_size=10,
        n_workers=2,
        chunk_size_days=10.0,
        positions=positions_view,
    )

    pairs_inmem = {(r[0], r[1]) for r in result_inmem}
    pairs_zarr = {(r[0], r[1]) for r in result_zarr}
    assert pairs_inmem == pairs_zarr

    by_pair_inmem = {(r[0], r[1]): r[3] for r in result_inmem}
    by_pair_zarr = {(r[0], r[1]): r[3] for r in result_zarr}
    for pair, d in by_pair_inmem.items():
        assert abs(d - by_pair_zarr[pair]) < 1e-6
