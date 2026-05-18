"""Test that refine_candidates uses the trajectory cache when supplied.

When ``positions`` and ``time_grid`` are passed, the refiner must read
distances from the N-body cache rather than re-propagating with Kepler.
This prevents the bug where rebound runs were silently downgraded to
2-body distances by the refinement step.
"""

from __future__ import annotations

import numpy as np
import polars as pl

from src.detect.refine import refine_candidates


def _build_elements() -> pl.DataFrame:
    """Two synthetic asteroids — orbital elements aren't actually used in the
    cache path, but the schema must match."""
    return pl.DataFrame(
        {
            "number": np.array([10, 4803], dtype=np.int32),
            "designation": ["(10) Hygiea", "(4803) Birkle"],
            "a_au": [3.142, 2.903],
            "e": [0.115, 0.037],
            "i_deg": [3.84, 2.92],
            "Omega_deg": [283.4, 46.2],
            "omega_deg": [312.1, 182.5],
            "M_deg": [264.5, 260.0],
            "epoch_jd": [2457200.5, 2457200.5],
        },
        schema_overrides={"number": pl.Int32, "designation": pl.Utf8},
    )


def test_refine_uses_cache_when_provided() -> None:
    """The cache-based refinement must report the cache's distance, not Kepler's."""
    elements = _build_elements()
    # Synthetic time grid: 24-hour spacing × 5 steps.
    time_grid = np.array([0.0, 1.0, 2.0, 3.0, 4.0])

    # Build a (T, N, 3) cache where the minimum distance occurs at t = 2.0
    # and equals exactly 0.011924 AU (the JPL value for Hygiea-Birkle).
    positions = np.zeros((5, 2, 3), dtype=np.float32)
    positions[:, 0, 0] = [10.0, 5.0, 0.0, 5.0, 10.0]  # asteroid 0 along x
    # Asteroid 1 sits at (0, 0.011924, 0); distance(t=2) = 0.011924 exactly.
    positions[:, 1, 0] = 0.0
    positions[:, 1, 1] = 0.011924
    positions[:, 1, 2] = 0.0

    # Coarse candidate around t=2 (the true min)
    candidates = [(0, 1, 2.0, 5.0)]

    out = refine_candidates(
        elements,
        candidates,
        threshold_au=0.05,
        positions=positions,
        time_grid=time_grid,
    )

    assert len(out) == 1
    row = out.row(0, named=True)
    # The refined distance must come from the cache (0.011924), not from the
    # Kepler propagation of the synthetic elements (which would give a totally
    # different number because the elements were chosen at random).
    assert abs(row["dist_au"] - 0.011924) < 1e-5
    assert abs(row["jd_tdb"] - 2.0) < 1.5  # parabola fit is near the cache min


def test_refine_kepler_path_unaffected() -> None:
    """Without positions/time_grid, the refiner must still use Kepler analytics."""
    elements = _build_elements()
    # A candidate at the (synthetic) epoch — the Kepler path will integrate.
    candidates = [(0, 1, 2457201.0, 0.04)]

    out = refine_candidates(
        elements,
        candidates,
        threshold_au=0.05,
        fine_step_seconds=60.0,
        window_hours=2.0,
    )
    # We don't assert exact distance — just that the Kepler path still runs
    # without errors and either keeps or drops the candidate.
    assert len(out) <= 1
