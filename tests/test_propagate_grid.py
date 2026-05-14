"""Tests for src/propagate/grid.py."""

from __future__ import annotations

import numpy as np
import polars as pl
import pytest

from src.propagate.grid import make_time_grid, propagate_grid

# ---------------------------------------------------------------------------
# make_time_grid
# ---------------------------------------------------------------------------


def test_grid_length_and_step() -> None:
    grid = make_time_grid(2451545.0, 2451546.0, 6.0)
    assert len(grid) == 5  # 0, 6h, 12h, 18h, 24h
    np.testing.assert_allclose(np.diff(grid), 0.25, atol=1e-12)


def test_grid_includes_start() -> None:
    grid = make_time_grid(2451545.0, 2451546.0, 6.0)
    assert grid[0] == pytest.approx(2451545.0)


def test_grid_includes_end() -> None:
    grid = make_time_grid(2451545.0, 2451546.0, 6.0)
    assert grid[-1] == pytest.approx(2451546.0)


def test_grid_single_step() -> None:
    grid = make_time_grid(2451545.0, 2451545.0, 1.0)
    assert len(grid) == 1
    assert grid[0] == pytest.approx(2451545.0)


def test_grid_fractional_step() -> None:
    # 3 days at 36 h step → 0, 1.5, 3.0 days → 3 points
    grid = make_time_grid(2451545.0, 2451548.0, 36.0)
    assert len(grid) == 3


# ---------------------------------------------------------------------------
# propagate_grid
# ---------------------------------------------------------------------------

_ELEMENTS = pl.DataFrame(
    {
        "a_au": [2.769, 2.362],
        "e": [0.0758, 0.0889],
        "i_deg": [10.59, 7.14],
        "Omega_deg": [80.3, 103.8],
        "omega_deg": [73.6, 149.6],
        "M_deg": [95.99, 173.1],
        "epoch_jd": [2451545.0, 2451545.0],
    }
)


def test_propagate_grid_step_count() -> None:
    grid = make_time_grid(2451545.0, 2451547.0, 12.0)
    steps = list(propagate_grid(_ELEMENTS, grid))
    assert len(steps) == len(grid)


def test_propagate_grid_position_shape() -> None:
    grid = make_time_grid(2451545.0, 2451546.0, 6.0)
    for _t, pos in propagate_grid(_ELEMENTS, grid):
        assert pos.shape == (2, 3)


def test_propagate_grid_time_values_match() -> None:
    grid = make_time_grid(2451545.0, 2451546.0, 6.0)
    times = [t for t, _ in propagate_grid(_ELEMENTS, grid)]
    np.testing.assert_allclose(times, grid, atol=1e-12)


def test_propagate_grid_radii_physical() -> None:
    grid = make_time_grid(2451545.0, 2451547.0, 24.0)
    for _t, pos in propagate_grid(_ELEMENTS, grid):
        r = np.linalg.norm(pos, axis=1)
        for idx, (a, e) in enumerate([(2.769, 0.0758), (2.362, 0.0889)]):
            assert a * (1 - e) - 1e-8 <= r[idx] <= a * (1 + e) + 1e-8
