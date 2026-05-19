"""Tests for src/propagate/nbody_perturber.py."""

from __future__ import annotations

import numpy as np
import pytest

from src.propagate.nbody_perturber import propagate_target_with_perturber


def _make_elements(
    a_au: float = 2.5,
    e: float = 0.1,
    i_deg: float = 5.0,
    Omega_deg: float = 30.0,  # noqa: N803
    omega_deg: float = 60.0,
    M_deg: float = 100.0,  # noqa: N803
    epoch_jd: float = 2456200.5,
) -> dict:
    return {
        "a_au": a_au,
        "e": e,
        "i_deg": i_deg,
        "Omega_deg": Omega_deg,
        "omega_deg": omega_deg,
        "M_deg": M_deg,
        "epoch_jd": epoch_jd,
    }


def test_propagate_shape() -> None:
    """Output shape matches the time grid."""
    target = _make_elements()
    perturber = _make_elements(a_au=3.0)
    grid = np.array([2456200.5, 2456250.5, 2456300.5, 2456400.5])
    out = propagate_target_with_perturber(
        target_elements=target,
        perturber_elements=perturber,
        perturber_mass_kg=0.0,
        time_grid_jd_tdb=grid,
    )
    assert out.shape == (4, 3)


def test_zero_mass_perturber_matches_no_perturber() -> None:
    """Two runs with mass=0 give same trajectory regardless of perturber elements."""
    target = _make_elements()
    pert_a = _make_elements(a_au=3.0)
    pert_b = _make_elements(a_au=4.0, Omega_deg=180.0)
    grid = np.array([2456200.5 + d for d in (1.0, 30.0, 100.0)])
    out_a = propagate_target_with_perturber(
        target, pert_a, 0.0, grid, include_planets=("sun",)
    )
    out_b = propagate_target_with_perturber(
        target, pert_b, 0.0, grid, include_planets=("sun",)
    )
    # Massless perturbers contribute nothing → identical target trajectories
    assert np.allclose(out_a, out_b, atol=1e-10)


def test_mass_changes_trajectory() -> None:
    """Non-zero perturber mass produces a different trajectory."""
    target = _make_elements(a_au=2.5, M_deg=0.0)
    # Perturber on similar orbit so they interact (close pass over the window)
    perturber = _make_elements(a_au=2.5, M_deg=2.0)  # small offset
    grid = np.array([2456200.5 + d for d in (1.0, 100.0, 300.0)])
    out_zero = propagate_target_with_perturber(
        target, perturber, 0.0, grid, include_planets=("sun",)
    )
    # Massive perturber: ~Ceres mass (9e20 kg)
    out_mass = propagate_target_with_perturber(
        target, perturber, 9.4e20, grid, include_planets=("sun",)
    )
    # Trajectories should differ
    diff = np.linalg.norm(out_mass - out_zero, axis=1).max()
    assert diff > 1e-9, f"Massless and massive runs identical (diff={diff:.2e} AU)"


def test_epoch_mismatch_raises() -> None:
    target = _make_elements(epoch_jd=2456200.5)
    perturber = _make_elements(epoch_jd=2456500.5)
    grid = np.array([2456300.5])
    with pytest.raises(ValueError, match="osculating epoch"):
        propagate_target_with_perturber(target, perturber, 0.0, grid)


def test_negative_mass_raises() -> None:
    target = _make_elements()
    perturber = _make_elements()
    grid = np.array([2456200.5])
    with pytest.raises(ValueError, match="non-negative"):
        propagate_target_with_perturber(target, perturber, -1.0, grid)
