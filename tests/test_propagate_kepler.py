"""Tests for src/propagate/kepler.py."""

from __future__ import annotations

import numpy as np
import polars as pl

from src.propagate.kepler import _DEG, _K, kepler_to_cartesian, propagate_df, solve_kepler

# ---------------------------------------------------------------------------
# solve_kepler
# ---------------------------------------------------------------------------


def test_solve_kepler_m_zero_any_e() -> None:
    result = solve_kepler(0.0, 0.5)
    assert abs(result[0]) < 1e-12


def test_solve_kepler_m_pi_any_e() -> None:
    result = solve_kepler(np.pi, 0.3)
    np.testing.assert_allclose(result, [np.pi], atol=1e-10)


def test_solve_kepler_residual_low_e() -> None:
    e = 0.2
    m_vals = np.linspace(0.0, 2.0 * np.pi, 300)
    ecc_vals = solve_kepler(m_vals, e)
    np.testing.assert_allclose(m_vals - (ecc_vals - e * np.sin(ecc_vals)), 0.0, atol=1e-10)


def test_solve_kepler_residual_high_e() -> None:
    e = 0.85
    m_vals = np.linspace(0.05, 2.0 * np.pi - 0.05, 100)
    ecc_vals = solve_kepler(m_vals, e)
    np.testing.assert_allclose(m_vals - (ecc_vals - e * np.sin(ecc_vals)), 0.0, atol=1e-8)


def test_solve_kepler_vectorized_e() -> None:
    m_vals = np.array([np.pi / 4, np.pi / 2, np.pi])
    e_vals = np.array([0.1, 0.3, 0.5])
    ecc_vals = solve_kepler(m_vals, e_vals)
    np.testing.assert_allclose(
        m_vals - (ecc_vals - e_vals * np.sin(ecc_vals)), 0.0, atol=1e-10
    )


# ---------------------------------------------------------------------------
# kepler_to_cartesian — scalar / simple cases
# ---------------------------------------------------------------------------


def test_circular_equatorial_at_epoch() -> None:
    # e=0, i=0, Ω=0, ω=0, M0=0, t=epoch  → (a, 0, 0)
    pos = kepler_to_cartesian(2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2451545.0, 2451545.0)
    np.testing.assert_allclose(pos, [[2.0, 0.0, 0.0]], atol=1e-12)


def test_circular_orbit_quarter_period() -> None:
    # After T/4, M = π/2, argument of latitude = π/2 → (0, a, 0)
    a = 2.0
    n = _K / a**1.5
    period = 2.0 * np.pi / n
    t0 = 2451545.0
    pos = kepler_to_cartesian(a, 0.0, 0.0, 0.0, 0.0, 0.0, t0, t0 + period / 4.0)
    np.testing.assert_allclose(pos[0, 0], 0.0, atol=1e-10)
    np.testing.assert_allclose(pos[0, 1], a, atol=1e-10)
    np.testing.assert_allclose(pos[0, 2], 0.0, atol=1e-12)


def test_period_conservation() -> None:
    # After exactly one orbital period the position must repeat.
    a, e, inc = 2.769, 0.0758, 10.59 * _DEG
    node, peri, m0 = 80.3 * _DEG, 73.6 * _DEG, 95.99 * _DEG
    t0 = 2451545.0
    period = 2.0 * np.pi / (_K / a**1.5)
    pos0 = kepler_to_cartesian(a, e, inc, node, peri, m0, t0, t0)
    pos1 = kepler_to_cartesian(a, e, inc, node, peri, m0, t0, t0 + period)
    np.testing.assert_allclose(pos0, pos1, atol=1e-8)


def test_radius_within_periapsis_apoapsis() -> None:
    a, e = 2.769, 0.0758
    t0 = 2451545.0
    for dt in (0, 100, 300, 600, 900):
        pos = kepler_to_cartesian(a, e, 0.3, 0.8, 1.2, 0.5, t0, t0 + dt)
        r = float(np.linalg.norm(pos))
        assert a * (1 - e) - 1e-10 <= r <= a * (1 + e) + 1e-10


def test_inclination_produces_nonzero_z() -> None:
    # i=90°, ω=0, Ω=0, M0=π/2 → argument of latitude = π/2 → |z| ≈ a
    a = 2.0
    pos = kepler_to_cartesian(
        a, 0.0, np.pi / 2.0, 0.0, 0.0, np.pi / 2.0, 2451545.0, 2451545.0
    )
    np.testing.assert_allclose(abs(pos[0, 2]), a, atol=1e-10)


def test_equatorial_orbit_z_is_zero() -> None:
    pos = kepler_to_cartesian(2.5, 0.1, 0.0, 1.0, 0.5, 0.3, 2451545.0, 2451600.0)
    np.testing.assert_allclose(pos[0, 2], 0.0, atol=1e-12)


# ---------------------------------------------------------------------------
# kepler_to_cartesian — vectorized / many asteroids
# ---------------------------------------------------------------------------


def test_vectorized_many_radii_within_bounds() -> None:
    rng = np.random.default_rng(42)
    n_bodies = 200
    a = rng.uniform(1.5, 3.5, n_bodies)
    e = rng.uniform(0.0, 0.3, n_bodies)
    pos = kepler_to_cartesian(
        a, e,
        rng.uniform(0, 0.5, n_bodies),
        rng.uniform(0, 2 * np.pi, n_bodies),
        rng.uniform(0, 2 * np.pi, n_bodies),
        rng.uniform(0, 2 * np.pi, n_bodies),
        2451545.0,
        2451600.0,
    )
    assert pos.shape == (n_bodies, 3)
    r = np.linalg.norm(pos, axis=1)
    np.testing.assert_array_less(a * (1 - e) - 1e-8, r)
    np.testing.assert_array_less(r, a * (1 + e) + 1e-8)


# ---------------------------------------------------------------------------
# propagate_df
# ---------------------------------------------------------------------------

_CERES_ELEMENTS = {
    "a_au": [2.7691652],
    "e": [0.0758458],
    "i_deg": [10.5935],
    "Omega_deg": [80.3099],
    "omega_deg": [73.5975],
    "M_deg": [95.9892],
    "epoch_jd": [2454200.5],
}

_THREE_ELEMENTS = {
    "a_au": [2.769, 2.362, 2.773],
    "e": [0.0758, 0.0889, 0.2316],
    "i_deg": [10.59, 7.14, 34.84],
    "Omega_deg": [80.3, 103.8, 173.0],
    "omega_deg": [73.6, 149.6, 310.1],
    "M_deg": [95.99, 173.1, 232.6],
    "epoch_jd": [2451545.0, 2451545.0, 2451545.0],
}


def test_propagate_df_output_shape() -> None:
    df = pl.DataFrame(_THREE_ELEMENTS)
    pos = propagate_df(df, 2451600.0)
    assert pos.shape == (3, 3)


def test_propagate_df_radius_physical() -> None:
    df = pl.DataFrame(_CERES_ELEMENTS)
    pos = propagate_df(df, 2456863.5)  # within Gaia observation window
    r = float(np.linalg.norm(pos[0]))
    a, e = 2.7691652, 0.0758458
    assert a * (1 - e) <= r <= a * (1 + e) + 1e-8


def test_propagate_df_matches_scalar_call() -> None:
    df = pl.DataFrame(_CERES_ELEMENTS)
    t = 2456863.5
    pos_df = propagate_df(df, t)
    pos_scalar = kepler_to_cartesian(
        2.7691652, 0.0758458,
        10.5935 * _DEG, 80.3099 * _DEG, 73.5975 * _DEG, 95.9892 * _DEG,
        2454200.5, t,
    )
    np.testing.assert_allclose(pos_df[0], pos_scalar[0], atol=1e-12)
