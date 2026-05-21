"""Unit tests for src/astrometry/transforms.py."""

from __future__ import annotations

import numpy as np
import pytest

from src.astrometry.transforms import (
    ecliptic_to_equatorial,
    equatorial_to_ecliptic,
    heliocentric_to_barycentric_icrs,
    light_time_iterate,
    radec_to_unit_vec,
    stellar_aberration,
    sun_barycentric_au,
    xyz_to_radec,
)

# ---------------------------------------------------------------------------
# Frame rotation tests
# ---------------------------------------------------------------------------


def test_ecliptic_x_axis_unchanged() -> None:
    """The X axis is the rotation axis — should be invariant."""
    v = np.array([1.0, 0.0, 0.0])
    assert np.allclose(ecliptic_to_equatorial(v), v)


def test_ecliptic_y_axis_rotates_correctly() -> None:
    """+Y_ecl should map to (0, cos ε, sin ε)."""
    v = np.array([0.0, 1.0, 0.0])
    out = ecliptic_to_equatorial(v)
    eps = np.radians(23.43928083)
    expected = np.array([0.0, np.cos(eps), np.sin(eps)])
    assert np.allclose(out, expected)


def test_rotation_inverse() -> None:
    """ecliptic→equatorial→ecliptic should be identity."""
    rng = np.random.default_rng(42)
    v = rng.standard_normal((10, 3))
    back = equatorial_to_ecliptic(ecliptic_to_equatorial(v))
    assert np.allclose(back, v, atol=1e-14)


# ---------------------------------------------------------------------------
# Sun barycentric position
# ---------------------------------------------------------------------------


def test_sun_barycentric_small() -> None:
    """The Sun is displaced from the barycenter by ≲ 0.01 AU (Jupiter)."""
    jd = 2457000.5  # 2014-12-09 UTC ish
    sun = sun_barycentric_au(jd)
    assert sun.shape == (3,)
    assert np.linalg.norm(sun) < 0.01


def test_sun_barycentric_vectorised() -> None:
    """Vectorised input returns matching shape."""
    jds = np.array([2457000.5, 2457500.5, 2458000.5])
    sun = sun_barycentric_au(jds)
    assert sun.shape == (3, 3)


# ---------------------------------------------------------------------------
# Heliocentric → barycentric
# ---------------------------------------------------------------------------


def test_helio_to_bary_close_to_helio_eq() -> None:
    """Barycentric position should differ from rotated heliocentric by ≲ 0.01 AU."""
    pos_ecl = np.array([2.5, 1.0, 0.1])
    jd = 2457000.5
    pos_bary = heliocentric_to_barycentric_icrs(pos_ecl, jd)
    pos_eq = ecliptic_to_equatorial(pos_ecl)
    assert np.linalg.norm(pos_bary - pos_eq) < 0.01


# ---------------------------------------------------------------------------
# Cartesian ↔ spherical
# ---------------------------------------------------------------------------


def test_xyz_to_radec_roundtrip() -> None:
    rng = np.random.default_rng(1)
    n = 10
    ra = rng.uniform(0, 360, n)
    dec = rng.uniform(-89, 89, n)
    vec = radec_to_unit_vec(ra, dec)
    ra2, dec2 = xyz_to_radec(vec)
    assert np.allclose(ra2, ra, atol=1e-10)
    assert np.allclose(dec2, dec, atol=1e-10)


def test_radec_unit_vec_norm() -> None:
    vec = radec_to_unit_vec(np.array([45.0, 180.0]), np.array([30.0, -45.0]))
    norms = np.linalg.norm(vec, axis=-1)
    assert np.allclose(norms, 1.0)


# ---------------------------------------------------------------------------
# Light-time correction
# ---------------------------------------------------------------------------


def test_light_time_constant_target() -> None:
    """A stationary target at known distance gives expected light-time."""
    distance_au = 2.0
    target = np.array([distance_au, 0.0, 0.0])

    def target_pos(jd: float) -> np.ndarray:  # noqa: ARG001
        return target

    gaia = np.array([0.0, 0.0, 0.0])
    pos, tau = light_time_iterate(target_pos, jd_tdb_obs=2457000.5, gaia_xyz_bary=gaia)
    # τ = d / c, with c = 173.144 AU/day
    expected_tau = distance_au / 173.144632674
    assert tau == pytest.approx(expected_tau, rel=1e-6)
    assert np.allclose(pos, target)


def test_light_time_moving_target() -> None:
    """A target with linear motion: retarded position should be earlier."""
    velocity = np.array([0.0, 0.01, 0.0])  # 0.01 AU/day

    def target_pos(jd: float) -> np.ndarray:
        return np.array([2.0, 0.0, 0.0]) + velocity * (jd - 2457000.5)

    gaia = np.array([0.0, 0.0, 0.0])
    pos, tau = light_time_iterate(target_pos, jd_tdb_obs=2457000.5, gaia_xyz_bary=gaia)
    # Sanity: tau ≈ 2/c ≈ 0.012 days; position retarded by velocity*tau
    expected_tau = 2.0 / 173.144632674
    assert tau == pytest.approx(expected_tau, rel=1e-3)
    # The retarded position is slightly earlier than at t_obs
    assert pos[1] < 0.0
    assert abs(pos[1] - (-velocity[1] * tau)) < 1e-6


# ---------------------------------------------------------------------------
# Stellar aberration
# ---------------------------------------------------------------------------


def test_stellar_aberration_zero_velocity() -> None:
    """No velocity → no aberration."""
    los = np.array([1.0, 0.0, 0.0])
    v = np.zeros(3)
    out = stellar_aberration(los, v)
    assert np.allclose(out, los)


def test_stellar_aberration_magnitude() -> None:
    """Earth velocity ≈ 30 km/s → annual aberration ~20 arcsec."""
    # Velocity in AU/day: 30 km/s = 0.01731 AU/day (perpendicular to LOS)
    v = np.array([0.0, 0.01731, 0.0])
    los = np.array([1.0, 0.0, 0.0])
    apparent = stellar_aberration(los, v)
    # Aberration angle ≈ |v|/c ≈ 1e-4 rad ≈ 20 arcsec
    angle = np.arccos(np.dot(apparent, los) / np.linalg.norm(apparent))
    arcsec = np.degrees(angle) * 3600.0
    assert 15.0 < arcsec < 25.0
