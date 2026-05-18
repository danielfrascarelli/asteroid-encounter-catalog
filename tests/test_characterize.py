"""Tests for Phase 5 characterization modules."""

from __future__ import annotations

import numpy as np
import pytest

from src.characterize.geometry import (
    dist_au_to_km,
    vel_au_per_day_to_km_s,
    vel_au_per_day_to_m_s,
)
from src.characterize.observability import (
    apparent_mag_hg,
    is_gaia_observable,
    solar_elongation_deg,
)
from src.characterize.physical import classify_orbit, diameter_km

# ------------------------------------------------------------------ #
# physical.py                                                          #
# ------------------------------------------------------------------ #


class TestDiameterKm:
    def test_ceres(self) -> None:
        # (1) Ceres: H=3.34, actual albedo=0.09 → ~940 km
        d = float(diameter_km(3.34, albedo=0.09))
        assert abs(d - 940.0) / 940.0 < 0.05, f"Ceres diameter {d:.1f} km not within 5% of 940"

    def test_vesta(self) -> None:
        # (4) Vesta: H=3.25, actual albedo=0.34 → ~525 km
        d = float(diameter_km(3.25, albedo=0.34))
        assert abs(d - 525.0) / 525.0 < 0.05, f"Vesta diameter {d:.1f} km not within 5% of 525"

    def test_array_input(self) -> None:
        h_vals = np.array([3.34, 3.25, 10.0])
        p = np.array([0.09, 0.34, 0.14])
        d = diameter_km(h_vals, p)
        assert d.shape == (3,)
        assert np.all(d > 0)

    def test_brighter_is_larger(self) -> None:
        # Smaller H = brighter = larger object
        assert diameter_km(5.0) > diameter_km(10.0)

    def test_higher_albedo_is_smaller(self) -> None:
        # Higher albedo → smaller object for same H
        assert diameter_km(10.0, albedo=0.05) > diameter_km(10.0, albedo=0.30)

    def test_nan_h_returns_nan(self) -> None:
        assert np.isnan(float(diameter_km(np.nan)))

    def test_nan_array_propagates(self) -> None:
        h = np.array([3.34, np.nan, 10.0])
        result = diameter_km(h)
        assert not np.isnan(result[0])
        assert np.isnan(result[1])
        assert not np.isnan(result[2])


class TestClassifyOrbit:
    def test_mba(self) -> None:
        assert classify_orbit(2.5, 0.1) == "MBA"

    def test_nea(self) -> None:
        # q = 1.0*(1-0.7) = 0.3 < 1.3
        assert classify_orbit(1.0, 0.7) == "NEA"

    def test_trojan(self) -> None:
        assert classify_orbit(5.2, 0.05) == "Trojan"

    def test_centaur(self) -> None:
        assert classify_orbit(15.0, 0.2) == "Centaur"

    def test_tno(self) -> None:
        assert classify_orbit(45.0, 0.05) == "TNO"

    def test_array_input(self) -> None:
        a = np.array([2.5, 1.0, 5.2, 15.0, 45.0])
        e = np.array([0.1, 0.7, 0.05, 0.2, 0.05])
        cls = classify_orbit(a, e)
        assert cls[0] == "MBA"
        assert cls[1] == "NEA"
        assert cls[2] == "Trojan"
        assert cls[3] == "Centaur"
        assert cls[4] == "TNO"


# ------------------------------------------------------------------ #
# geometry.py                                                          #
# ------------------------------------------------------------------ #


class TestGeometry:
    def test_au_to_km(self) -> None:
        assert abs(float(dist_au_to_km(1.0)) - 149_597_870.7) < 1.0

    def test_earth_orbital_speed(self) -> None:
        # Earth orbits at ~0.01720 AU/day ≈ 29.78 km/s
        v_km = float(vel_au_per_day_to_km_s(0.01720))
        assert abs(v_km - 29.78) < 0.5

    def test_m_per_s_vs_km_s(self) -> None:
        v = 0.005
        assert abs(float(vel_au_per_day_to_m_s(v)) - float(vel_au_per_day_to_km_s(v)) * 1000) < 1e-6

    def test_array_input(self) -> None:
        d = dist_au_to_km(np.array([1.0, 2.0]))
        assert d[1] == pytest.approx(2 * d[0])


# ------------------------------------------------------------------ #
# observability.py                                                     #
# ------------------------------------------------------------------ #


class TestSolarElongation:
    def test_opposition(self) -> None:
        # Asteroid directly opposite the Sun as seen from Earth
        # Earth at (1, 0, 0), asteroid at (2, 0, 0): elongation = 180°
        earth = np.array([[1.0, 0.0, 0.0]])
        enc = np.array([[2.0, 0.0, 0.0]])
        elong = solar_elongation_deg(enc, earth)
        assert abs(float(elong[0]) - 180.0) < 1.0

    def test_conjunction(self) -> None:
        # Asteroid at (-2, 0, 0), Earth at (1, 0, 0): elongation ≈ 0°
        earth = np.array([[1.0, 0.0, 0.0]])
        enc = np.array([[-2.0, 0.0, 0.0]])
        elong = solar_elongation_deg(enc, earth)
        assert float(elong[0]) < 5.0

    def test_quadrature(self) -> None:
        # For 90° elongation: Earth at (1,0,0), Sun at origin → Sun dir from Earth = (-1,0,0).
        # Perpendicular direction from Earth is (0,1,0), so asteroid at (1,2,0) gives elong=90°.
        earth = np.array([[1.0, 0.0, 0.0]])
        enc = np.array([[1.0, 2.0, 0.0]])
        elong = solar_elongation_deg(enc, earth)
        assert abs(float(elong[0]) - 90.0) < 1.0

    def test_array_input(self) -> None:
        earth = np.array([[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        enc = np.array([[2.0, 0.0, 0.0], [-2.0, 0.0, 0.0]])
        elong = solar_elongation_deg(enc, earth)
        assert len(elong) == 2


class TestApparentMag:
    def test_increases_with_distance(self) -> None:
        # Farther away → fainter (larger magnitude)
        m1 = float(apparent_mag_hg(10.0, 2.0, 1.0))
        m2 = float(apparent_mag_hg(10.0, 3.0, 2.0))
        assert m2 > m1

    def test_reasonable_mba(self) -> None:
        # MBA at 2.5 AU from Sun, 1.5 AU from Earth, H=12 → roughly 17–19 mag
        m = float(apparent_mag_hg(12.0, 2.5, 1.5))
        assert 14.0 < m < 22.0


class TestGaiaObservable:
    def test_observable(self) -> None:
        elong = np.array([90.0])
        mag = np.array([18.0])
        assert is_gaia_observable(elong, mag)[0]

    def test_too_close_to_sun(self) -> None:
        elong = np.array([30.0])
        mag = np.array([18.0])
        assert not is_gaia_observable(elong, mag)[0]

    def test_too_faint(self) -> None:
        elong = np.array([90.0])
        mag = np.array([22.0])
        assert not is_gaia_observable(elong, mag)[0]
