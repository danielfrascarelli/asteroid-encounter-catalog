"""Tests de src/orbdet/constants.py."""

from __future__ import annotations

import math

import pytest

from src.orbdet import constants as const


def test_gm_sun_matches_gauss_k_squared() -> None:
    assert const.GM_SUN == pytest.approx(const.GAUSS_K**2)
    # Valor clásico ≈ 2.9591220828e-4 AU^3/día^2
    assert const.GM_SUN == pytest.approx(2.959122082855911e-4, rel=1e-12)


def test_gm_sun_si_consistent_with_au_day_value() -> None:
    """GM_sun en AU^3/día^2 derivado del valor SI coincide con k^2 a ~1e-3."""
    gm_from_si = const.GM_SUN_SI * (const.DAY_S**2) / (const.AU_M**3)
    # k^2 (sistema antiguo) y GM_sun SI (IAU 2015) difieren ~1e-3 relativo.
    assert gm_from_si == pytest.approx(const.GM_SUN, rel=2e-3)


def test_speed_of_light_au_per_day() -> None:
    # Valor estándar ≈ 173.1446 AU/día
    assert const.C_AU_PER_DAY == pytest.approx(173.144632674, rel=1e-7)


def test_obliquity_value() -> None:
    assert const.OBLIQUITY_J2000_RAD == pytest.approx(math.radians(23.4392946), rel=1e-6)


def test_solar_mass_kg() -> None:
    # M_sun ≈ 1.988e30 kg
    assert const.M_SUN_KG == pytest.approx(1.988e30, rel=1e-3)


def test_gm_mass_roundtrip() -> None:
    mass = 9.4e20  # ~Ceres, kg
    gm = const.gm_from_mass_kg(mass)
    assert const.mass_kg_from_gm(gm) == pytest.approx(mass, rel=1e-12)


def test_ceres_gm_is_tiny_fraction_of_sun() -> None:
    """GM de Ceres debe ser ~4.7e-10 del GM solar (sanity de escala)."""
    gm_ceres = const.gm_from_mass_kg(9.38e20)
    assert gm_ceres / const.GM_SUN == pytest.approx(4.7e-10, rel=0.1)
