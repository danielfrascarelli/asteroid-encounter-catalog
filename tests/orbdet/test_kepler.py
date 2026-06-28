"""Tests de src/orbdet/kepler.py — invariantes físicos y round-trips."""

from __future__ import annotations

import math

import numpy as np
import pytest

from src.orbdet.constants import GM_SUN
from src.orbdet.kepler import (
    KeplerElements,
    dstate_delements,
    elements_to_state,
    mean_motion,
    period,
    propagate,
    solve_kepler,
    state_to_elements,
)

# Órbita genérica inclinada y excéntrica (evita singularidades e≈0, i≈0).
_EL = KeplerElements(
    a=2.7,  # AU, cinturón principal
    e=0.15,
    i=math.radians(10.0),
    Omega=math.radians(80.0),
    omega=math.radians(60.0),
    M=math.radians(45.0),
)


# --- Ecuación de Kepler -----------------------------------------------------


@pytest.mark.parametrize("e", [0.0, 0.1, 0.5, 0.9, 0.99])
@pytest.mark.parametrize("M_deg", [0.0, 30.0, 90.0, 180.0, 270.0, 359.0])
def test_solve_kepler_residual(e: float, M_deg: float) -> None:
    M = math.radians(M_deg)
    E = solve_kepler(M, e)
    # Verifica M = E - e·sinE (módulo 2π)
    residual = (E - e * math.sin(E)) - math.fmod(M, 2.0 * math.pi)
    assert abs(residual) < 1e-11


def test_solve_kepler_rejects_hyperbolic() -> None:
    with pytest.raises(ValueError):
        solve_kepler(1.0, 1.0)


# --- Round-trip elementos ↔ estado ------------------------------------------


def test_elements_state_roundtrip() -> None:
    r, v = elements_to_state(_EL)
    el2 = state_to_elements(r, v)
    assert el2.a == pytest.approx(_EL.a, rel=1e-11)
    assert el2.e == pytest.approx(_EL.e, rel=1e-10)
    assert el2.i == pytest.approx(_EL.i, rel=1e-10)
    assert el2.Omega == pytest.approx(_EL.Omega, rel=1e-10)
    assert el2.omega == pytest.approx(_EL.omega, rel=1e-10)
    assert el2.M == pytest.approx(_EL.M, rel=1e-10)


# --- Invariantes físicos ----------------------------------------------------


def test_specific_energy_matches_semimajor_axis() -> None:
    r, v = elements_to_state(_EL)
    rmag = float(np.linalg.norm(r))
    vmag = float(np.linalg.norm(v))
    energy = vmag**2 / 2.0 - GM_SUN / rmag
    assert energy == pytest.approx(-GM_SUN / (2.0 * _EL.a), rel=1e-12)


def test_angular_momentum_matches_elements() -> None:
    r, v = elements_to_state(_EL)
    h = np.cross(r, v)
    hmag = float(np.linalg.norm(h))
    expected = math.sqrt(GM_SUN * _EL.a * (1.0 - _EL.e**2))
    assert hmag == pytest.approx(expected, rel=1e-12)
    # Inclinación: cos i = h_z / |h|
    assert math.acos(h[2] / hmag) == pytest.approx(_EL.i, rel=1e-11)


def test_perihelion_distance() -> None:
    # En M=0 (perihelio) r = a(1-e).
    el_peri = KeplerElements(a=_EL.a, e=_EL.e, i=_EL.i, Omega=_EL.Omega, omega=_EL.omega, M=0.0)
    r, _ = elements_to_state(el_peri)
    assert float(np.linalg.norm(r)) == pytest.approx(_EL.a * (1.0 - _EL.e), rel=1e-12)


# --- Propagación ------------------------------------------------------------


def test_propagate_full_period_returns_to_start() -> None:
    P = period(_EL.a)
    r0, v0 = elements_to_state(_EL)
    el_T = propagate(_EL, P)
    rT, vT = elements_to_state(el_T)
    assert np.allclose(rT, r0, atol=1e-9)
    assert np.allclose(vT, v0, atol=1e-9)


def test_propagate_conserves_orbit_shape() -> None:
    el2 = propagate(_EL, 123.456)
    # a, e, i, Ω, ω invariantes en dos cuerpos; solo cambia M.
    assert el2.a == pytest.approx(_EL.a)
    assert el2.e == pytest.approx(_EL.e)
    assert el2.i == pytest.approx(_EL.i)
    assert el2.Omega == pytest.approx(_EL.Omega)
    assert el2.omega == pytest.approx(_EL.omega)
    assert el2.M != pytest.approx(_EL.M)


def test_mean_motion_period_consistency() -> None:
    assert mean_motion(_EL.a) * period(_EL.a) == pytest.approx(2.0 * math.pi, rel=1e-13)


def test_period_of_one_au_is_one_year() -> None:
    # a=1 AU → P ≈ 365.25 días (k define el año gaussiano ≈ 365.2568 d).
    assert period(1.0) == pytest.approx(365.2568983, rel=1e-6)


# --- Jacobiano analítico ∂[r,v]/∂elementos ----------------------------------


@pytest.mark.parametrize(
    "el",
    [
        _EL,
        KeplerElements(a=1.3, e=0.4, i=math.radians(25.0), Omega=1.1, omega=2.3, M=0.7),
        KeplerElements(a=3.2, e=0.05, i=math.radians(3.0), Omega=5.0, omega=0.2, M=3.0),
    ],
)
def test_dstate_delements_matches_finite_difference(el: KeplerElements) -> None:
    """El Jacobiano analítico coincide con la diferencia finita central del mapa
    elementos→estado (mapa estático kepleriano, sin dinámica)."""
    jac = dstate_delements(el)
    base = np.array(el.as_array(), dtype=float)
    # Pasos absolutos por elemento: relativo para a, absolutos chicos para el resto.
    steps = np.array([1e-7 * el.a, 1e-7, 1e-7, 1e-7, 1e-7, 1e-7])
    jac_fd = np.zeros((6, 6))
    for j in range(6):
        plus, minus = base.copy(), base.copy()
        plus[j] += steps[j]
        minus[j] -= steps[j]
        rp, vp = elements_to_state(KeplerElements(*plus))
        rm, vm = elements_to_state(KeplerElements(*minus))
        jac_fd[0:3, j] = (rp - rm) / (2.0 * steps[j])
        jac_fd[3:6, j] = (vp - vm) / (2.0 * steps[j])
    scale = np.maximum(np.abs(jac), np.abs(jac_fd))
    rel = np.abs(jac - jac_fd) / np.where(scale > 0.0, scale, 1.0)
    assert rel.max() < 1e-6, f"max rel error {rel.max():.2e}\n{rel}"
