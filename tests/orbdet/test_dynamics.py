"""Tests de src/orbdet/dynamics.py.

- Correctness sin red: el modelo N-cuerpos restringido al Sol debe reproducir
  el propagador kepleriano analítico a precisión del integrador (prueba
  integrador + frame + GM compartido).
- Sanity de efemérides planetarias (astropy builtin, sin red).
- Validación contra JPL Horizons: marcada ``horizons`` (deseleccionada en CI;
  requiere red).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from src.orbdet.dynamics import (
    PLANET_MASSES_MSUN,
    planet_state_ecliptic,
    propagate,
)
from src.orbdet.kepler import KeplerElements, elements_to_state
from src.orbdet.kepler import propagate as kepler_propagate

_EPOCH = 2_457_000.5  # JD TDB, ~2014-12 (era Gaia)
_EL = KeplerElements(
    a=2.7,
    e=0.15,
    i=math.radians(10.0),
    Omega=math.radians(80.0),
    omega=math.radians(60.0),
    M=math.radians(45.0),
)


# --- Efemérides planetarias (sin red, ephemeris builtin) --------------------


def test_planet_state_distances_sane() -> None:
    sun_p, _ = planet_state_ecliptic("sun", _EPOCH)
    earth_p, earth_v = planet_state_ecliptic("earth", _EPOCH)
    jup_p, _ = planet_state_ecliptic("jupiter", _EPOCH)
    # Tierra ~1 AU del Sol, Júpiter ~5 AU.
    assert np.linalg.norm(earth_p - sun_p) == pytest.approx(1.0, abs=0.05)
    assert 4.8 < np.linalg.norm(jup_p - sun_p) < 5.5
    # Velocidad orbital de la Tierra ~2π AU/año ≈ 0.0172 AU/día.
    assert np.linalg.norm(earth_v) == pytest.approx(0.0172, abs=0.003)


def test_planet_masses_present() -> None:
    assert PLANET_MASSES_MSUN["sun"] == 1.0
    assert PLANET_MASSES_MSUN["jupiter"] == pytest.approx(9.5479e-4, rel=1e-3)


# --- Consistencia con dos cuerpos (rebound, sin red) ------------------------


@pytest.mark.slow
def test_two_body_limit_matches_kepler() -> None:
    """Con solo el Sol, el N-cuerpos = Kepler analítico.

    En este caso el Sol no siente fuerzas (la partícula de prueba es no masiva),
    así que se mueve en línea recta: pos_sun(dt) = sun_p0 + sun_v0·dt. La
    posición baricéntrica esperada es la heliocéntrica kepleriana + el Sol.
    """
    sun_p0, sun_v0 = planet_state_ecliptic("sun", _EPOCH)
    dts = np.array([-200.0, -50.0, 0.0, 100.0, 365.0])
    out_epochs = _EPOCH + dts

    got = propagate(_EL, _EPOCH, out_epochs, perturbers=("sun",), integrator="ias15")

    for i, dt in enumerate(dts):
        helio = elements_to_state(kepler_propagate(_EL, float(dt)))[0]
        expected = helio + sun_p0 + sun_v0 * dt
        assert np.allclose(got[i], expected, atol=1e-8), f"dt={dt}: {got[i]} vs {expected}"


@pytest.mark.slow
def test_dt_zero_returns_initial_state() -> None:
    sun_p0, _ = planet_state_ecliptic("sun", _EPOCH)
    helio0 = elements_to_state(_EL)[0]
    got = propagate(_EL, _EPOCH, np.array([_EPOCH]), perturbers=("sun", "jupiter"))
    assert np.allclose(got[0], helio0 + sun_p0, atol=1e-9)


@pytest.mark.slow
def test_planets_perturb_the_orbit() -> None:
    """Agregar planetas cambia la posición a un año vista (señal no nula y chica)."""
    out = np.array([_EPOCH + 365.0])
    two_body = propagate(_EL, _EPOCH, out, perturbers=("sun",))[0]
    with_planets = propagate(_EL, _EPOCH, out, perturbers=("sun", "jupiter", "saturn"))[0]
    diff = float(np.linalg.norm(with_planets - two_body))
    # No nulo (los planetas importan) pero pequeño (< 0.01 AU en 1 año).
    assert 1e-7 < diff < 1e-2


@pytest.mark.slow
def test_velocity_return_shape_and_speed() -> None:
    out = _EPOCH + np.array([10.0, 50.0])
    pos, vel = propagate(_EL, _EPOCH, out, perturbers=("sun",), return_velocity=True)
    assert pos.shape == (2, 3)
    assert vel.shape == (2, 3)
    # Velocidad baricéntrica de un MBA: ~0.01 AU/día (orden de magnitud).
    assert 0.005 < float(np.linalg.norm(vel[0])) < 0.03


# --- Validación contra JPL Horizons (red; deseleccionada en CI) -------------


@pytest.mark.horizons
def test_matches_horizons_for_ceres() -> None:
    """Propaga (1) Ceres desde elementos osculadores de Horizons y compara la
    posición baricéntrica eclíptica contra los vectores de Horizons.
    """
    from astroquery.jplhorizons import Horizons

    # Elementos osculadores heliocéntricos eclípticos en la época.
    el_tab = Horizons(id="1", location="@sun", epochs=_EPOCH).elements(refplane="ecliptic")
    el = KeplerElements(
        a=float(el_tab["a"][0]),
        e=float(el_tab["e"][0]),
        i=math.radians(float(el_tab["incl"][0])),
        Omega=math.radians(float(el_tab["Omega"][0])),
        omega=math.radians(float(el_tab["w"][0])),
        M=math.radians(float(el_tab["M"][0])),
    )
    out_epochs = _EPOCH + np.array([30.0, 180.0, 365.0])
    got = propagate(
        el,
        _EPOCH,
        out_epochs,
        perturbers=(
            "sun",
            "mercury",
            "venus",
            "earth",
            "mars",
            "jupiter",
            "saturn",
            "uranus",
            "neptune",
        ),
        integrator="ias15",
    )
    vec = Horizons(id="1", location="@0", epochs=list(out_epochs)).vectors(refplane="ecliptic")
    ref = np.column_stack([vec["x"], vec["y"], vec["z"]]).astype(float)

    resid_au = np.linalg.norm(got - ref, axis=1)
    # Sanity dinámico: acuerdo sub-arcsec a ~3 AU (~1e-5 AU). El nivel mas exige
    # igualar DE441 + modelo de fuerzas completo (T8).
    assert np.max(resid_au) < 1e-4, f"residuos vs Horizons (AU): {resid_au}"
