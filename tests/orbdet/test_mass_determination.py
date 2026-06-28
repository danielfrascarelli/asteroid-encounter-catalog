"""Tests de src/orbdet/mass_determination.py — ajuste conjunto órbita+masa.

Gate de T6 (closing-loop): inyectar una masa de perturbador en observaciones
sintéticas del objetivo y recuperarla con el ajuste conjunto a ratio ≈ 1.0 y σ
realista — lo que el LOO secuencial nunca logró.

El perturbador se coloca cerca del objetivo a mitad del arco (encuentro cercano)
para garantizar leverage de masa; sin leverage, ninguna metodología recupera la
masa (es el problema intrínseco de DR3, no del método).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from src.orbdet.dynamics import AsteroidPerturber, planet_state_ecliptic
from src.orbdet.frames import ecliptic_to_equatorial
from src.orbdet.kepler import KeplerElements, elements_to_state, state_to_elements
from src.orbdet.mass_determination import determine_mass_and_orbit
from src.orbdet.observation import predict_radec

_EPOCH = 2_457_000.5
_TARGET = KeplerElements(
    a=2.70,
    e=0.15,
    i=math.radians(10.0),
    Omega=math.radians(80.0),
    omega=math.radians(60.0),
    M=math.radians(45.0),
)
_MASS_TRUE = 3e-9  # M_sun (~6× Ceres): señal de deflexión clara
_PERTURBERS = ("sun", "jupiter")


def _perturber_elements() -> KeplerElements:
    """Perturbador situado ~0.004 AU del objetivo en la época (encuentro cercano).

    Se deriva del estado heliocéntrico del objetivo + un offset de posición y un
    pequeño delta de velocidad, de modo que pasan cerca a mitad del arco con
    velocidad relativa no nula → leverage de masa garantizado.
    """
    r_t, v_t = elements_to_state(_TARGET)
    r_p = r_t + np.array([0.004, 0.0, 0.0])
    v_p = v_t + np.array([0.0, 5e-4, 0.0])
    return state_to_elements(r_p, v_p)


def _setup(n_obs: int, span_days: float):
    obs = _EPOCH + np.linspace(-span_days / 2, span_days / 2, n_obs)
    gaia_ecl = np.array([planet_state_ecliptic("earth", float(t))[0] for t in obs])
    gaia_icrs = ecliptic_to_equatorial(gaia_ecl)
    pert_el = _perturber_elements()
    pert_true = AsteroidPerturber("pert", _MASS_TRUE, pert_el)
    ra_t, dec_t = predict_radec(
        _TARGET, _EPOCH, obs, gaia_icrs, perturbers=_PERTURBERS, asteroid_perturbers=(pert_true,)
    )
    return obs, gaia_icrs, pert_el, ra_t, dec_t


@pytest.mark.slow
def test_closing_loop_noiseless() -> None:
    """GATE: sin ruido, recupera la masa inyectada con ratio ≈ 1.0."""
    obs, gaia_icrs, pert_el, ra_t, dec_t = _setup(n_obs=30, span_days=500.0)
    n = obs.size
    pa = np.linspace(15.0, 165.0, n)
    sigma_al = np.full(n, 1.0)

    mass_fit, _el_fit, res = determine_mass_and_orbit(
        _TARGET,
        0.4 * _MASS_TRUE,  # semilla de masa errada (0.4×)
        pert_el,
        _EPOCH,
        obs,
        ra_t,
        dec_t,
        pa,
        sigma_al,
        gaia_icrs,
        perturbers=_PERTURBERS,
    )
    assert res.converged
    assert res.chi2 < 1e-4  # residuos blanqueados → 0 en la masa+órbita verdadera
    assert mass_fit / _MASS_TRUE == pytest.approx(1.0, abs=2e-3)


@pytest.mark.slow
def test_closing_loop_with_noise_sigma_realistic() -> None:
    """GATE: con ruido AL, recupera la masa dentro de pocos σ y con σ finita."""
    obs, gaia_icrs, pert_el, ra_t, dec_t = _setup(n_obs=50, span_days=600.0)
    n = obs.size
    rng = np.random.default_rng(2024)
    pa = rng.uniform(0.0, 180.0, n)
    sigma_al = np.full(n, 0.5)  # mas
    noise = rng.normal(0.0, sigma_al)
    ra_obs = ra_t + noise * np.sin(np.radians(pa)) / (np.cos(np.radians(dec_t)) * 3.6e6)
    dec_obs = dec_t + noise * np.cos(np.radians(pa)) / 3.6e6

    # Semilla: masa 0.5×, elementos levemente perturbados.
    seed_el = KeplerElements(
        a=_TARGET.a + 1e-5,
        e=_TARGET.e + 1e-5,
        i=_TARGET.i + 1e-5,
        Omega=_TARGET.Omega + 1e-5,
        omega=_TARGET.omega + 1e-5,
        M=_TARGET.M + 1e-5,
    )
    mass_fit, _el_fit, res = determine_mass_and_orbit(
        seed_el,
        0.5 * _MASS_TRUE,
        pert_el,
        _EPOCH,
        obs,
        ra_obs,
        dec_obs,
        pa,
        sigma_al,
        gaia_icrs,
        perturbers=_PERTURBERS,
    )
    sigma_mass = math.sqrt(res.covariance[0, 0])
    assert res.converged
    assert 0.4 < res.chi2_reduced < 2.0
    assert sigma_mass > 0.0 and np.isfinite(sigma_mass)
    # Masa recuperada dentro de 3σ de la verdad y σ informativa (<50% de la masa).
    assert abs(mass_fit - _MASS_TRUE) < 3.0 * sigma_mass
    assert sigma_mass < 0.5 * _MASS_TRUE
