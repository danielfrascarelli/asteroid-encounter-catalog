"""Tests de src/orbdet/orbit_determination.py — corrector diferencial.

Gate de T5: recuperar una órbita conocida a partir de observaciones sintéticas.
  - Sin ruido: el ajuste recupera los elementos verdaderos a alta precisión desde
    una semilla perturbada (residuos → 0).
  - Con ruido AL: el ajuste converge, χ²_red ≈ 1, y los elementos quedan dentro de
    pocos σ de la verdad (cubre la incertidumbre de la covarianza).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from src.orbdet.dynamics import planet_state_ecliptic
from src.orbdet.frames import ecliptic_to_equatorial
from src.orbdet.kepler import KeplerElements
from src.orbdet.observation import predict_radec
from src.orbdet.orbit_determination import determine_orbit

_EPOCH = 2_457_000.5
_TRUTH = KeplerElements(
    a=2.61,
    e=0.12,
    i=math.radians(8.0),
    Omega=math.radians(74.0),
    omega=math.radians(120.0),
    M=math.radians(200.0),
)
_PERTURBERS = ("sun", "jupiter")


def _synthetic_setup(n_obs: int = 20, span_days: float = 500.0):
    """Épocas, posiciones de Gaia (proxy Tierra) y RA/Dec verdaderos."""
    obs = _EPOCH + np.linspace(-span_days / 2, span_days / 2, n_obs)
    gaia_ecl = np.array([planet_state_ecliptic("earth", float(t))[0] for t in obs])
    gaia_icrs = ecliptic_to_equatorial(gaia_ecl)
    ra_t, dec_t = predict_radec(_TRUTH, _EPOCH, obs, gaia_icrs, perturbers=_PERTURBERS)
    return obs, gaia_icrs, ra_t, dec_t


def _perturbed_seed() -> KeplerElements:
    return KeplerElements(
        a=_TRUTH.a + 5e-4,
        e=_TRUTH.e + 1e-3,
        i=_TRUTH.i + math.radians(0.05),
        Omega=_TRUTH.Omega - math.radians(0.05),
        omega=_TRUTH.omega + math.radians(0.05),
        M=_TRUTH.M - math.radians(0.05),
    )


@pytest.mark.slow
def test_recover_known_orbit_noiseless() -> None:
    """GATE: sin ruido, el corrector recupera los elementos verdaderos."""
    obs, gaia_icrs, ra_t, dec_t = _synthetic_setup()
    n = obs.size
    pa = np.linspace(10.0, 170.0, n)  # ángulos de barrido variados
    sigma_al = np.full(n, 1.0)

    fitted, res = determine_orbit(
        _perturbed_seed(),
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
    assert res.chi2 < 1e-6  # residuos blanqueados → 0
    assert fitted.a == pytest.approx(_TRUTH.a, abs=1e-7)
    assert fitted.e == pytest.approx(_TRUTH.e, abs=1e-7)
    assert fitted.i == pytest.approx(_TRUTH.i, abs=1e-7)
    assert fitted.Omega == pytest.approx(_TRUTH.Omega, abs=1e-7)
    assert fitted.omega == pytest.approx(_TRUTH.omega, abs=1e-7)
    assert fitted.M == pytest.approx(_TRUTH.M, abs=1e-7)


@pytest.mark.slow
def test_recover_known_orbit_with_al_noise() -> None:
    """GATE: con ruido AL, converge a χ²_red≈1 y elementos dentro de pocos σ."""
    obs, gaia_icrs, ra_t, dec_t = _synthetic_setup(n_obs=40, span_days=600.0)
    n = obs.size
    rng = np.random.default_rng(7)
    pa = rng.uniform(0.0, 180.0, n)
    sigma_al = np.full(n, 1.0)  # mas
    noise = rng.normal(0.0, sigma_al)  # along-scan, mas
    ra_obs = ra_t + noise * np.sin(np.radians(pa)) / (np.cos(np.radians(dec_t)) * 3.6e6)
    dec_obs = dec_t + noise * np.cos(np.radians(pa)) / 3.6e6

    fitted, res = determine_orbit(
        _perturbed_seed(),
        _EPOCH,
        obs,
        ra_obs,
        dec_obs,
        pa,
        sigma_al,
        gaia_icrs,
        perturbers=_PERTURBERS,
    )
    assert res.converged
    assert 0.5 < res.chi2_reduced < 1.8
    # Elementos dentro de ~5σ de la covarianza del ajuste.
    sigma = np.sqrt(np.diag(res.covariance))
    delta = np.abs(np.array(fitted.as_array()) - np.array(_TRUTH.as_array()))
    assert np.all(delta < 5.0 * sigma), f"delta={delta}, sigma={sigma}"
