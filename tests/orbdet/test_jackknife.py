"""Tests del jackknife dejar-un-objetivo-fuera para σ(masa) externa (F1).

``jackknife_mass_sigma`` re-ajusta la masa compartida ``N`` veces, cada una
excluyendo un objetivo, y reporta la dispersión de las réplicas como σ externa —
la que captura el error de regresión masa↔órbita que la σ formal (Fisher) no ve.

Gate de comportamiento:
- ``N < 3`` → σ_jack ``nan`` (jackknife no informativo) sin ajustar nada.
- ``N ≥ 3`` con objetivos ruidosos → σ_jack finita y positiva; tantas réplicas
  como objetivos convergidos.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from src.orbdet.dynamics import AsteroidPerturber, planet_state_ecliptic
from src.orbdet.frames import ecliptic_to_equatorial
from src.orbdet.kepler import KeplerElements, elements_to_state, state_to_elements
from src.orbdet.mass_determination import (
    TargetObservations,
    determine_shared_mass,
    jackknife_mass_sigma,
)
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
_MASS_TRUE = 3e-9
_PERTURBERS = ("sun",)


def _perturber_elements() -> KeplerElements:
    r_t, v_t = elements_to_state(_TARGET)
    r_p = r_t + np.array([0.004, 0.0, 0.0])
    v_p = v_t + np.array([0.0, 5e-4, 0.0])
    return state_to_elements(r_p, v_p)


def _make_target(
    pert_el: KeplerElements, noise_mas: float, rng: np.random.Generator
) -> TargetObservations:
    obs = _EPOCH + np.linspace(-200.0, 200.0, 18)
    gaia_ecl = np.array([planet_state_ecliptic("earth", float(t))[0] for t in obs])
    gaia_icrs = ecliptic_to_equatorial(gaia_ecl)
    pert_true = AsteroidPerturber("pert", _MASS_TRUE, pert_el)
    ra_t, dec_t = predict_radec(
        _TARGET, _EPOCH, obs, gaia_icrs, perturbers=_PERTURBERS, asteroid_perturbers=(pert_true,)
    )
    n = obs.size
    # Ruido independiente por objetivo → cada uno constriñe una masa algo distinta,
    # así quitar uno mueve la media y σ_jack > 0 (réplicas idénticas darían 0 exacto).
    deg_per_mas = 1.0 / 3.6e6
    ra_t = ra_t + rng.normal(0.0, noise_mas * deg_per_mas, n) / np.cos(np.radians(dec_t))
    dec_t = dec_t + rng.normal(0.0, noise_mas * deg_per_mas, n)
    pa = np.linspace(20.0, 160.0, n)
    return TargetObservations(
        initial_elements=_TARGET,
        obs_jd_tdb=obs,
        ra_obs_deg=ra_t,
        dec_obs_deg=dec_t,
        pa_scan_deg=pa,
        sigma_al_mas=np.full(n, noise_mas),
        gaia_bary_icrs=gaia_icrs,
    )


def test_jackknife_returns_nan_for_small_n() -> None:
    """Con N < 3 el jackknife no es informativo → σ_jack nan, sin ajustar."""
    pert_el = _perturber_elements()
    rng = np.random.default_rng(0)
    targets = [_make_target(pert_el, 1.0, rng) for _ in range(2)]
    res = jackknife_mass_sigma(
        targets,
        _MASS_TRUE,
        [_TARGET, _TARGET],
        pert_el,
        _EPOCH,
        perturbers=_PERTURBERS,
        backend="rebound",
    )
    assert math.isnan(res.sigma_jack_msun)
    assert res.masses_msun.size == 0


def test_jackknife_length_mismatch_raises() -> None:
    """targets y fitted_elements deben tener la misma longitud."""
    pert_el = _perturber_elements()
    rng = np.random.default_rng(0)
    targets = [_make_target(pert_el, 1.0, rng) for _ in range(3)]
    with pytest.raises(ValueError):
        jackknife_mass_sigma(
            targets, _MASS_TRUE, [_TARGET], pert_el, _EPOCH, perturbers=_PERTURBERS
        )


@pytest.mark.slow
def test_jackknife_sigma_positive_and_finite() -> None:
    """GATE: con N≥3 objetivos ruidosos, σ_jack es finita y > 0 y hay N réplicas."""
    pert_el = _perturber_elements()
    rng = np.random.default_rng(42)
    n = 4
    targets = [_make_target(pert_el, 2.0, rng) for _ in range(n)]
    mass_fit, fitted, result = determine_shared_mass(
        targets, 0.6 * _MASS_TRUE, pert_el, _EPOCH, perturbers=_PERTURBERS
    )
    assert result.converged
    # La masa se recupera dentro de unos pocos σ formales pese al ruido.
    assert mass_fit / _MASS_TRUE == pytest.approx(1.0, abs=0.3)

    jk = jackknife_mass_sigma(
        targets,
        mass_fit,
        fitted,
        pert_el,
        _EPOCH,
        perturbers=_PERTURBERS,
        backend="rebound",
    )
    assert jk.masses_msun.size + jk.n_failed == n
    assert jk.masses_msun.size >= 3
    assert math.isfinite(jk.sigma_jack_msun)
    assert jk.sigma_jack_msun > 0.0
