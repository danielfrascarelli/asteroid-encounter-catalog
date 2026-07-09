"""Tests del bootstrap no paramétrico para σ(masa) externa (B6).

``bootstrap_mass_sigma`` resamplea los ``N`` objetivos con reemplazo ``B`` veces y
re-ajusta la masa compartida en cada muestra. A diferencia del jackknife (sensible
a que una sola réplica de alto leverage domine la varianza), el bootstrap promedia
esa influencia, dando una σ y un intervalo percentil robustos.

Gate de comportamiento:
- ``N < 3`` → σ_boot ``nan`` sin ajustar nada.
- longitudes desalineadas → ``ValueError``.
- ``N ≥ 3`` con objetivos ruidosos → σ_boot finita y positiva, CI95 ordenado, y
  el resultado es reproducible con la misma semilla.
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
    bootstrap_mass_sigma,
    determine_shared_mass,
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


def test_bootstrap_returns_nan_for_small_n() -> None:
    """Con N < 3 el bootstrap no es informativo → σ_boot nan, sin ajustar."""
    pert_el = _perturber_elements()
    rng = np.random.default_rng(0)
    targets = [_make_target(pert_el, 1.0, rng) for _ in range(2)]
    res = bootstrap_mass_sigma(
        targets,
        _MASS_TRUE,
        [_TARGET, _TARGET],
        pert_el,
        _EPOCH,
        n_boot=10,
        perturbers=_PERTURBERS,
        backend="rebound",
    )
    assert math.isnan(res.sigma_boot_msun)
    assert res.masses_msun.size == 0
    assert all(math.isnan(x) for x in res.ci95_msun)


def test_bootstrap_length_mismatch_raises() -> None:
    """targets y fitted_elements deben tener la misma longitud."""
    pert_el = _perturber_elements()
    rng = np.random.default_rng(0)
    targets = [_make_target(pert_el, 1.0, rng) for _ in range(3)]
    with pytest.raises(ValueError):
        bootstrap_mass_sigma(
            targets, _MASS_TRUE, [_TARGET], pert_el, _EPOCH, n_boot=10, perturbers=_PERTURBERS
        )


@pytest.mark.slow
def test_bootstrap_sigma_positive_and_reproducible() -> None:
    """GATE: con N≥3 ruidosos, σ_boot finita > 0, CI95 ordenado, y reproducible."""
    pert_el = _perturber_elements()
    rng = np.random.default_rng(42)
    n = 4
    targets = [_make_target(pert_el, 2.0, rng) for _ in range(n)]
    mass_fit, fitted, result = determine_shared_mass(
        targets, 0.6 * _MASS_TRUE, pert_el, _EPOCH, perturbers=_PERTURBERS
    )
    assert result.converged

    kwargs = dict(perturbers=_PERTURBERS, backend="rebound", n_boot=16, seed=7)
    b1 = bootstrap_mass_sigma(targets, mass_fit, fitted, pert_el, _EPOCH, **kwargs)
    assert math.isfinite(b1.sigma_boot_msun)
    assert b1.sigma_boot_msun > 0.0
    assert b1.ci95_msun[0] <= b1.median_msun <= b1.ci95_msun[1]
    assert b1.masses_msun.size + b1.n_failed == 16

    # Misma semilla → mismas masas bootstrap (reproducibilidad).
    b2 = bootstrap_mass_sigma(targets, mass_fit, fitted, pert_el, _EPOCH, **kwargs)
    np.testing.assert_allclose(b1.masses_msun, b2.masses_msun)
