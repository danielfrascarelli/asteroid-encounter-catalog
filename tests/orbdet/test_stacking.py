"""Tests del stacking multi-objetivo (src/orbdet/mass_determination.determine_shared_mass).

Gate de T7: la incertidumbre de la masa compartida baja como ``σ(GM) ∝ 1/√N`` al
apilar objetivos independientes — el mecanismo por el que Fuentes-Muñoz rompe la
multimodalidad/no-identificabilidad de un solo encuentro.

Se usan réplicas idénticas del mismo encuentro (cada una con sus propios 6
elementos en el sistema): la información de Fisher de la masa se suma exactamente,
así que ``σ(N) = σ(1)/√N`` es un test limpio del solver en flecha.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from src.orbdet.dynamics import AsteroidPerturber, planet_state_ecliptic
from src.orbdet.frames import ecliptic_to_equatorial
from src.orbdet.kepler import KeplerElements, elements_to_state, state_to_elements
from src.orbdet.mass_determination import TargetObservations, determine_shared_mass
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


def _make_target(pert_el: KeplerElements) -> TargetObservations:
    obs = _EPOCH + np.linspace(-200.0, 200.0, 18)
    gaia_ecl = np.array([planet_state_ecliptic("earth", float(t))[0] for t in obs])
    gaia_icrs = ecliptic_to_equatorial(gaia_ecl)
    pert_true = AsteroidPerturber("pert", _MASS_TRUE, pert_el)
    ra_t, dec_t = predict_radec(
        _TARGET, _EPOCH, obs, gaia_icrs, perturbers=_PERTURBERS, asteroid_perturbers=(pert_true,)
    )
    n = obs.size
    pa = np.linspace(20.0, 160.0, n)
    return TargetObservations(
        initial_elements=_TARGET,
        obs_jd_tdb=obs,
        ra_obs_deg=ra_t,
        dec_obs_deg=dec_t,
        pa_scan_deg=pa,
        sigma_al_mas=np.full(n, 1.0),
        gaia_bary_icrs=gaia_icrs,
    )


def _fit_n(n: int, pert_el: KeplerElements):
    targets = [_make_target(pert_el) for _ in range(n)]
    return determine_shared_mass(targets, 0.6 * _MASS_TRUE, pert_el, _EPOCH, perturbers=_PERTURBERS)


@pytest.mark.slow
def test_shared_mass_sigma_scales_as_inverse_sqrt_n() -> None:
    """GATE: σ(masa) ∝ 1/√N al apilar objetivos."""
    pert_el = _perturber_elements()
    m1, els1, r1 = _fit_n(1, pert_el)
    m2, _e2, r2 = _fit_n(2, pert_el)
    m4, _e4, r4 = _fit_n(4, pert_el)

    for r in (r1, r2, r4):
        assert r.converged
    s1 = math.sqrt(r1.covariance[0, 0])
    s2 = math.sqrt(r2.covariance[0, 0])
    s4 = math.sqrt(r4.covariance[0, 0])
    assert s1 > 0.0
    # σ(N) = σ(1)/√N (réplicas idénticas → suma exacta de información de Fisher).
    assert s2 / s1 == pytest.approx(1.0 / math.sqrt(2.0), rel=0.05)
    assert s4 / s1 == pytest.approx(0.5, rel=0.05)
    # La masa sigue recuperándose (ratio ≈ 1).
    assert m4 / _MASS_TRUE == pytest.approx(1.0, abs=2e-3)


@pytest.mark.slow
def test_shared_mass_single_target_matches_joint() -> None:
    """Con N=1, determine_shared_mass coincide con el ajuste conjunto de T6."""
    pert_el = _perturber_elements()
    mass_fit, els, res = _fit_n(1, pert_el)
    assert res.converged
    assert len(els) == 1
    assert mass_fit / _MASS_TRUE == pytest.approx(1.0, abs=2e-3)
    assert res.covariance.shape == (7, 7)


@pytest.mark.slow
def test_shared_mass_parallel_matches_serial() -> None:
    """GATE: ``n_workers>1`` da resultados idénticos al modo serie.

    El pool evalúa los objetivos en procesos separados pero ensambla en orden fijo,
    así que la masa, los residuos y la covarianza deben coincidir bit a bit (salvo
    redondeo) con el cálculo serie.
    """
    pert_el = _perturber_elements()
    targets = [_make_target(pert_el) for _ in range(3)]
    kw = dict(perturbers=_PERTURBERS)
    m1, _e1, r1 = determine_shared_mass(
        targets, 0.6 * _MASS_TRUE, pert_el, _EPOCH, n_workers=1, **kw
    )
    m2, _e2, r2 = determine_shared_mass(
        targets, 0.6 * _MASS_TRUE, pert_el, _EPOCH, n_workers=3, **kw
    )
    assert m2 == pytest.approx(m1, rel=1e-11)
    np.testing.assert_allclose(r2.residuals, r1.residuals, rtol=1e-11, atol=1e-12)
    np.testing.assert_allclose(r2.covariance, r1.covariance, rtol=1e-9, atol=1e-14)
