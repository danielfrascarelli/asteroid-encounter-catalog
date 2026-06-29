"""Tests del adaptador de datos reales → motor (T9).

Cubre las cuatro conversiones sin atajos del adaptador:
- σ_AL = proyección de la covarianza (RA, Dec) sobre la dirección de barrido.
- elementos MPCORB (grados) → KeplerElements (radianes).
- propagación de elementos entre épocas con el modelo N-cuerpos (round-trip).
- ensamblado de TargetObservations consistente con el modelo de observación.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from src.orbdet.dynamics import planet_state_ecliptic
from src.orbdet.frames import ecliptic_to_equatorial
from src.orbdet.gaia_adapter import (
    build_target_observations,
    elements_from_mpcorb,
    propagate_elements,
    sigma_al_from_radec_covariance,
)
from src.orbdet.kepler import KeplerElements
from src.orbdet.observation import predict_radec
from src.orbdet.time_scales import J2010_TCB_JD, tdb_to_tcb
from tests.orbdet._ephem import requires_ephem

# --- σ_AL desde la covarianza (RA, Dec) -------------------------------------


def test_sigma_al_alignment_limits() -> None:
    """Con ρ=0 y familia aleatoria nula, PA=90° → σ_AL=σ_RA; PA=0° → σ_AL=σ_Dec."""
    zeros = np.zeros(2)
    s_ra = np.array([2.0, 2.0])
    s_dec = np.array([5.0, 5.0])
    pa = np.array([90.0, 0.0])  # AL a lo largo de RA, luego a lo largo de Dec
    sigma_al = sigma_al_from_radec_covariance(pa, s_ra, s_dec, zeros, zeros, zeros, zeros)
    assert sigma_al[0] == pytest.approx(2.0, rel=1e-12)
    assert sigma_al[1] == pytest.approx(5.0, rel=1e-12)


def test_sigma_al_matches_quadratic_form() -> None:
    """σ²_AL coincide con û_ALᵀ(Σ_sys+Σ_rand)û_AL evaluado explícitamente."""
    rng = np.random.default_rng(0)
    n = 7
    pa = rng.uniform(0, 180, n)
    s_ra_s, s_dec_s = rng.uniform(0.2, 2, n), rng.uniform(0.2, 2, n)
    rho_s = rng.uniform(-0.5, 0.5, n)
    s_ra_r, s_dec_r = rng.uniform(0.2, 2, n), rng.uniform(0.2, 2, n)
    rho_r = rng.uniform(-0.5, 0.5, n)
    got = sigma_al_from_radec_covariance(pa, s_ra_s, s_dec_s, rho_s, s_ra_r, s_dec_r, rho_r)

    expect = np.empty(n)
    for k in range(n):
        u = np.array([math.sin(math.radians(pa[k])), math.cos(math.radians(pa[k]))])
        cov_s = np.array(
            [
                [s_ra_s[k] ** 2, rho_s[k] * s_ra_s[k] * s_dec_s[k]],
                [rho_s[k] * s_ra_s[k] * s_dec_s[k], s_dec_s[k] ** 2],
            ]
        )
        cov_r = np.array(
            [
                [s_ra_r[k] ** 2, rho_r[k] * s_ra_r[k] * s_dec_r[k]],
                [rho_r[k] * s_ra_r[k] * s_dec_r[k], s_dec_r[k] ** 2],
            ]
        )
        expect[k] = math.sqrt(u @ (cov_s + cov_r) @ u)
    np.testing.assert_allclose(got, expect, rtol=1e-12)


# --- elementos MPCORB → KeplerElements --------------------------------------


def test_elements_from_mpcorb_degrees_to_radians() -> None:
    el = elements_from_mpcorb(2.5, 0.1, 10.0, 80.0, 60.0, 45.0)
    assert el.a == 2.5
    assert el.e == 0.1
    assert el.i == pytest.approx(math.radians(10.0))
    assert el.Omega == pytest.approx(math.radians(80.0))
    assert el.omega == pytest.approx(math.radians(60.0))
    assert el.M == pytest.approx(math.radians(45.0))


# --- propagación de elementos entre épocas ----------------------------------


@pytest.mark.slow
def test_propagate_elements_identity() -> None:
    el = elements_from_mpcorb(2.7, 0.12, 8.0, 100.0, 50.0, 30.0)
    out = propagate_elements(el, 2_457_000.5, 2_457_000.5)
    np.testing.assert_allclose(out.as_array(), el.as_array())


@requires_ephem
@pytest.mark.slow
def test_propagate_elements_roundtrip_assist() -> None:
    """Propagar adelante y volver recupera los elementos (reversibilidad ASSIST/IAS15).

    Con el backend ASSIST los planetas siguen la efeméride exactamente, así que la
    propagación es reversible a precisión de integrador (a diferencia de los planetas
    integrados libremente del backend rebound, que derivan ~1e-5 AU sobre 800 días).
    """
    el = elements_from_mpcorb(2.7, 0.12, 8.0, 100.0, 50.0, 30.0)
    epoch0 = 2_457_000.5
    epoch1 = epoch0 + 400.0
    fwd = propagate_elements(el, epoch0, epoch1, backend="assist")
    back = propagate_elements(fwd, epoch1, epoch0, backend="assist")
    # Piso ~1e-6 AU = paso adaptativo no-simétrico de IAS15 (fwd≠bwd), no del modelo;
    # ~0.1 mas a 2 AU, despreciable para una semilla y 10× mejor que rebound.
    np.testing.assert_allclose(back.a, el.a, rtol=1e-5)
    np.testing.assert_allclose(back.e, el.e, atol=1e-5)
    np.testing.assert_allclose(back.i, el.i, atol=1e-5)
    np.testing.assert_allclose(back.Omega, el.Omega, atol=1e-5)
    np.testing.assert_allclose(back.omega, el.omega, atol=1e-5)
    np.testing.assert_allclose(math.cos(back.M), math.cos(el.M), atol=1e-5)


# --- ensamblado de TargetObservations ----------------------------------------


@pytest.mark.slow
def test_build_target_observations_consistency() -> None:
    """TargetObservations armado desde columnas Gaia reproduce la geometría sintética.

    Genera observaciones sintéticas con :func:`predict_radec` (dinámica N-cuerpos),
    las empaqueta como columnas crudas de Gaia, y verifica que el adaptador
    reconstruye jd_tdb, σ_AL y posiciones de Gaia coherentes, con residuo ≈ 0 al
    re-predecir con los mismos elementos.
    """
    epoch = 2_457_000.5
    el = KeplerElements(
        a=2.70,
        e=0.15,
        i=math.radians(10.0),
        Omega=math.radians(80.0),
        omega=math.radians(60.0),
        M=math.radians(45.0),
    )
    n = 12
    obs_jd = epoch + np.linspace(-200.0, 200.0, n)
    gaia_ecl = np.array([planet_state_ecliptic("earth", float(t))[0] for t in obs_jd])
    gaia_icrs = ecliptic_to_equatorial(gaia_ecl)
    ra, dec = predict_radec(el, epoch, obs_jd, gaia_icrs, perturbers=("sun", "jupiter"))

    # Empaqueta como columnas crudas: invierte el mapeo de tiempos (JD TDB → época TCB).
    epoch_days = tdb_to_tcb(obs_jd) - J2010_TCB_JD
    pa = np.linspace(15.0, 165.0, n)
    s = np.full(n, 1.0)
    rho = np.zeros(n)

    tobs = build_target_observations(
        el,
        epoch,
        epoch_days_tcb=epoch_days,
        ra_deg=ra,
        dec_deg=dec,
        pa_scan_deg=pa,
        ra_err_sys=s,
        dec_err_sys=s * 100.0,
        corr_sys=rho,
        ra_err_rand=s,
        dec_err_rand=s * 100.0,
        corr_rand=rho,
        x_gaia=gaia_icrs[:, 0],
        y_gaia=gaia_icrs[:, 1],
        z_gaia=gaia_icrs[:, 2],
    )

    # Tiempos reconstruidos coinciden con los originales (round-trip TDB→TCB→TDB).
    np.testing.assert_allclose(tobs.obs_jd_tdb, obs_jd, atol=1e-6)
    # Posición de Gaia preservada.
    np.testing.assert_allclose(tobs.gaia_bary_icrs, gaia_icrs)
    # σ_AL > 0 y finita.
    assert np.all(np.isfinite(tobs.sigma_al_mas)) and np.all(tobs.sigma_al_mas > 0)

    # Re-predecir con los mismos elementos da residuo AL ≈ 0.
    ra_p, dec_p = predict_radec(
        tobs.initial_elements,
        epoch,
        tobs.obs_jd_tdb,
        tobs.gaia_bary_icrs,
        perturbers=("sun", "jupiter"),
    )
    np.testing.assert_allclose(ra_p, ra, atol=1e-9)
    np.testing.assert_allclose(dec_p, dec, atol=1e-9)
