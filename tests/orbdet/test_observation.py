"""Tests de src/orbdet/observation.py — modelo de observación + covarianza AL.

Gate de T4: los residuos de una órbita conocida quedan al nivel del ruido AL.
Se valida en dos partes:
  - El chain dinámica→ICRS→RA/Dec con light-time reproduce un cálculo analítico
    independiente (kepler dos-cuerpos) muy por debajo del ruido AL (<0.1 mas).
  - La covarianza anisotrópica blanquea ruido AL inyectado a χ²/obs ≈ 1.

La geometría y la covarianza son puras (rápidas, sin rebound). El gate end-to-end
con dinámica N-cuerpos va marcado ``slow``.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from src.orbdet.dynamics import planet_state_ecliptic
from src.orbdet.frames import ecliptic_to_equatorial
from src.orbdet.kepler import KeplerElements, elements_to_state
from src.orbdet.kepler import propagate as kepler_propagate
from src.orbdet.observation import (
    along_scan_residual,
    anisotropic_covariance,
    light_time_correct,
    predict_radec,
    radec_from_positions,
    radec_to_unit_vec,
    tangent_residuals_mas,
    whiten_residuals_2d,
    xyz_to_radec,
)

_EPOCH = 2_457_000.5
_EL = KeplerElements(
    a=2.7,
    e=0.15,
    i=math.radians(10.0),
    Omega=math.radians(80.0),
    omega=math.radians(60.0),
    M=math.radians(45.0),
)


# --- Geometría RA/Dec -------------------------------------------------------


def test_xyz_radec_roundtrip() -> None:
    ra = np.array([0.0, 90.0, 187.3, 359.9])
    dec = np.array([0.0, 45.0, -30.0, 12.5])
    vec = radec_to_unit_vec(ra, dec)
    ra2, dec2 = xyz_to_radec(vec)
    assert np.allclose(ra, ra2, atol=1e-10)
    assert np.allclose(dec, dec2, atol=1e-10)


def test_radec_from_positions_axes() -> None:
    gaia = np.array([[0.0, 0.0, 0.0]] * 3)
    ast = np.array([[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]])
    ra, dec = radec_from_positions(ast, gaia)
    assert np.allclose(ra[:2], [0.0, 90.0], atol=1e-9)
    assert dec[2] == pytest.approx(90.0, abs=1e-9)


def test_tangent_residuals_cos_dec_factor() -> None:
    # A 1e-5 deg offset en RA a Dec=60° → dra = 1e-5·cos60·3.6e6 = 18 mas.
    dra, ddec = tangent_residuals_mas(
        np.array([10.00001]), np.array([60.0]), np.array([10.0]), np.array([60.0])
    )
    # rel laxo: el wrap (x+540)%360−180 pierde ~5 dígitos de un offset diminuto por
    # cancelación (~4e-7 mas, despreciable), pero valida el factor cosδ y la escala.
    assert dra[0] == pytest.approx(1e-5 * math.cos(math.radians(60.0)) * 3.6e6, rel=1e-6)
    assert ddec[0] == pytest.approx(0.0, abs=1e-6)


def test_tangent_residuals_ra_wrap() -> None:
    # obs=0.0001°, pred=359.9999° → diferencia +0.0002° (no -359.9998°).
    dra, _ = tangent_residuals_mas(
        np.array([0.0001]), np.array([0.0]), np.array([359.9999]), np.array([0.0])
    )
    assert dra[0] == pytest.approx(0.0002 * 3.6e6, rel=1e-9)


# --- Covarianza anisotrópica ------------------------------------------------


def test_anisotropic_covariance_pa_zero() -> None:
    # PA=0 → along-scan a lo largo de Dec (y): Σyy=σ_AL², Σxx=σ_AC², sin correlación.
    sxx, sxy, syy = anisotropic_covariance(np.array([0.0]), np.array([0.5]), np.array([300.0]))
    assert syy[0] == pytest.approx(0.25)
    assert sxx[0] == pytest.approx(9e4)
    assert sxy[0] == pytest.approx(0.0, abs=1e-12)


def test_anisotropic_covariance_isotropic() -> None:
    sxx, sxy, syy = anisotropic_covariance(np.array([37.0]), np.array([2.0]), np.array([2.0]))
    assert sxx[0] == pytest.approx(4.0)
    assert syy[0] == pytest.approx(4.0)
    assert sxy[0] == pytest.approx(0.0, abs=1e-12)


def test_whiten_chi2_matches_mahalanobis() -> None:
    rng = np.random.default_rng(0)
    n = 500
    pa = rng.uniform(0.0, 180.0, n)
    sal = rng.uniform(0.2, 2.0, n)
    sac = rng.uniform(100.0, 400.0, n)
    dra = rng.normal(0.0, 5.0, n)
    ddec = rng.normal(0.0, 5.0, n)
    whitened, chi2 = whiten_residuals_2d(dra, ddec, pa, sal, sac)
    # r1²+r2² == chi2 por observación.
    assert np.allclose(whitened[0::2] ** 2 + whitened[1::2] ** 2, chi2, rtol=1e-9)
    # chi2 == δᵀ Σ⁻¹ δ por inversión explícita de la 2×2.
    sxx, sxy, syy = anisotropic_covariance(pa, sal, sac)
    chi2_ref = np.empty(n)
    for k in range(n):
        sigma = np.array([[sxx[k], sxy[k]], [sxy[k], syy[k]]])
        d = np.array([dra[k], ddec[k]])
        chi2_ref[k] = d @ np.linalg.solve(sigma, d)
    assert np.allclose(chi2, chi2_ref, rtol=1e-8)


def test_across_scan_error_barely_penalized() -> None:
    # Un residuo puramente across-scan de tamaño D con σ_AC grande → χ² ≈ (D/σ_AC)².
    pa, sal, sac, big = 30.0, 1.0, 1000.0, 50.0
    rad = math.radians(pa)
    dra = np.array([big * math.cos(rad)])  # D·û_AC = D·(cosPA, −sinPA)
    ddec = np.array([-big * math.sin(rad)])
    _w, chi2 = whiten_residuals_2d(dra, ddec, np.array([pa]), np.array([sal]), np.array([sac]))
    assert chi2[0] == pytest.approx((big / sac) ** 2, rel=1e-6)


def test_along_scan_unit_residual_chi2_one() -> None:
    # Residuo de tamaño σ_AL a lo largo de along-scan → χ² ≈ 1.
    pa, sal, sac = 30.0, 1.3, 500.0
    rad = math.radians(pa)
    dra = np.array([sal * math.sin(rad)])  # σ_AL·û_AL
    ddec = np.array([sal * math.cos(rad)])
    _w, chi2 = whiten_residuals_2d(dra, ddec, np.array([pa]), np.array([sal]), np.array([sac]))
    assert chi2[0] == pytest.approx(1.0, rel=1e-4)


def test_along_scan_residual_projection() -> None:
    dra = np.array([3.0])
    ddec = np.array([4.0])
    pa = np.array([90.0])  # sinPA=1, cosPA=0 → r_AL = dra
    r_al, sigma = along_scan_residual(dra, ddec, pa, np.array([0.7]))
    assert r_al[0] == pytest.approx(3.0)
    assert sigma[0] == pytest.approx(0.7)


# --- Light-time (analítico, sin rebound) ------------------------------------


def test_light_time_fixed_point() -> None:
    """El tiempo retardado satisface la ecuación implícita τ = ρ(t−τ)/c."""
    from src.orbdet.constants import C_AU_PER_DAY

    p0 = np.array([2.0, 0.5, 0.1])
    vel = np.array([0.0, 0.01, 0.0])  # AU/día
    t0 = _EPOCH

    def bary_ecl_at(jd: np.ndarray) -> np.ndarray:
        jd = np.atleast_1d(jd)
        return p0[None, :] + vel[None, :] * (jd[:, None] - t0)

    obs = _EPOCH + np.array([0.0, 10.0, -10.0])
    gaia_icrs = np.tile(ecliptic_to_equatorial(np.array([0.0, 0.0, 0.0])), (3, 1))
    jd_ret, ast_icrs = light_time_correct(bary_ecl_at, obs, gaia_icrs, n_iter=5)
    tau = obs - jd_ret
    rho = np.linalg.norm(ast_icrs - gaia_icrs, axis=1)
    assert np.allclose(tau, rho / C_AU_PER_DAY, rtol=1e-10)
    assert np.all(tau > 0.0)


def test_light_time_shifts_position() -> None:
    """Con light-time la posición usada es la retardada, no la del instante de obs."""
    p0 = np.array([2.0, 0.0, 0.0])
    vel = np.array([0.0, 0.02, 0.0])
    t0 = _EPOCH

    def bary_ecl_at(jd: np.ndarray) -> np.ndarray:
        jd = np.atleast_1d(jd)
        return p0[None, :] + vel[None, :] * (jd[:, None] - t0)

    obs = np.array([_EPOCH])
    gaia_icrs = np.tile(ecliptic_to_equatorial(np.array([0.0, 0.0, 0.0])), (1, 1))
    _jd_ret, ast_icrs = light_time_correct(bary_ecl_at, obs, gaia_icrs, n_iter=4)
    ast_no_lt = ecliptic_to_equatorial(bary_ecl_at(obs))
    # El desplazamiento es ~ v·τ, con τ = 2/c ≈ 0.01156 d → ~2.3e-4 AU en y.
    shift = np.linalg.norm(ast_icrs - ast_no_lt)
    assert shift == pytest.approx(0.02 * (2.0 / 173.1446), rel=1e-3)


# --- GATE T4: chain dinámica end-to-end vs oráculo analítico ----------------


@pytest.mark.slow
def test_predict_radec_matches_kepler_oracle() -> None:
    """GATE: predict_radec (N-cuerpos solo-Sol, con light-time) reproduce un
    cálculo analítico kepleriano independiente muy por debajo del ruido AL.

    Solo-Sol, el N-cuerpos = Kepler analítico (validado en T2 a 1e-8 AU), así que
    la diferencia angular debe ser ≪ 0.1 mas (el ruido AL es 0.2–2 mas)."""
    obs = _EPOCH + np.linspace(-250.0, 250.0, 12)
    gaia_ecl = np.array([planet_state_ecliptic("earth", float(t))[0] for t in obs])
    gaia_icrs = ecliptic_to_equatorial(gaia_ecl)

    # Oráculo analítico: Kepler heliocéntrico + Sol baricéntrico. En el modelo
    # solo-Sol, el Sol no siente fuerzas (la prueba es no masiva) → se mueve en
    # línea recta desde su estado en la época; el oráculo usa ese mismo Sol lineal,
    # no la efeméride curva, para coincidir con el N-cuerpos solo-Sol.
    sun_p0, sun_v0 = planet_state_ecliptic("sun", _EPOCH)

    def kepler_bary_ecl(jd: np.ndarray) -> np.ndarray:
        jd = np.atleast_1d(jd)
        out = np.empty((jd.size, 3))
        for k, t in enumerate(jd):
            el_t = kepler_propagate(_EL, float(t) - _EPOCH)
            r_h = elements_to_state(el_t)[0]
            out[k] = r_h + sun_p0 + sun_v0 * (float(t) - _EPOCH)
        return out

    _jd_ret, ast_ref = light_time_correct(kepler_bary_ecl, obs, gaia_icrs, n_iter=3)
    ra_ref, dec_ref = radec_from_positions(ast_ref, gaia_icrs)

    ra_pred, dec_pred = predict_radec(
        _EL, _EPOCH, obs, gaia_icrs, perturbers=("sun",), integrator="ias15"
    )
    dra, ddec = tangent_residuals_mas(ra_ref, dec_ref, ra_pred, dec_pred)
    max_resid_mas = float(np.max(np.hypot(dra, ddec)))
    assert max_resid_mas < 0.1, f"residuo máximo {max_resid_mas:.3e} mas"


@pytest.mark.slow
def test_predict_radec_al_noise_chi2_unity() -> None:
    """GATE: ruido AL gaussiano de σ_AL inyectado sobre la verdad, blanqueado con
    la covarianza anisotrópica, da χ²/obs ≈ 1 (residuos al nivel del ruido AL)."""
    obs = _EPOCH + np.linspace(-250.0, 250.0, 60)
    gaia_ecl = np.array([planet_state_ecliptic("earth", float(t))[0] for t in obs])
    gaia_icrs = ecliptic_to_equatorial(gaia_ecl)
    ra_t, dec_t = predict_radec(_EL, _EPOCH, obs, gaia_icrs, perturbers=("sun", "jupiter"))

    rng = np.random.default_rng(12345)
    n = obs.size
    pa = rng.uniform(0.0, 180.0, n)
    sigma_al, sigma_ac = 0.5, 400.0
    noise = rng.normal(0.0, sigma_al, n)  # mas a lo largo de along-scan
    dra_inj = noise * np.sin(np.radians(pa))
    ddec_inj = noise * np.cos(np.radians(pa))
    ra_obs = ra_t + dra_inj / (np.cos(np.radians(dec_t)) * 3.6e6)
    dec_obs = dec_t + ddec_inj / 3.6e6

    dra, ddec = tangent_residuals_mas(ra_obs, dec_obs, ra_t, dec_t)
    _w, chi2 = whiten_residuals_2d(dra, ddec, pa, np.full(n, sigma_al), np.full(n, sigma_ac))
    assert 0.6 < float(np.mean(chi2)) < 1.6, f"chi2/obs medio = {np.mean(chi2):.3f}"
