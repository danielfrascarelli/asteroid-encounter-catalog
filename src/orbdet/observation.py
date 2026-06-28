"""Modelo de observación de Gaia: estado → RA/Dec + covarianza along-scan.

Cierra el lazo entre la dinámica (``dynamics``/``variational``, que dan el estado
baricéntrico eclíptico del objetivo) y lo que Gaia mide (RA/Dec en el frame
astrométrico baricéntrico ICRS), con dos piezas que el método de masas necesita:

1. **Modelo geométrico de observación** — estado baricéntrico eclíptico → ICRS →
   línea de visión desde la posición de Gaia → RA/Dec, con **corrección de
   light-time iterativa** (la posición observada es la del asteroide cuando emitió
   la luz, no cuando Gaia la recibió). Para no acoplar este módulo al integrador,
   la light-time toma un *callable* ``bary_ecl_at(jd) -> (..., 3)`` que entrega la
   posición baricéntrica eclíptica a tiempos arbitrarios (lo provee la dinámica, o
   un propagador kepleriano analítico en los tests).
2. **Covarianza anisotrópica along-scan / across-scan** — Gaia mide con gran
   precisión a lo largo de la dirección de barrido (AL, σ ~0.2–2 mas) y muy mal en
   la perpendicular (AC, σ ~cientos de mas). La covarianza por observación es la
   elipse ``Σ = σ_AL² û_AL û_ALᵀ + σ_AC² û_AC û_ACᵀ`` en el plano tangente
   (RA·cosδ, δ), orientada por el ángulo de posición de barrido. Se provee tanto la
   proyección AL pura (línea base) como el blanqueo Mahalanobis 2D completo.

Convención del ángulo de barrido (idéntica a ``src/mass``): el versor along-scan
en el plano tangente es ``û_AL = (sin PA, cos PA)`` con PA medido desde el Norte
(eje +δ) hacia el Este (eje +RA·cosδ); el versor across-scan es su perpendicular
``û_AC = (cos PA, −sin PA)``.

Unidades: ángulos de salida en grados (interfaz); residuos tangenciales en mas;
posiciones en AU; tiempos en JD TDB.
"""

from __future__ import annotations

import numpy as np

from .constants import C_AU_PER_DAY, OBLIQUITY_J2000_RAD
from .frames import ecliptic_to_equatorial, rotation_x

_MAS_PER_DEG: float = 3_600_000.0
_MAS_PER_RAD: float = np.degrees(1.0) * _MAS_PER_DEG  # 1 rad en mas (≈2.0626e8)
_DEG: float = np.pi / 180.0
# Piso del determinante de la covarianza 2×2 (mas⁴) antes de caer a diagonal.
_DET_FLOOR: float = 1e-12


# --- Geometría: Cartesiano ICRS ↔ RA/Dec ------------------------------------


def xyz_to_radec(vec: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Cartesiano ICRS ``(..., 3)`` → ``(RA_deg, Dec_deg)`` con RA en ``[0, 360)``."""
    v = np.asarray(vec, dtype=float)
    x, y, z = v[..., 0], v[..., 1], v[..., 2]
    rho = np.sqrt(x * x + y * y)
    ra = np.degrees(np.arctan2(y, x)) % 360.0
    dec = np.degrees(np.arctan2(z, rho))
    return ra, dec


def radec_to_unit_vec(ra_deg: np.ndarray, dec_deg: np.ndarray) -> np.ndarray:
    """``(RA_deg, Dec_deg)`` → versores ICRS ``(..., 3)``."""
    ra = np.asarray(ra_deg, dtype=float) * _DEG
    dec = np.asarray(dec_deg, dtype=float) * _DEG
    cos_dec = np.cos(dec)
    return np.stack([cos_dec * np.cos(ra), cos_dec * np.sin(ra), np.sin(dec)], axis=-1)


def radec_from_positions(
    ast_bary_icrs: np.ndarray, gaia_bary_icrs: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """RA/Dec (grados) del asteroide visto desde Gaia (geométrico, sin light-time).

    Ambas posiciones en **ICRS baricéntrico**, AU. La línea de visión es
    ``r_ast − r_gaia`` (frame astrométrico baricéntrico: Gaia ya remueve la
    aberración, por eso se usa la posición baricéntrica de Gaia y no su velocidad).
    """
    return xyz_to_radec(
        np.asarray(ast_bary_icrs, dtype=float) - np.asarray(gaia_bary_icrs, dtype=float)
    )


# --- Corrección de light-time -----------------------------------------------


def light_time_correct(
    bary_ecl_at,
    obs_jd_tdb: np.ndarray,
    gaia_bary_icrs: np.ndarray,
    *,
    n_iter: int = 3,
    c_au_per_day: float = C_AU_PER_DAY,
) -> tuple[np.ndarray, np.ndarray]:
    """Itera el tiempo retardado de emisión y devuelve la posición del asteroide.

    La luz que Gaia recibe en ``t_obs`` fue emitida en ``t_emit = t_obs − τ`` con
    ``τ = |r_ast(t_emit) − r_gaia(t_obs)| / c``. Se itera (3 pasos convergen a
    ≪ 1 s a 3 AU). ``bary_ecl_at(jd_array)`` devuelve la posición baricéntrica
    **eclíptica** ``(N, 3)`` del asteroide; acá se rota a ICRS.

    Returns
    -------
    (jd_retarded, ast_bary_icrs)
        ``jd_retarded`` ``(N,)`` y la posición ICRS ``(N, 3)`` del asteroide en el
        tiempo de emisión, lista para :func:`radec_from_positions`.
    """
    obs = np.atleast_1d(np.asarray(obs_jd_tdb, dtype=float))
    gaia = np.asarray(gaia_bary_icrs, dtype=float)
    tau = np.zeros(obs.shape[0], dtype=float)
    ast_icrs = np.zeros((obs.shape[0], 3), dtype=float)
    for _ in range(max(1, n_iter)):
        jd_ret = obs - tau
        ast_ecl = np.atleast_2d(np.asarray(bary_ecl_at(jd_ret), dtype=float))
        ast_icrs = ecliptic_to_equatorial(ast_ecl)
        tau = np.linalg.norm(ast_icrs - gaia, axis=1) / c_au_per_day
    return obs - tau, ast_icrs


# --- Residuos tangenciales --------------------------------------------------


def tangent_residuals_mas(
    ra_obs_deg: np.ndarray,
    dec_obs_deg: np.ndarray,
    ra_pred_deg: np.ndarray,
    dec_pred_deg: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Residuos tangenciales ``(obs − pred)`` en mas alrededor de la predicción.

    ``dra = (RA_obs − RA_pred)·cos(Dec_pred)`` (envuelto a ±180°) y
    ``ddec = Dec_obs − Dec_pred``.
    """
    ra_obs = np.asarray(ra_obs_deg, dtype=float)
    dec_obs = np.asarray(dec_obs_deg, dtype=float)
    ra_pred = np.asarray(ra_pred_deg, dtype=float)
    dec_pred = np.asarray(dec_pred_deg, dtype=float)
    dra = ((ra_obs - ra_pred + 540.0) % 360.0 - 180.0) * np.cos(dec_pred * _DEG) * _MAS_PER_DEG
    ddec = (dec_obs - dec_pred) * _MAS_PER_DEG
    return dra, ddec


# --- Covarianza anisotrópica along-scan / across-scan -----------------------


def anisotropic_covariance(
    pa_scan_deg: np.ndarray, sigma_al_mas: np.ndarray, sigma_ac_mas: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Componentes ``(Σxx, Σxy, Σyy)`` de la covarianza 2×2 por observación (mas²).

    En el plano tangente ``(x, y) = (RA·cosδ, δ)``:
    ``Σ = σ_AL² û_AL û_ALᵀ + σ_AC² û_AC û_ACᵀ`` con ``û_AL = (sinPA, cosPA)`` y
    ``û_AC = (cosPA, −sinPA)``. El límite ``σ_AC → ∞`` recupera la proyección AL
    pura (la dirección across-scan deja de pesar).
    """
    pa = np.asarray(pa_scan_deg, dtype=float) * _DEG
    s, c = np.sin(pa), np.cos(pa)
    val = np.asarray(sigma_al_mas, dtype=float) ** 2
    vac = np.asarray(sigma_ac_mas, dtype=float) ** 2
    sxx = val * s * s + vac * c * c
    syy = val * c * c + vac * s * s
    sxy = (val - vac) * s * c
    return sxx, sxy, syy


def along_scan_residual(
    dra_mas: np.ndarray,
    ddec_mas: np.ndarray,
    pa_scan_deg: np.ndarray,
    sigma_al_mas: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Proyección along-scan: ``r_AL = dra·sinPA + ddec·cosPA`` y su σ_AL.

    Línea base de una observación de un solo tránsito: conserva solo la componente
    along-scan (la bien medida) y descarta la across-scan.
    """
    pa = np.asarray(pa_scan_deg, dtype=float) * _DEG
    r_al = np.asarray(dra_mas, dtype=float) * np.sin(pa) + np.asarray(
        ddec_mas, dtype=float
    ) * np.cos(pa)
    return r_al, np.asarray(sigma_al_mas, dtype=float)


def whiten_residuals_2d(
    dra_mas: np.ndarray,
    ddec_mas: np.ndarray,
    pa_scan_deg: np.ndarray,
    sigma_al_mas: np.ndarray,
    sigma_ac_mas: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Blanquea el residuo tangencial 2D con la covarianza anisotrópica AL/AC.

    Devuelve ``(whitened, chi2_per_obs)`` con ``whitened`` ``(2N,)`` intercalado
    ``(r1, r2)`` por observación tal que ``Σ r1²+r2² == Σ chi2_per_obs`` =
    Mahalanobis total ``δᵀ Σ⁻¹ δ``. Mantiene la componente across-scan con su σ_AC
    grande (peso ~0) en vez de descartarla, que es el modelo de error correcto.
    """
    dra = np.asarray(dra_mas, dtype=float)
    ddec = np.asarray(ddec_mas, dtype=float)
    sxx, sxy, syy = anisotropic_covariance(pa_scan_deg, sigma_al_mas, sigma_ac_mas)

    det = sxx * syy - sxy * sxy
    fallback = det < _DET_FLOOR
    safe_det = np.where(fallback, 1.0, det)
    inv_xx = np.where(fallback, 1.0 / np.maximum(sxx, _DET_FLOOR), syy / safe_det)
    inv_xy = np.where(fallback, 0.0, -sxy / safe_det)
    inv_yy = np.where(fallback, 1.0 / np.maximum(syy, _DET_FLOOR), sxx / safe_det)

    chi2 = dra * dra * inv_xx + 2.0 * dra * ddec * inv_xy + ddec * ddec * inv_yy

    # Cholesky superior de Σ⁻¹ = [[ixx, ixy],[ixy, iyy]]: r1²+r2² == chi2.
    l11 = np.sqrt(np.maximum(inv_xx, _DET_FLOOR))
    l12 = inv_xy / l11
    l22 = np.sqrt(np.maximum(inv_yy - l12 * l12, _DET_FLOOR))
    r1 = l11 * dra + l12 * ddec
    r2 = l22 * ddec
    whitened = np.empty(2 * dra.shape[0], dtype=float)
    whitened[0::2] = r1
    whitened[1::2] = r2
    return whitened, chi2


# --- Jacobiano de la observación respecto a parámetros del modelo -----------


def radec_jacobian_wrt_position(
    ast_bary_icrs: np.ndarray, gaia_bary_icrs: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """``∂RA/∂ρ`` y ``∂Dec/∂ρ`` (rad/AU) para la línea de visión ICRS ``ρ = r_ast − r_gaia``.

    Devuelve dos arrays ``(N, 3)``: las parciales de RA y Dec (radianes) respecto a
    la posición ICRS del asteroide.
    """
    rho = np.asarray(ast_bary_icrs, dtype=float) - np.asarray(gaia_bary_icrs, dtype=float)
    x, y, z = rho[:, 0], rho[:, 1], rho[:, 2]
    rxy2 = x * x + y * y
    rxy = np.sqrt(rxy2)
    r2 = rxy2 + z * z
    zeros = np.zeros_like(x)
    dra = np.stack([-y / rxy2, x / rxy2, zeros], axis=1)
    ddec = np.stack([-x * z / (r2 * rxy), -y * z / (r2 * rxy), rxy / r2], axis=1)
    return dra, ddec


def tangent_partials_mas(
    ast_bary_icrs: np.ndarray,
    gaia_bary_icrs: np.ndarray,
    dpos_ecl_dparam: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Parciales del plano tangente **predicho** (mas) respecto a K parámetros.

    Dada ``∂r_ast_ecl/∂param`` ``(N, 3, K)`` (parciales de la posición baricéntrica
    **eclíptica** del asteroide respecto a los parámetros — p. ej. de las
    ecuaciones variacionales), aplica la cadena
    ``∂(RA,Dec)/∂ρ · R_ecl→icrs · ∂r_ecl/∂param`` y escala a mas:

    Returns
    -------
    (dra_pred_mas, ddec_pred_mas)
        ``(N, K)`` cada uno: ``∂(RA·cosδ)/∂param`` y ``∂Dec/∂param`` de la posición
        **predicha**, en mas. El residuo es ``obs − pred``, así que su Jacobiano es
        el **negativo** de estos (lo aplica quien arma el ajuste).
    """
    rot_e2i = rotation_x(OBLIQUITY_J2000_RAD)
    dra_drho, ddec_drho = radec_jacobian_wrt_position(ast_bary_icrs, gaia_bary_icrs)
    dparam = np.asarray(dpos_ecl_dparam, dtype=float)
    drho_dparam = np.einsum("ij,njk->nik", rot_e2i, dparam)  # (N, 3, K)
    dra_rad = np.einsum("ni,nik->nk", dra_drho, drho_dparam)  # (N, K)
    ddec_rad = np.einsum("ni,nik->nk", ddec_drho, drho_dparam)
    _ra, dec = radec_from_positions(ast_bary_icrs, gaia_bary_icrs)
    cos_dec = np.cos(np.radians(dec))[:, None]
    return _MAS_PER_RAD * cos_dec * dra_rad, _MAS_PER_RAD * ddec_rad


def along_scan_jacobian(
    ast_bary_icrs: np.ndarray,
    gaia_bary_icrs: np.ndarray,
    dpos_ecl_dparam: np.ndarray,
    pa_scan_deg: np.ndarray,
    sigma_al_mas: np.ndarray,
) -> np.ndarray:
    """Jacobiano del residuo along-scan **blanqueado** respecto a K parámetros.

    ``∂(r_AL/σ_AL)/∂param = −(sinPA·∂dra_pred + cosPA·∂ddec_pred) / σ_AL`` ``(N, K)``
    (signo negativo porque el residuo es ``obs − pred``).
    """
    dra_pred, ddec_pred = tangent_partials_mas(ast_bary_icrs, gaia_bary_icrs, dpos_ecl_dparam)
    pa = np.asarray(pa_scan_deg, dtype=float) * _DEG
    sig = np.asarray(sigma_al_mas, dtype=float)[:, None]
    dr_al = np.sin(pa)[:, None] * dra_pred + np.cos(pa)[:, None] * ddec_pred
    return np.asarray(-dr_al / sig, dtype=float)


# --- Conveniencia: predicción RA/Dec con dinámica N-cuerpos ------------------


def predict_radec(
    test_elements,
    epoch_jd_tdb: float,
    obs_jd_tdb: np.ndarray,
    gaia_bary_icrs: np.ndarray,
    *,
    perturbers=None,
    asteroid_perturbers=(),
    integrator: str = "ias15",
    dt_days: float = 1.0,
    n_lighttime_iter: int = 3,
) -> tuple[np.ndarray, np.ndarray]:
    """Predice ``(RA_deg, Dec_deg)`` en cada época de observación con la dinámica.

    Envuelve :func:`orbdet.dynamics.propagate` como ``bary_ecl_at`` para la
    light-time. ``gaia_bary_icrs`` ``(N, 3)`` es la posición baricéntrica ICRS de
    Gaia en cada época. **Depende de:** ``dynamics`` (rebound).
    """
    from .dynamics import DEFAULT_PERTURBERS, propagate

    pert = DEFAULT_PERTURBERS if perturbers is None else perturbers

    def bary_ecl_at(jd_array: np.ndarray) -> np.ndarray:
        # return_velocity=False → propagate devuelve solo posiciones (ndarray).
        pos = propagate(
            test_elements,
            epoch_jd_tdb,
            np.atleast_1d(jd_array),
            perturbers=pert,
            integrator=integrator,
            dt_days=dt_days,
            asteroid_perturbers=asteroid_perturbers,
        )
        return np.asarray(pos, dtype=float)

    _jd_ret, ast_icrs = light_time_correct(
        bary_ecl_at, obs_jd_tdb, gaia_bary_icrs, n_iter=n_lighttime_iter
    )
    return radec_from_positions(ast_icrs, gaia_bary_icrs)


__all__ = [
    "xyz_to_radec",
    "radec_to_unit_vec",
    "radec_from_positions",
    "light_time_correct",
    "tangent_residuals_mas",
    "anisotropic_covariance",
    "along_scan_residual",
    "whiten_residuals_2d",
    "radec_jacobian_wrt_position",
    "tangent_partials_mas",
    "along_scan_jacobian",
    "predict_radec",
]
