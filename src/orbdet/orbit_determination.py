"""Determinación de órbitas por mínimos cuadrados sobre el arco completo.

Ensambla las piezas del motor en un corrector diferencial que ajusta los **6
elementos keplerianos** del objetivo a las observaciones astrométricas de Gaia,
sobre todo el arco a la vez (sin split LOO):

    elementos → dinámica N-cuerpos (``dynamics``) + light-time (``observation``)
              → RA/Dec predichos
    Jacobiano ∂(RA,Dec)/∂elementos  =  ∂(RA,Dec)/∂ρ · R_ecl→icrs · ∂r/∂elementos
              con ∂r/∂elementos de las ecuaciones variacionales (``variational``)
    residuos along-scan blanqueados (``observation``) → Levenberg-Marquardt
              (``least_squares``)

La dependencia de la light-time respecto a los elementos es de segundo orden
(τ ~ 0.01 d) y se omite en el Jacobiano, como es estándar en OD; sí se aplica en
el residuo (la posición se evalúa en el tiempo de emisión retardado).

Este es el gate de T5: recuperar una órbita conocida a partir de observaciones
sintéticas. Con el parámetro de masa agregado al vector de estado, es la base del
ajuste conjunto órbita+masa (T6).
"""

from __future__ import annotations

import numpy as np

from .dynamics import DEFAULT_PERTURBERS, AsteroidPerturber, propagate
from .frames import ecliptic_to_equatorial
from .kepler import KeplerElements
from .least_squares import LeastSquaresResult, levenberg_marquardt
from .observation import (
    along_scan_jacobian,
    along_scan_residual,
    light_time_correct,
    radec_from_positions,
    tangent_residuals_mas,
)
from .variational import partials_wrt_elements


def determine_orbit(
    initial_elements: KeplerElements,
    epoch_jd_tdb: float,
    obs_jd_tdb: np.ndarray,
    ra_obs_deg: np.ndarray,
    dec_obs_deg: np.ndarray,
    pa_scan_deg: np.ndarray,
    sigma_al_mas: np.ndarray,
    gaia_bary_icrs: np.ndarray,
    *,
    perturbers: tuple[str, ...] = DEFAULT_PERTURBERS,
    asteroid_perturbers: tuple[AsteroidPerturber, ...] = (),
    integrator: str = "ias15",
    dt_days: float = 1.0,
    n_lighttime_iter: int = 3,
    max_iter: int = 60,
    **lm_kwargs,
) -> tuple[KeplerElements, LeastSquaresResult]:
    """Ajusta los 6 elementos del objetivo a las observaciones Gaia (along-scan).

    Parameters
    ----------
    initial_elements:
        Semilla del ajuste (elementos heliocéntricos eclípticos en la época).
    epoch_jd_tdb:
        Época de los elementos.
    obs_jd_tdb:
        Épocas de observación ``(N,)`` (JD TDB).
    ra_obs_deg, dec_obs_deg:
        RA/Dec observados ``(N,)`` (grados, frame astrométrico baricéntrico ICRS).
    pa_scan_deg:
        Ángulo de posición de barrido por observación ``(N,)`` (grados).
    sigma_al_mas:
        σ along-scan por observación ``(N,)`` (mas).
    gaia_bary_icrs:
        Posición baricéntrica ICRS de Gaia por observación ``(N, 3)`` (AU).
    perturbers, asteroid_perturbers, integrator, dt_days, n_lighttime_iter:
        Configuración del modelo dinámico/observacional.
    max_iter, **lm_kwargs:
        Pasados a :func:`orbdet.least_squares.levenberg_marquardt`.

    Returns
    -------
    (fitted_elements, result)
        Los elementos ajustados y el :class:`LeastSquaresResult` (covarianza,
        χ²_red, convergencia).
    """
    obs_jd = np.atleast_1d(np.asarray(obs_jd_tdb, dtype=float))
    ra_obs = np.asarray(ra_obs_deg, dtype=float)
    dec_obs = np.asarray(dec_obs_deg, dtype=float)
    pa = np.asarray(pa_scan_deg, dtype=float)
    sigma_al = np.asarray(sigma_al_mas, dtype=float)
    gaia = np.asarray(gaia_bary_icrs, dtype=float)

    def residual_and_jac(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        el = KeplerElements(*x)

        def bary_ecl_at(jd: np.ndarray) -> np.ndarray:
            return np.asarray(
                propagate(
                    el,
                    epoch_jd_tdb,
                    np.atleast_1d(jd),
                    perturbers=perturbers,
                    integrator=integrator,
                    dt_days=dt_days,
                    asteroid_perturbers=asteroid_perturbers,
                ),
                dtype=float,
            )

        jd_ret, _ = light_time_correct(bary_ecl_at, obs_jd, gaia, n_iter=n_lighttime_iter)
        pos, _vel, dstate = partials_wrt_elements(
            el,
            epoch_jd_tdb,
            jd_ret,
            perturbers=perturbers,
            integrator=integrator,
            dt_days=dt_days,
            asteroid_perturbers=asteroid_perturbers,
        )
        ast_icrs = ecliptic_to_equatorial(pos)
        ra_pred, dec_pred = radec_from_positions(ast_icrs, gaia)
        dra, ddec = tangent_residuals_mas(ra_obs, dec_obs, ra_pred, dec_pred)
        r_al, sig = along_scan_residual(dra, ddec, pa, sigma_al)
        resid = r_al / sig

        # Jacobiano del residuo AL blanqueado respecto a los 6 elementos, vía la
        # cadena ∂(RA,Dec)/∂ρ · R_ecl→icrs · ∂r_ecl/∂elementos.
        jac = along_scan_jacobian(ast_icrs, gaia, dstate[:, 0:3, :], pa, sigma_al)
        return resid, jac

    result = levenberg_marquardt(
        residual_and_jac,
        np.asarray(initial_elements.as_array(), dtype=float),
        max_iter=max_iter,
        **lm_kwargs,
    )
    return KeplerElements(*result.x), result


__all__ = ["determine_orbit"]
