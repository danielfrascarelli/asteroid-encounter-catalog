"""Ajuste conjunto órbita + masa de un perturbador (Fase 1, T6).

Agrega el ``GM`` (vía la masa en M_sun) del perturbador al **mismo** sistema de
mínimos cuadrados que los 6 elementos del objetivo, resuelto sobre el arco
completo. Esto es lo que rompe la degeneración masa↔drift que hundió al LOO
secuencial: la información que distingue una deflexión real de un error de órbita
creciente se conserva dentro de la covarianza conjunta en vez de descartarse al
separar los pasos.

Vector de parámetros: ``x = [mass_msun, a, e, i, Ω, ω, M]`` (7). El perturbador se
modela con su órbita **conocida** (fija) y masa libre; el objetivo con sus 6
elementos libres. El Jacobiano combina:

- ``∂r/∂elementos`` de las ecuaciones variacionales bajo la dinámica que **incluye
  al perturbador** (``variational.partials_wrt_elements``), y
- ``∂r/∂GM`` por diferencias finitas centrales (``variational.partial_wrt_gm``),
  escalado a ``∂r/∂mass = GM_SUN · ∂r/∂GM``.

Gate (closing-loop): inyectar una masa en datos sintéticos y recuperarla con
ratio ≈ 1.0 y σ realista.
"""

from __future__ import annotations

import numpy as np

from .constants import GM_SUN
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
from .variational import partial_wrt_gm, partials_wrt_elements


def determine_mass_and_orbit(
    initial_elements: KeplerElements,
    initial_mass_msun: float,
    perturber_elements: KeplerElements,
    epoch_jd_tdb: float,
    obs_jd_tdb: np.ndarray,
    ra_obs_deg: np.ndarray,
    dec_obs_deg: np.ndarray,
    pa_scan_deg: np.ndarray,
    sigma_al_mas: np.ndarray,
    gaia_bary_icrs: np.ndarray,
    *,
    perturber_name: str = "perturber",
    background_perturbers: tuple[AsteroidPerturber, ...] = (),
    perturbers: tuple[str, ...] = DEFAULT_PERTURBERS,
    integrator: str = "ias15",
    dt_days: float = 1.0,
    n_lighttime_iter: int = 3,
    gm_rel_delta: float = 1e-3,
    max_iter: int = 80,
    **lm_kwargs,
) -> tuple[float, KeplerElements, LeastSquaresResult]:
    """Ajuste conjunto de la masa del perturbador y los 6 elementos del objetivo.

    Parameters
    ----------
    initial_elements, initial_mass_msun:
        Semillas del objetivo (elementos) y del perturbador (masa en M_sun).
    perturber_elements:
        Órbita **conocida** del perturbador (fija durante el ajuste).
    epoch_jd_tdb, obs_jd_tdb, ra_obs_deg, dec_obs_deg, pa_scan_deg, sigma_al_mas, gaia_bary_icrs:
        Como en :func:`orbdet.orbit_determination.determine_orbit` — observaciones
        del **objetivo** (el cuerpo deflectado).
    background_perturbers:
        Otros asteroides masivos de fondo con masa fija (no se ajustan).
    gm_rel_delta:
        Paso relativo de masa de la diferencia finita de ``∂r/∂GM``.
    max_iter, **lm_kwargs:
        Pasados a :func:`orbdet.least_squares.levenberg_marquardt`.

    Returns
    -------
    (mass_msun, fitted_elements, result)
        Masa ajustada (M_sun), elementos del objetivo y el ``LeastSquaresResult``.
        ``result.covariance`` es la covarianza conjunta 7×7 (índice 0 = masa).
    """
    obs_jd = np.atleast_1d(np.asarray(obs_jd_tdb, dtype=float))
    ra_obs = np.asarray(ra_obs_deg, dtype=float)
    dec_obs = np.asarray(dec_obs_deg, dtype=float)
    pa = np.asarray(pa_scan_deg, dtype=float)
    sigma_al = np.asarray(sigma_al_mas, dtype=float)
    gaia = np.asarray(gaia_bary_icrs, dtype=float)

    def residual_and_jac(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        mass = float(x[0])
        el = KeplerElements(*x[1:])
        studied = AsteroidPerturber(perturber_name, mass, perturber_elements)
        ast_perts = (studied, *background_perturbers)

        def bary_ecl_at(jd: np.ndarray) -> np.ndarray:
            return np.asarray(
                propagate(
                    el,
                    epoch_jd_tdb,
                    np.atleast_1d(jd),
                    perturbers=perturbers,
                    integrator=integrator,
                    dt_days=dt_days,
                    asteroid_perturbers=ast_perts,
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
            asteroid_perturbers=ast_perts,
        )
        dgm = partial_wrt_gm(
            el,
            epoch_jd_tdb,
            jd_ret,
            perturber_index=0,
            perturbers=perturbers,
            integrator=integrator,
            dt_days=dt_days,
            asteroid_perturbers=ast_perts,
            rel_delta=gm_rel_delta,
        )
        # ∂r_ecl/∂mass = GM_SUN · ∂r_ecl/∂GM (posición = primeras 3 componentes).
        dpos_dmass = dgm[:, 0:3] * GM_SUN  # (N, 3)

        ast_icrs = ecliptic_to_equatorial(pos)
        ra_pred, dec_pred = radec_from_positions(ast_icrs, gaia)
        dra, ddec = tangent_residuals_mas(ra_obs, dec_obs, ra_pred, dec_pred)
        r_al, sig = along_scan_residual(dra, ddec, pa, sigma_al)
        resid = r_al / sig

        # ∂r_ecl/∂param para los 7 parámetros: [masa | 6 elementos].
        dpos_dparam = np.concatenate([dpos_dmass[:, :, None], dstate[:, 0:3, :]], axis=2)
        jac = along_scan_jacobian(ast_icrs, gaia, dpos_dparam, pa, sigma_al)
        return resid, jac

    x0 = np.concatenate([[float(initial_mass_msun)], np.asarray(initial_elements.as_array())])
    result = levenberg_marquardt(residual_and_jac, x0, max_iter=max_iter, **lm_kwargs)
    return float(result.x[0]), KeplerElements(*result.x[1:]), result


__all__ = ["determine_mass_and_orbit"]
