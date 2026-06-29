"""Ajuste conjunto órbita + masa de un perturbador (Fase 1, T6/T7).

Agrega el ``GM`` (vía la masa en M_sun) del perturbador al **mismo** sistema de
mínimos cuadrados que los elementos del/los objetivo(s), resuelto sobre el arco
completo. Esto es lo que rompe la degeneración masa↔drift que hundió al LOO
secuencial: la información que distingue una deflexión real de un error de órbita
creciente se conserva dentro de la covarianza conjunta en vez de descartarse al
separar los pasos.

- **T6 — un objetivo** (:func:`determine_mass_and_orbit`): vector
  ``x = [mass, a, e, i, Ω, ω, M]`` (7). Gate: closing-loop (masa inyectada
  recuperada a ratio ≈ 1.0).
- **T7 — stacking multi-objetivo** (:func:`determine_shared_mass`): ``GM``
  compartido entre ``N`` objetivos, 6 elementos por objetivo, en un único sistema
  de ``1 + 6N`` parámetros con Jacobiano en flecha (la columna de masa es densa,
  los bloques de elementos son diagonales por objetivo). Gate: ``σ(GM) ∝ 1/√N``,
  rompiendo la multimodalidad de un solo encuentro.

El Jacobiano combina ``∂r/∂elementos`` (ecuaciones variacionales bajo la dinámica
que incluye al perturbador) y ``∂r/∂mass = GM_SUN · ∂r/∂GM`` (diferencias finitas).
"""

from __future__ import annotations

from dataclasses import dataclass

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


@dataclass(frozen=True)
class TargetObservations:
    """Observaciones de un objetivo (cuerpo deflectado) para el ajuste de masa.

    Todas las arrays tienen longitud ``N`` (número de observaciones del objetivo);
    ``gaia_bary_icrs`` es ``(N, 3)``. ``initial_elements`` es la semilla orbital
    del objetivo en la época común del ajuste.
    """

    initial_elements: KeplerElements
    obs_jd_tdb: np.ndarray
    ra_obs_deg: np.ndarray
    dec_obs_deg: np.ndarray
    pa_scan_deg: np.ndarray
    sigma_al_mas: np.ndarray
    gaia_bary_icrs: np.ndarray


@dataclass(frozen=True)
class _ModelConfig:
    """Configuración de dinámica/observación compartida por el ajuste."""

    epoch_jd_tdb: float
    perturber_elements: KeplerElements
    perturber_name: str
    background_perturbers: tuple[AsteroidPerturber, ...]
    perturbers: tuple[str, ...]
    integrator: str
    dt_days: float
    n_lighttime_iter: int
    gm_rel_delta: float
    backend: str = "rebound"
    gr: bool = True
    elem_fd_steps: tuple[float, ...] = (1e-7, 1e-7, 1e-7, 1e-7, 1e-7, 1e-7)


# Pasos de diferencias finitas para ∂pos/∂elementos en el backend ASSIST
# (donde la STM analítica de rebound no captura las fuerzas de la efeméride).
# Para ``a`` el paso es relativo (·a); para e/i/Ω/ω/M es absoluto.


def _assist_positions(
    el: KeplerElements, mass: float, jd: np.ndarray, cfg: _ModelConfig
) -> np.ndarray:
    """Posición baricéntrica eclíptica ``(N, 3)`` del objetivo con el backend ASSIST.

    ``mass`` es la masa (M_sun) del perturbador bajo estudio (índice 0); el fondo va
    en ``cfg.background_perturbers``.
    """
    from .dynamics_assist import propagate_assist

    studied = AsteroidPerturber(cfg.perturber_name, mass, cfg.perturber_elements)
    ast_perts = (studied, *cfg.background_perturbers)
    return np.atleast_2d(
        np.asarray(
            propagate_assist(
                el, cfg.epoch_jd_tdb, np.atleast_1d(jd), asteroid_perturbers=ast_perts, gr=cfg.gr
            ),
            dtype=float,
        )
    )


def _assist_pos_and_partials(
    el: KeplerElements, mass: float, jd_ret: np.ndarray, cfg: _ModelConfig
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(pos, ∂pos/∂mass, ∂pos/∂elementos)`` por diferencias finitas centrales (ASSIST).

    ``pos`` ``(N,3)``, ``∂pos/∂mass`` ``(N,3)``, ``∂pos/∂elementos`` ``(N,3,6)``. Bajo
    ASSIST las fuerzas (Sol/planetas/GR) vienen de la efeméride y rebound no propaga
    las partículas variacionales a través de ellas, por lo que TODAS las parciales se
    obtienen por FD central sobre :func:`propagate_assist` (exactas a O(δ²); el paso de
    masa usa ``gm_rel_delta`` igual que el backend rebound).
    """
    pos0 = _assist_positions(el, mass, jd_ret, cfg)

    dm = cfg.gm_rel_delta * mass
    dpos_dmass = (
        _assist_positions(el, mass + dm, jd_ret, cfg)
        - _assist_positions(el, mass - dm, jd_ret, cfg)
    ) / (2.0 * dm)

    x = np.asarray(el.as_array(), dtype=float)
    steps = np.asarray(cfg.elem_fd_steps, dtype=float).copy()
    steps[0] = steps[0] * max(abs(x[0]), 1.0)  # paso relativo para a
    dstate = np.zeros((pos0.shape[0], 3, 6), dtype=float)
    for k in range(6):
        xp = x.copy()
        xp[k] += steps[k]
        xm = x.copy()
        xm[k] -= steps[k]
        pp = _assist_positions(KeplerElements(*xp), mass, jd_ret, cfg)
        pm = _assist_positions(KeplerElements(*xm), mass, jd_ret, cfg)
        dstate[:, :, k] = (pp - pm) / (2.0 * steps[k])
    return pos0, dpos_dmass, dstate


def _target_resid_and_blocks(
    mass: float, el: KeplerElements, tobs: TargetObservations, cfg: _ModelConfig
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Residuo blanqueado y bloques del Jacobiano de un objetivo.

    Returns
    -------
    (resid, jac_mass, jac_elem)
        ``resid`` ``(N,)``; ``jac_mass`` ``(N,)`` = ∂resid/∂mass; ``jac_elem``
        ``(N, 6)`` = ∂resid/∂elementos.
    """
    obs_jd = np.atleast_1d(np.asarray(tobs.obs_jd_tdb, dtype=float))
    ra_obs = np.asarray(tobs.ra_obs_deg, dtype=float)
    dec_obs = np.asarray(tobs.dec_obs_deg, dtype=float)
    pa = np.asarray(tobs.pa_scan_deg, dtype=float)
    sigma_al = np.asarray(tobs.sigma_al_mas, dtype=float)
    gaia = np.asarray(tobs.gaia_bary_icrs, dtype=float)

    studied = AsteroidPerturber(cfg.perturber_name, mass, cfg.perturber_elements)
    ast_perts = (studied, *cfg.background_perturbers)

    if cfg.backend == "assist":

        def bary_ecl_at(jd: np.ndarray) -> np.ndarray:
            return _assist_positions(el, mass, np.atleast_1d(jd), cfg)

        jd_ret, _ = light_time_correct(bary_ecl_at, obs_jd, gaia, n_iter=cfg.n_lighttime_iter)
        pos, dpos_dmass, dstate = _assist_pos_and_partials(el, mass, jd_ret, cfg)
    else:

        def bary_ecl_at(jd: np.ndarray) -> np.ndarray:
            return np.asarray(
                propagate(
                    el,
                    cfg.epoch_jd_tdb,
                    np.atleast_1d(jd),
                    perturbers=cfg.perturbers,
                    integrator=cfg.integrator,
                    dt_days=cfg.dt_days,
                    asteroid_perturbers=ast_perts,
                ),
                dtype=float,
            )

        jd_ret, _ = light_time_correct(bary_ecl_at, obs_jd, gaia, n_iter=cfg.n_lighttime_iter)
        pos, _vel, dstate = partials_wrt_elements(
            el,
            cfg.epoch_jd_tdb,
            jd_ret,
            perturbers=cfg.perturbers,
            integrator=cfg.integrator,
            dt_days=cfg.dt_days,
            asteroid_perturbers=ast_perts,
        )
        dgm = partial_wrt_gm(
            el,
            cfg.epoch_jd_tdb,
            jd_ret,
            perturber_index=0,
            perturbers=cfg.perturbers,
            integrator=cfg.integrator,
            dt_days=cfg.dt_days,
            asteroid_perturbers=ast_perts,
            rel_delta=cfg.gm_rel_delta,
        )
        dpos_dmass = dgm[:, 0:3] * GM_SUN  # ∂r_ecl/∂mass = GM_SUN·∂r_ecl/∂GM

    ast_icrs = ecliptic_to_equatorial(pos)
    ra_pred, dec_pred = radec_from_positions(ast_icrs, gaia)
    dra, ddec = tangent_residuals_mas(ra_obs, dec_obs, ra_pred, dec_pred)
    r_al, sig = along_scan_residual(dra, ddec, pa, sigma_al)
    resid = r_al / sig

    dpos_dparam = np.concatenate([dpos_dmass[:, :, None], dstate[:, 0:3, :]], axis=2)  # (N,3,7)
    jac_all = along_scan_jacobian(ast_icrs, gaia, dpos_dparam, pa, sigma_al)  # (N, 7)
    return resid, jac_all[:, 0], jac_all[:, 1:]


def _make_config(
    epoch_jd_tdb: float,
    perturber_elements: KeplerElements,
    perturber_name: str,
    background_perturbers: tuple[AsteroidPerturber, ...],
    perturbers: tuple[str, ...],
    integrator: str,
    dt_days: float,
    n_lighttime_iter: int,
    gm_rel_delta: float,
    backend: str = "rebound",
    gr: bool = True,
) -> _ModelConfig:
    return _ModelConfig(
        epoch_jd_tdb=epoch_jd_tdb,
        perturber_elements=perturber_elements,
        perturber_name=perturber_name,
        background_perturbers=background_perturbers,
        perturbers=perturbers,
        integrator=integrator,
        dt_days=dt_days,
        n_lighttime_iter=n_lighttime_iter,
        gm_rel_delta=gm_rel_delta,
        backend=backend,
        gr=gr,
    )


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
    backend: str = "rebound",
    gr: bool = True,
    max_iter: int = 80,
    **lm_kwargs,
) -> tuple[float, KeplerElements, LeastSquaresResult]:
    """Ajuste conjunto de la masa del perturbador y los 6 elementos de un objetivo.

    ``backend="assist"`` usa el modelo de fuerzas DE440 + GR + perturbadores (T8);
    en ese caso ``background_perturbers`` debe traer los 15 asteroides grandes
    restantes (ver :func:`orbdet.dynamics_assist.big_asteroid_perturbers`).

    Returns
    -------
    (mass_msun, fitted_elements, result)
        ``result.covariance`` es la covarianza conjunta 7×7 (índice 0 = masa).
    """
    cfg = _make_config(
        epoch_jd_tdb,
        perturber_elements,
        perturber_name,
        background_perturbers,
        perturbers,
        integrator,
        dt_days,
        n_lighttime_iter,
        gm_rel_delta,
        backend=backend,
        gr=gr,
    )
    tobs = TargetObservations(
        initial_elements=initial_elements,
        obs_jd_tdb=obs_jd_tdb,
        ra_obs_deg=ra_obs_deg,
        dec_obs_deg=dec_obs_deg,
        pa_scan_deg=pa_scan_deg,
        sigma_al_mas=sigma_al_mas,
        gaia_bary_icrs=gaia_bary_icrs,
    )

    def residual_and_jac(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        resid, jac_mass, jac_elem = _target_resid_and_blocks(
            float(x[0]), KeplerElements(*x[1:]), tobs, cfg
        )
        return resid, np.column_stack([jac_mass, jac_elem])

    x0 = np.concatenate([[float(initial_mass_msun)], np.asarray(initial_elements.as_array())])
    result = levenberg_marquardt(residual_and_jac, x0, max_iter=max_iter, **lm_kwargs)
    return float(result.x[0]), KeplerElements(*result.x[1:]), result


def determine_shared_mass(
    targets: list[TargetObservations],
    initial_mass_msun: float,
    perturber_elements: KeplerElements,
    epoch_jd_tdb: float,
    *,
    perturber_name: str = "perturber",
    background_perturbers: tuple[AsteroidPerturber, ...] = (),
    perturbers: tuple[str, ...] = DEFAULT_PERTURBERS,
    integrator: str = "ias15",
    dt_days: float = 1.0,
    n_lighttime_iter: int = 3,
    gm_rel_delta: float = 1e-3,
    backend: str = "rebound",
    gr: bool = True,
    max_iter: int = 80,
    **lm_kwargs,
) -> tuple[float, list[KeplerElements], LeastSquaresResult]:
    """Stacking multi-objetivo: una masa de perturbador compartida por ``N`` objetivos.

    Vector de ``1 + 6N`` parámetros ``[mass, elementos_1, …, elementos_N]``. El
    Jacobiano es en flecha: la columna de masa es densa y cada bloque de 6
    elementos solo afecta a las observaciones de su objetivo. Al apilar objetivos
    independientes, la información de Fisher de la masa se suma y ``σ(GM) ∝ 1/√N``.

    Returns
    -------
    (mass_msun, fitted_elements_per_target, result)
        ``result.covariance`` es ``(1+6N, 1+6N)`` (índice 0 = masa compartida).
    """
    if not targets:
        raise ValueError("determine_shared_mass requiere al menos un objetivo")
    cfg = _make_config(
        epoch_jd_tdb,
        perturber_elements,
        perturber_name,
        background_perturbers,
        perturbers,
        integrator,
        dt_days,
        n_lighttime_iter,
        gm_rel_delta,
        backend=backend,
        gr=gr,
    )
    n_t = len(targets)
    n_par = 1 + 6 * n_t

    def residual_and_jac(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        mass = float(x[0])
        resids: list[np.ndarray] = []
        blocks: list[np.ndarray] = []
        for k, tobs in enumerate(targets):
            el_k = KeplerElements(*x[1 + 6 * k : 1 + 6 * k + 6])
            resid_k, jac_mass_k, jac_elem_k = _target_resid_and_blocks(mass, el_k, tobs, cfg)
            block = np.zeros((resid_k.size, n_par), dtype=float)
            block[:, 0] = jac_mass_k
            block[:, 1 + 6 * k : 1 + 6 * k + 6] = jac_elem_k
            resids.append(resid_k)
            blocks.append(block)
        return np.concatenate(resids), np.vstack(blocks)

    x0 = np.concatenate(
        [[float(initial_mass_msun)]] + [np.asarray(t.initial_elements.as_array()) for t in targets]
    )
    result = levenberg_marquardt(residual_and_jac, x0, max_iter=max_iter, **lm_kwargs)
    fitted = [KeplerElements(*result.x[1 + 6 * k : 1 + 6 * k + 6]) for k in range(n_t)]
    return float(result.x[0]), fitted, result


__all__ = ["TargetObservations", "determine_mass_and_orbit", "determine_shared_mass"]
