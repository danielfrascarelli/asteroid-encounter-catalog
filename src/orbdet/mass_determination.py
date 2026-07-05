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

import dataclasses
import os
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
from .variational import (
    partial_wrt_gm,
    partial_wrt_gm_variational,
    partials_wrt_elements,
)


def _pool_worker_init() -> None:
    """Inicializador de los procesos worker: evita la sobre-suscripción de hilos.

    Cada worker hace integraciones N-cuerpos seriales; si rebound/BLAS abrieran sus
    propios hilos, ``N_workers × N_hilos`` saturaría la CPU. Fijar a 1 hilo por
    worker mantiene el paralelismo a nivel de objetivo (uno por core)."""
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[var] = "1"


def _make_pool(n_workers: int, n_tasks: int):
    """Crea un ``multiprocessing.Pool`` (contexto ``fork``) o ``None`` si serie.

    ``fork`` hace que los workers hereden el proceso padre por copy-on-write, así la
    efeméride ASSIST (~750 MB, mapeada en memoria) se comparte sin recargarse. Se usan
    ``min(n_workers, n_tasks)`` procesos (no tiene sentido más workers que objetivos).
    """
    if not n_workers or n_workers <= 1 or n_tasks <= 1:
        return None
    import multiprocessing as mp

    n = min(int(n_workers), int(n_tasks))
    return mp.get_context("fork").Pool(processes=n, initializer=_pool_worker_init)


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
    fov_group: np.ndarray | None = None
    """Etiqueta entera del cruce de plano focal (FOV transit) de cada observación.

    Gaia mide cada objeto con hasta ~9 CCDs por cruce, separados ~segundos y que
    comparten el mismo error de actitud/centroide → residuos correlacionados. Si se
    provee, el ajuste blanquea con una **covarianza en bloques** (un bloque por
    valor distinto) en vez de tratar cada CCD como independiente, lo que de otro
    modo sobre-cuenta la información y subestima σ(masa). ``None`` → blanqueo
    diagonal clásico (cada observación independiente)."""


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
    gm_variational: bool = True
    """Usar la parcial ∂x/∂GM **analítica** (partícula variacional de masa, F6) en el
    backend ``rebound``, en vez de diferencias finitas. Ahorra dos propagaciones por
    Jacobiano y es exacta a mejor que 1e-6. Ignorado en el backend ``assist`` (las
    fuerzas de la efeméride no propagan variacionales → siempre FD)."""
    sys_floor_mas: float = 0.0
    """Piso de error sistemático along-scan **correlacionado dentro de cada FOV
    transit** (mas). Modela el error de actitud/centroide común a los CCDs de un
    mismo cruce: ``C_bloque = diag(σ_AL²) + sys_floor² · 11ᵀ``. Se calibra para que
    χ²_red ≈ 1. ``0`` → sin piso (covarianza diagonal)."""
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

    # Guardia C6 (tribunal 2026-07-04): si LM llevó la masa a ≤ 0 (p. ej. Davida),
    # `gm_rel_delta * mass` da un paso nulo o negativo → FD degenerada. El paso se
    # toma sobre |masa| con un piso absoluto (~2×10¹⁶ kg); solo afecta el paso de
    # diferenciación, no el modelo.
    dm = cfg.gm_rel_delta * max(abs(mass), 1e-14)
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
        if k == 1 and xm[1] < 0.0:
            # Guardia C6: el paso central cruzaría e < 0 (solve_kepler lanza).
            # FD adelantada de un lado — O(δ) en vez de O(δ²), solo para e ≈ 0.
            dstate[:, :, k] = (pp - pos0) / steps[k]
            continue
        pm = _assist_positions(KeplerElements(*xm), mass, jd_ret, cfg)
        dstate[:, :, k] = (pp - pm) / (2.0 * steps[k])
    return pos0, dpos_dmass, dstate


def _block_whiten(
    r_al: np.ndarray,
    jac_raw: np.ndarray,
    sigma_al: np.ndarray,
    fov_group: np.ndarray | None,
    sys_floor_mas: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Blanquea residuo y Jacobiano along-scan con covarianza en bloques por FOV.

    Modelo de error: dentro de cada cruce FOV (mismo ``fov_group``) los CCDs
    comparten un error sistemático along-scan de varianza ``sys_floor²`` (actitud/
    centroide) además de su ruido independiente ``σ_AL,i``:
    ``C_bloque = diag(σ_AL²) + sys_floor² · 11ᵀ``. Devuelve ``(r_w, jac_w)`` con
    ``r_wᵀr_w = r_alᵀ C⁻¹ r_al`` (whitening de Cholesky por bloque), de modo que la
    covarianza ``(JᵀC⁻¹J)⁻¹`` del ajuste ya es la honesta sin reescalado posterior.

    Sin ``fov_group`` o con ``sys_floor ≤ 0`` colapsa al blanqueo diagonal clásico
    ``r/σ`` (preserva exactamente los gates sintéticos T6/T7).
    """
    sigma_al = np.asarray(sigma_al, dtype=float)
    if sys_floor_mas <= 0.0 or fov_group is None:
        return r_al / sigma_al, jac_raw / sigma_al[:, None]

    s2 = float(sys_floor_mas) ** 2
    fov = np.asarray(fov_group)
    r_w = np.empty_like(r_al)
    jac_w = np.empty_like(jac_raw)
    for g in np.unique(fov):
        idx = np.where(fov == g)[0]
        if idx.size == 1:
            d = np.sqrt(sigma_al[idx[0]] ** 2 + s2)
            r_w[idx[0]] = r_al[idx[0]] / d
            jac_w[idx[0]] = jac_raw[idx[0]] / d
            continue
        cov = np.diag(sigma_al[idx] ** 2) + s2  # diag + s² 11ᵀ
        chol = np.linalg.cholesky(cov)
        r_w[idx] = np.linalg.solve(chol, r_al[idx])
        jac_w[idx] = np.linalg.solve(chol, jac_raw[idx])
    return r_w, jac_w


def _forward_al(
    mass: float, el: KeplerElements, tobs: TargetObservations, cfg: _ModelConfig
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Modelo directo: residuo y Jacobiano along-scan **crudos** (sin blanquear).

    Returns
    -------
    (r_al, jac_raw, sigma_al)
        ``r_al`` ``(N,)`` residuo along-scan en mas; ``jac_raw`` ``(N, 7)`` =
        ∂r_al/∂[mass, elementos] en mas/parámetro; ``sigma_al`` ``(N,)`` en mas.
        Aquí está toda la integración N-cuerpos; el blanqueo (diagonal o en bloques)
        se aplica aparte, de modo que el calibrador del piso pueda reponderar sin
        re-integrar.
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
        if cfg.gm_variational:
            # F6: parcial ∂x/∂GM analítica (una propagación por sentido) en vez de FD.
            dgm = partial_wrt_gm_variational(
                el,
                cfg.epoch_jd_tdb,
                jd_ret,
                perturber_index=0,
                perturbers=cfg.perturbers,
                integrator=cfg.integrator,
                dt_days=cfg.dt_days,
                asteroid_perturbers=ast_perts,
            )
        else:
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
    r_al, _sig = along_scan_residual(dra, ddec, pa, sigma_al)  # residuo crudo (mas)

    dpos_dparam = np.concatenate([dpos_dmass[:, :, None], dstate[:, 0:3, :]], axis=2)  # (N,3,7)
    # Jacobiano crudo (mas/parámetro): σ=1 evita el blanqueo diagonal interno.
    jac_raw = along_scan_jacobian(ast_icrs, gaia, dpos_dparam, pa, np.ones_like(sigma_al))  # (N,7)
    return r_al, jac_raw, sigma_al


def _forward_al_rsig(
    mass: float, el: KeplerElements, tobs: TargetObservations, cfg: _ModelConfig
) -> tuple[np.ndarray, np.ndarray]:
    """``(r_al, sigma_al)`` crudos de un objetivo (worker picklable para el pool de
    :func:`calibrate_sys_floor`)."""
    r_al, _jac_raw, sigma_al = _forward_al(mass, el, tobs, cfg)
    return r_al, sigma_al


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
    r_al, jac_raw, sigma_al = _forward_al(mass, el, tobs, cfg)
    # Blanqueo con covarianza en bloques por FOV (o diagonal si no hay grupos/piso).
    resid, jac_all = _block_whiten(r_al, jac_raw, sigma_al, tobs.fov_group, cfg.sys_floor_mas)
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
    sys_floor_mas: float = 0.0,
    gm_variational: bool = True,
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
        sys_floor_mas=sys_floor_mas,
        gm_variational=gm_variational,
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
    sys_floor_mas: float = 0.0,
    gm_variational: bool = True,
    max_iter: int = 80,
    **lm_kwargs,
) -> tuple[float, KeplerElements, LeastSquaresResult]:
    """Ajuste conjunto de la masa del perturbador y los 6 elementos de un objetivo.

    .. warning::
        El default ``backend="rebound"`` usa las efemérides ``builtin`` de astropy
        (series analíticas, ~km de error planetario — NO DE440). Para resultados de
        producción usar ``backend="assist"`` explícitamente: modelo de fuerzas
        DE440 + GR + perturbadores (T8); en ese caso ``background_perturbers`` debe
        traer los 15 asteroides grandes restantes (ver
        :func:`orbdet.dynamics_assist.big_asteroid_perturbers`).

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
        sys_floor_mas=sys_floor_mas,
        gm_variational=gm_variational,
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
    sys_floor_mas: float = 0.0,
    gm_variational: bool = True,
    max_iter: int = 80,
    n_workers: int = 1,
    **lm_kwargs,
) -> tuple[float, list[KeplerElements], LeastSquaresResult]:
    """Stacking multi-objetivo: una masa de perturbador compartida por ``N`` objetivos.

    .. warning::
        El default ``backend="rebound"`` usa las efemérides ``builtin`` de astropy
        (~km de error planetario, NO DE440); producción usa ``backend="assist"``.

    Vector de ``1 + 6N`` parámetros ``[mass, elementos_1, …, elementos_N]``. El
    Jacobiano es en flecha: la columna de masa es densa y cada bloque de 6
    elementos solo afecta a las observaciones de su objetivo. Al apilar objetivos
    independientes, la información de Fisher de la masa se suma y ``σ(GM) ∝ 1/√N``.

    ``n_workers>1`` evalúa los ``N`` objetivos (cada uno una integración N-cuerpos +
    su Jacobiano por diferencias finitas, ~15 propagaciones) **en paralelo** sobre un
    ``multiprocessing.Pool`` persistente (un proceso por core; la efeméride ASSIST se
    mapea en memoria compartida por COW de ``fork``). Los resultados se ensamblan en
    orden fijo ⇒ idénticos al modo serie. Es el cuello de botella del ajuste.

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
        sys_floor_mas=sys_floor_mas,
        gm_variational=gm_variational,
    )
    n_t = len(targets)
    n_par = 1 + 6 * n_t

    def _assemble(results) -> tuple[np.ndarray, np.ndarray]:
        resids: list[np.ndarray] = []
        blocks: list[np.ndarray] = []
        for k, (resid_k, jac_mass_k, jac_elem_k) in enumerate(results):
            block = np.zeros((resid_k.size, n_par), dtype=float)
            block[:, 0] = jac_mass_k
            block[:, 1 + 6 * k : 1 + 6 * k + 6] = jac_elem_k
            resids.append(resid_k)
            blocks.append(block)
        return np.concatenate(resids), np.vstack(blocks)

    pool = _make_pool(n_workers, n_t)
    try:

        def residual_and_jac(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
            mass = float(x[0])
            args = [
                (mass, KeplerElements(*x[1 + 6 * k : 1 + 6 * k + 6]), targets[k], cfg)
                for k in range(n_t)
            ]
            if pool is not None:
                results = pool.starmap(_target_resid_and_blocks, args)
            else:
                results = [_target_resid_and_blocks(*a) for a in args]
            return _assemble(results)

        x0 = np.concatenate(
            [[float(initial_mass_msun)]]
            + [np.asarray(t.initial_elements.as_array()) for t in targets]
        )
        result = levenberg_marquardt(residual_and_jac, x0, max_iter=max_iter, **lm_kwargs)
    finally:
        if pool is not None:
            pool.close()
            pool.join()
    fitted = [KeplerElements(*result.x[1 + 6 * k : 1 + 6 * k + 6]) for k in range(n_t)]
    return float(result.x[0]), fitted, result


def calibrate_sys_floor(
    targets: list[TargetObservations],
    mass_msun: float,
    elements_list: list[KeplerElements],
    epoch_jd_tdb: float,
    *,
    perturber_elements: KeplerElements,
    perturber_name: str = "perturber",
    background_perturbers: tuple[AsteroidPerturber, ...] = (),
    backend: str = "assist",
    gr: bool = True,
    n_params: int | None = None,
    max_floor_mas: float = 30.0,
    n_bisect: int = 40,
    n_workers: int = 1,
) -> tuple[float, float]:
    """Calibra el piso sistemático ``s_c`` para que ``χ²_red ≈ 1`` con bloques por FOV.

    Evalúa el modelo directo **una sola vez** por objetivo en la solución convergida
    (``mass_msun``, ``elements_list``), cachea los residuos along-scan crudos y, sobre
    ese caché, busca por bisección el ``s_c`` tal que la χ² blanqueada en bloques
    iguale los grados de libertad (``N_obs − n_params``). Como ``s_c`` apenas mueve el
    punto de mínimos cuadrados (sobre todo reescala σ), calibrar en la solución fija es
    una aproximación excelente; conviene re-ajustar una vez con el ``s_c`` resultante.

    Returns
    -------
    (sys_floor_mas, chi2_red_at_floor)
        Si los datos ya tienen ``χ²_red ≤ 1`` sin piso, devuelve ``(0.0, χ²_red₀)``.
    """
    cfg0 = _make_config(
        epoch_jd_tdb,
        perturber_elements,
        perturber_name,
        background_perturbers,
        DEFAULT_PERTURBERS,
        "ias15",
        1.0,
        3,
        1e-3,
        backend=backend,
        gr=gr,
        sys_floor_mas=0.0,
    )
    pool = _make_pool(n_workers, len(targets))
    try:
        args = [(mass_msun, el, t, cfg0) for t, el in zip(targets, elements_list)]
        if pool is not None:
            rsig = pool.starmap(_forward_al_rsig, args)
        else:
            rsig = [_forward_al_rsig(*a) for a in args]
    finally:
        if pool is not None:
            pool.close()
            pool.join()
    comps = [
        (r_al, sig, np.asarray(t.fov_group)) for (r_al, sig), t in zip(rsig, targets)
    ]  # lista de (r_al, sigma_al, fov_group)
    n_obs = int(sum(r.size for r, _, _ in comps))
    n_par = n_params if n_params is not None else 1 + 6 * len(targets)
    dof = max(n_obs - n_par, 1)

    def chi2_red(s: float) -> float:
        tot = 0.0
        for r_al, sig, fov in comps:
            rw, _ = _block_whiten(r_al, np.zeros((r_al.size, 1)), sig, fov, s)
            tot += float(rw @ rw)
        return tot / dof

    if chi2_red(0.0) <= 1.0:
        return 0.0, chi2_red(0.0)
    lo, hi = 0.0, max_floor_mas
    for _ in range(n_bisect):
        mid = 0.5 * (lo + hi)
        if chi2_red(mid) > 1.0:
            lo = mid
        else:
            hi = mid
    s_c = 0.5 * (lo + hi)
    return s_c, chi2_red(s_c)


@dataclass(frozen=True)
class JackknifeResult:
    """Resultado del jackknife dejar-un-objetivo-fuera sobre σ(masa).

    ``sigma_jack_msun`` es la σ externa del estimador; ``masses_msun`` son las ``N``
    masas de las réplicas (una por objetivo excluido); ``mean_msun`` su media;
    ``bias_msun`` la corrección de sesgo jackknife ``(N−1)(m̄ − m_full)`` (informativa,
    no se aplica). ``n_failed`` cuenta réplicas con masa no finita (excluidas del
    cómputo); las marcadas ``converged=False`` por LM pero con masa finita SÍ se
    incluyen, igual que el ajuste de producción.
    """

    sigma_jack_msun: float
    masses_msun: np.ndarray
    mean_msun: float
    bias_msun: float
    n_failed: int


def jackknife_mass_sigma(
    targets: list[TargetObservations],
    mass_msun: float,
    fitted_elements: list[KeplerElements],
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
    backend: str = "assist",
    gr: bool = True,
    sys_floor_mas: float = 0.0,
    max_iter: int = 40,
    n_workers: int = 1,
    **lm_kwargs,
) -> JackknifeResult:
    """σ(masa) externa por jackknife dejar-un-objetivo-fuera (leave-one-target-out).

    Re-ajusta la masa compartida ``N`` veces, cada una excluyendo un objetivo distinto,
    partiendo (**warm-start**) de la solución convergida del ajuste completo: la masa
    ``mass_msun`` como semilla y, para cada objetivo retenido, sus elementos ya
    ajustados ``fitted_elements``. Como cada réplica arranca en (o muy cerca de) su
    óptimo, converge en pocas iteraciones LM, así que el coste es ``~N`` ajustes tibios,
    no ``N`` fríos.

    La σ jackknife
    ``σ_jack = sqrt((N−1)/N · Σ_i (m_(i) − m̄)²)``, con ``m̄ = (1/N) Σ_i m_(i)``,
    mide cuánto se mueve la masa al quitar objetivos individuales. Captura el error de
    **regresión masa↔órbita** que la σ formal (Fisher de ``(JᵀC⁻¹J)⁻¹``) no ve: en
    perturbadores cuya deflexión queda bajo el ruido por-encuentro la masa es sensible a
    qué objetivos entran, y σ_jack ≫ σ_formal. Se reporta ``max(σ_formal, σ_jack)``.

    Parameters
    ----------
    targets, fitted_elements
        Los ``N`` objetivos del ajuste completo y sus elementos ajustados (mismo orden
        y longitud). ``targets[k].initial_elements`` se reemplaza por
        ``fitted_elements[k]`` antes de cada réplica.
    mass_msun
        Masa compartida convergida del ajuste completo (semilla de cada réplica).

    Returns
    -------
    JackknifeResult
        Requiere ``N ≥ 3``; con menos, ``sigma_jack_msun`` es ``nan`` (jackknife no
        informativo).
    """
    n = len(targets)
    if n != len(fitted_elements):
        raise ValueError("targets y fitted_elements deben tener la misma longitud")
    warm = [dataclasses.replace(t, initial_elements=el) for t, el in zip(targets, fitted_elements)]
    if n < 3:
        return JackknifeResult(
            sigma_jack_msun=float("nan"),
            masses_msun=np.array([], dtype=float),
            mean_msun=float("nan"),
            bias_msun=float("nan"),
            n_failed=0,
        )

    replicate_masses: list[float] = []
    n_failed = 0
    for i in range(n):
        subset = warm[:i] + warm[i + 1 :]
        m_i, _fitted_i, result_i = determine_shared_mass(
            subset,
            float(mass_msun),
            perturber_elements,
            epoch_jd_tdb,
            perturber_name=perturber_name,
            background_perturbers=background_perturbers,
            perturbers=perturbers,
            integrator=integrator,
            dt_days=dt_days,
            n_lighttime_iter=n_lighttime_iter,
            gm_rel_delta=gm_rel_delta,
            backend=backend,
            gr=gr,
            sys_floor_mas=sys_floor_mas,
            max_iter=max_iter,
            n_workers=n_workers,
            **lm_kwargs,
        )
        # Aceptar réplicas con masa finita aunque LM marque ``converged=False``: el
        # ajuste de producción reporta sus fits igual (alcanzan χ²_red≈1 al tope de
        # iteraciones sin cumplir la tolerancia estricta de paso/gradiente de LM).
        # Exigir convergencia estricta descartaría réplicas válidas y ruidizaría σ_jack.
        if not np.isfinite(m_i):
            n_failed += 1
            continue
        replicate_masses.append(float(m_i))

    masses = np.asarray(replicate_masses, dtype=float)
    if masses.size < 3:
        return JackknifeResult(
            sigma_jack_msun=float("nan"),
            masses_msun=masses,
            mean_msun=float(masses.mean()) if masses.size else float("nan"),
            bias_msun=float("nan"),
            n_failed=n_failed,
        )
    n_eff = masses.size
    mean = float(masses.mean())
    sigma_jack = float(np.sqrt((n_eff - 1) / n_eff * np.sum((masses - mean) ** 2)))
    bias = float((n_eff - 1) * (mean - float(mass_msun)))
    return JackknifeResult(
        sigma_jack_msun=sigma_jack,
        masses_msun=masses,
        mean_msun=mean,
        bias_msun=bias,
        n_failed=n_failed,
    )


__all__ = [
    "TargetObservations",
    "JackknifeResult",
    "determine_mass_and_orbit",
    "determine_shared_mass",
    "calibrate_sys_floor",
    "jackknife_mass_sigma",
]
