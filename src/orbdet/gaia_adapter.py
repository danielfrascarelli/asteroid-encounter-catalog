"""Adaptador de datos reales (Gaia + MPCORB) → motor de determinación (T9).

Convierte las observaciones por-tránsito de Gaia y los elementos orbitales del MPC
a la convención interna del motor (:class:`orbdet.mass_determination.TargetObservations`),
sin atajos:

- **Tiempos:** época Gaia (días desde J2010.0 TCB) → JD TDB exacto vía astropy
  (:func:`orbdet.time_scales.gaia_epoch_to_jd_tdb`).
- **Pesos:** σ_along-scan por tránsito proyectando la **elipse de covarianza
  (RA, Dec) completa** de Gaia (familias sistemática + aleatoria, con su
  correlación) sobre la dirección de barrido. Es la incertidumbre real de Gaia en
  la dirección bien medida; la across-scan (~cientos de mas) se descarta porque no
  aporta información — exactamente el modelo de error que usa la determinación de
  masas del campo (Fuentes-Muñoz).
- **Elementos:** fila MPCORB (grados, J2000 eclíptico, heliocéntrico) →
  :class:`orbdet.kepler.KeplerElements` (radianes). Armonización de épocas entre el
  MPCORB (época estándar) y la época común del ajuste **propagando con el propio
  modelo N-cuerpos** y re-extrayendo elementos osculadores — no Kepler de dos
  cuerpos — para no introducir error de propagación en la semilla ni en el
  perturbador.

Contrato de aislamiento: este módulo NO carga parquet ni habla con el TAP (eso usa
polars/astroquery, fuera de la whitelist de ``orbdet``). Recibe arrays numpy ya
cargados; el IO vive en ``src/ingest`` o en los scripts de ``scripts/mass``.

Convención del ángulo de barrido idéntica al resto del motor: el versor along-scan
en el plano tangente ``(RA·cosδ, δ)`` es ``û_AL = (sin PA, cos PA)`` con PA medido
desde el Norte (+δ) hacia el Este (+RA·cosδ).
"""

from __future__ import annotations

import numpy as np

from .constants import GM_SUN
from .dynamics import DEFAULT_PERTURBERS, AsteroidPerturber, planet_state_ecliptic, propagate
from .kepler import KeplerElements, state_to_elements
from .mass_determination import TargetObservations
from .time_scales import J2010_TCB_JD, gaia_epoch_to_jd_tdb

# Piso de varianza AL (mas²) para evitar pesos infinitos por covarianzas degeneradas.
_VAR_AL_FLOOR: float = 1e-6
# Clamp del coeficiente de correlación: el archivo garantiza |ρ|≤1 pero exports
# defectuosos lo superan ocasionalmente, lo que rompería la positividad.
_RHO_CLAMP: float = 0.9999

_DEG: float = np.pi / 180.0


def sigma_al_from_radec_covariance(
    pa_scan_deg: np.ndarray,
    ra_err_sys: np.ndarray,
    dec_err_sys: np.ndarray,
    corr_sys: np.ndarray,
    ra_err_rand: np.ndarray,
    dec_err_rand: np.ndarray,
    corr_rand: np.ndarray,
) -> np.ndarray:
    """σ along-scan por tránsito proyectando la covarianza (RA, Dec) de Gaia.

    La covarianza total es la suma de las familias sistemática y aleatoria,
    ``Σ = Σ_sys + Σ_rand`` con ``Σ = [[σ_RA², ρ σ_RA σ_Dec], [·, σ_Dec²]]`` en el
    plano tangente ``(RA·cosδ, δ)`` (los ``*_error_*`` de Gaia ya vienen en RA·cosδ,
    en mas). La varianza along-scan es ``σ²_AL = û_ALᵀ Σ û_AL`` con
    ``û_AL = (sin PA, cos PA)``:

    ``σ²_AL = sin²PA·σ_RA² + 2 sinPA cosPA ρ σ_RA σ_Dec + cos²PA·σ_Dec²``.

    Parameters
    ----------
    pa_scan_deg:
        Ángulo de posición de barrido (grados, Norte→Este), ``(N,)``.
    ra_err_sys, dec_err_sys, corr_sys:
        σ y correlación de la familia **sistemática** (mas), ``(N,)``.
    ra_err_rand, dec_err_rand, corr_rand:
        σ y correlación de la familia **aleatoria** (ruido de fotones), ``(N,)``.

    Returns
    -------
    np.ndarray
        σ_AL por tránsito en mas, ``(N,)``.
    """
    pa = np.asarray(pa_scan_deg, dtype=float) * _DEG
    e_ra, e_dec = np.sin(pa), np.cos(pa)

    def _var_al(s_ra: np.ndarray, s_dec: np.ndarray, rho: np.ndarray) -> np.ndarray:
        s_ra = np.asarray(s_ra, dtype=float)
        s_dec = np.asarray(s_dec, dtype=float)
        rho = np.clip(np.asarray(rho, dtype=float), -_RHO_CLAMP, _RHO_CLAMP)
        var: np.ndarray = (
            e_ra * e_ra * s_ra * s_ra
            + 2.0 * e_ra * e_dec * rho * s_ra * s_dec
            + e_dec * e_dec * s_dec * s_dec
        )
        return var

    var_al = _var_al(ra_err_sys, dec_err_sys, corr_sys) + _var_al(
        ra_err_rand, dec_err_rand, corr_rand
    )
    sigma_al: np.ndarray = np.sqrt(np.maximum(var_al, _VAR_AL_FLOOR))
    return sigma_al


def elements_from_mpcorb(
    a_au: float,
    e: float,
    i_deg: float,
    Omega_deg: float,
    omega_deg: float,
    M_deg: float,
) -> KeplerElements:
    """Fila MPCORB (grados, J2000 eclíptico heliocéntrico) → :class:`KeplerElements`.

    MPCORB entrega los ángulos en grados: ``M`` (anomalía media), ``ω`` (argumento
    del perihelio), ``Ω`` (nodo ascendente) e ``i`` (inclinación a la eclíptica
    J2000). El motor los usa en radianes y en el orden ``(a, e, i, Ω, ω, M)``.
    """
    return KeplerElements(
        a=float(a_au),
        e=float(e),
        i=float(i_deg) * _DEG,
        Omega=float(Omega_deg) * _DEG,
        omega=float(omega_deg) * _DEG,
        M=float(M_deg) * _DEG,
    )


def propagate_elements(
    el: KeplerElements,
    epoch_from_jd_tdb: float,
    epoch_to_jd_tdb: float,
    *,
    backend: str = "assist",
    perturbers: tuple[str, ...] = DEFAULT_PERTURBERS,
    asteroid_perturbers: tuple[AsteroidPerturber, ...] = (),
    integrator: str = "ias15",
    dt_days: float = 1.0,
) -> KeplerElements:
    """Propaga elementos heliocéntricos de una época a otra con el modelo N-cuerpos.

    Integra el estado bajo la dinámica completa (no Kepler de dos cuerpos) y
    re-extrae los **elementos osculadores heliocéntricos** en ``epoch_to``. Se usa
    para llevar todos los cuerpos (objetivos, perturbador, fondo) a una época común
    del ajuste cercana al arco de datos —mejor condicionamiento del sistema de
    mínimos cuadrados— sin introducir error de dos cuerpos en la semilla.

    ``backend="assist"`` (por defecto) usa el modelo DE440 + GR (exacto a ~mas sobre
    años); ``"rebound"`` usa los planetas integrados libremente del backend previo.
    Si ``epoch_to == epoch_from`` devuelve los elementos sin tocar.
    """
    if epoch_to_jd_tdb == epoch_from_jd_tdb:
        return el
    if backend == "assist":
        from .dynamics_assist import propagate_assist

        pos_bary, vel_bary = propagate_assist(
            el,
            epoch_from_jd_tdb,
            np.array([epoch_to_jd_tdb], dtype=float),
            asteroid_perturbers=asteroid_perturbers,
            return_velocity=True,
        )
    else:
        pos_bary, vel_bary = propagate(
            el,
            epoch_from_jd_tdb,
            np.array([epoch_to_jd_tdb], dtype=float),
            perturbers=perturbers,
            integrator=integrator,
            dt_days=dt_days,
            asteroid_perturbers=asteroid_perturbers,
            return_velocity=True,
        )
    # Nota (tribunal 2026-07-04, menor 3): el estado del Sol viene de la efeméride
    # builtin de astropy (erfa), mientras pos_bary puede venir de ASSIST/DE440.
    # La mezcla introduce ~km en r_helio (fraccionalmente ~1e-5 en los elementos
    # heliocéntricos del perturbador en modo perturber-orbit=mpcorb). Acotado por
    # el gate Horizons (0.17 mas); unificar efeméride si se requiere mejor.
    sun_p, sun_v = planet_state_ecliptic("sun", epoch_to_jd_tdb)
    r_helio = np.asarray(pos_bary[0], dtype=float) - sun_p
    v_helio = np.asarray(vel_bary[0], dtype=float) - sun_v
    return state_to_elements(r_helio, v_helio, GM_SUN)


def gaia_positions_icrs(x_gaia: np.ndarray, y_gaia: np.ndarray, z_gaia: np.ndarray) -> np.ndarray:
    """Apila ``(x_gaia, y_gaia, z_gaia)`` (AU, ICRS baricéntrico) en ``(N, 3)``.

    Las columnas ``*_gaia`` del archivo Gaia ya están en el BCRS (orientación ICRS,
    AU), que es justo lo que espera :func:`orbdet.observation.radec_from_positions`.
    """
    return np.column_stack(
        [
            np.asarray(x_gaia, dtype=float),
            np.asarray(y_gaia, dtype=float),
            np.asarray(z_gaia, dtype=float),
        ]
    )


def fov_groups_from_epochs(obs_jd: np.ndarray, gap_days: float = 0.01) -> np.ndarray:
    """Etiqueta cada observación con el cruce de plano focal (FOV transit) al que pertenece.

    Gaia mide cada objeto con hasta ~9 CCDs por cruce, separados ~segundos
    (~1e-4 d), mientras que cruces distintos distan ≳0.05 d. Un corte en
    ``gap_days`` (~14 min por defecto) separa limpiamente los CCDs de un mismo cruce
    de los de cruces distintos. Devuelve un array de enteros ``(N,)`` en el **orden
    original** de ``obs_jd`` (no reordena). Las observaciones de un mismo grupo
    comparten error sistemático → se blanquean con covarianza en bloques.
    """
    obs = np.asarray(obs_jd, dtype=float)
    order = np.argsort(obs, kind="stable")
    labels = np.empty(obs.size, dtype=int)
    g = 0
    prev = None
    for rank, i in enumerate(order):
        if prev is not None and (obs[i] - prev) > gap_days:
            g += 1
        labels[i] = g
        prev = obs[i]
    return labels


def build_target_observations(
    initial_elements: KeplerElements,
    epoch_jd_tdb: float,
    *,
    epoch_days_tcb: np.ndarray,
    ra_deg: np.ndarray,
    dec_deg: np.ndarray,
    pa_scan_deg: np.ndarray,
    ra_err_sys: np.ndarray,
    dec_err_sys: np.ndarray,
    corr_sys: np.ndarray,
    ra_err_rand: np.ndarray,
    dec_err_rand: np.ndarray,
    corr_rand: np.ndarray,
    x_gaia: np.ndarray,
    y_gaia: np.ndarray,
    z_gaia: np.ndarray,
    epoch_ref_jd_tcb: float = J2010_TCB_JD,
    fov_gap_days: float = 0.01,
) -> TargetObservations:
    """Ensambla un :class:`TargetObservations` a partir de columnas crudas de Gaia.

    Convierte tiempos (TCB→TDB), deriva σ_AL de la covarianza (RA, Dec) y arma la
    posición de Gaia. ``initial_elements`` deben estar ya referidos a
    ``epoch_jd_tdb`` (usar :func:`propagate_elements` si vienen de otra época).
    """
    obs_jd_tdb = gaia_epoch_to_jd_tdb(epoch_days_tcb, epoch_ref_jd_tcb)
    sigma_al = sigma_al_from_radec_covariance(
        pa_scan_deg,
        ra_err_sys,
        dec_err_sys,
        corr_sys,
        ra_err_rand,
        dec_err_rand,
        corr_rand,
    )
    obs_jd_tdb = np.asarray(obs_jd_tdb, dtype=float)
    return TargetObservations(
        initial_elements=initial_elements,
        obs_jd_tdb=obs_jd_tdb,
        ra_obs_deg=np.asarray(ra_deg, dtype=float),
        dec_obs_deg=np.asarray(dec_deg, dtype=float),
        pa_scan_deg=np.asarray(pa_scan_deg, dtype=float),
        sigma_al_mas=np.asarray(sigma_al, dtype=float),
        gaia_bary_icrs=gaia_positions_icrs(x_gaia, y_gaia, z_gaia),
        fov_group=fov_groups_from_epochs(obs_jd_tdb, fov_gap_days),
    )


__all__ = [
    "sigma_al_from_radec_covariance",
    "elements_from_mpcorb",
    "propagate_elements",
    "gaia_positions_icrs",
    "fov_groups_from_epochs",
    "build_target_observations",
]
