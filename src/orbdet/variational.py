"""Ecuaciones variacionales — el corazón de la determinación de órbitas+masa.

Provee las parciales que alimentan el Jacobiano del corrector diferencial (T5) y
del ajuste conjunto órbita+masa (T6):

- **Matriz de transición de estado** ``Φ(t) = ∂[r, v](t) / ∂[r, v]₀`` (6×6),
  integrada analíticamente junto a la trayectoria con las partículas
  variacionales de primer orden de REBOUND (``sim.add_variation``). No usa
  diferencias finitas: las ecuaciones variacionales se integran con el mismo
  integrador que la órbita.
- **Parciales respecto a los elementos** ``∂[r, v](t) / ∂elementos`` (6×6), por
  regla de la cadena ``Φ(t) · J_elem``, donde ``J_elem = ∂[r, v]₀ / ∂elementos``
  es el Jacobiano analítico del mapa kepleriano en la época
  (:func:`orbdet.kepler.dstate_delements`).
- **Parcial respecto a GM del perturbador** ``∂[r, v](t) / ∂GM``, en dos variantes:

  - :func:`partial_wrt_gm` — **diferencias finitas centrales** con verificación de
    convergencia (Richardson), según la decisión T3 (2026-06-01). Es la variante
    usada por el backend ASSIST, donde las fuerzas de la efeméride no propagan
    partículas variacionales.
  - :func:`partial_wrt_gm_variational` — **partícula variacional analítica** de la
    masa (F6): integra ``∂[r, v]/∂GM`` junto a la trayectoria en una sola
    propagación por sentido (``Variation.vary(index, "m")``). Válida en el backend
    :mod:`orbdet.dynamics` (Sol + planetas como partículas masivas de REBOUND);
    coincide con la FD a mejor que 1e-6 relativo y ahorra dos propagaciones por
    Jacobiano.

Convenciones de marco/unidades idénticas a :mod:`orbdet.dynamics` (baricéntrico
eclíptico J2000; AU, día, M_sun; ``GM`` en AU³/día²).
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

from .constants import GM_SUN
from .dynamics import (
    DEFAULT_PERTURBERS,
    AsteroidPerturber,
    _build_simulation,
    planet_state_ecliptic,
    propagate,
)
from .kepler import KeplerElements, dstate_delements, elements_to_state

# Perturbación unitaria que siembra cada partícula variacional (columna de Φ),
# en el orden de fase (x, y, z, vx, vy, vz).
_PHASE_COORDS: tuple[str, ...] = ("x", "y", "z", "vx", "vy", "vz")


def _sorted_fwd_bwd(dts: np.ndarray) -> tuple[list[int], list[int]]:
    """Índices ordenados para integrar hacia adelante (dt≥0) y atrás (dt<0).

    Misma partición que :func:`orbdet.dynamics.propagate`: hacia adelante en orden
    ascendente de ``dt``, hacia atrás en orden descendente.
    """
    order = np.argsort(dts)
    fwd = [int(i) for i in order if dts[i] >= 0.0]
    bwd = [int(i) for i in order[::-1] if dts[i] < 0.0]
    return fwd, bwd


@dataclass(frozen=True)
class StateTransition:
    """Resultado de :func:`propagate_with_stm`.

    Attributes
    ----------
    positions, velocities:
        Estado baricéntrico eclíptico ``(N, 3)`` en AU y AU/día (idéntico a
        :func:`orbdet.dynamics.propagate`).
    stm:
        Matriz de transición de estado ``(N, 6, 6)``: ``stm[k]`` = ``∂[r, v](tₖ) /
        ∂[r, v]₀`` con filas ``(rx, ry, rz, vx, vy, vz)`` y columnas en el mismo
        orden de fase respecto al estado inicial.
    """

    positions: np.ndarray
    velocities: np.ndarray
    stm: np.ndarray


def propagate_with_stm(
    test_elements: KeplerElements,
    epoch_jd_tdb: float,
    out_epochs_jd_tdb: np.ndarray,
    *,
    perturbers: tuple[str, ...] = DEFAULT_PERTURBERS,
    integrator: str = "ias15",
    dt_days: float = 1.0,
    asteroid_perturbers: tuple[AsteroidPerturber, ...] = (),
) -> StateTransition:
    """Propaga el objetivo y su matriz de transición de estado ``Φ(t)``.

    La trayectoria y las seis partículas variacionales (una por coordenada de
    fase inicial) se integran simultáneamente con el mismo integrador, de modo que
    ``Φ`` es analítica (ecuaciones variacionales de primer orden), no por
    diferencias finitas.

    Parameters
    ----------
    test_elements, epoch_jd_tdb, out_epochs_jd_tdb, perturbers, integrator, dt_days, asteroid_perturbers:
        Idénticos a :func:`orbdet.dynamics.propagate`.

    Returns
    -------
    StateTransition
        Estado y ``stm`` ``(N, 6, 6)`` en las épocas pedidas.
    """
    out_epochs = np.atleast_1d(np.asarray(out_epochs_jd_tdb, dtype=float))
    n = out_epochs.size
    pos_out = np.full((n, 3), np.nan, dtype=float)
    vel_out = np.full((n, 3), np.nan, dtype=float)
    stm_out = np.full((n, 6, 6), np.nan, dtype=float)

    r_h, v_h = elements_to_state(test_elements, GM_SUN)
    sun_p, sun_v = planet_state_ecliptic("sun", epoch_jd_tdb)
    test_state = (r_h + sun_p, v_h + sun_v)

    dts = out_epochs - epoch_jd_tdb
    fwd, bwd = _sorted_fwd_bwd(dts)

    for idx_list in (fwd, bwd):
        if not idx_list:
            continue
        sim = _build_simulation(
            epoch_jd_tdb, test_state, perturbers, integrator, dt_days, asteroid_perturbers
        )
        test_idx = sim.N - 1

        # Seis variaciones de primer orden, cada una con una perturbación unitaria
        # en una coordenada de fase de la partícula de prueba. Las demás partículas
        # variacionales quedan en cero (la prueba es no masiva → no perturba a los
        # cuerpos activos, y derivamos solo respecto a su estado inicial).
        variations = []
        for coord in _PHASE_COORDS:
            var = sim.add_variation(order=1)
            setattr(var.particles[test_idx], coord, 1.0)
            variations.append(var)

        for i in idx_list:
            sim.integrate(float(dts[i]))
            tp = sim.particles[test_idx]
            pos_out[i] = (tp.x, tp.y, tp.z)
            vel_out[i] = (tp.vx, tp.vy, tp.vz)
            for col, var in enumerate(variations):
                vp = var.particles[test_idx]
                stm_out[i, :, col] = (vp.x, vp.y, vp.z, vp.vx, vp.vy, vp.vz)

    return StateTransition(positions=pos_out, velocities=vel_out, stm=stm_out)


def partials_wrt_elements(
    test_elements: KeplerElements,
    epoch_jd_tdb: float,
    out_epochs_jd_tdb: np.ndarray,
    *,
    perturbers: tuple[str, ...] = DEFAULT_PERTURBERS,
    integrator: str = "ias15",
    dt_days: float = 1.0,
    asteroid_perturbers: tuple[AsteroidPerturber, ...] = (),
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Parciales ``∂[r, v](t) / ∂elementos`` por regla de la cadena.

    ``∂[r, v](t)/∂elementos = Φ(t) · J_elem`` con ``J_elem = ∂[r, v]₀/∂elementos``
    (analítico, :func:`orbdet.kepler.dstate_delements`). El offset baricéntrico del
    Sol en la época no depende de los elementos, así que no entra en ``J_elem``.

    Returns
    -------
    (positions, velocities, dstate)
        ``positions``/``velocities`` ``(N, 3)``; ``dstate`` ``(N, 6, 6)`` con
        columnas en el orden ``(a, e, i, Ω, ω, M)``.
    """
    st = propagate_with_stm(
        test_elements,
        epoch_jd_tdb,
        out_epochs_jd_tdb,
        perturbers=perturbers,
        integrator=integrator,
        dt_days=dt_days,
        asteroid_perturbers=asteroid_perturbers,
    )
    j_elem = dstate_delements(test_elements, GM_SUN)  # (6, 6)
    dstate = st.stm @ j_elem  # (N, 6, 6) por broadcasting
    return st.positions, st.velocities, dstate


def _propagate_state(
    test_elements: KeplerElements,
    epoch_jd_tdb: float,
    out_epochs: np.ndarray,
    perturbers: tuple[str, ...],
    integrator: str,
    dt_days: float,
    asteroid_perturbers: tuple[AsteroidPerturber, ...],
) -> np.ndarray:
    """Estado ``(N, 6)`` = ``[r | v]`` apilado (helper para diferencias finitas)."""
    pos, vel = propagate(
        test_elements,
        epoch_jd_tdb,
        out_epochs,
        perturbers=perturbers,
        integrator=integrator,
        dt_days=dt_days,
        asteroid_perturbers=asteroid_perturbers,
        return_velocity=True,
    )
    return np.hstack([pos, vel])


def _scaled_perturber_mass(
    asteroid_perturbers: tuple[AsteroidPerturber, ...], index: int, dmass: float
) -> tuple[AsteroidPerturber, ...]:
    """Copia de *asteroid_perturbers* con la masa del *index* desplazada *dmass*."""
    ap = asteroid_perturbers[index]
    bumped = replace(ap, mass_msun=ap.mass_msun + dmass)
    return tuple(bumped if k == index else p for k, p in enumerate(asteroid_perturbers))


def partial_wrt_gm(
    test_elements: KeplerElements,
    epoch_jd_tdb: float,
    out_epochs_jd_tdb: np.ndarray,
    *,
    perturber_index: int,
    perturbers: tuple[str, ...] = DEFAULT_PERTURBERS,
    integrator: str = "ias15",
    dt_days: float = 1.0,
    asteroid_perturbers: tuple[AsteroidPerturber, ...],
    rel_delta: float = 1e-3,
) -> np.ndarray:
    """Parcial ``∂[r, v](t) / ∂GM_perturbador`` por diferencias finitas centrales.

    Se desplaza la masa del perturbador ``asteroid_perturbers[perturber_index]``
    en ``±δm = rel_delta · m`` y se toma la diferencia central del estado
    propagado. Como ``GM = GM_SUN · m`` (el ``G`` de la simulación es ``GM_SUN`` y
    las masas van en M_sun), la parcial respecto a ``GM`` es la parcial respecto a
    ``m`` dividida por ``GM_SUN``:

    ``∂x/∂GM = [x(m+δm) − x(m−δm)] / (2 · δm · GM_SUN)``.

    Cambiar la masa del perturbador solo altera las fuerzas (su estado inicial se
    fija desde sus elementos, independiente de la masa), que es justo el efecto
    físico buscado.

    Parameters
    ----------
    perturber_index:
        Índice del perturbador en *asteroid_perturbers* cuyo ``GM`` se deriva.
    asteroid_perturbers:
        Debe contener al menos al perturbador indexado con masa nominal > 0.
    rel_delta:
        Paso relativo de masa. ``1e-3`` es la elección por defecto de la decisión
        T3; usar :func:`richardson_convergence_dgm` para verificar la meseta.

    Returns
    -------
    np.ndarray
        ``(N, 6)`` con columnas ``(rx, ry, rz, vx, vy, vz)``: la derivada del
        estado respecto a ``GM`` (unidades AU/(AU³·día⁻²) y AU·día⁻¹/(AU³·día⁻²)).
    """
    m0 = asteroid_perturbers[perturber_index].mass_msun
    if m0 <= 0.0:
        raise ValueError(
            "partial_wrt_gm requiere masa nominal > 0 para la diferencia central "
            f"relativa; mass_msun={m0}"
        )
    dmass = rel_delta * m0

    out_epochs = np.atleast_1d(np.asarray(out_epochs_jd_tdb, dtype=float))
    common = (perturbers, integrator, dt_days)
    state_plus = _propagate_state(
        test_elements,
        epoch_jd_tdb,
        out_epochs,
        *common,
        _scaled_perturber_mass(asteroid_perturbers, perturber_index, +dmass),
    )
    state_minus = _propagate_state(
        test_elements,
        epoch_jd_tdb,
        out_epochs,
        *common,
        _scaled_perturber_mass(asteroid_perturbers, perturber_index, -dmass),
    )
    return np.asarray((state_plus - state_minus) / (2.0 * dmass * GM_SUN), dtype=float)


def _perturber_particle_index(perturbers: tuple[str, ...], perturber_index: int) -> int:
    """Índice rebound de ``asteroid_perturbers[perturber_index]`` en la simulación.

    El orden de :func:`orbdet.dynamics._build_simulation` es: el Sol y los planetas
    (con el Sol antepuesto si falta), luego los ``asteroid_perturbers``, luego el
    objetivo. Por tanto el asteroide ``perturber_index`` cae en
    ``len(nombres_planetarios) + perturber_index``.
    """
    names = [p.lower() for p in perturbers]
    if "sun" not in names:
        names = ["sun", *names]
    return len(names) + int(perturber_index)


def partial_wrt_gm_variational(
    test_elements: KeplerElements,
    epoch_jd_tdb: float,
    out_epochs_jd_tdb: np.ndarray,
    *,
    perturber_index: int,
    perturbers: tuple[str, ...] = DEFAULT_PERTURBERS,
    integrator: str = "ias15",
    dt_days: float = 1.0,
    asteroid_perturbers: tuple[AsteroidPerturber, ...],
) -> np.ndarray:
    """Parcial ``∂[r, v](t) / ∂GM_perturbador`` **analítica** (F6).

    Integra la ecuación variacional de primer orden respecto a la masa del
    perturbador junto a la trayectoria, con la partícula variacional de masa de
    REBOUND (``Variation.vary(index, "m")``). A diferencia de
    :func:`partial_wrt_gm` (diferencias finitas, dos propagaciones completas), aquí
    la sensibilidad se propaga con el mismo integrador que la órbita en **una sola**
    integración por sentido, y captura exactamente el acoplamiento perturbador↔resto
    del sistema (el objetivo no masivo y todos los cuerpos activos son partículas de
    REBOUND, así que ``sim.G`` y la gravedad mutua entran en la variacional).

    .. note::
       Sólo es válida en el backend :mod:`orbdet.dynamics` (Sol + planetas como
       partículas masivas de REBOUND). Bajo el backend ASSIST las fuerzas del
       Sol/planetas/GR son ``additional_forces`` de la efeméride que **no** propagan
       partículas variacionales, por lo que allí la parcial sigue siendo por
       diferencias finitas (:func:`partial_wrt_gm`).

    REBOUND deriva respecto a la *masa* de la partícula (``m`` = ``mass_msun``, pues
    ``sim.G = GM_SUN``). Como ``GM = GM_SUN · mass_msun``, la parcial respecto a
    ``GM`` es la parcial respecto a la masa dividida por ``GM_SUN`` — idéntica
    convención y unidades que :func:`partial_wrt_gm`.

    Parameters
    ----------
    perturber_index:
        Índice del perturbador en *asteroid_perturbers* cuyo ``GM`` se deriva.
    asteroid_perturbers:
        Debe contener al menos al perturbador indexado (su masa puede ser cualquier
        valor ≥ 0: la variacional es lineal y no la usa como paso).

    Returns
    -------
    np.ndarray
        ``(N, 6)`` con columnas ``(rx, ry, rz, vx, vy, vz)``: ``∂[r, v]/∂GM`` en el
        marco baricéntrico eclíptico, misma convención que :func:`partial_wrt_gm`.
    """
    out_epochs = np.atleast_1d(np.asarray(out_epochs_jd_tdb, dtype=float))
    n = out_epochs.size
    dgm_out = np.full((n, 6), np.nan, dtype=float)

    r_h, v_h = elements_to_state(test_elements, GM_SUN)
    sun_p, sun_v = planet_state_ecliptic("sun", epoch_jd_tdb)
    test_state = (r_h + sun_p, v_h + sun_v)

    pert_idx = _perturber_particle_index(perturbers, perturber_index)

    dts = out_epochs - epoch_jd_tdb
    fwd, bwd = _sorted_fwd_bwd(dts)

    for idx_list in (fwd, bwd):
        if not idx_list:
            continue
        sim = _build_simulation(
            epoch_jd_tdb, test_state, perturbers, integrator, dt_days, asteroid_perturbers
        )
        test_idx = sim.N - 1
        # Una partícula variacional de primer orden respecto a la MASA del
        # perturbador estudiado; el resto de las semillas quedan en cero.
        var = sim.add_variation(order=1)
        var.vary(pert_idx, "m")

        for i in idx_list:
            sim.integrate(float(dts[i]))
            vp = var.particles[test_idx]
            # ∂[r,v]_test/∂mass_msun → ∂/∂GM dividiendo por GM_SUN.
            dgm_out[i] = (vp.x, vp.y, vp.z, vp.vx, vp.vy, vp.vz)

    return np.asarray(dgm_out / GM_SUN, dtype=float)


def richardson_convergence_dgm(
    test_elements: KeplerElements,
    epoch_jd_tdb: float,
    out_epochs_jd_tdb: np.ndarray,
    *,
    perturber_index: int,
    perturbers: tuple[str, ...] = DEFAULT_PERTURBERS,
    integrator: str = "ias15",
    dt_days: float = 1.0,
    asteroid_perturbers: tuple[AsteroidPerturber, ...],
    rel_deltas: tuple[float, ...] = (4e-3, 2e-3, 1e-3, 5e-4),
) -> dict[str, np.ndarray]:
    """Barre el paso ``rel_delta`` para verificar la meseta de Richardson de ∂x/∂GM.

    Devuelve un diccionario con:

    - ``rel_deltas``: los pasos evaluados ``(K,)``.
    - ``partials``: ``(K, N, 6)`` la parcial a cada paso.
    - ``max_rel_change``: ``(K-1,)`` máximo cambio relativo entre pasos sucesivos
      (norma del cambio / norma de la parcial), agregando sobre épocas y
      componentes. En la meseta numérica debe decrecer y estabilizarse muy por
      debajo de 1 (las diferencias centrales convergen como ``O(δ²)``).
    """
    partials = np.stack(
        [
            partial_wrt_gm(
                test_elements,
                epoch_jd_tdb,
                out_epochs_jd_tdb,
                perturber_index=perturber_index,
                perturbers=perturbers,
                integrator=integrator,
                dt_days=dt_days,
                asteroid_perturbers=asteroid_perturbers,
                rel_delta=d,
            )
            for d in rel_deltas
        ]
    )
    diffs = np.linalg.norm(partials[1:] - partials[:-1], axis=(1, 2))
    base = np.linalg.norm(partials[:-1], axis=(1, 2))
    max_rel_change = diffs / np.where(base > 0.0, base, 1.0)
    return {
        "rel_deltas": np.asarray(rel_deltas, dtype=float),
        "partials": partials,
        "max_rel_change": max_rel_change,
    }
