"""Backend de dinámica state-of-the-art vía ASSIST (efemérides JPL DE440 + GR).

Reemplaza la integración de planetas "libres" de :mod:`orbdet.dynamics` (sembrados
desde la efeméride en la época y dejados gravitar en rebound, que derivan de la
efeméride sobre arcos largos y omiten relatividad) por el modelo de fuerzas que
usan JPL y los trabajos modernos de determinación de masas (Fuentes-Muñoz, Holman):

- **Sol + 8 planetas + Luna + Plutón** leídos directamente de la efeméride DE440
  (posiciones exactas en cada paso, no integradas), vía ASSIST.
- **Relatividad general** (términos Einstein–Infeld–Hoffmann, ``GR_EIH``): el
  avance relativista del perihelio es significativo a nivel mas sobre años.
- **Perturbadores asteroidales masivos** (los 16 grandes de DE441): se agregan como
  partículas masivas propias de rebound —no vía la fuerza ``ASTEROIDS`` de la
  efeméride— para poder **variar la masa** del perturbador bajo estudio (la fuerza
  ``ASTEROIDS`` queda excluida para no doblar el conteo). Tras adjuntar ASSIST se
  reactiva ``sim.gravity`` para que rebound compute la gravedad mutua de esas
  partículas y su tirón sobre el objetivo.

Contrato de salida idéntico a :func:`orbdet.dynamics.propagate`: posiciones (y
opcionalmente velocidades) **baricéntricas eclípticas J2000** ``(N, 3)`` en AU
(AU/día). Internamente ASSIST trabaja en **ICRF ecuatorial baricéntrico** con
``G = 1`` y masa = GM en AU³/día²; este módulo hace las conversiones de marco
(eclíptica↔ecuatorial) y de unidades (M_sun↔GM) en los bordes.

Convenciones de tiempo: ``sim.t = JD_TDB − ephem.jd_ref`` (``jd_ref = 2451545.0``).
"""

from __future__ import annotations

import os
from functools import lru_cache

import numpy as np

from .constants import GM_SUN
from .dynamics import AsteroidPerturber
from .frames import ecliptic_to_equatorial, equatorial_to_ecliptic
from .kepler import KeplerElements, elements_to_state, state_to_elements

# Los 16 perturbadores asteroidales del modelo dinámico DE441 (sb441-n16.bsp),
# por nombre tal como los expone ASSIST. Son los que dominan las perturbaciones
# mutuas del cinturón principal.
BIG_ASTEROIDS: tuple[str, ...] = (
    "Ceres", "Vesta", "Pallas", "Hygiea", "Euphrosyne", "Interamnia", "Davida",
    "Eunomia", "Juno", "Psyche", "Cybele", "Thisbe", "Europa", "Sylvia",
    "Camilla", "Iris",
)

_DEFAULT_PLANETS_FILE = "linux_p1550p2650.440"
_DEFAULT_ASTEROIDS_FILE = "sb441-n16.bsp"


def _ephem_dir() -> str:
    return os.environ.get("ORBDET_EPHEM_DIR", os.path.join("data", "raw", "ephem"))


@lru_cache(maxsize=4)
def load_ephem(planets_path: str | None = None, asteroids_path: str | None = None):
    """Carga (cacheada) de ``assist.Ephem`` a partir de los archivos DE440.

    Rutas por defecto: ``$ORBDET_EPHEM_DIR`` (o ``data/raw/ephem``) +
    ``linux_p1550p2650.440`` y ``sb441-n16.bsp``. El cache evita re-mapear ~750 MB
    en cada propagación.
    """
    import assist

    base = _ephem_dir()
    pp = planets_path or os.path.join(base, _DEFAULT_PLANETS_FILE)
    ap = asteroids_path or os.path.join(base, _DEFAULT_ASTEROIDS_FILE)
    return assist.Ephem(pp, ap)


def _ephem_pos_vel(ephem, name: str, t: float) -> tuple[np.ndarray, np.ndarray]:
    """Estado baricéntrico ecuatorial ``(r, v)`` (AU, AU/día) de *name* en ``t`` sim.

    Los segmentos SPK de los asteroides (sb441-n16) sólo entregan posición vía
    ``get_particle`` (la velocidad sale ``nan``); se reconstruye por **diferencia
    central** de la posición, exacta a O(dt²) dada la alta precisión posicional del
    SPK. Para Sol/planetas (efeméride .440) la velocidad es nativa y se usa directo.
    """
    p = ephem.get_particle(name, t)
    r = np.array([p.x, p.y, p.z], dtype=float)
    v = np.array([p.vx, p.vy, p.vz], dtype=float)
    if not np.all(np.isfinite(v)):
        h = 1e-2  # días
        pp = ephem.get_particle(name, t + h)
        pm = ephem.get_particle(name, t - h)
        v = (
            np.array([pp.x, pp.y, pp.z], dtype=float)
            - np.array([pm.x, pm.y, pm.z], dtype=float)
        ) / (2.0 * h)
    return r, v


def _bary_eq_state_from_elements(
    el: KeplerElements, sun_xyz: np.ndarray, sun_vxyz: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Elementos heliocéntricos eclípticos → estado baricéntrico **ecuatorial**.

    Convierte (helio, ecl) → (helio, eq) por rotación y suma el estado baricéntrico
    del Sol (ya en ecuatorial, de la efeméride). Unidades AU, AU/día.
    """
    r_h_ecl, v_h_ecl = elements_to_state(el, GM_SUN)
    r_h_eq = ecliptic_to_equatorial(r_h_ecl)
    v_h_eq = ecliptic_to_equatorial(v_h_ecl)
    return r_h_eq + sun_xyz, v_h_eq + sun_vxyz


def _build_sim(
    epoch_jd_tdb: float,
    test_elements: KeplerElements,
    asteroid_perturbers: tuple[AsteroidPerturber, ...],
    ephem,
    gr: bool,
):
    """Arma la rebound.Simulation + assist.Extras para la propagación.

    Orden de partículas: primero los ``K`` perturbadores masivos (``N_active = K``),
    luego el objetivo de prueba (no masivo). ``sim.gravity`` se reactiva para que
    rebound calcule la gravedad de los perturbadores sobre el objetivo.
    """
    import assist
    import rebound

    t0 = epoch_jd_tdb - ephem.jd_ref
    sim = rebound.Simulation()
    sim.t = t0
    sun = ephem.get_particle("Sun", t0)
    sun_xyz = np.array([sun.x, sun.y, sun.z])
    sun_vxyz = np.array([sun.vx, sun.vy, sun.vz])

    for ap in asteroid_perturbers:
        r, v = _bary_eq_state_from_elements(ap.elements, sun_xyz, sun_vxyz)
        # G = 1 en ASSIST → masa rebound = GM del cuerpo = mass_msun · GM_SUN.
        sim.add(m=ap.mass_msun * GM_SUN, x=r[0], y=r[1], z=r[2], vx=v[0], vy=v[1], vz=v[2])

    r_t, v_t = _bary_eq_state_from_elements(test_elements, sun_xyz, sun_vxyz)
    sim.add(m=0.0, x=r_t[0], y=r_t[1], z=r_t[2], vx=v_t[0], vy=v_t[1], vz=v_t[2])

    extras = assist.Extras(sim, ephem)
    extras.forces = ["SUN", "PLANETS", "GR_EIH"] if gr else ["SUN", "PLANETS"]
    # ASSIST pone sim.gravity='none'; reactivamos la gravedad directa de rebound
    # para los perturbadores masivos (las fuerzas del Sol/planetas/GR siguen
    # viniendo de ASSIST como fuerzas adicionales).
    n_active = len(asteroid_perturbers)
    sim.gravity = "basic"
    sim.N_active = n_active if n_active > 0 else sim.N
    return sim, extras


def propagate_assist(
    test_elements: KeplerElements,
    epoch_jd_tdb: float,
    out_epochs_jd_tdb: np.ndarray,
    *,
    asteroid_perturbers: tuple[AsteroidPerturber, ...] = (),
    gr: bool = True,
    planets_path: str | None = None,
    asteroids_path: str | None = None,
    return_velocity: bool = False,
) -> np.ndarray | tuple[np.ndarray, np.ndarray]:
    """Propaga un objetivo bajo el modelo ASSIST (DE440 + GR + perturbadores).

    Misma firma/contrato que :func:`orbdet.dynamics.propagate` salvo que el modelo
    de fuerzas es la efeméride JPL en vez de planetas integrados libremente, y que
    el conjunto de planetas es el completo de DE440. Devuelve posiciones
    **baricéntricas eclípticas** ``(N, 3)`` (y velocidades si ``return_velocity``).
    """
    ephem = load_ephem(planets_path, asteroids_path)
    out_epochs = np.atleast_1d(np.asarray(out_epochs_jd_tdb, dtype=float))
    n = out_epochs.size
    pos_eq = np.full((n, 3), np.nan, dtype=float)
    vel_eq = np.full((n, 3), np.nan, dtype=float)

    dts = out_epochs - epoch_jd_tdb
    order = np.argsort(dts)
    fwd = [int(i) for i in order if dts[i] >= 0.0]
    bwd = [int(i) for i in order[::-1] if dts[i] < 0.0]

    for idx_list in (fwd, bwd):
        if not idx_list:
            continue
        sim, _extras = _build_sim(epoch_jd_tdb, test_elements, asteroid_perturbers, ephem, gr)
        test_idx = sim.N - 1
        for i in idx_list:
            sim.integrate(float(out_epochs[i] - ephem.jd_ref))
            p = sim.particles[test_idx]
            pos_eq[i] = (p.x, p.y, p.z)
            vel_eq[i] = (p.vx, p.vy, p.vz)

    pos_ecl = equatorial_to_ecliptic(pos_eq)
    if return_velocity:
        return pos_ecl, equatorial_to_ecliptic(vel_eq)
    return pos_ecl


def big_asteroid_perturbers(
    epoch_jd_tdb: float,
    *,
    exclude: tuple[str, ...] = (),
    names: tuple[str, ...] = BIG_ASTEROIDS,
    planets_path: str | None = None,
    asteroids_path: str | None = None,
) -> tuple[AsteroidPerturber, ...]:
    """Construye los perturbadores asteroidales de fondo desde la efeméride DE441.

    Para cada asteroide grande (salvo los de ``exclude``), lee su estado y GM de la
    efeméride en ``epoch_jd_tdb`` y arma un :class:`AsteroidPerturber` con elementos
    osculadores heliocéntricos eclípticos y masa en M_sun (``= GM / GM_SUN``). Se usa
    como ``background_perturbers`` del ajuste; el perturbador **bajo estudio** se
    pasa por separado con su masa como parámetro libre y debe estar en ``exclude``.
    """
    ephem = load_ephem(planets_path, asteroids_path)
    t = epoch_jd_tdb - ephem.jd_ref
    sun = ephem.get_particle("Sun", t)
    sun_xyz = np.array([sun.x, sun.y, sun.z])
    sun_vxyz = np.array([sun.vx, sun.vy, sun.vz])
    excl = {e.lower() for e in exclude}

    out: list[AsteroidPerturber] = []
    for name in names:
        if name.lower() in excl:
            continue
        p = ephem.get_particle(name, t)
        r_abs, v_abs = _ephem_pos_vel(ephem, name, t)
        r_eq = r_abs - sun_xyz
        v_eq = v_abs - sun_vxyz
        r_ecl = equatorial_to_ecliptic(r_eq)
        v_ecl = equatorial_to_ecliptic(v_eq)
        el = state_to_elements(r_ecl, v_ecl, GM_SUN)
        out.append(AsteroidPerturber(name=name, mass_msun=float(p.m) / GM_SUN, elements=el))
    return tuple(out)


def ephem_asteroid_mass_msun(
    name: str,
    epoch_jd_tdb: float = 2_451_545.0,
    *,
    planets_path: str | None = None,
    asteroids_path: str | None = None,
) -> float:
    """Masa (M_sun) del asteroide *name* según la efeméride DE441 (GM/GM_SUN)."""
    ephem = load_ephem(planets_path, asteroids_path)
    p = ephem.get_particle(name, epoch_jd_tdb - ephem.jd_ref)
    return float(p.m) / GM_SUN


__all__ = [
    "BIG_ASTEROIDS",
    "load_ephem",
    "propagate_assist",
    "big_asteroid_perturbers",
    "ephem_asteroid_mass_msun",
]
