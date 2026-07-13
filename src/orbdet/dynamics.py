"""Modelo dinámico N-cuerpos (fuerzas) vía REBOUND.

Integra una partícula de prueba (asteroide objetivo) bajo el Sol + planetas
mayores (+ opcionalmente asteroides perturbadores masivos, incluido el perturber
cuya masa se quiere determinar). Es la base sobre la que van las ecuaciones
variacionales (T3) y el corrector diferencial (T5).

Convenciones:
- Marco: eclíptico J2000, baricéntrico.
- Unidades: AU, días, M_sun. ``sim.G`` se fija explícitamente a ``GM_SUN`` (= k²)
  para que el GM solar coincida **exactamente** con el del propagador kepleriano
  (``kepler.py``) — así el límite de dos cuerpos es consistente bit a bit.
- Estados planetarios: ``astropy`` con efemérides ``builtin`` (series analíticas
  erfa ``epv00``/``plan94`` — precisión ~km, NO una versión de DE440; el backend
  ASSIST usa DE440 real y es el de producción), rotados ICRS→eclíptica con la
  oblicuidad de ``constants``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from astropy import units as u
from astropy.coordinates import get_body_barycentric_posvel, solar_system_ephemeris
from astropy.time import Time

from .constants import GM_SUN
from .frames import equatorial_to_ecliptic
from .kepler import KeplerElements, elements_to_state

# Masas en M_sun (IAU 2015 / DE440 nominal). 'earth' es el baricentro Tierra-Luna.
PLANET_MASSES_MSUN: dict[str, float] = {
    "sun": 1.0,
    "mercury": 1.660_113_797e-7,
    "venus": 2.447_835_335e-6,
    "earth": 3.040_432_648e-6,
    "mars": 3.227_156_038e-7,
    "jupiter": 9.547_919_152e-4,
    "saturn": 2.858_859_807e-4,
    "uranus": 4.366_244_043e-5,
    "neptune": 5.151_389_021e-5,
}

DEFAULT_PERTURBERS: tuple[str, ...] = ("sun", "jupiter", "saturn")


@dataclass(frozen=True)
class AsteroidPerturber:
    """Asteroide masivo perturbador (p. ej. el cuerpo cuya masa se determina)."""

    name: str
    mass_msun: float
    elements: KeplerElements  # heliocéntricos eclípticos, en la época del modelo


def planet_state_ecliptic(body: str, jd_tdb: float) -> tuple[np.ndarray, np.ndarray]:
    """Estado baricéntrico eclíptico J2000 de *body* en *jd_tdb*.

    Devuelve ``(pos_AU, vel_AU_por_día)`` como arrays ``(3,)``.
    """
    t = Time(jd_tdb, format="jd", scale="tdb")
    with solar_system_ephemeris.set("builtin"):
        pos, vel = get_body_barycentric_posvel(body, t)
    p_icrs = np.asarray(pos.xyz.to(u.AU).value, dtype=float)
    v_icrs = np.asarray(vel.xyz.to(u.AU / u.day).value, dtype=float)
    return equatorial_to_ecliptic(p_icrs), equatorial_to_ecliptic(v_icrs)


def _build_simulation(
    epoch_jd_tdb: float,
    test_state: tuple[np.ndarray, np.ndarray],
    perturbers: tuple[str, ...],
    integrator: str,
    dt_days: float,
    asteroid_perturbers: tuple[AsteroidPerturber, ...],
):
    """Arma una rebound.Simulation con los cuerpos masivos + la partícula de prueba.

    La partícula de prueba (índice ``N-1``) es no masiva. El reloj ``sim.t``
    arranca en 0 = *epoch_jd_tdb*.
    """
    import rebound

    sim = rebound.Simulation()
    sim.G = GM_SUN  # AU^3 / (M_sun · día^2): GM_sun = GM_SUN exactamente
    sim.integrator = integrator
    if integrator == "whfast":
        sim.dt = dt_days

    names = [p.lower() for p in perturbers]
    if "sun" not in names:
        names = ["sun", *names]
    unknown = [n for n in names if n not in PLANET_MASSES_MSUN]
    if unknown:
        raise ValueError(f"perturbadores desconocidos: {unknown}")

    for name in names:
        p, v = planet_state_ecliptic(name, epoch_jd_tdb)
        sim.add(m=PLANET_MASSES_MSUN[name], x=p[0], y=p[1], z=p[2], vx=v[0], vy=v[1], vz=v[2])

    sun_p, sun_v = planet_state_ecliptic("sun", epoch_jd_tdb)
    for ap in asteroid_perturbers:
        r_h, v_h = elements_to_state(ap.elements, GM_SUN)
        pb, vb = r_h + sun_p, v_h + sun_v
        sim.add(m=ap.mass_msun, x=pb[0], y=pb[1], z=pb[2], vx=vb[0], vy=vb[1], vz=vb[2])

    sim.N_active = sim.N  # todos los anteriores son masivos/activos

    r0, v0 = test_state
    sim.add(m=0.0, x=r0[0], y=r0[1], z=r0[2], vx=v0[0], vy=v0[1], vz=v0[2])
    return sim


def propagate(
    test_elements: KeplerElements,
    epoch_jd_tdb: float,
    out_epochs_jd_tdb: np.ndarray,
    *,
    perturbers: tuple[str, ...] = DEFAULT_PERTURBERS,
    integrator: str = "ias15",
    dt_days: float = 1.0,
    asteroid_perturbers: tuple[AsteroidPerturber, ...] = (),
    return_velocity: bool = False,
) -> np.ndarray | tuple[np.ndarray, np.ndarray]:
    """Propaga un asteroide objetivo bajo el modelo N-cuerpos.

    Parameters
    ----------
    test_elements:
        Elementos keplerianos **heliocéntricos eclípticos** del objetivo en
        *epoch_jd_tdb*.
    epoch_jd_tdb:
        Época de los elementos (JD TDB). El modelo construye los cuerpos masivos
        en esta época.
    out_epochs_jd_tdb:
        Épocas (JD TDB) en las que se quiere la posición. Pueden estar antes y
        después de la época (se integra hacia ambos lados).

    Returns
    -------
    Posiciones **baricéntricas eclípticas** ``(N, 3)`` en AU (y velocidades
    ``(N, 3)`` en AU/día si ``return_velocity``).
    """
    out_epochs = np.atleast_1d(np.asarray(out_epochs_jd_tdb, dtype=float))
    n = out_epochs.size
    pos_out = np.full((n, 3), np.nan, dtype=float)
    vel_out = np.full((n, 3), np.nan, dtype=float)

    r_h, v_h = elements_to_state(test_elements, GM_SUN)
    sun_p, sun_v = planet_state_ecliptic("sun", epoch_jd_tdb)
    test_state = (r_h + sun_p, v_h + sun_v)

    dts = out_epochs - epoch_jd_tdb
    fwd = [int(i) for i in np.argsort(dts) if dts[i] >= 0.0]
    bwd = [int(i) for i in np.argsort(dts)[::-1] if dts[i] < 0.0]

    for idx_list in (fwd, bwd):
        if not idx_list:
            continue
        sim = _build_simulation(
            epoch_jd_tdb, test_state, perturbers, integrator, dt_days, asteroid_perturbers
        )
        test_idx = sim.N - 1
        for i in idx_list:
            sim.integrate(float(dts[i]))
            p = sim.particles[test_idx]
            pos_out[i] = (p.x, p.y, p.z)
            vel_out[i] = (p.vx, p.vy, p.vz)

    if return_velocity:
        return pos_out, vel_out
    return pos_out
