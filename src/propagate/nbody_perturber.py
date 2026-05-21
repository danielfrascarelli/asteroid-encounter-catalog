"""N-body propagator for a single (target, perturber) system.

Specialised wrapper around REBOUND for the mass-determination fit, where:

  - ONE target asteroid is propagated as a massless test particle
  - ONE perturber is included as a massive body with a configurable mass
  - The full set of "background" perturbers (Sun, Jupiter, Saturn, optionally
    the big-4 asteroids) is included so we can isolate the perturber's effect

The function is designed to be CHEAP to call inside an optimiser inner loop:
each call sets up its own ``rebound.Simulation``, integrates for the requested
time range, and returns the target's positions at the time grid.

Units & frame
-------------
- Positions returned in **heliocentric ecliptic J2000** (AU)
- Time grid in **JD TDB**
- Perturber mass in **kg** (converted internally to M_sun for REBOUND)

Typical run time for a ±180 day window: < 1 s with WHFast + 1-day step.
"""

from __future__ import annotations

import logging

import numpy as np

from src.propagate.nbody import (
    _MAJOR_ASTEROIDS,
    _PLANET_GMS,
    _heliocentric_kepler_state,
    _planet_state_at,
)

logger = logging.getLogger(__name__)


_KG_PER_MSUN = 1.989e30


def propagate_target_with_perturber(
    target_elements: dict,
    perturber_elements: dict,
    perturber_mass_kg: float,
    time_grid_jd_tdb: np.ndarray,
    *,
    include_planets: tuple[str, ...] = (
        "sun",
        "mercury",
        "venus",
        "earth",
        "mars",
        "jupiter",
        "saturn",
        "uranus",
        "neptune",
    ),
    include_background: bool = False,
    background_elements: dict[str, dict] | None = None,
    # Legacy aliases kept for backward compatibility
    include_big4: bool = False,
    big4_elements: dict[str, dict] | None = None,
    integrator: str = "whfast",
    dt_days: float = 1.0,
) -> np.ndarray:
    """Propagate *target_elements* under the gravity of *perturber* + planets.

    Parameters
    ----------
    target_elements:
        Dict with MPCORB columns: ``a_au``, ``e``, ``i_deg``, ``Omega_deg``,
        ``omega_deg``, ``M_deg``, ``epoch_jd``.
    perturber_elements:
        Same dict for the perturber.  Its mass is taken from
        *perturber_mass_kg*, NOT from any column.
    perturber_mass_kg:
        Mass of the perturber in kg. **This is the fit parameter.**
    time_grid_jd_tdb:
        (T,) array of epochs at which to record the target's position.
        Must be monotonically increasing (or decreasing).  Need not contain
        the integration start epoch.
    include_planets:
        Tuple of planet names to include as massive perturbers.  Sun is always
        added regardless.  Defaults to all 8 planets.
    include_background:
        If True, add massive background asteroids from *background_elements*
        with masses from ``_MAJOR_ASTEROIDS``.  Supersedes ``include_big4``.
    background_elements:
        Dict keyed by lowercase asteroid name with their MPCORB element dicts.
        Required when ``include_background=True``.
    include_big4:
        Deprecated alias for ``include_background``.
    big4_elements:
        Deprecated alias for ``background_elements``.
    integrator:
        REBOUND integrator (``"whfast"`` or ``"ias15"``).
    dt_days:
        Time step for symplectic integrators (days).

    Returns
    -------
    np.ndarray of shape (T, 3)
        Target's heliocentric ecliptic J2000 position (AU) at each grid point.

    Notes
    -----
    The perturber and target must share the same epoch_jd (i.e. come from the
    same MPCORB snapshot), so their orbital elements are osculating at the
    same instant.
    """
    import rebound

    # Resolve legacy aliases
    if include_big4 and not include_background:
        include_background = True
    if big4_elements is not None and background_elements is None:
        background_elements = big4_elements

    if perturber_mass_kg < 0:
        raise ValueError("perturber_mass_kg must be non-negative")

    epoch_jd = float(target_elements["epoch_jd"])
    if abs(float(perturber_elements["epoch_jd"]) - epoch_jd) > 1e-3:
        raise ValueError(
            "target and perturber must share the same osculating epoch; "
            f"got {target_elements['epoch_jd']} vs {perturber_elements['epoch_jd']}"
        )

    t_grid = np.asarray(time_grid_jd_tdb, dtype=float)
    if t_grid.ndim != 1:
        raise ValueError("time_grid_jd_tdb must be 1-D")

    # Build the simulation
    sim = rebound.Simulation()
    sim.units = ("AU", "days", "Msun")
    sim.integrator = integrator
    if integrator == "whfast":
        sim.dt = dt_days

    # Add planets (massive)
    planet_names = list(include_planets)
    if "sun" not in planet_names:
        planet_names = ["sun"] + planet_names
    for name in planet_names:
        pos, vel = _planet_state_at(name, epoch_jd)
        sim.add(
            m=_PLANET_GMS[name],
            x=float(pos[0]),
            y=float(pos[1]),
            z=float(pos[2]),
            vx=float(vel[0]),
            vy=float(vel[1]),
            vz=float(vel[2]),
        )

    sun_pos_bary, sun_vel_bary = _planet_state_at("sun", epoch_jd)

    # Add background asteroids if requested
    if include_background:
        if background_elements is None:
            raise ValueError("include_background=True requires background_elements dict")
        for name, (_number, gm_msun) in _MAJOR_ASTEROIDS.items():
            if name not in background_elements:
                continue  # silently skip bodies not supplied by caller
            p_helio, v_helio = _heliocentric_kepler_state(background_elements[name], epoch_jd)
            sim.add(
                m=gm_msun,
                x=float(p_helio[0] + sun_pos_bary[0]),
                y=float(p_helio[1] + sun_pos_bary[1]),
                z=float(p_helio[2] + sun_pos_bary[2]),
                vx=float(v_helio[0] + sun_vel_bary[0]),
                vy=float(v_helio[1] + sun_vel_bary[1]),
                vz=float(v_helio[2] + sun_vel_bary[2]),
            )

    # Add the perturber (massive, mass = fit parameter)
    p_pert, v_pert = _heliocentric_kepler_state(perturber_elements, epoch_jd)
    mass_msun = perturber_mass_kg / _KG_PER_MSUN
    sim.add(
        m=mass_msun,
        x=float(p_pert[0] + sun_pos_bary[0]),
        y=float(p_pert[1] + sun_pos_bary[1]),
        z=float(p_pert[2] + sun_pos_bary[2]),
        vx=float(v_pert[0] + sun_vel_bary[0]),
        vy=float(v_pert[1] + sun_vel_bary[1]),
        vz=float(v_pert[2] + sun_vel_bary[2]),
    )

    # Add the target (massless test particle)
    p_t, v_t = _heliocentric_kepler_state(target_elements, epoch_jd)
    # Use the live Sun particle to be consistent with whatever COM shift
    # REBOUND applied when adding the perturber
    p_sun_now = np.array([sim.particles[0].x, sim.particles[0].y, sim.particles[0].z])
    v_sun_now = np.array([sim.particles[0].vx, sim.particles[0].vy, sim.particles[0].vz])
    sim.add(
        m=0.0,
        x=float(p_t[0] + p_sun_now[0]),
        y=float(p_t[1] + p_sun_now[1]),
        z=float(p_t[2] + p_sun_now[2]),
        vx=float(v_t[0] + v_sun_now[0]),
        vy=float(v_t[1] + v_sun_now[1]),
        vz=float(v_t[2] + v_sun_now[2]),
    )
    sim.N_active = len(sim.particles) - 1  # all but the target are massive

    # Index of the target particle
    target_idx = len(sim.particles) - 1

    # Integrate to each grid point and record the target's heliocentric position
    out = np.empty((len(t_grid), 3), dtype=float)
    for i, t_target in enumerate(t_grid):
        sim.integrate(float(t_target) - epoch_jd)
        target_p = sim.particles[target_idx]
        sun_p = sim.particles[0]
        # Heliocentric ecliptic position
        helio_eq = np.array([target_p.x - sun_p.x, target_p.y - sun_p.y, target_p.z - sun_p.z])
        # State vectors in REBOUND are in ecliptic frame (we set up planets in
        # ecliptic from _planet_state_at), so no rotation needed.
        out[i] = helio_eq

    return out


def positions_to_barycentric_icrs(
    helio_ecl_positions: np.ndarray,
    time_grid_jd_tdb: np.ndarray,
) -> np.ndarray:
    """Convert (T, 3) heliocentric ecliptic positions to barycentric ICRS.

    Adds Sun's barycentric position at each epoch and rotates ecliptic →
    equatorial.
    """
    from src.astrometry.transforms import heliocentric_to_barycentric_icrs

    out = np.empty_like(helio_ecl_positions)
    for i, (pos, jd) in enumerate(zip(helio_ecl_positions, time_grid_jd_tdb, strict=True)):
        out[i] = np.asarray(heliocentric_to_barycentric_icrs(pos, float(jd))).reshape(3)
    return out
