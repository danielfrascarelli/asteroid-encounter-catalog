"""N-body propagator using REBOUND (Sun + planets ± major asteroids as perturbers).

The asteroids of interest are added as massless test particles, so they do not
perturb each other or the massive bodies — this keeps the integration
embarrassingly cheap *per asteroid* once the heavy bodies' trajectory is
established by the symplectic integrator.

Units & frame
-------------
- Distances : AU
- Times     : JD TDB internally; REBOUND's clock is days since the snapshot epoch.
- Frame     : heliocentric mean ecliptic of J2000.  Massive-body initial state
              vectors are obtained from the JPL DE built-in ephemeris bundled
              with astropy (``solar_system_ephemeris.set('builtin')``) — fully
              offline, low-precision but more than adequate for the perturber
              role on a 3-year window.  ICRS → ecliptic via the J2000 obliquity
              ε = 23.4392911° (IAU 2006).

Integrator
----------
WHFast (symplectic Wisdom-Holman) with a 1-day time step is used by default.
The configuration's ``rebound.integrator`` key is honoured if provided
(``ias15`` is also a valid choice for the asteroid-only sanity tests).

Public API
----------
- :func:`propagate_grid_nbody` — full-grid integration, returns ``(T, N, 3)``.
- :func:`propagate_grid_iter`  — iterator wrapper that yields ``(t, pos)`` per
  step, matching :func:`src.propagate.grid.propagate_grid` semantics.

Phase-3 hook
------------
When ``include_major_asteroids`` is true, (1) Ceres, (2) Pallas, (4) Vesta and
(10) Hygiea are promoted from massless test particles to massive members of the
integration using GMs from the JPL Small-Body Database (units of M_sun).
"""

from __future__ import annotations

import logging
from collections.abc import Iterator

import numpy as np
import polars as pl
from astropy import units as u
from astropy.coordinates import get_body_barycentric_posvel, solar_system_ephemeris
from astropy.time import Time

from src.propagate.kepler import kepler_to_cartesian

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_DEG = np.pi / 180.0

# Obliquity of the ecliptic at J2000 (IAU 2006), used to rotate ICRS → ecliptic.
_EPS_J2000 = 23.4392911 * _DEG
_COS_EPS = float(np.cos(_EPS_J2000))
_SIN_EPS = float(np.sin(_EPS_J2000))

# ICRS → ecliptic rotation matrix (3x3).  Pre-computed to avoid per-call cost.
# Acts as: p_ecl = R @ p_icrs.
_ICRS_TO_ECL: np.ndarray = np.array(
    [
        [1.0, 0.0, 0.0],
        [0.0, _COS_EPS, _SIN_EPS],
        [0.0, -_SIN_EPS, _COS_EPS],
    ]
)

# Planet masses in M_sun (IAU 2015 / DE440 nominal values).  Used when the
# corresponding planet is requested in ``include_planets``.  Names are the
# lowercase strings expected from config.
_PLANET_GMS: dict[str, float] = {
    "sun": 1.0,
    "mercury": 1.660_113_797e-7,
    "venus": 2.447_835_335e-6,
    "earth": 3.040_432_648e-6,  # Earth+Moon (Earth-Moon barycenter)
    "mars": 3.227_156_038e-7,
    "jupiter": 9.547_919_152e-4,
    "saturn": 2.858_859_807e-4,
    "uranus": 4.366_244_043e-5,
    "neptune": 5.151_389_021e-5,
}

# Major asteroid GMs in M_sun (JPL SBDB, M_sun = 1.989e30 kg).
# Ceres GM = 62.628 km^3/s^2 → 4.72e-10 M_sun, etc.  Source comment: SBDB
# snapshot 2024.  Used only when include_major_asteroids=True.
_MAJOR_ASTEROIDS: dict[str, tuple[int, float]] = {
    # name: (MPC number, mass in M_sun)
    "ceres": (1, 4.72e-10),
    "pallas": (2, 1.06e-10),
    "vesta": (4, 1.30e-10),
    "hygiea": (10, 4.40e-11),
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _icrs_to_ecliptic(vec: np.ndarray) -> np.ndarray:
    """Rotate an ICRS-equatorial 3-vector to ecliptic J2000."""
    return np.asarray(_ICRS_TO_ECL @ vec)


def _planet_state_at(body: str, epoch_jd_tdb: float) -> tuple[np.ndarray, np.ndarray]:
    """Return (pos_AU, vel_AU_per_day) of *body* in barycentric ecliptic J2000.

    Uses astropy's built-in low-precision JPL ephemeris (no network, no files
    beyond the wheels).  Inputs are in JD TDB.

    Parameters
    ----------
    body:
        Body name accepted by ``astropy.coordinates.get_body_barycentric_posvel``
        ("sun", "mercury", ..., "neptune").
    epoch_jd_tdb:
        Epoch in JD TDB.

    Returns
    -------
    (pos_au, vel_au_per_day)
        Both arrays of shape (3,).
    """
    t = Time(epoch_jd_tdb, format="jd", scale="tdb")
    with solar_system_ephemeris.set("builtin"):
        pos, vel = get_body_barycentric_posvel(body, t)
    p_icrs = pos.xyz.to(u.AU).value
    v_icrs = vel.xyz.to(u.AU / u.day).value
    return _icrs_to_ecliptic(p_icrs), _icrs_to_ecliptic(v_icrs)


def _heliocentric_kepler_state(
    row: dict,
    epoch_jd: float,
    *,
    delta_days: float = 1.0e-3,
) -> tuple[np.ndarray, np.ndarray]:
    """Initial heliocentric ecliptic state vector for a Keplerian orbit.

    The position is taken from :func:`kepler_to_cartesian` at *epoch_jd*; the
    velocity is approximated by a centred finite difference over ±*delta_days*.
    This avoids re-deriving the analytic velocity expression while keeping
    accuracy near machine precision for ~mean-belt orbits.

    Parameters
    ----------
    row:
        Dict with keys ``a_au, e, i_deg, Omega_deg, omega_deg, M_deg, epoch_jd``.
    epoch_jd:
        Target epoch (JD TDB) at which to evaluate position and velocity.
    delta_days:
        Half-step (days) for the centred finite-difference velocity.

    Returns
    -------
    (pos, vel)
        Heliocentric ecliptic J2000 position (AU) and velocity (AU/day).
    """
    args = dict(
        a_au=row["a_au"],
        e=row["e"],
        i_rad=float(row["i_deg"]) * _DEG,
        Omega_rad=float(row["Omega_deg"]) * _DEG,
        omega_rad=float(row["omega_deg"]) * _DEG,
        M0_rad=float(row["M_deg"]) * _DEG,
        epoch_jd=row["epoch_jd"],
    )
    p_minus = kepler_to_cartesian(t_jd=epoch_jd - delta_days, **args).reshape(3)
    p_plus = kepler_to_cartesian(t_jd=epoch_jd + delta_days, **args).reshape(3)
    p_now = kepler_to_cartesian(t_jd=epoch_jd, **args).reshape(3)
    vel = (p_plus - p_minus) / (2.0 * delta_days)
    return p_now, vel


def _vectorised_kepler_state(
    elements: pl.DataFrame,
    epoch_jd: float,
    *,
    delta_days: float = 1.0e-3,
) -> tuple[np.ndarray, np.ndarray]:
    """Vectorised analogue of :func:`_heliocentric_kepler_state` over a DataFrame.

    Returns
    -------
    (pos, vel)
        Arrays of shape (N, 3).
    """
    a = elements["a_au"].to_numpy()
    e = elements["e"].to_numpy()
    i = elements["i_deg"].to_numpy() * _DEG
    cap_o = elements["Omega_deg"].to_numpy() * _DEG
    o = elements["omega_deg"].to_numpy() * _DEG
    m0 = elements["M_deg"].to_numpy() * _DEG
    ep = elements["epoch_jd"].to_numpy()

    def _pos(t: float) -> np.ndarray:
        return kepler_to_cartesian(
            a_au=a,
            e=e,
            i_rad=i,
            Omega_rad=cap_o,
            omega_rad=o,
            M0_rad=m0,
            epoch_jd=ep,
            t_jd=t,
        )

    p_now = _pos(epoch_jd)
    p_plus = _pos(epoch_jd + delta_days)
    p_minus = _pos(epoch_jd - delta_days)
    vel = (p_plus - p_minus) / (2.0 * delta_days)
    return p_now, vel


# ---------------------------------------------------------------------------
# Public propagator
# ---------------------------------------------------------------------------


def propagate_grid_nbody(
    elements: pl.DataFrame,
    time_grid: np.ndarray,
    *,
    include_planets: list[str] | None = None,
    include_major_asteroids: bool = False,
    integrator: str = "whfast",
    dt_days: float = 1.0,
    epoch_jd: float | None = None,
    out: np.ndarray | None = None,
) -> np.ndarray:
    """Integrate *elements* across *time_grid* with planetary perturbers.

    Parameters
    ----------
    elements:
        Orbital elements DataFrame (same schema as
        :func:`src.propagate.kepler.propagate_df`).
    time_grid:
        JD TDB epochs at which to record heliocentric positions.  Must be
        monotonically increasing.
    include_planets:
        Names of massive bodies to include from ``_PLANET_GMS``.  ``"sun"`` is
        always added regardless of whether it appears in the list.  Defaults to
        ``["sun", "jupiter", "saturn"]`` (Phase 1 minimum).
    include_major_asteroids:
        If true, promote (1) Ceres, (2) Pallas, (4) Vesta and (10) Hygiea from
        test particles to massive perturbers.  When they also appear in
        *elements* they are still recorded in the output.
    integrator:
        REBOUND integrator name.  Defaults to ``"whfast"`` (symplectic, fast).
        ``"ias15"`` is also supported for higher accuracy at extra cost.
    dt_days:
        Time step for the symplectic integrator (days).  Ignored by
        adaptive integrators like IAS15.
    epoch_jd:
        Initial epoch for the integration in JD TDB.  Defaults to the
        ``epoch_jd`` column of *elements* (which is uniform within a single
        MPCORB snapshot).

    Returns
    -------
    np.ndarray of shape (T, N, 3), float32
        Heliocentric ecliptic J2000 positions in AU, one slab per time step.

    Notes
    -----
    Memory: ``T × N × 3 × 4 B``.  At T=25 000, N=100 000 this is ~30 GB —
    use :mod:`src.propagate.cache` to spill to disk, or restrict the
    subset (config ``subset.max_asteroids``) for in-memory runs.
    """
    import rebound

    if include_planets is None:
        include_planets = ["sun", "jupiter", "saturn"]

    n_ast = len(elements)
    n_steps = len(time_grid)
    if n_ast == 0:
        return np.empty((n_steps, 0, 3), dtype=np.float32)
    if n_steps == 0:
        return np.empty((0, n_ast, 3), dtype=np.float32)

    if epoch_jd is None:
        # All asteroids in a single MPCORB snapshot share an epoch.  Validate.
        eps = elements["epoch_jd"].unique().to_list()
        if len(eps) != 1:
            raise ValueError(
                f"N-body propagation requires a single epoch; got {len(eps)} distinct values"
            )
        epoch_jd = float(eps[0])

    if not np.all(np.diff(time_grid) > 0):
        raise ValueError("time_grid must be strictly monotonically increasing")

    # ---------------------------------------------------------------------
    # Assemble massive bodies and initial barycentric state vectors
    # ---------------------------------------------------------------------
    planet_names = [p.lower() for p in include_planets]
    if "sun" not in planet_names:
        planet_names = ["sun"] + planet_names

    unknown = [p for p in planet_names if p not in _PLANET_GMS]
    if unknown:
        raise ValueError(f"Unknown planet names in include_planets: {unknown}")

    logger.info(
        "N-body propagation: %d asteroids | %d steps | dt=%.3f d | integrator=%s | planets=%s | major_asteroids=%s",
        n_ast,
        n_steps,
        dt_days,
        integrator,
        ",".join(planet_names),
        include_major_asteroids,
    )

    sim = rebound.Simulation()
    sim.units = ("AU", "days", "Msun")
    sim.integrator = integrator
    if integrator == "whfast":
        sim.dt = dt_days
    # Add planets with barycentric state vectors from astropy at epoch_jd
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
    n_planets = len(planet_names)

    # Major asteroids as massive perturbers (Phase 3).  Initial state from
    # MPCORB Keplerian elements converted to heliocentric Cartesian, then
    # boosted to barycentric using the Sun's state.
    sun_pos_bary, sun_vel_bary = _planet_state_at("sun", epoch_jd)
    major_ast_indices: list[int] = []  # row indices in *elements* (if present)
    if include_major_asteroids:
        for name, (number, gm) in _MAJOR_ASTEROIDS.items():
            # Find this body in elements DF (by number)
            mask = (elements["number"] == number).to_numpy()
            if not mask.any():
                logger.warning(
                    "Major asteroid (%d) %s requested but not in elements DF; skipping",
                    number,
                    name,
                )
                continue
            row_idx = int(np.argmax(mask))
            row = {col: elements[col][row_idx] for col in elements.columns}
            p_helio, v_helio = _heliocentric_kepler_state(row, epoch_jd)
            p_bary = p_helio + sun_pos_bary
            v_bary = v_helio + sun_vel_bary
            sim.add(
                m=gm,
                x=float(p_bary[0]),
                y=float(p_bary[1]),
                z=float(p_bary[2]),
                vx=float(v_bary[0]),
                vy=float(v_bary[1]),
                vz=float(v_bary[2]),
            )
            major_ast_indices.append(row_idx)
            logger.info("  Added massive perturber (%d) %s  GM=%.2e Msun", number, name, gm)

    # ---------------------------------------------------------------------
    # Add asteroids as massless test particles, in their elements-DF order
    # ---------------------------------------------------------------------
    p_helio_all, v_helio_all = _vectorised_kepler_state(elements, epoch_jd)
    # Boost to barycentric using the integrator's current Sun (index 0).  We
    # *cannot* reuse sun_pos_bary directly because rebound shifts to the
    # centre-of-momentum frame after adding the first massive particle in some
    # configurations; using the live particle state guarantees consistency.
    p_sun_now = np.array([sim.particles[0].x, sim.particles[0].y, sim.particles[0].z])
    v_sun_now = np.array([sim.particles[0].vx, sim.particles[0].vy, sim.particles[0].vz])

    for k in range(n_ast):
        if include_major_asteroids and k in major_ast_indices:
            # Already added as a massive body; we'll fill its trajectory from
            # the corresponding massive-particle slot at output time.
            continue
        p = p_helio_all[k] + p_sun_now
        v = v_helio_all[k] + v_sun_now
        sim.add(
            m=0.0,
            x=float(p[0]),
            y=float(p[1]),
            z=float(p[2]),
            vx=float(v[0]),
            vy=float(v[1]),
            vz=float(v[2]),
        )

    # Massive bodies = planets + (major asteroids if any).  Test particles
    # follow them.  REBOUND uses N_active to mark how many particles act as
    # gravitational sources.
    n_active = n_planets + (len(major_ast_indices) if include_major_asteroids else 0)
    sim.N_active = n_active

    # Map elements-DF index → REBOUND particle index
    #   massive perturbers occupy [0 .. n_active)
    #   test particles occupy [n_active .. n_active + (n_ast - len(major_ast_indices)))
    # We reconstruct the original elements ordering via this lookup table.
    rebound_index = np.empty(n_ast, dtype=np.int64)
    # Map major asteroids → their massive slot (planets first, then majors in order)
    if include_major_asteroids and major_ast_indices:
        for slot_offset, df_idx in enumerate(major_ast_indices):
            rebound_index[df_idx] = n_planets + slot_offset
    # Fill the rest in DF order, skipping major asteroids that are already mapped
    next_slot = n_active
    for k in range(n_ast):
        if include_major_asteroids and k in set(major_ast_indices):
            continue
        rebound_index[k] = next_slot
        next_slot += 1

    # ---------------------------------------------------------------------
    # Integrate and snapshot heliocentric positions
    # ---------------------------------------------------------------------
    # REBOUND's clock starts at sim.t = 0 corresponding to epoch_jd.
    sim.t = 0.0
    # WHFast safe-mode synchronisation must be on for accurate intermediate
    # output between sim.integrate(t) calls; mandatory because the KD-tree
    # scan samples the trajectory at every grid step (not just the endpoint).
    if integrator == "whfast":
        sim.ri_whfast.safe_mode = 1

    if out is None:
        out = np.empty((n_steps, n_ast, 3), dtype=np.float32)
    else:
        expected = (n_steps, n_ast, 3)
        if out.shape != expected:
            raise ValueError(f"out has shape {out.shape}, expected {expected}")
        if out.dtype != np.float32:
            raise ValueError(f"out has dtype {out.dtype}, expected float32")

    # Snapshot the initial state so we can re-run the simulation backward
    # without rebuilding it from scratch.  REBOUND stores particle state by
    # reference; freezing into ndarrays is what we need.
    init_xyz = np.empty((sim.N, 3), dtype=np.float64)
    init_vxvyvz = np.empty((sim.N, 3), dtype=np.float64)
    sim.serialize_particle_data(xyz=init_xyz, vxvyvz=init_vxvyvz)

    def _snapshot(idx: int) -> None:
        sun_p = np.array([sim.particles[0].x, sim.particles[0].y, sim.particles[0].z])
        all_pos = np.empty((sim.N, 3), dtype=np.float64)
        sim.serialize_particle_data(xyz=all_pos)
        helio = all_pos - sun_p
        out[idx, :, :] = helio[rebound_index].astype(np.float32)

    from tqdm import tqdm

    # --- Forward portion (t_jd >= epoch_jd) ---
    fwd_mask = time_grid >= epoch_jd - 1e-9
    bwd_mask = ~fwd_mask
    fwd_steps = np.where(fwd_mask)[0]
    bwd_steps = np.where(bwd_mask)[0]

    if fwd_steps.size > 0:
        if integrator == "whfast":
            sim.dt = dt_days
        for step_idx in tqdm(fwd_steps, desc="N-body integrate (fwd)", unit="step", leave=False):
            t_sim = float(time_grid[step_idx]) - epoch_jd
            if t_sim > sim.t:
                sim.integrate(t_sim)
            _snapshot(int(step_idx))

    # --- Backward portion (t_jd < epoch_jd): re-init and integrate backward ---
    if bwd_steps.size > 0:
        # Reset to initial state
        for k in range(sim.N):
            p = sim.particles[k]
            p.x, p.y, p.z = float(init_xyz[k, 0]), float(init_xyz[k, 1]), float(init_xyz[k, 2])
            p.vx, p.vy, p.vz = (
                float(init_vxvyvz[k, 0]),
                float(init_vxvyvz[k, 1]),
                float(init_vxvyvz[k, 2]),
            )
        sim.t = 0.0
        if integrator == "whfast":
            sim.dt = -abs(dt_days)
            sim.ri_whfast.recalculate_coordinates_this_timestep = 1
        # Walk backward steps in decreasing time order (most-negative last)
        for step_idx in tqdm(
            bwd_steps[::-1], desc="N-body integrate (bwd)", unit="step", leave=False
        ):
            t_sim = float(time_grid[step_idx]) - epoch_jd
            if t_sim < sim.t:
                sim.integrate(t_sim)
            _snapshot(int(step_idx))

    logger.info("N-body propagation finished: output shape %s", out.shape)
    return out


def propagate_grid_iter(
    elements: pl.DataFrame,
    time_grid: np.ndarray,
    *,
    include_planets: list[str] | None = None,
    include_major_asteroids: bool = False,
    integrator: str = "whfast",
    dt_days: float = 1.0,
    epoch_jd: float | None = None,
) -> Iterator[tuple[float, np.ndarray]]:
    """Iterator-style adapter over :func:`propagate_grid_nbody`.

    Yields ``(t_jd, positions_NxA)`` for each grid step, matching the contract
    of :func:`src.propagate.grid.propagate_grid` so the KD-tree scan code can
    consume either propagator.

    Note
    ----
    Memory-greedy: the full ``(T, N, 3)`` array is materialised before
    iteration begins.  For large runs use the cache module to memmap from
    disk instead.
    """
    positions = propagate_grid_nbody(
        elements,
        time_grid,
        include_planets=include_planets,
        include_major_asteroids=include_major_asteroids,
        integrator=integrator,
        dt_days=dt_days,
        epoch_jd=epoch_jd,
    )
    for k, t in enumerate(time_grid):
        yield float(t), positions[k]
