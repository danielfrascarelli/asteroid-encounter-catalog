"""Per-pair N-body refinement of a Kepler-detected close encounter.

Public entry point :func:`refine_pair_nbody` takes the orbital elements of
two asteroids (typically as parsed from MPCORB) and a Kepler-refined
encounter epoch, and returns the *N-body* minimum-distance, epoch and
relative velocity inside a ±``window_hours`` window around that epoch.

Setup
-----
A REBOUND simulation is initialised at the MPCORB epoch with massive
bodies = ``[Sun, Jupiter, Saturn]`` (optionally Ceres/Pallas/Vesta/Hygiea
if ``include_major_asteroids=True``).  Both target asteroids are added as
massless test particles using heliocentric ecliptic state vectors derived
from ``kepler_to_cartesian`` and boosted to the integrator's frame.

Two-phase integration
---------------------
1. **Warmup**:  WHFast with ``dt = warmup_dt_seconds`` from the MPCORB
   epoch to ``t_min − window_hours``.  Symplectic, no energy drift, fast
   over multi-year spans.
2. **Refine window**:  Switch to IAS15 (high-order adaptive) and
   integrate forward through ±``window_hours`` with positions snapshot
   every ``sample_dt_seconds``.  IAS15 chooses its own internal step;
   the sampler just calls ``sim.integrate(t)`` at requested epochs.

The minimum of ``dist(t)`` over the sampled grid is then refined sub-step
by fitting a parabola to its three neighbours.

Diagnostics returned
--------------------
``RefinementResult`` carries the refined distance, epoch, relative-
velocity, the integration's energy drift, a ``converged`` boolean and
the residual flagged by IAS15.

Notes
-----
- Frame: heliocentric ecliptic J2000.  All inputs and outputs follow
  this convention.
- Units: AU, days (JD), AU/day.
- The implementation is intentionally per-pair (one REBOUND simulation
  per call); it is meant to be driven from a multiprocessing pool by
  ``compare_kepler_vs_nbody.py``.

Usage (CLI smoke test)
----------------------
    docker compose run --rm pipeline python -m scripts.validate.refine_pair_nbody \\
        --num1 1 --num2 4 --jd 2457800.5
"""

from __future__ import annotations

import argparse
import logging
from dataclasses import dataclass

import numpy as np

from src.propagate.kepler import kepler_to_cartesian
from src.propagate.nbody import (
    _MAJOR_ASTEROIDS,
    _PLANET_GMS,
    _planet_state_at,
)

logger = logging.getLogger(__name__)

_DEG = np.pi / 180.0


@dataclass(frozen=True)
class RefinementResult:
    """Result of an N-body refinement on a single pair."""

    dist_au_nbody: float
    t_min_nbody_jd: float
    rel_vel_au_day: float
    energy_drift: float
    converged: bool
    n_samples: int


def _heliocentric_state(
    *,
    a_au: float,
    e: float,
    i_deg: float,
    Omega_deg: float,
    omega_deg: float,
    M_deg: float,
    epoch_jd: float,
    delta_days: float = 1.0e-3,
) -> tuple[np.ndarray, np.ndarray]:
    """Return (position, velocity) at ``epoch_jd`` (heliocentric ecliptic, AU & AU/day)."""
    kwargs = dict(
        a_au=a_au,
        e=e,
        i_rad=float(i_deg) * _DEG,
        Omega_rad=float(Omega_deg) * _DEG,
        omega_rad=float(omega_deg) * _DEG,
        M0_rad=float(M_deg) * _DEG,
        epoch_jd=epoch_jd,
    )
    p_now = kepler_to_cartesian(t_jd=epoch_jd, **kwargs).reshape(3)
    p_plus = kepler_to_cartesian(t_jd=epoch_jd + delta_days, **kwargs).reshape(3)
    p_minus = kepler_to_cartesian(t_jd=epoch_jd - delta_days, **kwargs).reshape(3)
    vel = (p_plus - p_minus) / (2.0 * delta_days)
    return p_now, vel


def _parabolic_min(t: np.ndarray, d: np.ndarray) -> tuple[float, float]:
    """Locate the parabolic minimum of (t, d).  Returns (t_min, d_min).

    Falls back to the discrete argmin if the discrete minimum is at the
    endpoints (no valid 3-point neighbourhood).
    """
    k = int(np.argmin(d))
    if k == 0 or k == len(d) - 1:
        return float(t[k]), float(d[k])
    x1, x2, x3 = t[k - 1], t[k], t[k + 1]
    y1, y2, y3 = d[k - 1], d[k], d[k + 1]
    denom = (x1 - x2) * (x1 - x3) * (x2 - x3)
    if denom == 0.0:
        return float(x2), float(y2)
    a = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom
    b = (x3 * x3 * (y1 - y2) + x2 * x2 * (y3 - y1) + x1 * x1 * (y2 - y3)) / denom
    if a <= 0.0:
        # Not a minimum (concave-down) — fall back to the discrete value.
        return float(x2), float(y2)
    t_star = -b / (2.0 * a)
    if t_star < x1 or t_star > x3:
        return float(x2), float(y2)
    c = y2 - a * x2 * x2 - b * x2
    d_star = a * t_star * t_star + b * t_star + c
    return float(t_star), float(d_star)


def refine_pair_nbody(
    *,
    elements_1: dict,
    elements_2: dict,
    t_center_jd: float,
    window_hours: float = 6.0,
    sample_dt_seconds: float = 60.0,
    warmup_dt_seconds: float = 600.0,
    include_planets: tuple[str, ...] = ("sun", "jupiter", "saturn"),
    include_major_asteroids: bool = True,
) -> RefinementResult:
    """Refine the minimum-distance, epoch and relative velocity of one pair under N-body.

    Parameters
    ----------
    elements_1, elements_2:
        Dicts with keys ``a_au, e, i_deg, Omega_deg, omega_deg, M_deg,
        epoch_jd`` (MPCORB row).  Both must share the same ``epoch_jd``.
        An optional ``number`` key (MPC catalog number) is honoured: if
        present, and it matches one of the major asteroids that would be
        added as a perturber, that major is excluded from the perturber
        list so the target is not double-counted as both a massive body
        and a test particle.
    t_center_jd:
        Centre of the refinement window (typically the Kepler-refined
        encounter epoch, JD TDB).
    window_hours:
        Half-width of the IAS15 refinement window (hours, default 6 h).
    sample_dt_seconds:
        Snapshot cadence inside the window (default 60 s).
    warmup_dt_seconds:
        WHFast step for the warmup phase from MPCORB epoch to window
        start (default 600 s).
    include_planets:
        Massive perturbers (lowercase names from ``_PLANET_GMS``).
    include_major_asteroids:
        If true, promote Ceres/Pallas/Vesta/Hygiea to massive bodies.

    Returns
    -------
    RefinementResult
    """
    import rebound

    if abs(elements_1["epoch_jd"] - elements_2["epoch_jd"]) > 1e-6:
        raise ValueError(
            f"Asteroids must share MPCORB epoch; got {elements_1['epoch_jd']} vs "
            f"{elements_2['epoch_jd']}"
        )
    epoch_jd = float(elements_1["epoch_jd"])

    window_days = window_hours / 24.0
    sample_dt_days = sample_dt_seconds / 86400.0
    warmup_dt_days = warmup_dt_seconds / 86400.0

    t_lo = float(t_center_jd) - window_days
    t_hi = float(t_center_jd) + window_days

    # --- Build simulation at MPCORB epoch ---
    sim = rebound.Simulation()
    sim.units = ("AU", "days", "Msun")
    sim.integrator = "whfast"
    sim.dt = warmup_dt_days if t_lo >= epoch_jd else -warmup_dt_days

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

    # If either target is itself one of the major asteroids, exclude it from
    # the perturber loop — otherwise it would appear twice (once as a massive
    # body, once as a test particle) and the refined distance is garbage.
    target_numbers: set[int] = set()
    for ele in (elements_1, elements_2):
        n = ele.get("number")
        if n is not None:
            target_numbers.add(int(n))

    if include_major_asteroids:
        _load_major_elements()
        for name, (number, gm) in _MAJOR_ASTEROIDS.items():
            if name not in {"ceres", "pallas", "vesta", "hygiea"}:
                continue
            if number in target_numbers:
                continue
            elem = _MAJOR_ELEMENTS[name]
            p_helio, v_helio = _heliocentric_state(**elem, epoch_jd=epoch_jd)
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

    n_active = sim.N
    sim.N_active = n_active

    # Re-grab Sun's current state (rebound may have shifted to the COM frame).
    p_sun_now = np.array([sim.particles[0].x, sim.particles[0].y, sim.particles[0].z])
    v_sun_now = np.array([sim.particles[0].vx, sim.particles[0].vy, sim.particles[0].vz])

    for ele in (elements_1, elements_2):
        p_h, v_h = _heliocentric_state(
            a_au=ele["a_au"],
            e=ele["e"],
            i_deg=ele["i_deg"],
            Omega_deg=ele["Omega_deg"],
            omega_deg=ele["omega_deg"],
            M_deg=ele["M_deg"],
            epoch_jd=epoch_jd,
        )
        p = p_h + p_sun_now
        v = v_h + v_sun_now
        sim.add(
            m=0.0,
            x=float(p[0]),
            y=float(p[1]),
            z=float(p[2]),
            vx=float(v[0]),
            vy=float(v[1]),
            vz=float(v[2]),
        )
    idx_1 = n_active
    idx_2 = n_active + 1

    e_init = sim.energy()

    # --- Warmup phase: WHFast from epoch to t_lo (in either direction) ---
    sim.t = 0.0
    sim.did_modify_particles = 1
    target_warmup_t = t_lo - epoch_jd
    if abs(target_warmup_t) > 1e-9:
        sim.integrate(target_warmup_t)

    # --- Refine phase: switch to IAS15, sample at 60 s cadence ---
    sim.integrator = "ias15"
    # IAS15 is adaptive, dt is just an initial guess.
    sim.dt = sample_dt_days

    sample_times_jd = np.arange(t_lo, t_hi + 0.5 * sample_dt_days, sample_dt_days)
    n_samples = len(sample_times_jd)
    dists = np.empty(n_samples, dtype=np.float64)
    rel_vel = np.empty(n_samples, dtype=np.float64)

    for k, t_jd in enumerate(sample_times_jd):
        sim.integrate(t_jd - epoch_jd)
        p1 = sim.particles[idx_1]
        p2 = sim.particles[idx_2]
        dx = p1.x - p2.x
        dy = p1.y - p2.y
        dz = p1.z - p2.z
        dvx = p1.vx - p2.vx
        dvy = p1.vy - p2.vy
        dvz = p1.vz - p2.vz
        dists[k] = float(np.sqrt(dx * dx + dy * dy + dz * dz))
        rel_vel[k] = float(np.sqrt(dvx * dvx + dvy * dvy + dvz * dvz))

    t_star_jd, d_star = _parabolic_min(sample_times_jd, dists)

    # Use the rel-vel at the discrete-min sample (sub-step variation is
    # below the precision we care about for this characterisation).
    k_min = int(np.argmin(dists))
    rv_at_min = float(rel_vel[k_min])

    e_final = sim.energy()
    drift = float(abs(e_final - e_init) / abs(e_init)) if e_init != 0.0 else float("nan")
    converged = bool(np.isfinite(d_star) and drift < 1e-7)

    return RefinementResult(
        dist_au_nbody=float(d_star),
        t_min_nbody_jd=float(t_star_jd),
        rel_vel_au_day=rv_at_min,
        energy_drift=drift,
        converged=converged,
        n_samples=n_samples,
    )


# Static heliocentric Keplerian elements of the major asteroids at MPCORB
# 2016-02-17 epoch (snapshot used by the frozen catalog).  Read once at
# import time so the refiner doesn't need to parse MPCORB for every call.
#
# Source: data/raw/mpcorb_archive/MPCORB_20160217.DAT
_MAJOR_ELEMENTS: dict[str, dict] = {
    # Filled lazily on first import — see _load_major_elements() below.
}


def _load_major_elements() -> None:
    """Populate ``_MAJOR_ELEMENTS`` from the MPCORB snapshot of the frozen catalog."""
    if _MAJOR_ELEMENTS:
        return
    from pathlib import Path

    from src.ingest.mpcorb import parse_mpcorb

    snapshot = Path("data/raw/mpcorb_archive/MPCORB_20160217.DAT")
    elements = parse_mpcorb(snapshot, only_numbered=True)
    numbers = {"ceres": 1, "pallas": 2, "vesta": 4, "hygiea": 10}
    for name, n in numbers.items():
        row = elements.filter(elements["number"] == n)
        if len(row) != 1:
            raise RuntimeError(f"Major asteroid '{name}' (#{n}) not found in {snapshot}")
        d = row.to_dicts()[0]
        _MAJOR_ELEMENTS[name] = {
            "a_au": d["a_au"],
            "e": d["e"],
            "i_deg": d["i_deg"],
            "Omega_deg": d["Omega_deg"],
            "omega_deg": d["omega_deg"],
            "M_deg": d["M_deg"],
        }


def _cli() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--num1", type=int, required=True, help="Asteroid #1 MPC number")
    parser.add_argument("--num2", type=int, required=True, help="Asteroid #2 MPC number")
    parser.add_argument("--jd", type=float, required=True, help="Kepler-refined encounter JD TDB")
    parser.add_argument("--window-hours", type=float, default=6.0)
    parser.add_argument("--sample-dt-seconds", type=float, default=60.0)
    parser.add_argument("--no-major-asteroids", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s"
    )

    from pathlib import Path

    from src.ingest.mpcorb import parse_mpcorb

    elements = parse_mpcorb(Path("data/raw/mpcorb_archive/MPCORB_20160217.DAT"))
    rows = {int(r["number"]): r for r in elements.to_dicts()}
    if args.num1 not in rows:
        raise SystemExit(f"#{args.num1} not in MPCORB")
    if args.num2 not in rows:
        raise SystemExit(f"#{args.num2} not in MPCORB")

    result = refine_pair_nbody(
        elements_1=rows[args.num1],
        elements_2=rows[args.num2],
        t_center_jd=args.jd,
        window_hours=args.window_hours,
        sample_dt_seconds=args.sample_dt_seconds,
        include_major_asteroids=not args.no_major_asteroids,
    )
    logger.info("Result: %s", result)
    return 0


if __name__ == "__main__":
    raise SystemExit(_cli())
