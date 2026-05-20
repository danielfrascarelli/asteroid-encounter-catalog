"""Forward astrometric model: orbital elements + perturber mass → Gaia RA/Dec.

This is the function that gets called inside the optimisation loop. Given the
target's osculating Kepler elements (6 parameters), the perturber's mass
(1 parameter), and the perturber's known elements, it returns the predicted
RA/Dec at every Gaia observation epoch.

Chain:
    elements + mass
    → N-body propagation (src/propagate/nbody_perturber)
    → heliocentric ecliptic positions at observation epochs
    → barycentric ICRS via heliocentric_to_barycentric_icrs
    → line of sight from Gaia (using x_gaia, y_gaia, z_gaia from SSO table)
    → light-time correction (iterative)
    → RA/Dec (in barycentric astrometric frame — Gaia removes aberration)

Convention: epochs in JD TDB, positions in AU, mass in kg.
"""

from __future__ import annotations

import numpy as np
from scipy.interpolate import CubicSpline

from src.astrometry.transforms import (
    heliocentric_to_barycentric_icrs,
    light_time_iterate,
    xyz_to_radec,
)
from src.propagate.nbody_perturber import propagate_target_with_perturber


def forward_model(
    target_elements: dict,
    perturber_elements: dict,
    perturber_mass_kg: float,
    obs_jd_tdb: np.ndarray,
    gaia_xyz_bary: np.ndarray,
    *,
    include_planets: tuple[str, ...] = (
        "sun", "mercury", "venus", "earth", "mars",
        "jupiter", "saturn", "uranus", "neptune",
    ),
    include_big4: bool = False,
    big4_elements: dict[str, dict] | None = None,
    dt_days: float = 1.0,
    integrator: str = "whfast",
    bracket_days: float = 2.0,
    n_bracket_points: int = 9,
) -> tuple[np.ndarray, np.ndarray]:
    """Predict (RA, Dec) at every Gaia observation in *obs_jd_tdb*.

    Internally integrates the target's orbit over a single contiguous time
    grid that covers all observation epochs plus a ``±bracket_days`` bracket
    on each end (needed by the light-time iteration, which evaluates the
    target's position at retarded times slightly before each obs epoch).

    Parameters
    ----------
    target_elements:
        MPCORB-style dict for the target.
    perturber_elements:
        MPCORB-style dict for the perturber.
    perturber_mass_kg:
        Free parameter: perturber's mass in kg.
    obs_jd_tdb:
        (N,) array of observation epochs in JD TDB.
    gaia_xyz_bary:
        (N, 3) Gaia barycentric ICRS positions at the same epochs (AU).
    include_planets:
        Planets to include in the N-body integration.
    include_big4:
        If True, include (1) Ceres, (2) Pallas, (4) Vesta, (10) Hygiea as
        massive perturbers.  Requires *big4_elements*.
    big4_elements:
        Dict keyed by lowercase name (``"ceres"``, ``"pallas"``, ``"vesta"``,
        ``"hygiea"``) with their MPCORB element dicts.  Required when
        ``include_big4=True``.
    dt_days:
        WHFast time step (days).
    bracket_days:
        Half-width of the time-grid extension on either side of the obs
        range, used for the light-time interpolation. ~25 minutes is the
        worst-case light-time at 3 AU; 2 days is plenty.
    n_bracket_points:
        Number of interior sample points spread over the obs range plus
        bracket, used to build the position spline. More points = finer
        interpolation but more N-body steps.

    Returns
    -------
    (ra_deg, dec_deg)
        Two (N,) arrays of predicted RA, Dec in degrees in the barycentric
        astrometric ICRS frame (= Gaia DR3 SSO frame).
    """
    obs_jd_tdb = np.asarray(obs_jd_tdb, dtype=float)
    gaia_xyz_bary = np.asarray(gaia_xyz_bary, dtype=float)
    if gaia_xyz_bary.shape != (len(obs_jd_tdb), 3):
        raise ValueError(
            f"gaia_xyz_bary shape {gaia_xyz_bary.shape} != ({len(obs_jd_tdb)}, 3)"
        )

    # Build the integration grid: a uniform set of nodes covering
    # [min(obs) − bracket, max(obs) + bracket] with at least n_bracket_points
    # nodes plus all observation epochs themselves (so we always have direct
    # samples at the obs epochs).
    t_lo = float(obs_jd_tdb.min()) - bracket_days
    t_hi = float(obs_jd_tdb.max()) + bracket_days

    n_uniform = max(n_bracket_points, int(np.ceil((t_hi - t_lo) / 5.0)))
    uniform = np.linspace(t_lo, t_hi, n_uniform)
    grid = np.unique(np.concatenate([uniform, obs_jd_tdb]))
    grid.sort()

    # Run N-body and get heliocentric ecliptic positions at every grid node.
    helio_ecl = propagate_target_with_perturber(
        target_elements=target_elements,
        perturber_elements=perturber_elements,
        perturber_mass_kg=perturber_mass_kg,
        time_grid_jd_tdb=grid,
        include_planets=include_planets,
        include_big4=include_big4,
        big4_elements=big4_elements,
        integrator=integrator,
        dt_days=dt_days,
    )

    # Convert to barycentric ICRS — heliocentric_to_barycentric_icrs handles
    # the rotation + Sun offset internally.  Done at grid nodes.
    bary_icrs = np.empty_like(helio_ecl)
    for i in range(len(grid)):
        bary_icrs[i] = np.asarray(
            heliocentric_to_barycentric_icrs(helio_ecl[i], float(grid[i]))
        ).reshape(3)

    # Spline the (T, 3) barycentric trajectory so we can evaluate at any time
    # for the light-time iteration.
    splines = [CubicSpline(grid, bary_icrs[:, k]) for k in range(3)]

    def target_pos_at(jd: float) -> np.ndarray:
        return np.array([s(jd) for s in splines])

    # For each observation, apply light-time and project to RA/Dec.
    ra = np.empty(len(obs_jd_tdb), dtype=float)
    dec = np.empty(len(obs_jd_tdb), dtype=float)
    for i, jd_obs in enumerate(obs_jd_tdb):
        pos, _tau = light_time_iterate(target_pos_at, float(jd_obs), gaia_xyz_bary[i])
        los = pos - gaia_xyz_bary[i]
        r, d = xyz_to_radec(los)
        ra[i] = r
        dec[i] = d
    return ra, dec


def residuals_mas(
    ra_obs: np.ndarray,
    dec_obs: np.ndarray,
    ra_pred: np.ndarray,
    dec_pred: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Signed residuals (obs − pred) in mas, on the tangent plane around pred."""
    deg = np.pi / 180.0
    dra = ((ra_obs - ra_pred + 540.0) % 360.0 - 180.0) * np.cos(dec_pred * deg) * 3_600_000.0
    ddec = (dec_obs - dec_pred) * 3_600_000.0
    return dra, ddec
