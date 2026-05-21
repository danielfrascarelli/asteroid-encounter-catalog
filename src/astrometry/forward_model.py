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
    xyz_to_radec,
)
from src.propagate.nbody_perturber import propagate_target_with_perturber

# Speed of light in AU/day — must match the value in transforms.py
_C_AU_DAY = 173.144632674


def forward_model(
    target_elements: dict,
    perturber_elements: dict,
    perturber_mass_kg: float,
    obs_jd_tdb: np.ndarray,
    gaia_xyz_bary: np.ndarray,
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
    # Legacy aliases
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
    include_background:
        If True, include massive background asteroids from *background_elements*
        with masses from ``_MAJOR_ASTEROIDS``.  Supersedes ``include_big4``.
    background_elements:
        Dict keyed by lowercase name with MPCORB element dicts.
    include_big4:
        Deprecated alias for ``include_background``.
    big4_elements:
        Deprecated alias for ``background_elements``.
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
        raise ValueError(f"gaia_xyz_bary shape {gaia_xyz_bary.shape} != ({len(obs_jd_tdb)}, 3)")

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
        include_background=include_background or include_big4,
        background_elements=background_elements or big4_elements,
        integrator=integrator,
        dt_days=dt_days,
    )

    # Vectorised: convert entire grid in one call (ecliptic→equatorial + Sun offset).
    bary_icrs = heliocentric_to_barycentric_icrs(helio_ecl, grid)

    # Spline the (T, 3) barycentric trajectory so we can evaluate at any time
    # for the light-time iteration.
    splines = [CubicSpline(grid, bary_icrs[:, k]) for k in range(3)]

    # Vectorised light-time correction: 3 fixed iterations converge to < 1 s
    # precision at 3 AU (same as the scalar light_time_iterate with max_iter=3).
    tau = np.zeros(len(obs_jd_tdb))
    for _ in range(3):
        jd_ret = obs_jd_tdb - tau
        pos_bary = np.stack([s(jd_ret) for s in splines], axis=-1)  # (N, 3)
        tau = np.linalg.norm(pos_bary - gaia_xyz_bary, axis=1) / _C_AU_DAY

    jd_ret = obs_jd_tdb - tau
    pos_bary = np.stack([s(jd_ret) for s in splines], axis=-1)
    return xyz_to_radec(pos_bary - gaia_xyz_bary)


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
