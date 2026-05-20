"""Astrometric coordinate transformations.

Conventions
-----------
- Positions are 3-vectors (..., 3) in AU.
- Velocities are 3-vectors (..., 3) in AU/day.
- "ICRS" = equatorial J2000 (the Gaia frame).
- "ecliptic" = ecliptic J2000 (the heliocentric Kepler frame).
- Times are JD TDB unless noted.
- Angles internal: radians.

All transforms use the IAU 2006 obliquity ε = 23.4392911° (consistent with
the value used in nbody.py to rotate planet ephemerides into ecliptic frame).
Using a mismatched obliquity introduces a ~37 mas systematic in the ecliptic
→ equatorial rotation.
"""

from __future__ import annotations

import numpy as np
from astropy.coordinates import get_body_barycentric
from astropy.time import Time

# IAU 2006 obliquity at J2000 — must match nbody.py's _EPS_J2000
_OBLIQUITY_DEG = 23.4392911
_EPS = np.radians(_OBLIQUITY_DEG)
_COS_EPS = np.cos(_EPS)
_SIN_EPS = np.sin(_EPS)

# Speed of light in AU/day
_C_AU_DAY = 173.144632674  # = 2.99792458e8 m/s × 86400 / 1.495978707e11

# Conversion factors
_RAD_TO_MAS = 180.0 / np.pi * 3_600_000.0
_DEG_TO_RAD = np.pi / 180.0


# ---------------------------------------------------------------------------
# Frame rotations
# ---------------------------------------------------------------------------


def ecliptic_to_equatorial(xyz_ecl: np.ndarray) -> np.ndarray:
    """Rotate (..., 3) vectors from ecliptic J2000 to equatorial J2000 (ICRS).

    Rotation about the X axis by +ε.
    """
    xyz_ecl = np.asarray(xyz_ecl, dtype=float)
    x = xyz_ecl[..., 0]
    y = _COS_EPS * xyz_ecl[..., 1] - _SIN_EPS * xyz_ecl[..., 2]
    z = _SIN_EPS * xyz_ecl[..., 1] + _COS_EPS * xyz_ecl[..., 2]
    return np.stack([x, y, z], axis=-1)


def equatorial_to_ecliptic(xyz_eq: np.ndarray) -> np.ndarray:
    """Rotate (..., 3) vectors from equatorial J2000 (ICRS) to ecliptic J2000.

    Rotation about the X axis by -ε.
    """
    xyz_eq = np.asarray(xyz_eq, dtype=float)
    x = xyz_eq[..., 0]
    y = _COS_EPS * xyz_eq[..., 1] + _SIN_EPS * xyz_eq[..., 2]
    z = -_SIN_EPS * xyz_eq[..., 1] + _COS_EPS * xyz_eq[..., 2]
    return np.stack([x, y, z], axis=-1)


# ---------------------------------------------------------------------------
# Heliocentric ↔ barycentric
# ---------------------------------------------------------------------------


def sun_barycentric_au(jd_tdb: np.ndarray | float) -> np.ndarray:
    """Sun's barycentric (ICRS) position at *jd_tdb* in AU.

    Uses astropy's default solar-system ephemeris.
    """
    jd_tdb_arr = np.atleast_1d(np.asarray(jd_tdb, dtype=float))
    out = np.empty((len(jd_tdb_arr), 3), dtype=float)
    for i, t in enumerate(jd_tdb_arr):
        sun = get_body_barycentric("sun", Time(t, format="jd", scale="tdb"))
        out[i, 0] = sun.x.to_value("AU")
        out[i, 1] = sun.y.to_value("AU")
        out[i, 2] = sun.z.to_value("AU")
    if np.ndim(jd_tdb) == 0:
        return out[0]
    return out


def heliocentric_to_barycentric_icrs(
    pos_helio_ecl: np.ndarray, jd_tdb: np.ndarray | float
) -> np.ndarray:
    """Convert heliocentric ecliptic position(s) to barycentric ICRS.

    Steps: rotate ecliptic → equatorial; add Sun's barycentric ICRS position.

    Parameters
    ----------
    pos_helio_ecl:
        (..., 3) heliocentric ecliptic position in AU.
    jd_tdb:
        Matching epoch(s) in JD TDB. Either scalar (broadcast over all rows)
        or matching shape with pos_helio_ecl's leading dimension.
    """
    pos_eq = ecliptic_to_equatorial(pos_helio_ecl)
    sun_bary = sun_barycentric_au(jd_tdb)
    # Broadcast sun_bary against pos_eq if needed
    return pos_eq + sun_bary


# ---------------------------------------------------------------------------
# Cartesian ↔ spherical
# ---------------------------------------------------------------------------


def xyz_to_radec(vec: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Convert Cartesian (..., 3) to (RA_deg, Dec_deg)."""
    vec = np.asarray(vec, dtype=float)
    x = vec[..., 0]
    y = vec[..., 1]
    z = vec[..., 2]
    rho = np.sqrt(x * x + y * y)
    ra = np.degrees(np.arctan2(y, x)) % 360.0
    dec = np.degrees(np.arctan2(z, rho))
    return ra, dec


def radec_to_unit_vec(ra_deg: np.ndarray, dec_deg: np.ndarray) -> np.ndarray:
    """Convert (RA_deg, Dec_deg) to (..., 3) unit vectors."""
    ra = np.asarray(ra_deg, dtype=float) * _DEG_TO_RAD
    dec = np.asarray(dec_deg, dtype=float) * _DEG_TO_RAD
    cos_dec = np.cos(dec)
    return np.stack(
        [cos_dec * np.cos(ra), cos_dec * np.sin(ra), np.sin(dec)], axis=-1
    )


# ---------------------------------------------------------------------------
# Light-time correction
# ---------------------------------------------------------------------------


def light_time_iterate(
    target_pos_func,
    jd_tdb_obs: float,
    gaia_xyz_bary: np.ndarray,
    max_iter: int = 3,
    tol_seconds: float = 1.0,
) -> tuple[np.ndarray, float]:
    """Apply iterative light-time correction.

    For a Gaia observation at JD t_obs, the asteroid was actually at its
    position at the retarded time t_ret = t_obs − |r_target(t_ret) − r_gaia| / c.

    Parameters
    ----------
    target_pos_func:
        Callable ``f(jd_tdb) → (3,) barycentric ICRS position in AU``.
    jd_tdb_obs:
        Observation epoch in JD TDB.
    gaia_xyz_bary:
        (3,) Gaia barycentric ICRS position at *jd_tdb_obs* in AU.
    max_iter:
        Maximum number of iterations (typically converges in 2).
    tol_seconds:
        Convergence tolerance (seconds of light-travel time).

    Returns
    -------
    (pos_target_at_retarded, tau_days)
        ``pos_target_at_retarded``: target position (3,) AU at the retarded
        time; ``tau_days``: light-travel time in days.
    """
    tau = 0.0  # initial guess: no light-time correction
    for _ in range(max_iter):
        jd_ret = jd_tdb_obs - tau
        r_t = target_pos_func(jd_ret)
        d = float(np.linalg.norm(r_t - gaia_xyz_bary))
        new_tau = d / _C_AU_DAY
        if abs(new_tau - tau) * 86400.0 < tol_seconds:
            tau = new_tau
            break
        tau = new_tau
    jd_ret = jd_tdb_obs - tau
    return target_pos_func(jd_ret), tau


# ---------------------------------------------------------------------------
# Stellar aberration
# ---------------------------------------------------------------------------


def stellar_aberration(
    los_vec: np.ndarray, observer_vel_au_day: np.ndarray
) -> np.ndarray:
    """Apply first-order stellar aberration to a line-of-sight unit vector.

    Given the true direction-of-arrival ``los_vec`` (unit vector from observer
    to source) and the observer's velocity ``observer_vel_au_day``, returns
    the apparent direction modified by aberration.

    Formula (first-order in v/c):
        n_app ≈ (n + β) / |n + β|     where β = v/c

    For Gaia at ~30 km/s the effect peaks at ~20 arcsec annually.
    """
    los_vec = np.asarray(los_vec, dtype=float)
    v = np.asarray(observer_vel_au_day, dtype=float)
    beta = v / _C_AU_DAY
    apparent = los_vec + beta
    return apparent / np.linalg.norm(apparent, axis=-1, keepdims=True)
