"""Keplerian two-body propagator (heliocentric ecliptic J2000).

Units
-----
- Distances : AU
- Times     : Julian Date TDB
- Angles    : radians (convert degrees before calling)

The Gaussian gravitational constant k = 0.01720209895 AU^(3/2) day^-1 gives
the correct mean motion n = k / a^(3/2) rad day^-1 for heliocentric orbits.
"""

from __future__ import annotations

import numpy as np
import polars as pl

_K = 0.01720209895  # Gaussian gravitational constant [AU^(3/2) day^-1]
_DEG = np.pi / 180.0


def solve_kepler(
    M: float | np.ndarray,  # noqa: N803
    e: float | np.ndarray,
    tol: float = 1e-12,
    max_iter: int = 50,
) -> np.ndarray:
    """Solve Kepler's equation M = E - e·sin(E) by Newton-Raphson.

    Parameters
    ----------
    M:
        Mean anomaly in radians.
    e:
        Eccentricity (0 ≤ e < 1).
    tol:
        Convergence tolerance on |ΔE|.
    max_iter:
        Maximum iterations.

    Returns
    -------
    np.ndarray
        Eccentric anomaly E in radians, broadcast shape of M and e.
    """
    M = np.atleast_1d(np.asarray(M, dtype=float))  # noqa: N806
    e = np.atleast_1d(np.asarray(e, dtype=float))
    M, e = np.broadcast_arrays(M, e)  # noqa: N806
    M = M.copy()  # noqa: N806
    e = e.copy()
    E = M.copy()  # noqa: N806
    for _ in range(max_iter):
        dE = (M - E + e * np.sin(E)) / (1.0 - e * np.cos(E))  # noqa: N806
        E += dE  # noqa: N806
        if np.all(np.abs(dE) < tol):
            break
    return E  # noqa: N806


def kepler_to_cartesian(
    a_au: float | np.ndarray,
    e: float | np.ndarray,
    i_rad: float | np.ndarray,
    Omega_rad: float | np.ndarray,  # noqa: N803
    omega_rad: float | np.ndarray,
    M0_rad: float | np.ndarray,  # noqa: N803
    epoch_jd: float | np.ndarray,
    t_jd: float | np.ndarray,
) -> np.ndarray:
    """Propagate Keplerian elements to heliocentric ecliptic J2000 Cartesian coordinates.

    All array inputs are broadcast together; scalar inputs are supported.

    Parameters
    ----------
    a_au:
        Semi-major axis in AU.
    e:
        Eccentricity.
    i_rad:
        Inclination in radians.
    Omega_rad:
        Longitude of ascending node in radians.
    omega_rad:
        Argument of perihelion in radians.
    M0_rad:
        Mean anomaly at epoch in radians.
    epoch_jd:
        Reference epoch in JD TDB.
    t_jd:
        Target time in JD TDB.

    Returns
    -------
    np.ndarray
        Heliocentric ecliptic J2000 positions in AU, shape (..., 3).
    """
    a = np.asarray(a_au, dtype=float)
    e_ = np.asarray(e, dtype=float)
    i = np.asarray(i_rad, dtype=float)
    cap_o = np.asarray(Omega_rad, dtype=float)
    o = np.asarray(omega_rad, dtype=float)
    m0 = np.asarray(M0_rad, dtype=float)
    t0 = np.asarray(epoch_jd, dtype=float)
    t = np.asarray(t_jd, dtype=float)

    n = _K / a**1.5  # mean motion [rad/day]
    m = (m0 + n * (t - t0)) % (2.0 * np.pi)  # noqa: N806
    ecc_anom = solve_kepler(m, e_)

    nu = 2.0 * np.arctan2(
        np.sqrt(1.0 + e_) * np.sin(ecc_anom / 2.0),
        np.sqrt(1.0 - e_) * np.cos(ecc_anom / 2.0),
    )
    r = a * (1.0 - e_ * np.cos(ecc_anom))
    u = o + nu  # argument of latitude

    cos_cap_o = np.cos(cap_o)
    sin_cap_o = np.sin(cap_o)
    cos_i = np.cos(i)
    sin_i = np.sin(i)
    cos_u = np.cos(u)
    sin_u = np.sin(u)

    x = r * (cos_cap_o * cos_u - sin_cap_o * sin_u * cos_i)
    y = r * (sin_cap_o * cos_u + cos_cap_o * sin_u * cos_i)
    z = r * sin_u * sin_i

    return np.stack([x, y, z], axis=-1)


def propagate_df(elements: pl.DataFrame, t_jd: float) -> np.ndarray:
    """Propagate all asteroids in *elements* to *t_jd*.

    Parameters
    ----------
    elements:
        DataFrame with columns: a_au, e, i_deg, Omega_deg, omega_deg, M_deg, epoch_jd.
    t_jd:
        Target time in JD TDB.

    Returns
    -------
    np.ndarray of shape (N, 3)
        Heliocentric ecliptic J2000 positions in AU.
    """
    return kepler_to_cartesian(
        a_au=elements["a_au"].to_numpy(),
        e=elements["e"].to_numpy(),
        i_rad=elements["i_deg"].to_numpy() * _DEG,
        Omega_rad=elements["Omega_deg"].to_numpy() * _DEG,
        omega_rad=elements["omega_deg"].to_numpy() * _DEG,
        M0_rad=elements["M_deg"].to_numpy() * _DEG,
        epoch_jd=elements["epoch_jd"].to_numpy(),
        t_jd=t_jd,
    )
