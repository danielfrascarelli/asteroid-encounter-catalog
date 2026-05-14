"""Physical properties of asteroids derived from H magnitude and orbital elements."""

from __future__ import annotations

import numpy as np


def diameter_km(
    h: float | np.ndarray,
    albedo: float | np.ndarray = 0.14,
) -> np.ndarray:
    """Estimate diameter in km from absolute magnitude H and geometric albedo.

    Uses the standard relation: D = (1329 / sqrt(p)) * 10^(-H/5)

    Parameters
    ----------
    h:
        Absolute (V-band) magnitude.
    albedo:
        Geometric albedo.  Default 0.14 (average C-type).

    Returns
    -------
    np.ndarray
        Diameter in km.
    """
    return (1329.0 / np.sqrt(albedo)) * 10.0 ** (-np.asarray(h, dtype=float) / 5.0)


def classify_orbit(
    a: float | np.ndarray,
    e: float | np.ndarray,
) -> np.ndarray:
    """Classify orbit type from semi-major axis (AU) and eccentricity.

    Returns a string array with one of:
    ``"NEA"``, ``"MBA"``, ``"Trojan"``, ``"Centaur"``, ``"TNO"``, ``"Other"``.

    Parameters
    ----------
    a:
        Semi-major axis in AU.
    e:
        Eccentricity.
    """
    a = np.asarray(a, dtype=float)
    e = np.asarray(e, dtype=float)
    q = a * (1.0 - e)  # perihelion distance

    result = np.full(a.shape or (1,), "Other", dtype=object)
    result[a > 30.0] = "TNO"
    result[(a > 5.5) & (a <= 30.0)] = "Centaur"
    result[np.abs(a - 5.205) < 0.5] = "Trojan"
    result[(a >= 1.7) & (a < 5.5) & (np.abs(a - 5.205) >= 0.5)] = "MBA"
    result[q < 1.3] = "NEA"  # NEA overrides MBA if perihelion inside Mars

    if result.shape == (1,) and np.ndim(a) == 0:
        return result[0]  # type: ignore[no-any-return]
    return result
