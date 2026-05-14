"""Gaia observability estimates for asteroid close encounters.

Uses Earth's barycentric position as a proxy for Gaia (L2, ~0.01 AU offset —
negligible for elongation and magnitude purposes).
"""

from __future__ import annotations

import numpy as np
from astropy.coordinates import get_body_barycentric
from astropy.time import Time

_GAIA_EXCLUSION_DEG: float = 45.0  # minimum solar elongation for Gaia scanning
_GAIA_MAG_LIMIT: float = 21.0  # faint limit of Gaia (G band)


def get_earth_positions_au(jd_tdb: np.ndarray) -> np.ndarray:
    """Return barycentric positions of Earth for an array of JD TDB times.

    Parameters
    ----------
    jd_tdb:
        Array of Julian Dates in TDB scale.

    Returns
    -------
    np.ndarray of shape (N, 3)
        Barycentric positions in AU (ecliptic J2000 frame approximation;
        astropy returns ICRS, which is close enough for elongation purposes).
    """
    t = Time(jd_tdb, format="jd", scale="tdb")
    pos = get_body_barycentric("earth", t)
    x = pos.x.to("AU").value
    y = pos.y.to("AU").value
    z = pos.z.to("AU").value
    return np.stack([x, y, z], axis=-1)


def solar_elongation_deg(
    encounter_xyz: np.ndarray,
    earth_xyz: np.ndarray,
) -> np.ndarray:
    """Solar elongation of the encounter point from Earth (degrees).

    Parameters
    ----------
    encounter_xyz:
        Barycentric positions of the encounter, shape (N, 3) in AU.
    earth_xyz:
        Barycentric positions of Earth, shape (N, 3) in AU.

    Returns
    -------
    np.ndarray of shape (N,)
        Solar elongation in degrees (0° = at the Sun, 180° = opposition).
    """
    # The Sun is near the barycenter; its position ≈ -earth_xyz (small correction).
    # Vector Earth → encounter
    enc = encounter_xyz - earth_xyz
    # Vector Earth → Sun  (≈ -earth position, since Sun ≈ barycenter)
    sun = -earth_xyz

    enc_norm = np.linalg.norm(enc, axis=-1, keepdims=True)
    sun_norm = np.linalg.norm(sun, axis=-1, keepdims=True)

    cos_elong = np.sum(enc * sun, axis=-1) / (enc_norm[:, 0] * sun_norm[:, 0])
    cos_elong = np.clip(cos_elong, -1.0, 1.0)
    return np.degrees(np.arccos(cos_elong))  # type: ignore[no-any-return]


def apparent_mag_hg(
    h: np.ndarray,
    r_helio: np.ndarray,
    delta: np.ndarray,
    G: np.ndarray | float = 0.15,  # noqa: N803
    phase_deg: np.ndarray | float = 20.0,
) -> np.ndarray:
    """Apparent magnitude using the H-G photometric system.

    Parameters
    ----------
    h:
        Absolute magnitude.
    r_helio:
        Heliocentric distance in AU.
    delta:
        Observer-target distance in AU.
    G:
        Slope parameter (default 0.15).
    phase_deg:
        Phase angle in degrees.  Default 20° when unknown.

    Returns
    -------
    np.ndarray
        Apparent magnitude.
    """
    alpha = np.radians(np.asarray(phase_deg, dtype=float))
    tan_half = np.tan(alpha / 2.0)
    phi1 = np.exp(-3.33 * tan_half**0.63)
    phi2 = np.exp(-1.87 * tan_half**1.22)
    g = np.asarray(G, dtype=float)
    return h - 2.5 * np.log10((1.0 - g) * phi1 + g * phi2) + 5.0 * np.log10(r_helio * delta)  # type: ignore[no-any-return]


def is_gaia_observable(
    elongation: np.ndarray,
    app_mag: np.ndarray,
) -> np.ndarray:
    """Return boolean mask: True if within Gaia's scanning capability.

    Criteria:
    - Solar elongation > 45° (outside exclusion zone).
    - Apparent magnitude < 21 (Gaia faint limit).
    """
    return (elongation > _GAIA_EXCLUSION_DEG) & (app_mag < _GAIA_MAG_LIMIT)
