"""Gaia observability estimates for asteroid close encounters.

All position vectors used here are **heliocentric mean ecliptic J2000** in AU —
the same frame as :func:`src.propagate.kepler.kepler_to_cartesian`.  Earth's
position is obtained from Astropy (barycentric ICRS) and converted to that
frame so that downstream subtractions and dot products are physically valid.

Earth's position is used as a proxy for Gaia (at L2, ~0.01 AU offset —
negligible for elongation and magnitude purposes).
"""

from __future__ import annotations

import numpy as np
from astropy.coordinates import get_body_barycentric
from astropy.time import Time

_GAIA_EXCLUSION_DEG: float = 45.0  # minimum solar elongation for Gaia scanning
_GAIA_MAG_LIMIT: float = 21.0  # faint limit of Gaia (G band)

# IAU 2006 mean obliquity of the ecliptic at J2000.0: ε₀ = 84381.406 arcsec.
# Sufficient for the Gaia DR3 epoch range; the secular drift over ±3 years
# is ~1.5 arcsec and irrelevant at the precision needed for elongation cuts.
_OBLIQ_J2000_RAD: float = np.radians(84381.406 / 3600.0)


def get_earth_positions_au(jd_tdb: np.ndarray) -> np.ndarray:
    """Return heliocentric mean ecliptic J2000 positions of Earth.

    The returned frame matches :func:`src.propagate.kepler.kepler_to_cartesian`
    so the two can be combined without a frame transformation.

    Astropy returns Earth in barycentric ICRS (equatorial).  This function:
    1. Subtracts the Sun's barycentric position to obtain heliocentric ICRS.
    2. Rotates by the J2000 obliquity about the X axis to obtain heliocentric
       mean ecliptic J2000.

    Parameters
    ----------
    jd_tdb:
        Array of Julian Dates in TDB scale.

    Returns
    -------
    np.ndarray of shape (N, 3)
        Heliocentric ecliptic J2000 positions in AU.
    """
    t = Time(jd_tdb, format="jd", scale="tdb")
    earth_bary = get_body_barycentric("earth", t)
    sun_bary = get_body_barycentric("sun", t)
    dx = (earth_bary.x - sun_bary.x).to("AU").value
    dy = (earth_bary.y - sun_bary.y).to("AU").value
    dz = (earth_bary.z - sun_bary.z).to("AU").value
    # ICRS (equatorial J2000) → ecliptic J2000: rotation about +X by ε.
    cos_e = np.cos(_OBLIQ_J2000_RAD)
    sin_e = np.sin(_OBLIQ_J2000_RAD)
    x_ec = dx
    y_ec = cos_e * dy + sin_e * dz
    z_ec = -sin_e * dy + cos_e * dz
    return np.stack([x_ec, y_ec, z_ec], axis=-1)


def solar_elongation_deg(
    encounter_xyz: np.ndarray,
    earth_xyz: np.ndarray,
) -> np.ndarray:
    """Solar elongation of the encounter point from Earth (degrees).

    Both inputs must be in the same frame.  In this codebase that frame is
    heliocentric mean ecliptic J2000 (see module docstring), so the Sun sits
    at the origin and the vector Earth → Sun is exactly ``-earth_xyz``.

    Parameters
    ----------
    encounter_xyz:
        Heliocentric ecliptic positions of the encounter, shape (N, 3) in AU.
    earth_xyz:
        Heliocentric ecliptic positions of Earth, shape (N, 3) in AU.

    Returns
    -------
    np.ndarray of shape (N,)
        Solar elongation in degrees (0° = at the Sun, 180° = opposition).
    """
    # Vector Earth → encounter
    enc = encounter_xyz - earth_xyz
    # Vector Earth → Sun (Sun is at the heliocentric origin)
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
    phi1 = np.exp(-3.33 * tan_half**0.63)  # Bowell et al. 1989, Asteroids II, eq. A4
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
