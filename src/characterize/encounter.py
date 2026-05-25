"""Encounter characterization — enriches the raw detection catalog.

Takes the output of the detection pipeline and adds physical properties,
unit conversions, orbit classification, and Gaia observability estimates.

Public entry point: :func:`characterize_catalog`.
"""

from __future__ import annotations

import logging

import numpy as np
import polars as pl
from astropy.time import Time

from src.characterize.geometry import (
    dist_au_to_km,
    vel_au_per_day_to_km_s,
    vel_au_per_day_to_m_s,
)
from src.characterize.observability import (
    apparent_mag_hg,
    get_earth_positions_au,
    is_gaia_observable,
    solar_elongation_deg,
)
from src.characterize.physical import classify_orbit, diameter_km
from src.propagate.kepler import kepler_to_cartesian

logger = logging.getLogger(__name__)

_DEG = np.pi / 180.0

# Columns from elements DataFrame needed for position propagation
_ELEM_COLS = ["number", "a_au", "e", "i_deg", "Omega_deg", "omega_deg", "M_deg", "epoch_jd"]


def _join_elements(
    encounters: pl.DataFrame,
    elements: pl.DataFrame,
    suffix: str,
) -> pl.DataFrame:
    """Left-join orbital elements onto encounters for one body (1 or 2).

    Parameters
    ----------
    encounters:
        Encounters DataFrame with ``number_{suffix}`` column.
    elements:
        Orbital elements DataFrame (must include ``number`` and ``_ELEM_COLS``).
    suffix:
        ``"1"`` or ``"2"``.

    Returns
    -------
    pl.DataFrame
        encounters with added columns ``a_au_{suffix}``, ``e_{suffix}``, etc.
    """
    elem_cols = [c for c in _ELEM_COLS if c != "number"]
    renamed = elements.select(["number"] + elem_cols).rename(
        {c: f"{c}_{suffix}" for c in elem_cols}
    )
    return encounters.join(renamed, left_on=f"number_{suffix}", right_on="number", how="left")


def _positions_from_elements(df: pl.DataFrame, suffix: str) -> np.ndarray:
    """Vectorised Kepler propagation for body *suffix* at each encounter JD."""
    return kepler_to_cartesian(
        a_au=df[f"a_au_{suffix}"].to_numpy(),
        e=df[f"e_{suffix}"].to_numpy(),
        i_rad=df[f"i_deg_{suffix}"].to_numpy() * _DEG,
        Omega_rad=df[f"Omega_deg_{suffix}"].to_numpy() * _DEG,
        omega_rad=df[f"omega_deg_{suffix}"].to_numpy() * _DEG,
        M0_rad=df[f"M_deg_{suffix}"].to_numpy() * _DEG,
        epoch_jd=df[f"epoch_jd_{suffix}"].to_numpy(),
        t_jd=df["jd_tdb"].to_numpy(),
    )


def characterize_catalog(
    encounters: pl.DataFrame,
    elements: pl.DataFrame,
    mpcorb: pl.DataFrame,
    albedo: float = 0.14,
) -> pl.DataFrame:
    """Enrich the detection catalog with physical and observational properties.

    Parameters
    ----------
    encounters:
        Output of :func:`src.detect.pipeline.detect_encounters`.  Required
        columns: ``number_1``, ``number_2``, ``designation_1``,
        ``designation_2``, ``jd_tdb``, ``dist_au``, ``rel_vel_au_day``.
    elements:
        Orbital elements (gaia_orbits + Horizons supplement).  Must include
        all columns in ``_ELEM_COLS``.
    mpcorb:
        Parsed MPCORB DataFrame with columns ``number``, ``H``, ``G``.
    albedo:
        Default geometric albedo for diameter estimation when an
        asteroid-specific value is unavailable.

    Returns
    -------
    pl.DataFrame
        Enriched catalog with columns:

        ``number_1``, ``number_2``, ``designation_1``, ``designation_2``,
        ``jd_tdb``, ``date_utc``,
        ``dist_au``, ``dist_km``,
        ``rel_vel_au_day``, ``rel_vel_km_s``, ``rel_vel_m_s``,
        ``H_1``, ``H_2``,
        ``diameter_1_km``, ``diameter_2_km``,
        ``class_1``, ``class_2``,
        ``solar_elongation_deg``, ``gaia_observable``.
    """
    n = len(encounters)
    logger.info("Characterizing %d encounters…", n)

    # ------------------------------------------------------------------ #
    # 1. Join orbital elements for both bodies                            #
    # ------------------------------------------------------------------ #
    df = _join_elements(encounters, elements, "1")
    df = _join_elements(df, elements, "2")

    # ------------------------------------------------------------------ #
    # 2. Join H and G from MPCORB                                         #
    # ------------------------------------------------------------------ #
    h_cols = mpcorb.select(["number", "H", "G"])
    df = df.join(
        h_cols.rename({"H": "H_1", "G": "G_1"}), left_on="number_1", right_on="number", how="left"
    )
    df = df.join(
        h_cols.rename({"H": "H_2", "G": "G_2"}), left_on="number_2", right_on="number", how="left"
    )

    # ------------------------------------------------------------------ #
    # 3. Unit conversions                                                  #
    # ------------------------------------------------------------------ #
    dist_km = dist_au_to_km(df["dist_au"].to_numpy())
    vel_km_s = vel_au_per_day_to_km_s(df["rel_vel_au_day"].to_numpy())
    vel_m_s = vel_au_per_day_to_m_s(df["rel_vel_au_day"].to_numpy())

    date_utc = [Time(jd, format="jd", scale="tdb").utc.iso[:10] for jd in df["jd_tdb"].to_list()]

    # ------------------------------------------------------------------ #
    # 4. Diameters                                                         #
    # ------------------------------------------------------------------ #
    h1 = df["H_1"].to_numpy()
    h2 = df["H_2"].to_numpy()
    diam_1 = diameter_km(h1, albedo)
    diam_2 = diameter_km(h2, albedo)

    # ------------------------------------------------------------------ #
    # 5. Orbit classification                                              #
    # ------------------------------------------------------------------ #
    a1 = df["a_au_1"].to_numpy()
    e1 = df["e_1"].to_numpy()
    a2 = df["a_au_2"].to_numpy()
    e2 = df["e_2"].to_numpy()
    cls1 = classify_orbit(a1, e1).astype(str)
    cls2 = classify_orbit(a2, e2).astype(str)

    # ------------------------------------------------------------------ #
    # 6. Heliocentric positions at encounter epoch (vectorised Kepler)    #
    # ------------------------------------------------------------------ #
    logger.info("Propagating asteroid positions for observability…")
    pos1 = _positions_from_elements(df, "1")  # (N, 3)
    pos2 = _positions_from_elements(df, "2")  # (N, 3)
    enc_pos = 0.5 * (pos1 + pos2)  # midpoint as encounter position

    r_helio_1 = np.linalg.norm(pos1, axis=1)  # heliocentric distance body 1

    # ------------------------------------------------------------------ #
    # 7. Earth positions and solar elongation                             #
    # ------------------------------------------------------------------ #
    # get_earth_positions_au returns heliocentric ecliptic J2000 to match
    # kepler_to_cartesian (which produced pos1, pos2 above) — both vectors
    # live in the same frame, so subtractions and dot products are valid.
    logger.info("Fetching Earth positions from astropy ephemeris…")
    jd_arr = df["jd_tdb"].to_numpy()
    earth_pos = get_earth_positions_au(jd_arr)  # (N, 3)

    elong = solar_elongation_deg(enc_pos, earth_pos)

    # ------------------------------------------------------------------ #
    # 8. Apparent magnitude and Gaia observability                        #
    # ------------------------------------------------------------------ #
    delta = np.linalg.norm(enc_pos - earth_pos, axis=1)  # Earth–encounter
    g1 = np.where(np.isnan(df["G_1"].to_numpy()), 0.15, df["G_1"].to_numpy())
    app_mag = np.where(
        np.isnan(h1),
        np.nan,
        apparent_mag_hg(np.where(np.isnan(h1), 20.0, h1), r_helio_1, delta, G=g1),
    )
    observable = is_gaia_observable(elong, app_mag)

    # ------------------------------------------------------------------ #
    # 9. Assemble output DataFrame                                         #
    # ------------------------------------------------------------------ #
    logger.info("Assembling enriched catalog…")
    return pl.DataFrame(
        {
            "number_1": df["number_1"],
            "number_2": df["number_2"],
            "designation_1": df["designation_1"],
            "designation_2": df["designation_2"],
            "jd_tdb": df["jd_tdb"],
            "date_utc": date_utc,
            "dist_au": df["dist_au"],
            "dist_km": dist_km,
            "rel_vel_au_day": df["rel_vel_au_day"],
            "rel_vel_km_s": vel_km_s,
            "rel_vel_m_s": vel_m_s,
            "H_1": h1,
            "H_2": h2,
            "diameter_1_km": diam_1,
            "diameter_2_km": diam_2,
            "class_1": cls1,
            "class_2": cls2,
            "solar_elongation_deg": elong,
            "gaia_observable": observable,
        }
    ).sort("dist_au")
