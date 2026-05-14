"""Typed schema definition for the close-encounter catalog."""

from __future__ import annotations

import polars as pl

CATALOG_SCHEMA: dict[str, type] = {
    "run_id": pl.Utf8,
    "number_1": pl.Int32,
    "number_2": pl.Int32,
    "designation_1": pl.Utf8,
    "designation_2": pl.Utf8,
    "jd_tdb": pl.Float64,
    "date_utc": pl.Utf8,
    "dist_au": pl.Float64,
    "dist_km": pl.Float64,
    "rel_vel_au_day": pl.Float64,
    "rel_vel_km_s": pl.Float64,
    "rel_vel_m_s": pl.Float64,
    "H_1": pl.Float64,
    "H_2": pl.Float64,
    "diameter_1_km": pl.Float64,
    "diameter_2_km": pl.Float64,
    "class_1": pl.Utf8,
    "class_2": pl.Utf8,
    "solar_elongation_deg": pl.Float64,
    "gaia_observable": pl.Boolean,
}

CATALOG_COLUMNS: list[str] = list(CATALOG_SCHEMA.keys())

COLUMN_DESCRIPTIONS: dict[str, str] = {
    "run_id": "Unique identifier for the pipeline run (ISO timestamp)",
    "number_1": "MPC number of the first asteroid (smaller number in pair)",
    "number_2": "MPC number of the second asteroid (larger number in pair)",
    "designation_1": "MPC packed designation of asteroid 1",
    "designation_2": "MPC packed designation of asteroid 2",
    "jd_tdb": "Julian Date of closest approach (TDB scale)",
    "date_utc": "Calendar date of closest approach (UTC, yyyy-mm-dd)",
    "dist_au": "Minimum 3D separation at closest approach (AU)",
    "dist_km": "Minimum 3D separation at closest approach (km)",
    "rel_vel_au_day": "Relative velocity at encounter (AU/day)",
    "rel_vel_km_s": "Relative velocity at encounter (km/s)",
    "rel_vel_m_s": "Relative velocity at encounter (m/s)",
    "H_1": "Absolute magnitude of asteroid 1 (V-band, from MPCORB)",
    "H_2": "Absolute magnitude of asteroid 2 (V-band, from MPCORB)",
    "diameter_1_km": "Estimated diameter of asteroid 1 (km) from H and default albedo",
    "diameter_2_km": "Estimated diameter of asteroid 2 (km) from H and default albedo",
    "class_1": "Dynamical class of asteroid 1: NEA/MBA/Trojan/Centaur/TNO/Other",
    "class_2": "Dynamical class of asteroid 2: NEA/MBA/Trojan/Centaur/TNO/Other",
    "solar_elongation_deg": "Solar elongation of encounter midpoint from Earth (degrees)",
    "gaia_observable": "True if elongation > 45° and estimated apparent mag < 21",
}
