"""Constantes físicas — única fuente de verdad para el motor orbdet.

Todas en el sistema de unidades interno (AU, día, radián) salvo las auxiliares
de conversión, que llevan el sufijo de unidad en el nombre.
"""

from __future__ import annotations

import math

# --- Tiempo -----------------------------------------------------------------
DAY_S: float = 86_400.0  # segundos por día

# --- Longitud ---------------------------------------------------------------
AU_KM: float = 149_597_870.7  # 1 AU en km (IAU 2012)
AU_M: float = AU_KM * 1_000.0

# --- Gravitación ------------------------------------------------------------
# Constante gravitacional de Gauss (AU^(3/2) / día), valor clásico exacto.
GAUSS_K: float = 0.017_202_098_95
# GM del Sol en AU^3/día^2 = k^2. Es la referencia canónica para el problema de
# dos cuerpos heliocéntrico (la masa del asteroide es despreciable frente al Sol).
GM_SUN: float = GAUSS_K * GAUSS_K  # ≈ 2.959122e-4 AU^3/día^2

# GM del Sol en SI (IAU 2015 nominal), para conversiones masa↔GM.
GM_SUN_SI: float = 1.327_124_400_18e20  # m^3 / s^2
# Masa del Sol en kg (M_sun = GM_sun_SI / G).
G_SI: float = 6.674_30e-11  # m^3 kg^-1 s^-2 (CODATA 2018)
M_SUN_KG: float = GM_SUN_SI / G_SI  # ≈ 1.98841e30 kg

# --- Luz --------------------------------------------------------------------
C_KM_S: float = 299_792.458  # velocidad de la luz, km/s (exacta por definición)
C_AU_PER_DAY: float = C_KM_S * DAY_S / AU_KM  # ≈ 173.1446 AU/día

# --- Orientación ------------------------------------------------------------
# Oblicuidad media de la eclíptica en J2000.0 (IAU 2006): 84381.406 arcsec.
OBLIQUITY_J2000_RAD: float = math.radians(84_381.406 / 3_600.0)  # ≈ 0.409092600600


def gm_from_mass_kg(mass_kg: float) -> float:
    """Convierte una masa en kg a GM en AU^3/día^2.

    GM[AU^3/día^2] = G·M[SI] · (día_s^2) / (AU_m^3).
    """
    gm_si = G_SI * mass_kg  # m^3/s^2
    return gm_si * (DAY_S * DAY_S) / (AU_M**3)


def mass_kg_from_gm(gm_au3_day2: float) -> float:
    """Inversa de :func:`gm_from_mass_kg`: GM (AU^3/día^2) → masa en kg."""
    gm_si = gm_au3_day2 * (AU_M**3) / (DAY_S * DAY_S)
    return gm_si / G_SI
