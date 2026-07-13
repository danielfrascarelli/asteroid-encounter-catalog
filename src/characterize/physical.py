"""Physical properties of asteroids derived from H magnitude and orbital elements.

Diameter priority chain (tribunal 2026-07-04, B3 — the old code applied a fixed
albedo 0.14 to every body, giving Ceres 763 km instead of 939):

1. **Measured diameter** (IRAS/AKARI/NEOWISE/occultations via JPL SBDB).
2. ``D(H, measured albedo)`` when only the albedo is measured.
3. ``D(H, zone albedo)`` — an *orbital-zone* average albedo (inner/mid/outer
   belt).  This is a dynamical proxy, NOT a taxonomic classification.
4. ``D(H, default_albedo)`` when the semi-major axis is unknown.
"""

from __future__ import annotations

import numpy as np

# Promedios de albedo geométrico por zona orbital (proxy dinámico, no taxonomía).
# Valores redondeados de las distribuciones NEOWISE por zona del cinturón
# (Masiero et al.): interior dominada por tipo S, exterior por tipo C.
ALBEDO_BY_ZONE: dict[str, float] = {
    "inner": 0.20,  # a < 2.5 AU
    "mid": 0.13,  # 2.5 ≤ a < 2.82 AU
    "outer": 0.06,  # a ≥ 2.82 AU
}

# Densidades bulk típicas por zona orbital (kg/m³): interior dominada por tipo S,
# exterior por tipo C (Carry 2012). Mismo proxy dinámico que ALBEDO_BY_ZONE.
DENSITY_BY_ZONE_KG_M3: dict[str, float] = {
    "inner": 2700.0,
    "mid": 2000.0,
    "outer": 1300.0,
}

_G_SI = 6.674e-11  # m³ kg⁻¹ s⁻²


def deflection_dv_m_s(
    diameter_km: np.ndarray,
    a_au: np.ndarray,
    dist_m: np.ndarray,
    v_rel_m_s: np.ndarray,
    density_by_zone: dict[str, float] | None = None,
) -> np.ndarray:
    """Kick de velocidad impartido por el perturbador sobre el otro cuerpo (m/s).

    Aproximación impulsiva de dos cuerpos: ``Δv ≈ 2·GM_pert/(b·v_rel)``, con
    ``GM = G·(π/6)·ρ·D³`` y ρ por zona orbital (proxy, no medición). Es la métrica
    de *ranking* de utilidad por par que piden Ivantsov 2018 / FM 2025 (tribunal
    M7): la señal astrométrica observable es ∝ Δv integrado sobre el arco, así
    que ordenar por Δv ordena por señal.

    Parameters
    ----------
    diameter_km:
        Diámetro del perturbador (cuerpo grande del par).
    a_au:
        Semieje del perturbador — selecciona la densidad de zona.
    dist_m:
        Distancia mínima del encuentro (m) — el parámetro de impacto ``b``.
    v_rel_m_s:
        Velocidad relativa en el encuentro (m/s).

    Returns
    -------
    np.ndarray
        Δv en m/s (NaN donde falte diámetro, distancia o velocidad).
    """
    zones = density_by_zone if density_by_zone is not None else DENSITY_BY_ZONE_KG_M3
    diameter_km = np.asarray(diameter_km, dtype=float)
    a_au = np.asarray(a_au, dtype=float)
    dist_m = np.asarray(dist_m, dtype=float)
    v_rel_m_s = np.asarray(v_rel_m_s, dtype=float)

    rho = np.full(diameter_km.shape, zones["mid"], dtype=float)
    rho[a_au < 2.5] = zones["inner"]
    rho[a_au >= 2.82] = zones["outer"]

    d_m = diameter_km * 1e3
    gm = _G_SI * (np.pi / 6.0) * rho * d_m**3  # m³/s²
    with np.errstate(divide="ignore", invalid="ignore"):
        dv = 2.0 * gm / (dist_m * v_rel_m_s)
    dv = np.where((dist_m > 0) & (v_rel_m_s > 0), dv, np.nan)
    return dv


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


def diameter_km_with_source(
    h: np.ndarray,
    a_au: np.ndarray,
    diameter_measured_km: np.ndarray | None = None,
    albedo_measured: np.ndarray | None = None,
    default_albedo: float = 0.14,
    albedo_by_zone: dict[str, float] | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Diameter with per-body provenance, applying the priority chain of B3.

    Parameters
    ----------
    h:
        Absolute magnitudes (NaN where unknown).
    a_au:
        Semi-major axes (NaN where unknown) — selects the zone albedo.
    diameter_measured_km:
        Measured diameters (NaN where unmeasured), e.g. from JPL SBDB.
    albedo_measured:
        Measured geometric albedos (NaN where unmeasured).
    default_albedo:
        Fallback albedo when neither measurement nor zone applies.
    albedo_by_zone:
        Zone-average albedos; defaults to :data:`ALBEDO_BY_ZONE`.

    Returns
    -------
    (diameter_km, source)
        ``diameter_km`` is Float64 (NaN when H and diameter are both unknown);
        ``source`` is an object array with one of ``"measured"``,
        ``"albedo_measured"``, ``"zone_albedo"``, ``"default_albedo"``,
        ``"unknown"``.
    """
    zones = albedo_by_zone if albedo_by_zone is not None else ALBEDO_BY_ZONE
    h = np.asarray(h, dtype=float)
    a_au = np.asarray(a_au, dtype=float)
    n = h.shape[0]
    d_meas = (
        np.asarray(diameter_measured_km, dtype=float)
        if diameter_measured_km is not None
        else np.full(n, np.nan)
    )
    p_meas = (
        np.asarray(albedo_measured, dtype=float)
        if albedo_measured is not None
        else np.full(n, np.nan)
    )

    # Zone albedo from semi-major axis; default where a is unknown.
    albedo_fallback = np.full(n, default_albedo, dtype=float)
    albedo_fallback[a_au < 2.5] = zones["inner"]
    albedo_fallback[(a_au >= 2.5) & (a_au < 2.82)] = zones["mid"]
    albedo_fallback[a_au >= 2.82] = zones["outer"]
    zone_known = ~np.isnan(a_au)

    diam = np.full(n, np.nan, dtype=float)
    source = np.full(n, "unknown", dtype=object)

    h_known = ~np.isnan(h)
    # 4 → 3 → 2 → 1: later (higher-priority) assignments overwrite earlier ones.
    diam[h_known] = diameter_km(h[h_known], default_albedo)
    source[h_known] = "default_albedo"

    m_zone = h_known & zone_known
    diam[m_zone] = diameter_km(h[m_zone], albedo_fallback[m_zone])
    source[m_zone] = "zone_albedo"

    m_alb = h_known & ~np.isnan(p_meas) & (p_meas > 0.0)
    diam[m_alb] = diameter_km(h[m_alb], p_meas[m_alb])
    source[m_alb] = "albedo_measured"

    m_diam = ~np.isnan(d_meas) & (d_meas > 0.0)
    diam[m_diam] = d_meas[m_diam]
    source[m_diam] = "measured"

    return diam, source


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
