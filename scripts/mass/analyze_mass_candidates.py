"""Rank Category B novel encounters as mass-determination candidates.

Category B = encounters whose perturber has *no* published mass
(``mass_unknown == true`` in ``data/output/relevant_novel_encounters.csv``).
These are the scientifically interesting candidates: each one is a potential
new mass measurement.

For each candidate this script:

1. Estimates the perturber mass from its catalog diameter assuming a typical
   main-belt density ρ = 1.5 g/cm³ (1500 kg/m³):

       M_est = ρ · (4/3) π (D/2)³

2. Computes the expected linear gravitational deflection of the target

       δ = 2 G M_est / (v² b)

   in micro-arcseconds, where ``b`` is the closest-approach distance and
   ``v`` is the relative velocity at the encounter.

3. Flags ``priority_by_impulse_score`` whenever δ ≥ 100 μas — the Gaia
   single-transit astrometric precision for a typical MBA.  This is a
   *priority* signal for follow-up mass-fitting, **not** a confirmation
   that the encounter is dynamically observable.

4. Optionally (``--with-jpl``) cross-checks the geometry against JPL Horizons
   with the same ±2-day / 30-min window used in the Cat A validator and
   writes ``jpl_dist_au`` / ``delta_au`` columns.

5. Checks whether the target was actually observed by Gaia (i.e. whether its
   MPC number appears in ``data/raw/gaia_sso.parquet`` under the ``number_mp``
   column).  This determines whether a mass measurement is realistically
   achievable from the existing DR3 transits.

The output is ``data/output/mass_candidates.csv`` plus a human-readable
ranking table in the logs.

Usage
-----
    docker compose run --rm pipeline python -m scripts.analyze_mass_candidates
    docker compose run --rm pipeline python -m scripts.analyze_mass_candidates \\
        --top-n 30 --with-jpl
"""

from __future__ import annotations

import argparse
import logging
import math
import time
from pathlib import Path

import numpy as np
import polars as pl
from astropy.time import Time

from src.utils.config import load_config

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Physics constants
# ---------------------------------------------------------------------------

_G = 6.674e-11  # N m² / kg²
_AU_M = 1.495978707e11  # meters per AU
_RAD_TO_MUAS = 2.06265e11  # μas per radian

# Average MBA density for mass-from-diameter estimates.
_RHO_KG_M3 = 1500.0  # 1.5 g/cm³

# Gaia single-transit astrometric precision for a typical MBA.
_GAIA_PRECISION_MUAS = 100.0

# Default paths.
_ENCOUNTERS_PATH = Path("data/output/relevant_novel_encounters.csv")
_OUTPUT_PATH = Path("data/output/mass_candidates.csv")
# Release-scoped by default (ingest writes gaia_sso_{release}.parquet). Pass
# --gaia_sso data/raw/gaia_sso_fpr.parquet to use the FPR bulk download.
_GAIA_SSO_PATH = Path("data/raw/gaia_sso_dr3.parquet")


# ---------------------------------------------------------------------------
# Horizons helpers (mirrored from validate_novel_a / validate_jpl_horizons)
# ---------------------------------------------------------------------------


def _horizons_vectors(
    asteroid_id: int,
    jd_center: float,
    half_window_days: float = 2.0,
    step: str = "30m",
) -> tuple[np.ndarray, np.ndarray]:
    """Fetch heliocentric (Sun-centred, @10) state vectors around *jd_center*."""
    from astroquery.jplhorizons import Horizons

    t_center = Time(jd_center, format="jd", scale="tdb")
    epochs = {
        "start": (t_center - half_window_days).utc.iso[:19],
        "stop": (t_center + half_window_days).utc.iso[:19],
        "step": step,
    }
    h = Horizons(id=str(asteroid_id), location="@10", epochs=epochs, id_type="smallbody")
    tbl = h.vectors(refplane="ecliptic")
    jd = np.array(tbl["datetime_jd"], dtype=float)
    xyz = np.column_stack(
        [
            np.array(tbl["x"], dtype=float),
            np.array(tbl["y"], dtype=float),
            np.array(tbl["z"], dtype=float),
        ]
    )
    return jd, xyz


def _jpl_min_distance(
    a: int,
    b: int,
    jd_center: float,
    half_window_days: float = 2.0,
    step: str = "30m",
) -> tuple[float, float]:
    """Return ``(min_distance_AU, jd_TDB_at_min)`` between *a* and *b* near *jd_center*."""
    jd_a, xyz_a = _horizons_vectors(a, jd_center, half_window_days, step)
    jd_b, xyz_b = _horizons_vectors(b, jd_center, half_window_days, step)
    if not np.allclose(jd_a, jd_b):
        common = np.intersect1d(jd_a, jd_b)
        idx_a = np.array([np.where(jd_a == t)[0][0] for t in common])
        idx_b = np.array([np.where(jd_b == t)[0][0] for t in common])
        xyz_a, xyz_b, jd = xyz_a[idx_a], xyz_b[idx_b], common
    else:
        jd = jd_a
    d = np.linalg.norm(xyz_a - xyz_b, axis=1)
    k = int(np.argmin(d))
    return float(d[k]), float(jd[k])


# ---------------------------------------------------------------------------
# Physics
# ---------------------------------------------------------------------------


def estimate_mass_kg(diameter_km: float, rho_kg_m3: float = _RHO_KG_M3) -> float:
    """Sphere-of-density mass estimate from a diameter in km."""
    if not math.isfinite(diameter_km) or diameter_km <= 0.0:
        return float("nan")
    d_m = diameter_km * 1000.0
    radius_m = 0.5 * d_m
    volume_m3 = (4.0 / 3.0) * math.pi * radius_m**3
    return rho_kg_m3 * volume_m3


def compute_deflection_muas(mass_kg: float, rel_vel_km_s: float, dist_au: float) -> float:
    """Linear gravitational deflection in μas; NaN on invalid inputs."""
    if not (math.isfinite(mass_kg) and math.isfinite(rel_vel_km_s) and math.isfinite(dist_au)):
        return float("nan")
    if rel_vel_km_s <= 0.0 or dist_au <= 0.0:
        return float("nan")
    b_m = dist_au * _AU_M
    v_ms = rel_vel_km_s * 1000.0
    delta_rad = 2.0 * _G * mass_kg / (v_ms * v_ms * b_m)
    return float(delta_rad * _RAD_TO_MUAS)


# ---------------------------------------------------------------------------
# Gaia SSO target lookup
# ---------------------------------------------------------------------------


def _load_gaia_target_numbers(path: Path) -> set[int] | None:
    """Return the set of MPC numbers present in the Gaia SSO parquet, or None."""
    if not path.exists():
        logger.warning("Gaia SSO parquet not found at %s; gaia_has_target will be None", path)
        return None
    try:
        df = pl.read_parquet(path, columns=["number_mp"])
    except Exception as exc:  # noqa: BLE001
        logger.warning("Failed to read %s (%s); gaia_has_target will be None", path, exc)
        return None
    numbers = (
        df.filter(pl.col("number_mp").is_not_null())
        .get_column("number_mp")
        .cast(pl.Int64, strict=False)
        .drop_nulls()
        .unique()
        .to_list()
    )
    logger.info("Loaded %d unique MPC numbers from Gaia SSO catalog", len(numbers))
    return set(numbers)


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------


def _date_utc_to_jd_tdb(date_utc: str) -> float:
    return float(Time(date_utc, format="iso", scale="utc").tdb.jd)


def _is_missing_number(value: object) -> bool:
    if value is None:
        return True
    if isinstance(value, float) and math.isnan(value):
        return True
    return False


def analyze_candidates(
    encounters_path: Path,
    output_path: Path,
    top_n: int,
    with_jpl: bool,
    rate_limit_s: float,
    half_window_days: float,
    gaia_sso_path: Path,
) -> pl.DataFrame:
    """Build the ranked mass-candidate catalog and write it to disk."""
    logger.info("Loading encounters from %s", encounters_path)
    df = pl.read_csv(encounters_path)
    logger.info("Total novel encounters: %d", df.height)

    cat_b = df.filter(pl.col("mass_unknown")).sort("deflection_score", descending=True).head(top_n)
    logger.info(
        "Category B (mass_unknown == true): selecting top %d by deflection_score",
        cat_b.height,
    )

    gaia_numbers = _load_gaia_target_numbers(gaia_sso_path)

    rows: list[dict] = []
    for rank, row in enumerate(cat_b.iter_rows(named=True), start=1):
        perturber = int(row["number_1"])
        perturber_name = row["designation_1"]
        perturber_diameter_km = float(row["diameter_1_km"])
        target_no_raw = row["number_2"]
        target_designation = row["designation_2"]
        target_diameter_km = (
            float(row["diameter_2_km"])
            if row["diameter_2_km"] is not None
            and not (isinstance(row["diameter_2_km"], float) and math.isnan(row["diameter_2_km"]))
            else float("nan")
        )
        date_utc = row["date_utc"]
        dist_au = float(row["dist_au"])
        rel_vel_km_s = float(row["rel_vel_km_s"])
        deflection_score = float(row["deflection_score"])

        target_no: int | None
        if _is_missing_number(target_no_raw):
            target_no = None
        else:
            target_no = int(target_no_raw)

        # Mass + deflection
        mass_est_kg = estimate_mass_kg(perturber_diameter_km)
        deflection_muas = compute_deflection_muas(mass_est_kg, rel_vel_km_s, dist_au)
        # NOTE: this flag indicates the impulse-formula deflection exceeds 100 μas.
        # It is a *priority score*, NOT a detection or publishability criterion.
        # A real detection requires an AL-weighted LOO orbit fit with N-body integration.
        priority_by_impulse_score = bool(
            math.isfinite(deflection_muas) and deflection_muas >= _GAIA_PRECISION_MUAS
        )

        # Gaia target lookup
        if gaia_numbers is None:
            gaia_has_target: bool | None = None
        elif target_no is None:
            gaia_has_target = False
        else:
            gaia_has_target = target_no in gaia_numbers

        logger.info(
            "[#%02d] (%d) %s + %s on %s :: d=%.6f AU, v=%.3f km/s, "
            "D₁=%.1f km, M_est=%.3e kg, δ=%.1f μas, priority=%s, gaia=%s",
            rank,
            perturber,
            perturber_name,
            target_designation,
            date_utc,
            dist_au,
            rel_vel_km_s,
            perturber_diameter_km,
            mass_est_kg,
            deflection_muas,
            priority_by_impulse_score,
            gaia_has_target,
        )

        # Special highlight for the closest known encounter in the catalog.
        if (
            perturber == 111
            and isinstance(target_designation, str)
            and "2000" in target_designation
            and "nt3" in target_designation.lower()
        ):
            logger.info(
                "  >>> NOTE: (111) Ate + 2000 NT3 is the closest encounter in the "
                "catalog (≈ 4.72e-4 AU). Prime mass-determination target."
            )

        # Optional JPL cross-check
        jpl_dist_au = float("nan")
        delta_au = float("nan")
        if with_jpl:
            if target_no is None:
                logger.warning("  target has no MPC number; skipping JPL query")
            else:
                try:
                    jd_center = _date_utc_to_jd_tdb(date_utc)
                    jpl_dist_au, _ = _jpl_min_distance(
                        perturber, target_no, jd_center, half_window_days, "30m"
                    )
                    delta_au = dist_au - jpl_dist_au
                    logger.info(
                        "  JPL min dist = %.6f AU (Δ = %+.2e AU)",
                        jpl_dist_au,
                        delta_au,
                    )
                except Exception as exc:  # noqa: BLE001
                    logger.warning("  JPL Horizons query failed: %s", exc)
                finally:
                    time.sleep(rate_limit_s)

        rows.append(
            {
                "rank": rank,
                "perturber_number": perturber,
                "perturber_name": perturber_name,
                "perturber_diameter_km": perturber_diameter_km,
                "target_number": target_no,
                "target_designation": target_designation,
                "target_diameter_km": target_diameter_km,
                "date_utc": date_utc,
                "dist_au": dist_au,
                "rel_vel_km_s": rel_vel_km_s,
                "deflection_score": deflection_score,
                "mass_est_kg": mass_est_kg,
                "deflection_muas": deflection_muas,
                "gaia_precision_muas": _GAIA_PRECISION_MUAS,
                "gaia_has_target": gaia_has_target,
                "priority_by_impulse_score": priority_by_impulse_score,
                "jpl_dist_au": jpl_dist_au,
                "delta_au": delta_au,
            }
        )

    result = pl.DataFrame(rows)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.write_csv(output_path)
    logger.info("Wrote %d rows to %s", result.height, output_path)

    # Summary
    n_prio = int(result["priority_by_impulse_score"].sum())
    logger.info(
        "Above-impulse-score-threshold candidates (δ ≥ %.0f μas): %d / %d",
        _GAIA_PRECISION_MUAS,
        n_prio,
        result.height,
    )
    if gaia_numbers is not None:
        n_observed = int(result.filter(pl.col("gaia_has_target") == True).height)  # noqa: E712
        logger.info("Targets present in Gaia SSO catalog: %d / %d", n_observed, result.height)
        n_both = int(
            result.filter(
                (pl.col("priority_by_impulse_score"))
                & (pl.col("gaia_has_target") == True)  # noqa: E712
            ).height
        )
        logger.info(
            "Both above-impulse-threshold AND observed by Gaia (prime targets): %d / %d",
            n_both,
            result.height,
        )

    # Pretty ranking table (top 10)
    logger.info("Top candidates (by deflection_score):")
    header = (
        f"{'rk':>3}  {'perturber':<22}  {'target':<18}  "
        f"{'date':<19}  {'d_AU':>10}  {'δ_μas':>10}  {'prio':>6}  {'gaia':>5}"
    )
    logger.info(header)
    logger.info("-" * len(header))
    for r in result.head(min(10, result.height)).iter_rows(named=True):
        pert_label = f"({r['perturber_number']}) {r['perturber_name']}"
        gaia_lbl = (
            "?" if r["gaia_has_target"] is None else ("yes" if r["gaia_has_target"] else "no")
        )
        logger.info(
            "%3d  %-22s  %-18s  %-19s  %10.6f  %10.1f  %6s  %5s",
            r["rank"],
            pert_label[:22],
            (r["target_designation"] or "")[:18],
            r["date_utc"][:19],
            r["dist_au"],
            r["deflection_muas"],
            "yes" if r["priority_by_impulse_score"] else "no",
            gaia_lbl,
        )

    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config.yaml", help="Path to config YAML.")
    parser.add_argument(
        "--top-n",
        type=int,
        default=20,
        help="Number of top-ranked Category B candidates to analyze.",
    )
    parser.add_argument(
        "--with-jpl",
        action="store_true",
        help="Also cross-check geometry against JPL Horizons (rate-limited).",
    )
    parser.add_argument(
        "--half-window-days",
        type=float,
        default=2.0,
        help="Half-window for JPL Horizons vector queries (days, only with --with-jpl).",
    )
    parser.add_argument(
        "--encounters",
        default=str(_ENCOUNTERS_PATH),
        help="Path to the filtered novel-encounters CSV.",
    )
    parser.add_argument(
        "--output",
        default=str(_OUTPUT_PATH),
        help="Output CSV path.",
    )
    parser.add_argument(
        "--gaia-sso",
        default=str(_GAIA_SSO_PATH),
        help="Path to the Gaia DR3 SSO parquet (for target lookup).",
    )
    args = parser.parse_args()

    cfg = load_config(args.config)
    rate_limit_s = float(cfg.sources.jpl_horizons.rate_limit_seconds)
    if args.with_jpl:
        logger.info("JPL cross-check ENABLED  (rate-limit: %.2f s)", rate_limit_s)
    else:
        logger.info("JPL cross-check disabled (pass --with-jpl to enable)")

    encounters_path = Path(args.encounters)
    output_path = Path(args.output)
    gaia_sso_path = Path(args.gaia_sso)

    if not encounters_path.exists():
        logger.error("Encounters file not found: %s", encounters_path)
        return 1

    analyze_candidates(
        encounters_path=encounters_path,
        output_path=output_path,
        top_n=args.top_n,
        with_jpl=args.with_jpl,
        rate_limit_s=rate_limit_s,
        half_window_days=args.half_window_days,
        gaia_sso_path=gaia_sso_path,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
