"""Validate Category A novel encounters against JPL Horizons.

Category A = encounters where the perturber has a *known* mass from in-situ
measurements (Dawn mission): (1) Ceres and (4) Vesta.  These rows in
``data/output/relevant_novel_encounters.csv`` are flagged with
``mass_unknown == false``.

For every Cat A row this script:

1. Queries JPL Horizons (DE440 ephemeris) for heliocentric state vectors of
   both asteroids on a dense grid (±2 days, 30-min step) around our predicted
   close-approach epoch and recomputes the true minimum-distance geometry.
2. Compares the JPL minimum distance to ours and writes ``delta_au`` and
   ``jpl_date_utc``.
3. Computes the *expected* astrometric deflection of the target by the
   known-mass perturber using the gravitational-deflection formula

       δ = 2 G M / (v² b)        [radians]

   where M is the Dawn-derived mass, v is the relative velocity at closest
   approach, and b is the impact parameter (the minimum distance itself for
   hyperbolic-grazing encounters at MBA velocities).
4. Flags ``viable_signal`` whenever the deflection ≥ 100 μas, i.e. the
   per-transit astrometric precision of Gaia for a typical MBA.

The output is ``data/output/cat_a_jpl_validation.csv`` plus a summary table
and the MAE of our minimum-distance estimate vs. JPL.

Usage
-----
    docker compose run --rm pipeline python -m scripts.validate_novel_a
    docker compose run --rm pipeline python -m scripts.validate_novel_a \\
        --config config.yaml --half-window-days 2.0
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
_RAD_TO_MUAS = 2.06265e11  # μas per radian  (= 180/π × 3600 × 1e6)

# Known perturber masses (Dawn mission gravity science).
_MASS_KG: dict[int, float] = {
    1: 9.384e20,  # (1) Ceres
    4: 2.591e20,  # (4) Vesta
}
_NAME: dict[int, str] = {1: "Ceres", 4: "Vesta"}

# Gaia per-transit astrometric precision for a typical MBA.
_GAIA_PRECISION_MUAS = 100.0


# ---------------------------------------------------------------------------
# Horizons helpers (adapted from scripts/validate_jpl_horizons.py)
# ---------------------------------------------------------------------------


def _horizons_vectors(
    asteroid_id: int,
    jd_center: float,
    half_window_days: float = 2.0,
    step: str = "30m",
) -> tuple[np.ndarray, np.ndarray]:
    """Fetch heliocentric (Sun-centred, @10) state vectors around *jd_center*.

    Returns
    -------
    (jd_array, xyz_array)
        ``jd_array`` is JD TDB; ``xyz_array`` is (N, 3) in AU.
    """
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
# Deflection model
# ---------------------------------------------------------------------------


def compute_deflection_muas(mass_kg: float, rel_vel_km_s: float, dist_au: float) -> float:
    """Linear gravitational deflection in micro-arcseconds.

    Uses ``δ = 2 G M / (v² b)`` with ``b`` = closest-approach distance.
    """
    if not (math.isfinite(mass_kg) and math.isfinite(rel_vel_km_s) and math.isfinite(dist_au)):
        return float("nan")
    if rel_vel_km_s <= 0.0 or dist_au <= 0.0:
        return float("nan")
    b_m = dist_au * _AU_M
    v_ms = rel_vel_km_s * 1000.0
    delta_rad = 2.0 * _G * mass_kg / (v_ms * v_ms * b_m)
    return float(delta_rad * _RAD_TO_MUAS)


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------


def _date_utc_to_jd_tdb(date_utc: str) -> float:
    """Convert an ISO UTC string to JD in TDB scale."""
    return float(Time(date_utc, format="iso", scale="utc").tdb.jd)


def _jd_tdb_to_utc_iso(jd_tdb: float) -> str:
    """Convert JD TDB back to ISO UTC string."""
    return Time(jd_tdb, format="jd", scale="tdb").utc.iso[:19]


def validate_cat_a(
    encounters_path: Path,
    output_path: Path,
    rate_limit_s: float,
    half_window_days: float,
) -> pl.DataFrame:
    """Run the full Cat A validation pipeline and return the results DataFrame."""
    logger.info("Loading encounters from %s", encounters_path)
    df = pl.read_csv(encounters_path)
    logger.info("Total novel encounters: %d", df.height)

    cat_a = df.filter(~pl.col("mass_unknown"))
    logger.info("Category A rows (mass_unknown == false): %d", cat_a.height)

    if cat_a.is_empty():
        logger.warning("No Category A encounters found; nothing to validate.")
        return pl.DataFrame()

    rows: list[dict] = []
    for i, row in enumerate(cat_a.iter_rows(named=True), start=1):
        perturber = int(row["number_1"])
        target_no = row["number_2"]
        target_designation = row["designation_2"]
        date_utc = row["date_utc"]
        our_dist_au = float(row["dist_au"])
        rel_vel_km_s = float(row["rel_vel_km_s"])

        perturber_name = _NAME.get(perturber, f"#{perturber}")
        mass_kg = _MASS_KG.get(perturber, float("nan"))

        logger.info(
            "[%d/%d] %s vs %s on %s (our dist=%.6f AU)",
            i,
            cat_a.height,
            perturber_name,
            target_designation,
            date_utc,
            our_dist_au,
        )

        jd_center = _date_utc_to_jd_tdb(date_utc)

        jpl_dist_au: float = float("nan")
        jpl_date_utc: str | None = None
        delta_au: float = float("nan")

        if target_no is None or (isinstance(target_no, float) and math.isnan(target_no)):
            logger.warning("  target has no MPC number; skipping JPL query")
        else:
            try:
                target_id = int(target_no)
                jpl_dist_au, jpl_jd = _jpl_min_distance(
                    perturber, target_id, jd_center, half_window_days, "30m"
                )
                jpl_date_utc = _jd_tdb_to_utc_iso(jpl_jd)
                delta_au = our_dist_au - jpl_dist_au
                logger.info(
                    "  JPL min dist = %.6f AU at %s  (Δ = %+.2e AU)",
                    jpl_dist_au,
                    jpl_date_utc,
                    delta_au,
                )
            except Exception as exc:  # noqa: BLE001
                logger.warning("  JPL Horizons query failed: %s", exc)
            finally:
                time.sleep(rate_limit_s)

        deflection_muas = compute_deflection_muas(mass_kg, rel_vel_km_s, our_dist_au)
        viable = bool(math.isfinite(deflection_muas) and deflection_muas >= _GAIA_PRECISION_MUAS)
        logger.info(
            "  expected deflection = %.1f μas  (viable signal: %s)",
            deflection_muas,
            viable,
        )

        rows.append(
            {
                "perturber_number": perturber,
                "perturber_name": perturber_name,
                "target_number": (
                    int(target_no)
                    if target_no is not None
                    and not (isinstance(target_no, float) and math.isnan(target_no))
                    else None
                ),
                "target_designation": target_designation,
                "date_utc": date_utc,
                "our_dist_au": our_dist_au,
                "jpl_dist_au": jpl_dist_au,
                "delta_au": delta_au,
                "jpl_date_utc": jpl_date_utc,
                "rel_vel_km_s": rel_vel_km_s,
                "expected_deflection_muas": deflection_muas,
                "gaia_precision_muas": _GAIA_PRECISION_MUAS,
                "viable_signal": viable,
            }
        )

    result = pl.DataFrame(rows)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.write_csv(output_path)
    logger.info("Wrote %d rows to %s", result.height, output_path)

    # Summary
    finite_mask = result["delta_au"].is_finite() & result["delta_au"].is_not_null()
    finite_deltas = result.filter(finite_mask)["delta_au"]
    if finite_deltas.len() > 0:
        mae_au = float(finite_deltas.abs().mean())
        logger.info("MAE(our − JPL) over %d rows: %.3e AU", finite_deltas.len(), mae_au)
    else:
        logger.warning("No successful JPL queries; cannot compute MAE.")

    n_viable = int(result["viable_signal"].sum())
    logger.info(
        "Viable signals (≥ %.0f μas): %d / %d",
        _GAIA_PRECISION_MUAS,
        n_viable,
        result.height,
    )

    # Pretty per-perturber breakdown
    for pert in sorted(result["perturber_number"].unique().to_list()):
        sub = result.filter(pl.col("perturber_number") == pert)
        logger.info(
            "  %s: %d encounters, %d viable",
            _NAME.get(pert, f"#{pert}"),
            sub.height,
            int(sub["viable_signal"].sum()),
        )

    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config.yaml", help="Path to config YAML.")
    parser.add_argument(
        "--half-window-days",
        type=float,
        default=2.0,
        help="Half-window for JPL Horizons vector queries (days).",
    )
    parser.add_argument(
        "--encounters",
        default="data/output/relevant_novel_encounters.csv",
        help="Path to the filtered novel-encounters CSV.",
    )
    parser.add_argument(
        "--output",
        default="data/output/cat_a_jpl_validation.csv",
        help="Output CSV path.",
    )
    args = parser.parse_args()

    cfg = load_config(args.config)
    rate_limit_s = float(cfg.sources.jpl_horizons.rate_limit_seconds)
    logger.info("JPL Horizons rate-limit: %.2f s between requests", rate_limit_s)

    encounters_path = Path(args.encounters)
    output_path = Path(args.output)

    if not encounters_path.exists():
        logger.error("Encounters file not found: %s", encounters_path)
        return 1

    validate_cat_a(
        encounters_path=encounters_path,
        output_path=output_path,
        rate_limit_s=rate_limit_s,
        half_window_days=args.half_window_days,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
