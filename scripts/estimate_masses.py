"""Post-hoc mass estimation from detected astrometric shifts.

Given the ``data/output/deflection_detections.csv`` produced by
``detect_deflections.py``, this script applies a simple kinematic model to
convert the measured RA/Dec shift into an implied perturber mass:

    Δv     = δ × v_pre = 2GM / (v_pre × b)             [m/s perpendicular]
    Δr(t)  = Δv × t                                    [m perpendicular]
    Δθ(t)  = Δr / d_observer                           [rad]
    ⇒ M    = Δθ × v_pre × b × d_observer / (2 G ⟨t⟩)

Here ⟨t⟩ ≈ (window_after / 2) ≈ 90 days is the mean elapsed time of the
post-encounter observations. ``d_observer`` is approximated as 2.5 AU
(typical target distance from Gaia).

This estimate is rough — it assumes:
- All of the measured shift comes from the perturber (no other contaminating
  perturbers, no orbital fit residual).
- The deflection is impulsive at the encounter epoch.
- The geometry factors (angle between Δv and line of sight) are unity.

Even with these limitations, the implied mass should be within an order
of magnitude of the true value, and can be cross-calibrated against known
masses (Ceres / Vesta in Category A).

Usage
-----
    docker compose run --rm pipeline python -m scripts.estimate_masses
"""

from __future__ import annotations

import argparse
import logging
import math
from pathlib import Path

import polars as pl

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)


_G = 6.674e-11  # N m² / kg²
_AU_M = 1.495978707e11  # meters per AU
_MAS_TO_RAD = math.pi / (180.0 * 3_600_000.0)

_T_EFFECTIVE_SECONDS = 90.0 * 86400.0  # mean post-encounter time
_D_OBSERVER_AU = 2.5  # typical target distance from Gaia
_D_OBSERVER_M = _D_OBSERVER_AU * _AU_M


# Known masses (Dawn mission) for self-calibration
_KNOWN_MASSES_KG = {1: 9.384e20, 4: 2.591e20}


def implied_mass_kg(shift_mas: float, rel_vel_km_s: float, dist_au: float) -> float:
    """Convert measured astrometric shift to implied perturber mass (kinematic)."""
    if not (math.isfinite(shift_mas) and math.isfinite(rel_vel_km_s) and math.isfinite(dist_au)):
        return float("nan")
    delta_theta_rad = abs(shift_mas) * _MAS_TO_RAD
    v_ms = rel_vel_km_s * 1000.0
    b_m = dist_au * _AU_M
    return delta_theta_rad * v_ms * b_m * _D_OBSERVER_M / (2.0 * _G * _T_EFFECTIVE_SECONDS)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument(
        "--detections",
        type=Path,
        default=Path("data/output/deflection_detections.csv"),
    )
    p.add_argument(
        "--candidates",
        type=Path,
        default=Path("data/output/publishable_mass_candidates.csv"),
    )
    p.add_argument(
        "--output",
        type=Path,
        default=Path("data/output/mass_estimates.csv"),
    )
    args = p.parse_args()

    if not args.detections.exists():
        logger.error("Detections file not found: %s", args.detections)
        return 1
    detections = pl.read_csv(args.detections)

    candidates = pl.read_csv(args.candidates) if args.candidates.exists() else None

    # Join distance / velocity from the candidates table
    if candidates is not None:
        cand_min = candidates.select(
            [
                "perturber_number",
                "target_number",
                "dist_au",
                pl.col("deflection_muas").alias("delta_expected_muas"),
            ]
        )
        # Rel velocity comes from relevant_novel_encounters via mass_candidates;
        # cat_a uses validate_novel_a output instead. Use candidates table where
        # it exists; fall back to deflection_detections (which has expected_muas
        # but not rel_vel directly).
        det = detections.join(cand_min, on=["perturber_number", "target_number"], how="left")
    else:
        det = detections

    # rel_vel_km_s must be looked up; relevant_novel_encounters has it
    rel_path = Path("data/output/relevant_novel_encounters.csv")
    if rel_path.exists():
        rel = (
            pl.read_csv(rel_path)
            .select(
                [
                    pl.col("number_1").alias("perturber_number"),
                    pl.col("number_2").alias("target_number"),
                    "rel_vel_km_s",
                ]
            )
            .unique(subset=["perturber_number", "target_number"])
        )
        det = det.join(rel, on=["perturber_number", "target_number"], how="left")

    rows: list[dict] = []
    for r in det.iter_rows(named=True):
        # Use the absolute total shift (root sum square of RA, Dec)
        shift_dra = r.get("shift_dra_mas") or 0.0
        shift_ddec = r.get("shift_ddec_mas") or 0.0
        total_shift_mas = math.sqrt(shift_dra * shift_dra + shift_ddec * shift_ddec)

        v = r.get("rel_vel_km_s")
        b = r.get("dist_au")
        if v is None or b is None or not math.isfinite(v) or not math.isfinite(b):
            m_est = float("nan")
        else:
            m_est = implied_mass_kg(total_shift_mas, float(v), float(b))

        known = _KNOWN_MASSES_KG.get(r["perturber_number"])
        ratio = (
            (m_est / known)
            if (known is not None and math.isfinite(m_est) and m_est > 0)
            else float("nan")
        )

        rows.append(
            {
                "perturber_number": r["perturber_number"],
                "perturber_name": r["perturber_name"],
                "target_number": r["target_number"],
                "target_designation": r["target_designation"],
                "date_utc": r["date_utc"],
                "shift_dra_mas": shift_dra,
                "shift_ddec_mas": shift_ddec,
                "total_shift_mas": total_shift_mas,
                "t_dra": r.get("t_dra"),
                "t_ddec": r.get("t_ddec"),
                "rel_vel_km_s": v,
                "dist_au": b,
                "implied_mass_kg": m_est,
                "known_mass_kg": known if known is not None else float("nan"),
                "ratio_implied_to_known": ratio,
                "detection": r["detection"],
            }
        )

    out = pl.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.write_csv(args.output)
    logger.info("Wrote %d rows to %s", out.height, args.output)

    # Pretty per-detection table
    detected = out.filter(pl.col("detection") == "yes").sort("total_shift_mas", descending=True)
    logger.info("")
    logger.info("Implied masses for the %d detected candidates:", detected.height)
    header = (
        f"  {'Perturber':<22} {'Target':<18} {'Δθ (mas)':>10} {'M_implied (kg)':>18} {'known':>12}"
    )
    logger.info(header)
    logger.info("-" * len(header))
    for r in detected.iter_rows(named=True):
        pert_label = f"({r['perturber_number']}) {r['perturber_name']}"
        known_str = f"{r['known_mass_kg']:.2e}" if math.isfinite(r["known_mass_kg"]) else "—"
        logger.info(
            "  %-22s %-18s %10.1f %18.3e %12s",
            pert_label[:22],
            (r["target_designation"] or "")[:18],
            r["total_shift_mas"],
            r["implied_mass_kg"],
            known_str,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
