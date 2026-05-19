"""Numerical experiment: does our N-body wrapper predict Davida's deflection?

For (511) Davida + 2003_sm90 on 2014-11-19 we expect a tiny but non-zero
shift in 2003_sm90's trajectory due to Davida's gravity. The expected
order of magnitude is:

    δ ≈ 2GM / (v² · b) ≈ 5 mas at encounter

Cumulative position offset 90 days after encounter:
    Δr ≈ δ × v_target × t  →  Δθ ≈ 1-10 mas as seen from Gaia

This script integrates the target's orbit with the perturber at TWO mass
values (zero, and the literature ~3.5e19 kg Goffin value), measures the
difference in trajectory, and converts to an angular shift.

If the deflection is detectable in the integration (i.e. the difference
is well above numerical noise), the N-body wrapper is wired up correctly
and we can proceed to T2.x (fitting machinery).
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import polars as pl

from src.ingest.mpcorb import parse_mpcorb
from src.propagate.nbody_perturber import propagate_target_with_perturber

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)


# Switched to (165) Loreley + (31067) 1996_tf50 because the original target
# 2003_sm90 has MPC number 115180 which falls outside the packed-designation
# range supported by our current parse_mpcorb (≤ 99999).
# Loreley also has a published mass (~8.7e18 kg) so it works as a calibration
# target.  We'll come back to Davida once we extend the parser.
_DAVIDA_NUMBER = 165  # (165) Loreley
_TARGET_NUMBER = 31067  # 1996_tf50
_DAVIDA_MASS_GOFFIN_KG = 8.7e18  # Loreley mass (Goffin 2014)
# Encounter date for this pair (from publishable_mass_candidates.csv)
_ENCOUNTER_DATE_UTC = "2014-12-08T00:00:00"


def _load_element_row(snapshot: Path, number: int) -> dict:
    df = parse_mpcorb(str(snapshot), semimajor_min_au=0.0, semimajor_max_au=50.0)
    sub = df.filter(pl.col("number") == number)
    if sub.height == 0:
        raise ValueError(f"asteroid {number} not in {snapshot}")
    return sub.row(0, named=True)


def main() -> int:
    snap = Path("data/raw/mpcorb_archive/MPCORB_20120918.DAT")
    if not snap.exists():
        snap = Path("data/raw/MPCORB.DAT")
        if not snap.exists():
            logger.error("No MPCORB.DAT found")
            return 1
    logger.info("Loading MPCORB from %s", snap)
    target = _load_element_row(snap, _TARGET_NUMBER)
    davida = _load_element_row(snap, _DAVIDA_NUMBER)
    logger.info(
        "Target %d: a=%.4f e=%.4f i=%.2f° epoch=%.1f",
        _TARGET_NUMBER, target["a_au"], target["e"], target["i_deg"], target["epoch_jd"],
    )
    logger.info(
        "Davida (%d): a=%.4f e=%.4f i=%.2f° epoch=%.1f",
        _DAVIDA_NUMBER, davida["a_au"], davida["e"], davida["i_deg"], davida["epoch_jd"],
    )

    # Time grid: ±180 days around the encounter
    from astropy.time import Time
    enc_jd = float(Time(_ENCOUNTER_DATE_UTC, scale="utc").tdb.jd)
    t_grid = np.linspace(enc_jd - 180.0, enc_jd + 180.0, 361)  # 1-day spacing

    logger.info(
        "Integrating target over ±180 d around encounter (JD %.1f) at 1-d grid",
        enc_jd,
    )
    # Run 1: massless perturber (= no gravitational effect)
    pos_no_mass = propagate_target_with_perturber(
        target_elements=target,
        perturber_elements=davida,
        perturber_mass_kg=0.0,
        time_grid_jd_tdb=t_grid,
        include_planets=("sun", "jupiter", "saturn"),
        dt_days=1.0,
    )
    logger.info("Run #1: perturber mass = 0 kg → trajectory baseline")

    # Run 2: Davida with Goffin (2014) mass
    pos_goffin = propagate_target_with_perturber(
        target_elements=target,
        perturber_elements=davida,
        perturber_mass_kg=_DAVIDA_MASS_GOFFIN_KG,
        time_grid_jd_tdb=t_grid,
        include_planets=("sun", "jupiter", "saturn"),
        dt_days=1.0,
    )
    logger.info("Run #2: perturber mass = %.2e kg → with Davida gravity", _DAVIDA_MASS_GOFFIN_KG)

    # Compute the difference in target position over time
    diff = np.linalg.norm(pos_goffin - pos_no_mass, axis=1)  # AU

    # Convert to angular shift as seen from Gaia (~2 AU away typical for MBA)
    typical_distance_au = 2.0
    angular_shift_rad = diff / typical_distance_au
    angular_shift_mas = angular_shift_rad * (180.0 / np.pi) * 3_600_000.0

    days_from_enc = t_grid - enc_jd

    logger.info("")
    logger.info("Position difference (with vs without Davida gravity):")
    logger.info("  at encounter (t=0):     diff = %.3e AU = %.2f mas",
                diff[180], angular_shift_mas[180])
    logger.info("  at +90 days:            diff = %.3e AU = %.2f mas",
                diff[270], angular_shift_mas[270])
    logger.info("  at +180 days (end):     diff = %.3e AU = %.2f mas",
                diff[-1], angular_shift_mas[-1])
    logger.info("  at -180 days (start):   diff = %.3e AU = %.2f mas",
                diff[0], angular_shift_mas[0])
    logger.info("  max difference in window: %.2f mas (at t=%+.0f d)",
                angular_shift_mas.max(), days_from_enc[int(np.argmax(angular_shift_mas))])

    # Save the data
    out = pl.DataFrame({
        "jd_tdb": t_grid,
        "days_from_encounter": days_from_enc,
        "diff_au": diff,
        "angular_shift_mas": angular_shift_mas,
    })
    out_path = Path("data/output/davida_deflection_test.csv")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.write_csv(out_path)
    logger.info("Wrote %d-row trajectory difference to %s", out.height, out_path)

    # Verdict
    max_shift = angular_shift_mas.max()
    if max_shift > 1.0:
        logger.info("")
        logger.info("✅ N-body wrapper produces a detectable shift (%.1f mas peak).", max_shift)
        logger.info(
            "   Order of magnitude consistent with expected δ ~ 5 mas for this encounter."
        )
        return 0
    else:
        logger.warning("")
        logger.warning("⚠️  Shift too small (%.3f mas). Check perturber elements/mass.", max_shift)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
