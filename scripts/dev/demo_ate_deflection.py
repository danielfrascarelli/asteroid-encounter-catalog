"""End-to-end demonstration: detect (111) Ate's perturbation on 2000_nt3.

Proof-of-concept that the encounter catalog can drive real asteroid-mass
determinations from Gaia DR3 data.

What this does
--------------
1. Downloads Gaia DR3 SSO observations of asteroid 18105 (2000_nt3) in a
   ±180-day window around its 2016-06-08 close encounter with (111) Ate.
2. Loads the MPCORB osculating orbital elements for asteroid 18105.
3. For each Gaia transit:
      • Propagates the MPCORB orbit as a pure 2-body Kepler problem to the
        observation epoch (no perturbation from (111) Ate).
      • Converts the heliocentric-ecliptic position to barycentric-ICRS by
        adding the Sun's barycentric position (from astropy) and rotating by
        the J2000 obliquity.
      • Subtracts Gaia's barycentric position (stored as x_gaia, y_gaia,
        z_gaia in the SSO table) to get the geometric line of sight.
      • Converts that to an RA/Dec prediction.
4. Computes residuals (observed − predicted) in milliarcseconds.
5. Splits residuals into pre- vs post-encounter and reports a one-sample
   t-statistic for the difference.

If Ate's perturbation is detectable in the existing data, the post-encounter
residuals should be biased relative to the pre-encounter ones.

CAVEATS
-------
- MPCORB osculating elements come from JPL's N-body fit which includes Ceres,
  Pallas, Vesta, Hygiea (the "big 4") and the major planets, but NOT (111)
  Ate. So the residuals SHOULD reveal Ate's signature.
- No light-time correction is applied; the 0.04% effect (~4 mas) is small
  compared to our 5 mas target signal but could matter.
- The propagation epoch in MPCORB is typically ~2026, so we propagate
  backward ~10 years; Kepler 2-body accumulates error from non-included
  perturbations during this window.
- The historical MPCORB snapshot from 2015 would give a better baseline,
  but the current MPCORB is used here for simplicity.

Usage
-----
    docker compose run --rm pipeline python -m scripts.demo_ate_deflection
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import polars as pl
from astropy.coordinates import get_body_barycentric
from astropy.time import Time
from astroquery.utils.tap.core import TapPlus

from src.ingest.mpcorb import parse_mpcorb
from src.propagate.kepler import kepler_to_cartesian
from src.utils.config import load_config

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_TARGET_NUMBER = 18105  # MPC number for 2000 NT3
_PERTURBER_NUMBER = 111  # (111) Ate
_ENCOUNTER_DATE = "2016-06-08T00:00:00"
_HALF_WINDOW_DAYS = 180.0
_BLACKOUT_DAYS = 7.0

# J2010 TCB epoch offset for Gaia SSO `epoch` column.
_J2010_TCB_JD = 2455197.5

# J2000 obliquity of the ecliptic (IAU 1980).
_OBLIQUITY_DEG = 23.43928083
_EPS = np.radians(_OBLIQUITY_DEG)
_COS_EPS = np.cos(_EPS)
_SIN_EPS = np.sin(_EPS)

# rad → mas
_RAD_TO_MAS = 180.0 / np.pi * 3_600_000.0


# ---------------------------------------------------------------------------
# Coordinate transforms
# ---------------------------------------------------------------------------


def ecliptic_to_equatorial(xyz_ecl: np.ndarray) -> np.ndarray:
    """Rotate heliocentric-ecliptic-J2000 to heliocentric-equatorial-J2000."""
    x = xyz_ecl[..., 0]
    y_eq = _COS_EPS * xyz_ecl[..., 1] - _SIN_EPS * xyz_ecl[..., 2]
    z_eq = _SIN_EPS * xyz_ecl[..., 1] + _COS_EPS * xyz_ecl[..., 2]
    return np.stack([x, y_eq, z_eq], axis=-1)


def xyz_to_radec(vec: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Convert an (..., 3) Cartesian vector to (RA, Dec) in degrees."""
    x = vec[..., 0]
    y = vec[..., 1]
    z = vec[..., 2]
    rho = np.sqrt(x * x + y * y)
    ra = np.degrees(np.arctan2(y, x)) % 360.0
    dec = np.degrees(np.arctan2(z, rho))
    return ra, dec


def sun_barycentric_au(jd_tdb: np.ndarray) -> np.ndarray:
    """Return the Sun's barycentric ICRS position (AU) at *jd_tdb*.

    Uses astropy's default ephemeris (built-in low-precision unless the user
    has configured a high-precision ephemeris).
    """
    out = np.empty((len(jd_tdb), 3), dtype=float)
    for i, t in enumerate(jd_tdb):
        sun = get_body_barycentric("sun", Time(t, format="jd", scale="tdb"))
        out[i, 0] = sun.x.to_value("AU")
        out[i, 1] = sun.y.to_value("AU")
        out[i, 2] = sun.z.to_value("AU")
    return out


# ---------------------------------------------------------------------------
# Gaia query
# ---------------------------------------------------------------------------


def fetch_target_observations(
    archive_url: str,
    target_number: int,
    days_min: float,
    days_max: float,
) -> pl.DataFrame:
    """Fetch transits of *target_number* in the days-since-J2010-TCB window."""
    adql = (
        "SELECT number_mp, epoch, epoch_utc, ra, dec, g_mag, "
        "x_gaia, y_gaia, z_gaia "
        "FROM gaiadr3.sso_observation "
        f"WHERE number_mp = {target_number} "
        f"AND epoch BETWEEN {days_min:.6f} AND {days_max:.6f} "
        "ORDER BY epoch"
    )
    logger.info("Querying Gaia TAP for asteroid %d…", target_number)
    tap = TapPlus(url=archive_url)
    job = tap.launch_job_async(adql)
    table = job.get_results()
    df = pl.from_pandas(table.to_pandas())
    df = df.rename({c: c.lower() for c in df.columns if c != c.lower()})
    logger.info("  → %d transits", df.height)
    return df


# ---------------------------------------------------------------------------
# MPCORB elements
# ---------------------------------------------------------------------------


def load_mpcorb_elements(target_number: int) -> pl.DataFrame:
    """Locate and load orbital elements for *target_number* from MPCORB.

    Prefers the historical 2015 snapshot if available (closer epoch reduces
    Kepler-vs-real drift over the 2014-2017 window).
    """
    raw_dir = Path("data/raw")
    archive_dir = raw_dir / "mpcorb_archive"
    candidates: list[Path] = []
    if archive_dir.exists():
        candidates.extend(sorted(archive_dir.glob("MPCORB_*.DAT")))
    main = raw_dir / "MPCORB.DAT"
    if main.exists():
        candidates.append(main)

    if not candidates:
        raise FileNotFoundError("No MPCORB.DAT found under data/raw/")

    for path in candidates:
        logger.info("Searching %s for asteroid %d…", path, target_number)
        df = parse_mpcorb(str(path), semimajor_min_au=0.0, semimajor_max_au=50.0)
        row = df.filter(pl.col("number") == target_number)
        if row.height > 0:
            logger.info("Found in %s", path)
            return row
    raise ValueError(f"Asteroid {target_number} not in any MPCORB snapshot")


# ---------------------------------------------------------------------------
# Prediction
# ---------------------------------------------------------------------------


def predict_radec_from_elements(
    elements_row: pl.DataFrame,
    jd_tdb_array: np.ndarray,
    gaia_xyz: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """For each observation epoch, predict (RA, Dec) of the target from Gaia.

    Parameters
    ----------
    elements_row:
        Single-row DataFrame with MPCORB columns.
    jd_tdb_array:
        (N,) array of observation epochs in JD TDB.
    gaia_xyz:
        (N, 3) Gaia barycentric ICRS position at each epoch in AU.

    Returns
    -------
    (ra_deg, dec_deg)
        (N,) arrays of predicted apparent RA, Dec in degrees.
    """
    row = elements_row.row(0, named=True)
    deg = np.pi / 180.0

    pos_helio_ecl = np.array(
        [
            kepler_to_cartesian(
                a_au=row["a_au"],
                e=row["e"],
                i_rad=row["i_deg"] * deg,
                Omega_rad=row["Omega_deg"] * deg,
                omega_rad=row["omega_deg"] * deg,
                M0_rad=row["M_deg"] * deg,
                epoch_jd=row["epoch_jd"],
                t_jd=t,
            )
            for t in jd_tdb_array
        ]
    )
    # kepler_to_cartesian returns (..., 3) — strip the singleton outer dim
    pos_helio_ecl = pos_helio_ecl.reshape(-1, 3)

    # heliocentric ecliptic → heliocentric equatorial
    pos_helio_eq = ecliptic_to_equatorial(pos_helio_ecl)

    # heliocentric → barycentric (add Sun's barycentric position)
    sun_bary = sun_barycentric_au(jd_tdb_array)
    pos_bary = pos_helio_eq + sun_bary

    # Vector from Gaia to the asteroid
    los = pos_bary - gaia_xyz

    return xyz_to_radec(los)


def _angular_separation_mas(
    ra1_deg: np.ndarray,
    dec1_deg: np.ndarray,
    ra2_deg: np.ndarray,
    dec2_deg: np.ndarray,
) -> np.ndarray:
    """Angular separation between (RA1, Dec1) and (RA2, Dec2) in mas."""
    deg = np.pi / 180.0
    a1 = ra1_deg * deg
    d1 = dec1_deg * deg
    a2 = ra2_deg * deg
    d2 = dec2_deg * deg
    # Haversine on the sphere
    sd2 = np.sin((d2 - d1) / 2.0) ** 2
    sa2 = np.sin((a2 - a1) / 2.0) ** 2
    h = sd2 + np.cos(d1) * np.cos(d2) * sa2
    return 2.0 * np.arcsin(np.minimum(1.0, np.sqrt(h))) * _RAD_TO_MAS


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/output/ate_2000nt3_residuals.csv"),
    )
    args = parser.parse_args()

    cfg = load_config(args.config)
    archive_url = cfg.sources.gaia_sso.archive_url

    # --- Time window ----------------------------------------------------------
    enc_jd_tdb = float(Time(_ENCOUNTER_DATE, scale="utc").tdb.jd)
    enc_days_j2010 = float(Time(_ENCOUNTER_DATE, scale="utc").tcb.jd) - _J2010_TCB_JD
    logger.info(
        "Encounter JD TDB = %.6f  (days since J2010 TCB = %.3f)", enc_jd_tdb, enc_days_j2010
    )

    # --- Fetch Gaia observations ---------------------------------------------
    obs = fetch_target_observations(
        archive_url,
        _TARGET_NUMBER,
        enc_days_j2010 - _HALF_WINDOW_DAYS,
        enc_days_j2010 + _HALF_WINDOW_DAYS,
    )
    if obs.height == 0:
        logger.error("No observations found")
        return 1

    # Convert Gaia epoch (days since J2010 TCB) → JD TDB
    epochs_days = obs["epoch"].to_numpy()
    jd_tcb = epochs_days + _J2010_TCB_JD
    # TCB → TDB conversion via astropy
    t_obs = Time(jd_tcb, format="jd", scale="tcb")
    jd_tdb = t_obs.tdb.jd.astype(float)

    gaia_xyz = np.column_stack(
        [
            obs["x_gaia"].to_numpy().astype(float),
            obs["y_gaia"].to_numpy().astype(float),
            obs["z_gaia"].to_numpy().astype(float),
        ]
    )

    # --- Load MPCORB orbit ----------------------------------------------------
    elements = load_mpcorb_elements(_TARGET_NUMBER)
    logger.info(
        "Loaded elements: a=%.4f AU  e=%.4f  i=%.2f°  epoch_jd=%.1f",
        elements.row(0, named=True)["a_au"],
        elements.row(0, named=True)["e"],
        elements.row(0, named=True)["i_deg"],
        elements.row(0, named=True)["epoch_jd"],
    )

    # --- Predict RA/Dec for every observation ---------------------------------
    logger.info("Propagating Kepler 2-body orbit to %d epochs…", len(jd_tdb))
    ra_pred, dec_pred = predict_radec_from_elements(elements, jd_tdb, gaia_xyz)

    # --- Compute residuals ---------------------------------------------------
    ra_obs = obs["ra"].to_numpy().astype(float)
    dec_obs = obs["dec"].to_numpy().astype(float)

    # Total angular separation
    sep_mas = _angular_separation_mas(ra_pred, dec_pred, ra_obs, dec_obs)

    # Component residuals (mas, on tangent plane around the predicted position)
    deg = np.pi / 180.0
    dra = (ra_obs - ra_pred + 540.0) % 360.0 - 180.0  # signed in degrees
    ddec = dec_obs - dec_pred
    dra_mas = dra * np.cos(dec_pred * deg) * 3_600_000.0  # cos(dec) factor
    ddec_mas = ddec * 3_600_000.0

    days_from_enc = jd_tdb - enc_jd_tdb

    result = pl.DataFrame(
        {
            "jd_tdb": jd_tdb,
            "days_from_encounter": days_from_enc,
            "ra_obs_deg": ra_obs,
            "dec_obs_deg": dec_obs,
            "ra_pred_deg": ra_pred,
            "dec_pred_deg": dec_pred,
            "dra_mas": dra_mas,
            "ddec_mas": ddec_mas,
            "sep_mas": sep_mas,
            "g_mag": obs["g_mag"].to_numpy().astype(float),
        }
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.write_csv(args.output)
    logger.info("Wrote %d rows to %s", result.height, args.output)

    # --- Summary statistics --------------------------------------------------
    before = result.filter(pl.col("days_from_encounter") < -_BLACKOUT_DAYS)
    after = result.filter(pl.col("days_from_encounter") > _BLACKOUT_DAYS)

    def _stats(df: pl.DataFrame, label: str) -> dict:
        if df.height == 0:
            return {}
        mu_dra = float(df["dra_mas"].mean())
        mu_ddec = float(df["ddec_mas"].mean())
        std_dra = float(df["dra_mas"].std())
        std_ddec = float(df["ddec_mas"].std())
        mu_sep = float(df["sep_mas"].mean())
        logger.info(
            "%s (N=%d)  ⟨ΔRA⟩=%+9.1f mas  ⟨ΔDec⟩=%+9.1f mas  ⟨sep⟩=%9.1f mas  "
            "σ(ΔRA)=%.1f  σ(ΔDec)=%.1f",
            label,
            df.height,
            mu_dra,
            mu_ddec,
            mu_sep,
            std_dra,
            std_ddec,
        )
        return {
            "n": df.height,
            "mu_dra": mu_dra,
            "mu_ddec": mu_ddec,
            "std_dra": std_dra,
            "std_ddec": std_ddec,
            "mu_sep": mu_sep,
        }

    stats_before = _stats(before, "BEFORE encounter")
    stats_after = _stats(after, "AFTER  encounter")

    if stats_before and stats_after:
        # Welch t-statistic on the difference of means (per axis)
        for axis, key in [("ΔRA", "dra"), ("ΔDec", "ddec")]:
            mb = stats_before[f"mu_{key}"]
            ma = stats_after[f"mu_{key}"]
            sb = stats_before[f"std_{key}"] / np.sqrt(stats_before["n"])
            sa = stats_after[f"std_{key}"] / np.sqrt(stats_after["n"])
            sigma = np.sqrt(sb * sb + sa * sa)
            t = (ma - mb) / sigma if sigma > 0 else float("inf")
            logger.info("Δ(after − before) on %s = %+.1f mas   t = %+.2fσ", axis, ma - mb, t)

    logger.info("")
    logger.info("Interpretation:")
    logger.info("  Expected Ate signal at encounter ≈ 5 mas (linear post-encounter).")
    logger.info("  Pre-encounter residuals reflect non-Ate perturbations + Kepler-fit drift.")
    logger.info("  Look at sep_mas as a function of days_from_encounter in the output CSV.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
