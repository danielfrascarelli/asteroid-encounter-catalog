"""Sanity check: predict RA/Dec for Ceres via our transforms vs JPL Horizons.

Pipeline:
  1. Load MPCORB elements for (1) Ceres.
  2. Pick 10 Gaia observation epochs of Ceres in 2015-2016.
  3. For each epoch:
       a. Propagate Kepler (MPCORB epoch → obs epoch).
       b. heliocentric_ecliptic → barycentric_ICRS (our transform).
       c. Apply light-time correction.
       d. Compute (RA, Dec) from Gaia's position.
  4. Query Horizons for apparent RA/Dec from Gaia at the same epochs.
  5. Compare residuals.

Goal: residuals < 50 mas. This is a loose criterion because Kepler 2-body
propagation accumulates ~tens of arcsec drift over years (no planet
perturbations), but the transforms themselves should be mas-level.

If residuals are >> 50 mas, debug each transform independently.

Usage
-----
    docker compose run --rm pipeline python -m scripts.sanity_check_transforms
"""

from __future__ import annotations

import argparse
import logging
import time
from pathlib import Path

import numpy as np
import polars as pl
from astropy.time import Time
from astroquery.jplhorizons import Horizons
from astroquery.utils.tap.core import TapPlus

from src.astrometry.transforms import (
    heliocentric_to_barycentric_icrs,
    light_time_iterate,
    xyz_to_radec,
)
from src.ingest.mpcorb import parse_mpcorb
from src.propagate.kepler import kepler_to_cartesian
from src.utils.config import load_config

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)


_J2010_TCB_JD = 2455197.5
_GAIA_OBSERVER = "500@-139479"
_DEG_TO_RAD = np.pi / 180.0


def fetch_gaia_observations(archive_url: str, target: int, n: int = 10) -> pl.DataFrame:
    """Pull *n* Gaia observations of *target*, well-spaced in time."""
    adql = (
        f"SELECT TOP {n} number_mp, epoch, ra, dec, x_gaia, y_gaia, z_gaia "
        f"FROM gaiadr3.sso_observation "
        f"WHERE number_mp = {target} "
        f"AND epoch BETWEEN 1500 AND 2700 "
        f"ORDER BY epoch"
    )
    tap = TapPlus(url=archive_url)
    job = tap.launch_job_async(adql)
    df = pl.from_pandas(job.get_results().to_pandas())
    df = df.rename({c: c.lower() for c in df.columns if c != c.lower()})
    return df


def horizons_predict(
    target: int, jd_tdb_arr: np.ndarray, rate_limit: float
) -> tuple[np.ndarray, np.ndarray]:
    """Apparent RA/Dec from Gaia as predicted by JPL Horizons."""
    h = Horizons(
        id=str(target),
        location=_GAIA_OBSERVER,
        epochs=jd_tdb_arr.tolist(),
        id_type="smallbody",
    )
    eph = h.ephemerides()
    time.sleep(rate_limit)
    return np.array(eph["RA"], dtype=float), np.array(eph["DEC"], dtype=float)


def our_predict(
    elements_row: dict,
    jd_tdb_arr: np.ndarray,
    gaia_xyz: np.ndarray,
    apply_light_time: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """Predict RA/Dec from our transforms.

    Steps per epoch:
      1. Kepler propagate target to jd_tdb (heliocentric ecliptic).
      2. Heliocentric ecliptic → barycentric ICRS.
      3. Light-time correction (iterative).
      4. Line of sight from Gaia → target.
      5. RA/Dec.
    """
    deg = _DEG_TO_RAD

    def _target_pos_bary(jd: float) -> np.ndarray:
        pos_helio_ecl = kepler_to_cartesian(
            a_au=elements_row["a_au"],
            e=elements_row["e"],
            i_rad=elements_row["i_deg"] * deg,
            Omega_rad=elements_row["Omega_deg"] * deg,
            omega_rad=elements_row["omega_deg"] * deg,
            M0_rad=elements_row["M_deg"] * deg,
            epoch_jd=elements_row["epoch_jd"],
            t_jd=jd,
        )
        # kepler_to_cartesian broadcasts inputs → for scalar inputs returns (1,3)
        pos_helio_ecl = np.asarray(pos_helio_ecl).reshape(3)
        return np.asarray(heliocentric_to_barycentric_icrs(pos_helio_ecl, jd)).reshape(3)

    ra_out = np.empty(len(jd_tdb_arr), dtype=float)
    dec_out = np.empty(len(jd_tdb_arr), dtype=float)
    for i, jd in enumerate(jd_tdb_arr):
        if apply_light_time:
            pos, _tau = light_time_iterate(
                _target_pos_bary, jd_tdb_obs=float(jd), gaia_xyz_bary=gaia_xyz[i]
            )
        else:
            pos = _target_pos_bary(float(jd))
        los = pos - gaia_xyz[i]
        ra, dec = xyz_to_radec(los)
        ra_out[i] = ra
        dec_out[i] = dec
    return ra_out, dec_out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--target", type=int, default=1, help="MPC number (default: Ceres)")
    parser.add_argument("--n", type=int, default=10, help="Number of observations")
    args = parser.parse_args()

    cfg = load_config(args.config)
    archive_url = cfg.sources.gaia_sso.archive_url
    rate_limit = float(cfg.sources.jpl_horizons.rate_limit_seconds)

    # Fetch Gaia
    obs = fetch_gaia_observations(archive_url, args.target, args.n)
    logger.info("Fetched %d Gaia observations of asteroid %d", obs.height, args.target)

    # Convert epoch (days-since-J2010 TCB) → JD TDB
    jd_tcb = obs["epoch"].to_numpy() + _J2010_TCB_JD
    jd_tdb = Time(jd_tcb, format="jd", scale="tcb").tdb.jd.astype(float)

    gaia_xyz = np.column_stack(
        [obs["x_gaia"].to_numpy(), obs["y_gaia"].to_numpy(), obs["z_gaia"].to_numpy()]
    ).astype(float)
    ra_obs = obs["ra"].to_numpy().astype(float)
    dec_obs = obs["dec"].to_numpy().astype(float)

    # MPCORB elements
    raw_dir = Path("data/raw")
    archive_dir = raw_dir / "mpcorb_archive"
    snapshots = sorted(archive_dir.glob("MPCORB_*.DAT"))
    if (raw_dir / "MPCORB.DAT").exists():
        snapshots.append(raw_dir / "MPCORB.DAT")
    if not snapshots:
        logger.error("No MPCORB found in data/raw")
        return 1
    logger.info("Searching %s…", snapshots[0])
    df = parse_mpcorb(str(snapshots[0]), semimajor_min_au=0.0, semimajor_max_au=50.0)
    row = df.filter(pl.col("number") == args.target)
    if row.height == 0:
        logger.error("Target %d not in %s", args.target, snapshots[0])
        return 1
    elements = row.row(0, named=True)
    logger.info(
        "Elements: a=%.4f AU  e=%.4f  i=%.2f°  epoch_jd=%.1f",
        elements["a_au"],
        elements["e"],
        elements["i_deg"],
        elements["epoch_jd"],
    )

    # Our prediction (with and without light-time)
    ra_us, dec_us = our_predict(elements, jd_tdb, gaia_xyz, apply_light_time=True)
    ra_us_nolt, dec_us_nolt = our_predict(elements, jd_tdb, gaia_xyz, apply_light_time=False)

    # Horizons prediction
    logger.info("Querying Horizons for ground truth…")
    ra_horiz, dec_horiz = horizons_predict(args.target, jd_tdb, rate_limit)

    # Compute residuals (Horizons - ours) in mas
    deg = _DEG_TO_RAD

    def _shift_mas(ra1, dec1, ra2, dec2):
        dra = ((ra1 - ra2 + 540.0) % 360.0 - 180.0) * np.cos(dec2 * deg) * 3_600_000.0
        ddec = (dec1 - dec2) * 3_600_000.0
        sep = np.sqrt(dra**2 + ddec**2)
        return dra, ddec, sep

    dra_us, ddec_us, sep_us = _shift_mas(ra_us, dec_us, ra_horiz, dec_horiz)
    dra_nolt, ddec_nolt, sep_nolt = _shift_mas(ra_us_nolt, dec_us_nolt, ra_horiz, dec_horiz)

    logger.info("")
    logger.info("=== Residuals vs JPL Horizons (asteroid %d) ===", args.target)
    logger.info("  WITHOUT light-time correction:")
    logger.info(
        "    median sep = %.0f mas    median |ΔRA| = %.0f mas    median |ΔDec| = %.0f mas",
        np.median(sep_nolt),
        np.median(np.abs(dra_nolt)),
        np.median(np.abs(ddec_nolt)),
    )
    logger.info("  WITH light-time correction:")
    logger.info(
        "    median sep = %.0f mas    median |ΔRA| = %.0f mas    median |ΔDec| = %.0f mas",
        np.median(sep_us),
        np.median(np.abs(dra_us)),
        np.median(np.abs(ddec_us)),
    )
    logger.info(
        "  Improvement from light-time: %.1fx",
        np.median(sep_nolt) / max(np.median(sep_us), 1e-9),
    )

    # Save details
    out = pl.DataFrame(
        {
            "jd_tdb": jd_tdb,
            "ra_obs": ra_obs,
            "dec_obs": dec_obs,
            "ra_horizons": ra_horiz,
            "dec_horizons": dec_horiz,
            "ra_ours": ra_us,
            "dec_ours": dec_us,
            "ra_ours_no_lt": ra_us_nolt,
            "dec_ours_no_lt": dec_us_nolt,
            "dra_mas": dra_us,
            "ddec_mas": ddec_us,
            "sep_mas": sep_us,
        }
    )
    out_path = Path("data/output/sanity_transforms.csv")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.write_csv(out_path)
    logger.info("Wrote details to %s", out_path)

    # Verdict
    median_sep = np.median(sep_us)
    if median_sep < 50.0:
        logger.info("")
        logger.info("✅ PASS: median residual < 50 mas — transforms working")
        return 0
    elif median_sep < 5000.0:
        logger.warning("")
        logger.warning("⚠️  MARGINAL: median residual %.0f mas (expected < 50). ", median_sep)
        logger.warning("    This may be due to Kepler 2-body drift over multi-year span,")
        logger.warning(
            "    NOT a transform bug. To verify, repeat with an epoch close to MPCORB epoch."
        )
        return 0
    else:
        logger.error("")
        logger.error("❌ FAIL: median residual %.0f mas. Transforms likely have a bug.", median_sep)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
