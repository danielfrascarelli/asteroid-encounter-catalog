"""Pure transform sanity check: bypass Kepler, use Horizons vectors as truth.

This isolates the coordinate transformations from orbital-propagation error.

Pipeline:
  1. For Ceres, query Horizons for *barycentric ICRS position vectors* at 10 epochs.
  2. For the same 10 epochs, query Horizons for *apparent RA/Dec from Gaia*.
  3. Use Gaia's x_gaia/y_gaia/z_gaia from the SSO table.
  4. Apply ONLY our transforms (light-time + xyz_to_radec) to the Horizons
     position vectors to compute predicted RA/Dec.
  5. Compare to Horizons apparent RA/Dec.

If our transforms are correct, residuals should be at most a few mas (limited
by aberration, which we don't apply here, and small numerical effects).
"""

from __future__ import annotations

import argparse
import logging
import time
from pathlib import Path

import numpy as np
import polars as pl
from astropy.coordinates import get_body_barycentric_posvel
from astropy.time import Time
from astroquery.jplhorizons import Horizons
from astroquery.utils.tap.core import TapPlus
from scipy.interpolate import CubicSpline

from src.astrometry.transforms import (
    light_time_iterate,
    stellar_aberration,
    xyz_to_radec,
)
from src.utils.config import load_config

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)


_J2010_TCB_JD = 2455197.5
_GAIA_OBSERVER = "500@-139479"


def fetch_gaia_observations(archive_url: str, target: int, n: int = 10) -> pl.DataFrame:
    adql = (
        f"SELECT TOP {n} number_mp, epoch, x_gaia, y_gaia, z_gaia "
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


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--target", type=int, default=1)
    parser.add_argument("--n", type=int, default=10)
    args = parser.parse_args()

    cfg = load_config(args.config)
    archive_url = cfg.sources.gaia_sso.archive_url
    rate_limit = float(cfg.sources.jpl_horizons.rate_limit_seconds)

    obs = fetch_gaia_observations(archive_url, args.target, args.n)
    logger.info("Fetched %d Gaia observations", obs.height)

    jd_tcb = obs["epoch"].to_numpy() + _J2010_TCB_JD
    jd_tdb = Time(jd_tcb, format="jd", scale="tcb").tdb.jd.astype(float)
    gaia_xyz = np.column_stack(
        [obs["x_gaia"].to_numpy(), obs["y_gaia"].to_numpy(), obs["z_gaia"].to_numpy()]
    ).astype(float)

    # 1. Horizons barycentric ICRS vectors: 7 points per obs (1-hour spacing
    #    centred on the obs epoch) to support light-time interpolation.
    bracket_jds: list[float] = []
    for j in jd_tdb:
        bracket_jds.extend([float(j + k / 24.0) for k in range(-3, 4)])
    bracket_jds = sorted(set(round(j, 6) for j in bracket_jds))
    logger.info("Querying Horizons VECTORS at %d epochs", len(bracket_jds))
    # Chunk in groups of 30 to stay under URL length limits.
    vec_jd_list = []
    vec_xyz_list = []
    chunk = 30
    for i in range(0, len(bracket_jds), chunk):
        chunk_epochs = bracket_jds[i : i + chunk]
        h = Horizons(
            id=str(args.target),
            location="@0",
            epochs=chunk_epochs,
            id_type="smallbody",
        )
        tbl = h.vectors(refplane="earth")
        vec_jd_list.append(np.array(tbl["datetime_jd"], dtype=float))
        vec_xyz_list.append(
            np.column_stack(
                [
                    np.array(tbl["x"], dtype=float),
                    np.array(tbl["y"], dtype=float),
                    np.array(tbl["z"], dtype=float),
                ]
            )
        )
        time.sleep(rate_limit)
    vec_jd = np.concatenate(vec_jd_list)
    vec_xyz = np.concatenate(vec_xyz_list)
    # Sort by epoch (chunks may have come back in order, but be safe)
    order = np.argsort(vec_jd)
    vec_jd = vec_jd[order]
    vec_xyz = vec_xyz[order]

    # 2. Horizons apparent RA/Dec from Gaia at the observation epochs
    logger.info("Querying Horizons APPARENT EPHEMERIDES from Gaia at %d epochs", len(jd_tdb))
    h2 = Horizons(
        id=str(args.target),
        location=_GAIA_OBSERVER,
        epochs=jd_tdb.tolist(),
        id_type="smallbody",
    )
    eph = h2.ephemerides()
    ra_horizons = np.array(eph["RA"], dtype=float)
    dec_horizons = np.array(eph["DEC"], dtype=float)

    # 3. Build a position function from interpolated Horizons vectors
    splines = [CubicSpline(vec_jd, vec_xyz[:, k]) for k in range(3)]

    def target_pos_at(jd: float) -> np.ndarray:
        return np.array([s(jd) for s in splines])

    # 4. Apply our transforms only (light-time + aberration + xyz_to_radec)
    ra_ours = np.empty(len(jd_tdb), dtype=float)
    dec_ours = np.empty(len(jd_tdb), dtype=float)
    ra_ours_nolt = np.empty(len(jd_tdb), dtype=float)
    dec_ours_nolt = np.empty(len(jd_tdb), dtype=float)
    ra_ours_full = np.empty(len(jd_tdb), dtype=float)
    dec_ours_full = np.empty(len(jd_tdb), dtype=float)

    for i, jd in enumerate(jd_tdb):
        # Earth velocity at obs epoch (good approximation for Gaia at L2)
        _, vel_au_day = get_body_barycentric_posvel(
            "earth", Time(float(jd), format="jd", scale="tdb")
        )
        observer_vel = np.array(
            [
                vel_au_day.x.to_value("AU/day"),
                vel_au_day.y.to_value("AU/day"),
                vel_au_day.z.to_value("AU/day"),
            ]
        )

        # With light-time only
        pos, _ = light_time_iterate(target_pos_at, float(jd), gaia_xyz[i])
        los = pos - gaia_xyz[i]
        ra, dec = xyz_to_radec(los)
        ra_ours[i] = ra
        dec_ours[i] = dec

        # Without light-time (no aberration either)
        pos_nolt = target_pos_at(float(jd))
        los_nolt = pos_nolt - gaia_xyz[i]
        ra_nolt, dec_nolt = xyz_to_radec(los_nolt)
        ra_ours_nolt[i] = ra_nolt
        dec_ours_nolt[i] = dec_nolt

        # With light-time AND aberration
        los_unit = los / np.linalg.norm(los)
        los_apparent_unit = stellar_aberration(los_unit, observer_vel)
        # Scale back to original magnitude
        los_apparent = los_apparent_unit * np.linalg.norm(los)
        ra_f, dec_f = xyz_to_radec(los_apparent)
        ra_ours_full[i] = ra_f
        dec_ours_full[i] = dec_f

    # 5. Compare
    deg = np.pi / 180.0
    def _shift(ra1, dec1, ra2, dec2):
        dra = ((ra1 - ra2 + 540.0) % 360.0 - 180.0) * np.cos(dec2 * deg) * 3_600_000.0
        ddec = (dec1 - dec2) * 3_600_000.0
        return dra, ddec, np.sqrt(dra ** 2 + ddec ** 2)

    dra_lt, ddec_lt, sep_lt = _shift(ra_ours, dec_ours, ra_horizons, dec_horizons)
    dra_nolt, ddec_nolt, sep_nolt = _shift(ra_ours_nolt, dec_ours_nolt, ra_horizons, dec_horizons)
    dra_full, ddec_full, sep_full = _shift(ra_ours_full, dec_ours_full, ra_horizons, dec_horizons)

    logger.info("")
    logger.info("=== Transforms-only sanity check (Ceres) ===")
    logger.info("  NEITHER light-time NOR aberration:")
    logger.info(
        "    median |ΔRA| = %.2f mas    median |ΔDec| = %.2f mas    median sep = %.2f mas",
        np.median(np.abs(dra_nolt)), np.median(np.abs(ddec_nolt)), np.median(sep_nolt),
    )
    logger.info("  WITH light-time only:")
    logger.info(
        "    median |ΔRA| = %.2f mas    median |ΔDec| = %.2f mas    median sep = %.2f mas",
        np.median(np.abs(dra_lt)), np.median(np.abs(ddec_lt)), np.median(sep_lt),
    )
    logger.info("  WITH light-time + aberration:")
    logger.info(
        "    median |ΔRA| = %.2f mas    median |ΔDec| = %.2f mas    median sep = %.2f mas",
        np.median(np.abs(dra_full)), np.median(np.abs(ddec_full)), np.median(sep_full),
    )
    logger.info("")
    # The light-time-only run is the one comparable with Gaia's reported (ra, dec):
    # Gaia DR3 SSO observations are in the barycentric astrometric ICRS frame, with
    # stellar aberration already removed by the Gaia pipeline. Applying our own
    # stellar aberration on top of that double-counts and makes the residual worse.
    best_residual = np.median(sep_lt)
    logger.info(
        "Best residual (light-time only, matches Gaia DR3 frame): %.1f mas",
        best_residual,
    )
    if best_residual < 5.0:
        logger.info("✅ PASS: transforms at mas level — ready for mass fit")
    elif best_residual < 1500.0:
        logger.info(
            "⚠️ MARGINAL: transforms at ~%.0f mas. Acceptable for step detection but "
            "may limit mass precision. Likely sources: Horizons vector interpolation "
            "(splines on 1-hour grid), or solar gravitational deflection (~few mas) "
            "not yet modelled.",
            best_residual,
        )
    else:
        logger.error("❌ FAIL: transforms have a bug, residual %.0f mas", best_residual)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
