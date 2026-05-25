"""Detect (111) Ate's perturbation on 2000_nt3 — Gaia residuals vs JPL Horizons.

Improved version of ``demo_ate_deflection.py``: instead of propagating MPCORB
elements with Kepler 2-body (which accumulates ~875 arcsec of error over 4
years and swamps the ~5 mas Ate signal), we query JPL Horizons for the
target's barycentric position at each Gaia observation epoch.

JPL Horizons uses DE440 + N-body integration with the major planets and the
"big 4" asteroids (Ceres, Pallas, Vesta, Hygiea), but NOT (111) Ate.
Therefore the residuals (Gaia_obs − Horizons_pred) should contain Ate's
gravitational perturbation as a coherent signal — most clearly as a
systematic offset that grows linearly after the 2016-06-08 encounter.

Output: data/output/ate_2000nt3_vs_horizons.csv

Usage
-----
    docker compose run --rm pipeline python -m scripts.demo_ate_vs_horizons
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

from src.utils.config import load_config

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)

_TARGET_NUMBER = 18105  # MPC number for 2000 NT3
_ENCOUNTER_DATE = "2016-06-08T00:00:00"
_HALF_WINDOW_DAYS = 180.0
_BLACKOUT_DAYS = 7.0
_J2010_TCB_JD = 2455197.5
_RAD_TO_MAS = 180.0 / np.pi * 3_600_000.0


def fetch_gaia(archive_url: str, target: int, d_min: float, d_max: float) -> pl.DataFrame:
    """Pull Gaia observations of *target* in the days-since-J2010-TCB window."""
    adql = (
        "SELECT number_mp, epoch, ra, dec, g_mag, x_gaia, y_gaia, z_gaia "
        "FROM gaiadr3.sso_observation "
        f"WHERE number_mp = {target} "
        f"AND epoch BETWEEN {d_min:.6f} AND {d_max:.6f} "
        "ORDER BY epoch"
    )
    tap = TapPlus(url=archive_url)
    job = tap.launch_job_async(adql)
    tbl = job.get_results()
    df = pl.from_pandas(tbl.to_pandas())
    df = df.rename({c: c.lower() for c in df.columns if c != c.lower()})
    logger.info("Gaia returned %d transits", df.height)
    return df


def horizons_barycentric_xyz(target: int, jd_tdb: np.ndarray, rate_limit_s: float) -> np.ndarray:
    """Query JPL Horizons for the asteroid's barycentric ICRS position at each epoch.

    Returns an (N, 3) array in AU.  Splits long requests into chunks of 100
    epochs to stay under the Horizons per-call limit.
    """
    out = np.empty((len(jd_tdb), 3), dtype=float)
    chunk = 100
    for start in range(0, len(jd_tdb), chunk):
        idx = np.arange(start, min(start + chunk, len(jd_tdb)))
        epochs = jd_tdb[idx].tolist()
        logger.info(
            "Querying Horizons for %d epochs (%d–%d)…",
            len(epochs),
            start,
            start + len(epochs) - 1,
        )
        # location='@0' = Solar System Barycenter
        h = Horizons(id=str(target), location="@0", epochs=epochs, id_type="smallbody")
        vec = h.vectors(refplane="earth")  # equatorial (ICRS) so it matches x_gaia
        out[idx, 0] = np.array(vec["x"], dtype=float)
        out[idx, 1] = np.array(vec["y"], dtype=float)
        out[idx, 2] = np.array(vec["z"], dtype=float)
        time.sleep(rate_limit_s)
    return out


def xyz_to_radec(vec: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Cartesian (..., 3) → (RA_deg, Dec_deg)."""
    x, y, z = vec[..., 0], vec[..., 1], vec[..., 2]
    rho = np.sqrt(x * x + y * y)
    ra = np.degrees(np.arctan2(y, x)) % 360.0
    dec = np.degrees(np.arctan2(z, rho))
    return ra, dec


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/output/ate_2000nt3_vs_horizons.csv"),
    )
    args = parser.parse_args()

    cfg = load_config(args.config)
    archive_url = cfg.sources.gaia_sso.archive_url
    rate_limit = float(cfg.sources.jpl_horizons.rate_limit_seconds)

    enc_jd_tdb = float(Time(_ENCOUNTER_DATE, scale="utc").tdb.jd)
    enc_days = float(Time(_ENCOUNTER_DATE, scale="utc").tcb.jd) - _J2010_TCB_JD

    obs = fetch_gaia(
        archive_url,
        _TARGET_NUMBER,
        enc_days - _HALF_WINDOW_DAYS,
        enc_days + _HALF_WINDOW_DAYS,
    )
    if obs.height == 0:
        logger.error("No Gaia observations found")
        return 1

    epochs_days = obs["epoch"].to_numpy()
    jd_tcb = epochs_days + _J2010_TCB_JD
    jd_tdb = Time(jd_tcb, format="jd", scale="tcb").tdb.jd.astype(float)

    gaia_xyz = np.column_stack(
        [obs["x_gaia"].to_numpy(), obs["y_gaia"].to_numpy(), obs["z_gaia"].to_numpy()]
    ).astype(float)

    # --- JPL Horizons predictions ---------------------------------------------
    target_bary = horizons_barycentric_xyz(_TARGET_NUMBER, jd_tdb, rate_limit)

    # Line of sight from Gaia to target
    los = target_bary - gaia_xyz
    ra_pred, dec_pred = xyz_to_radec(los)

    ra_obs = obs["ra"].to_numpy().astype(float)
    dec_obs = obs["dec"].to_numpy().astype(float)

    deg = np.pi / 180.0
    dra_mas = ((ra_obs - ra_pred + 540.0) % 360.0 - 180.0) * np.cos(dec_pred * deg) * 3_600_000.0
    ddec_mas = (dec_obs - dec_pred) * 3_600_000.0
    sep_mas = np.sqrt(dra_mas**2 + ddec_mas**2)

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

    # --- Statistics --------------------------------------------------------
    before = result.filter(pl.col("days_from_encounter") < -_BLACKOUT_DAYS)
    after = result.filter(pl.col("days_from_encounter") > _BLACKOUT_DAYS)

    def _stats(df: pl.DataFrame, label: str) -> dict[str, float]:
        if df.height == 0:
            return {}
        s = {
            "n": df.height,
            "mu_dra": float(df["dra_mas"].mean()),
            "mu_ddec": float(df["ddec_mas"].mean()),
            "std_dra": float(df["dra_mas"].std()),
            "std_ddec": float(df["ddec_mas"].std()),
            "median_sep": float(df["sep_mas"].median()),
        }
        logger.info(
            "%s (N=%d)  ⟨ΔRA⟩=%+8.2f mas  ⟨ΔDec⟩=%+8.2f mas  "
            "med(sep)=%8.2f mas  σ(ΔRA)=%.2f  σ(ΔDec)=%.2f",
            label,
            s["n"],
            s["mu_dra"],
            s["mu_ddec"],
            s["median_sep"],
            s["std_dra"],
            s["std_ddec"],
        )
        return s

    sb = _stats(before, "BEFORE encounter")
    sa = _stats(after, "AFTER  encounter")

    if sb and sa:
        for axis, key in [("ΔRA", "dra"), ("ΔDec", "ddec")]:
            mb = sb[f"mu_{key}"]
            ma = sa[f"mu_{key}"]
            se = np.sqrt((sb[f"std_{key}"] ** 2) / sb["n"] + (sa[f"std_{key}"] ** 2) / sa["n"])
            t = (ma - mb) / se if se > 0 else float("inf")
            logger.info(
                "Δ(after − before) on %s = %+8.2f ± %5.2f mas   t = %+5.2fσ",
                axis,
                ma - mb,
                se,
                t,
            )

    logger.info("")
    logger.info("Expected (111) Ate signal: ~5 mas at encounter, growing linearly.")
    logger.info("If post − pre offset is significant (|t| > 3) on either axis, signal detected.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
