"""Clean end-to-end deflection demo using JPL Horizons apparent ephemerides.

Best of both worlds: uses Horizons's ``ephemerides()`` API with Gaia as the
observer (location code ``500@-139479``). Horizons returns the *apparent*
RA/Dec from Gaia, fully corrected for:

  - light-time (the asteroid was where it was ~25 minutes earlier)
  - stellar aberration (Gaia's velocity adds ~20 arcsec annual modulation)
  - gravitational deflection by the Sun (~minor at 3 AU)

If the residuals (Gaia_obs − Horizons_pred) now hover at the few-mas level
and show a discontinuity at the encounter date, the (111) Ate perturbation
is detected (because Horizons does NOT include (111) Ate as a perturber).

Output: data/output/ate_2000nt3_clean.csv

Usage
-----
    docker compose run --rm pipeline python -m scripts.demo_ate_clean
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

_TARGET_NUMBER = 18105
_ENCOUNTER_DATE = "2016-06-08T00:00:00"
_HALF_WINDOW_DAYS = 180.0
_BLACKOUT_DAYS = 7.0
_J2010_TCB_JD = 2455197.5
_GAIA_OBSERVER = "500@-139479"  # Gaia spacecraft, MPC observatory code


def fetch_gaia(archive_url: str, target: int, d_min: float, d_max: float) -> pl.DataFrame:
    adql = (
        "SELECT number_mp, epoch, ra, dec, g_mag "
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


def horizons_apparent_radec(
    target: int, jd_tdb: np.ndarray, rate_limit_s: float
) -> tuple[np.ndarray, np.ndarray]:
    """Return (RA_apparent, Dec_apparent) from Gaia as predicted by Horizons."""
    ra_out = np.empty(len(jd_tdb), dtype=float)
    dec_out = np.empty(len(jd_tdb), dtype=float)
    chunk = 50
    for start in range(0, len(jd_tdb), chunk):
        idx = np.arange(start, min(start + chunk, len(jd_tdb)))
        epochs = jd_tdb[idx].tolist()
        logger.info(
            "Querying Horizons (apparent RA/Dec from Gaia) for epochs %d–%d…",
            start,
            start + len(epochs) - 1,
        )
        h = Horizons(id=str(target), location=_GAIA_OBSERVER, epochs=epochs, id_type="smallbody")
        eph = h.ephemerides()
        ra_out[idx] = np.array(eph["RA"], dtype=float)
        dec_out[idx] = np.array(eph["DEC"], dtype=float)
        time.sleep(rate_limit_s)
    return ra_out, dec_out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/output/ate_2000nt3_clean.csv"),
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

    # Horizons predictions: apparent RA/Dec from Gaia, fully corrected
    ra_pred, dec_pred = horizons_apparent_radec(_TARGET_NUMBER, jd_tdb, rate_limit)

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
            "med(sep)=%7.2f mas  σ(ΔRA)=%.2f  σ(ΔDec)=%.2f",
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
                "Δ(after − before) on %s = %+7.2f ± %5.2f mas   t = %+5.2fσ",
                axis,
                ma - mb,
                se,
                t,
            )

    logger.info("")
    logger.info("Expected (111) Ate signal: ~5 mas at encounter, growing post-encounter.")
    logger.info(
        "If residual scatter is at the mas level and Δ(after − before) > 3σ, signal is detected."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
