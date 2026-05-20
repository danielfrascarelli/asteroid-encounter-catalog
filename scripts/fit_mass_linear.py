"""Linear mass fit using Gaia−Horizons residuals over the full DR3 window.

Replaces the N-body orbit fit of fit_perturber_mass.py with a more robust
approach:

  1. Fetch ALL Gaia DR3 observations of the target (full mission window).
  2. Query JPL Horizons for the target's apparent RA/Dec at every epoch.
  3. Residuals R(t) = Gaia − Horizons contain:
       - A slowly-varying orbital drift (Horizons orbit ≠ true orbit)
       - The perturber signal (growing linearly after the encounter)
       - Measurement noise (~1–10 mas per transit)
  4. Model:
         R_RA(t)  = a + b·t + α · Δθ_RA(t)
         R_Dec(t) = c + d·t + α · Δθ_Dec(t)
     where Δθ(t) = N-body_prediction(M=M_est) − N-body_prediction(M=0)
     is the perturbation basis function (precomputed once), and α = M / M_est
     is the mass scale factor.
  5. Solve the (linear) least-squares system for [a, b, c, d, α].

Why this works better than fit_perturber_mass.py
-------------------------------------------------
* fit_perturber_mass used the N-body trajectory as its absolute orbit
  reference.  Our N-body misses unmodeled perturbers, giving ~97 mas
  systematic residuals even after orbit fitting.  With 9 pre-encounter
  obs and a 79 mas signal, Phase A can't converge.
* Here Horizons is the reference (few-mas accuracy).  The drift polynomial
  absorbs orbital-bias.  The N-body is only used to compute the DIFFERENTIAL
  perturbation Δθ, where forward-model systematics cancel.
* The signal grows linearly: at 895 days post-encounter Loreley's kick
  accumulates to ~40 mas — detectable above a detrended background.

Usage
-----
    docker compose run --rm pipeline python -m scripts.fit_mass_linear \\
        --perturber 165 --target 31067 --date 2014-12-08
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import time
from pathlib import Path

import numpy as np
import polars as pl
from astropy.time import Time
from astroquery.jplhorizons import Horizons
from astroquery.utils.tap.core import TapPlus

from src.astrometry.forward_model import forward_model
from src.ingest.mpcorb import parse_mpcorb
from src.utils.config import load_config

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)

_J2010_TCB_JD = 2455197.5
_GAIA_OBSERVER = "500@-139479"
_GAIA_START_JD_TCB = 2456863.5  # 2014-07-25 (Gaia mission start)
_GAIA_END_JD_TCB   = 2457910.5  # 2017-05-28 (Gaia DR3 end)
_AU_M = 1.495978707e11
_G = 6.674e-11


def fetch_gaia_full_window(archive_url: str, target: int) -> pl.DataFrame:
    """Fetch all Gaia DR3 observations of *target* in the mission window."""
    d_min = _GAIA_START_JD_TCB - _J2010_TCB_JD
    d_max = _GAIA_END_JD_TCB   - _J2010_TCB_JD
    adql = (
        "SELECT number_mp, epoch, ra, dec, g_mag, x_gaia, y_gaia, z_gaia "
        "FROM gaiadr3.sso_observation "
        f"WHERE number_mp = {target} "
        f"AND epoch BETWEEN {d_min:.6f} AND {d_max:.6f} "
        "ORDER BY epoch"
    )
    tap = TapPlus(url=archive_url)
    job = tap.launch_job_async(adql)
    df = pl.from_pandas(job.get_results().to_pandas())
    df = df.rename({c: c.lower() for c in df.columns if c != c.lower()})
    return df


def horizons_apparent_radec(
    target: int, jd_tdb: np.ndarray, rate_limit_s: float
) -> tuple[np.ndarray, np.ndarray]:
    """Query Horizons for apparent RA/Dec from Gaia at each epoch."""
    ra_out  = np.empty(len(jd_tdb), dtype=float)
    dec_out = np.empty(len(jd_tdb), dtype=float)
    chunk = 50
    for start in range(0, len(jd_tdb), chunk):
        idx = np.arange(start, min(start + chunk, len(jd_tdb)))
        h = Horizons(
            id=str(target), location=_GAIA_OBSERVER,
            epochs=jd_tdb[idx].tolist(), id_type="smallbody",
        )
        eph = h.ephemerides()
        ra_out[idx]  = np.array(eph["RA"],  dtype=float)
        dec_out[idx] = np.array(eph["DEC"], dtype=float)
        if start + chunk < len(jd_tdb):
            time.sleep(rate_limit_s)
    return ra_out, dec_out


def compute_perturbation_signal(
    target_elements: dict,
    perturber_elements: dict,
    m_est_kg: float,
    jd_tdb: np.ndarray,
    gaia_xyz: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute N-body perturbation basis Δθ_RA, Δθ_Dec (mas).

    Returns (Δθ_RA / M_est, Δθ_Dec / M_est) — normalized so that the fitted
    scale factor α = M_fitted / M_est directly gives the mass ratio.

    Systematics in the absolute trajectory cancel in the difference.
    """
    logger.info("Computing N-body perturbation basis (M=M_est)…")
    ra_with, dec_with = forward_model(
        target_elements, perturber_elements, m_est_kg, jd_tdb, gaia_xyz,
    )
    logger.info("Computing N-body perturbation basis (M=0)…")
    ra_zero, dec_zero = forward_model(
        target_elements, perturber_elements, 0.0, jd_tdb, gaia_xyz,
    )
    deg = np.pi / 180.0
    dra  = ((ra_with  - ra_zero  + 540.0) % 360.0 - 180.0) * np.cos(dec_zero * deg) * 3_600_000.0
    ddec = (dec_with - dec_zero) * 3_600_000.0
    # Normalized per M_est kg (so α = M/M_est)
    return dra, ddec


def fit_mass_linear(
    r_ra: np.ndarray,
    r_dec: np.ndarray,
    t_days: np.ndarray,
    sig_ra: np.ndarray,
    sig_dec: np.ndarray,
    dtheta_ra: np.ndarray,
    dtheta_dec: np.ndarray,
    enc_days: float,
) -> dict:
    """Solve the linear system for [a, b, c, d, α].

    Model:
        R_RA(t)  = a + b*(t - t_enc) + α * Δθ_RA(t)
        R_Dec(t) = c + d*(t - t_enc) + α * Δθ_Dec(t)

    Parameters
    ----------
    r_ra, r_dec  : residuals Gaia − Horizons in mas
    t_days       : time relative to mission start (days)
    sig_ra, sig_dec : per-obs noise estimate (mas). Set to 1 if unknown.
    dtheta_ra, dtheta_dec : perturbation basis functions (mas / M_est)
    enc_days     : encounter epoch in same t_days units

    Returns dict with mass scale factor α, its σ, and diagnostic columns.
    """
    dt = t_days - enc_days  # time from encounter

    n = len(t_days)
    # Build design matrix A (2n × 5): columns [a, b, c, d, α]
    # RA rows:  R_RA = a + b*dt + α * Δθ_RA
    # Dec rows: R_Dec = c + d*dt + α * Δθ_Dec
    A = np.zeros((2 * n, 5), dtype=float)
    A[:n, 0] = 1.0 / sig_ra           # a term (RA)
    A[:n, 1] = dt / sig_ra             # b term
    A[:n, 4] = dtheta_ra / sig_ra      # α term
    A[n:, 2] = 1.0 / sig_dec          # c term (Dec)
    A[n:, 3] = dt / sig_dec            # d term
    A[n:, 4] = dtheta_dec / sig_dec    # α term

    y = np.concatenate([r_ra / sig_ra, r_dec / sig_dec])

    # Solve via normal equations (equivalent to lstsq)
    result, residuals, rank, sv = np.linalg.lstsq(A, y, rcond=None)
    a_fit, b_fit, c_fit, d_fit, alpha_fit = result

    # Covariance matrix
    chi2 = float(np.sum((A @ result - y) ** 2))
    n_params = 5
    dof = 2 * n - n_params
    chi2_red = chi2 / max(1, dof)
    try:
        cov = np.linalg.inv(A.T @ A) * chi2_red
        sigma = np.sqrt(np.diag(cov))
    except np.linalg.LinAlgError:
        sigma = np.full(5, float("nan"))

    return {
        "a_mas": float(a_fit),
        "b_mas_per_day": float(b_fit),
        "c_mas": float(c_fit),
        "d_mas_per_day": float(d_fit),
        "alpha": float(alpha_fit),
        "alpha_sigma": float(sigma[4]),
        "chi2": float(chi2),
        "chi2_red": float(chi2_red),
        "dof": int(dof),
        "rank": int(rank),
    }


def load_element_row(snapshot: Path, number: int) -> dict:
    df = parse_mpcorb(str(snapshot), semimajor_min_au=0.0, semimajor_max_au=50.0)
    sub = df.filter(pl.col("number") == number)
    if sub.height == 0:
        raise ValueError(f"asteroid {number} not in {snapshot}")
    return sub.row(0, named=True)


def _mass_from_diameter(d_km: float, rho_kg_m3: float = 1500.0) -> float:
    if d_km is None or (isinstance(d_km, float) and math.isnan(d_km)) or d_km <= 0.0:
        return 1.0e18
    r = 0.5 * d_km * 1000.0
    return rho_kg_m3 * (4.0 / 3.0) * math.pi * r ** 3


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--config", default="config.yaml")
    p.add_argument("--perturber", type=int, required=True)
    p.add_argument("--target",    type=int, required=True)
    p.add_argument("--date",      required=True, help="Encounter date, ISO UTC")
    p.add_argument(
        "--mpcorb",
        type=Path,
        default=Path("data/raw/mpcorb_archive/MPCORB_20120918.DAT"),
    )
    p.add_argument(
        "--noise-mas", type=float, default=5.0,
        help="Assumed per-transit astrometric noise (mas). Used as uniform weight.",
    )
    p.add_argument(
        "--output", type=Path,
        default=None,
        help="Output JSON path. Defaults to data/output/fit_linear_<tag>.json",
    )
    args = p.parse_args()

    cfg = load_config(args.config)
    archive_url  = cfg.sources.gaia_sso.archive_url
    rate_limit   = float(cfg.sources.jpl_horizons.rate_limit_seconds)

    enc_jd_tdb = float(Time(args.date, scale="utc").tdb.jd)
    enc_jd_tcb = float(Time(args.date, scale="utc").tcb.jd)

    # Gaia observations (full DR3 window)
    logger.info("Fetching all Gaia DR3 observations of target %d…", args.target)
    obs = fetch_gaia_full_window(archive_url, args.target)
    if obs.height < 10:
        logger.error("Only %d Gaia transits — too few.", obs.height)
        return 1
    logger.info("  → %d transits over full DR3 window", obs.height)

    epochs_days = obs["epoch"].to_numpy()
    jd_tcb = epochs_days + _J2010_TCB_JD
    jd_tdb = Time(jd_tcb, format="jd", scale="tcb").tdb.jd.astype(float)
    gaia_xyz = np.column_stack([
        obs["x_gaia"].to_numpy(),
        obs["y_gaia"].to_numpy(),
        obs["z_gaia"].to_numpy(),
    ]).astype(float)
    ra_obs  = obs["ra"].to_numpy().astype(float)
    dec_obs = obs["dec"].to_numpy().astype(float)

    # Horizons predictions for full window
    logger.info("Querying Horizons for %d epochs…", len(jd_tdb))
    ra_pred, dec_pred = horizons_apparent_radec(args.target, jd_tdb, rate_limit)

    deg = np.pi / 180.0
    r_ra  = ((ra_obs  - ra_pred  + 540.0) % 360.0 - 180.0) * np.cos(dec_pred * deg) * 3_600_000.0
    r_dec = (dec_obs - dec_pred) * 3_600_000.0

    # Time axis (days from encounter, JD TDB scale)
    t_days = jd_tdb - enc_jd_tdb

    # Count pre/post
    n_pre  = int((t_days < -7.0).sum())
    n_post = int((t_days >  7.0).sum())
    logger.info("  %d pre-encounter + %d post-encounter observations", n_pre, n_post)

    # Load orbital elements
    logger.info("Loading MPCORB elements…")
    target_el    = load_element_row(args.mpcorb, args.target)
    perturber_el = load_element_row(args.mpcorb, args.perturber)

    h = perturber_el.get("H", None)
    if h is not None:
        d_km = (1329.0 / math.sqrt(0.14)) * 10.0 ** (-h / 5.0)
    else:
        d_km = 100.0
    m_est = _mass_from_diameter(d_km)
    logger.info("Perturber H=%.2f → D≈%.0f km → M_est=%.2e kg", h or 0.0, d_km, m_est)

    # N-body perturbation basis (only at observation epochs to save time)
    logger.info("Computing N-body perturbation signal (2 integrations)…")
    dtheta_ra, dtheta_dec = compute_perturbation_signal(
        target_el, perturber_el, m_est, jd_tdb, gaia_xyz,
    )
    max_signal = float(np.max(np.abs(dtheta_ra)))
    logger.info("  Max |Δθ_RA| over full window = %.1f mas  (at M=M_est=%.2e kg)",
                max_signal, m_est)

    # Uniform noise weights
    sig_ra  = np.full(len(jd_tdb), args.noise_mas)
    sig_dec = np.full(len(jd_tdb), args.noise_mas)

    # Linear fit
    logger.info("Running linear fit…")
    fit = fit_mass_linear(
        r_ra, r_dec, t_days, sig_ra, sig_dec,
        dtheta_ra, dtheta_dec,
        enc_days=0.0,  # t_days already relative to encounter
    )

    alpha     = fit["alpha"]
    alpha_sig = fit["alpha_sigma"]
    mass_kg   = alpha * m_est
    mass_sig  = alpha_sig * m_est

    logger.info("")
    logger.info("=== LINEAR FIT RESULT ===")
    logger.info("  α = %.4f ± %.4f  (M/M_est)", alpha, alpha_sig)
    logger.info("  Fitted mass = %.3e ± %.2e kg", mass_kg, mass_sig)
    logger.info("  chi²_red = %.2f  (dof = %d)", fit["chi2_red"], fit["dof"])
    logger.info("  Drift RA:  %.1f + %.4f × t mas/day", fit["a_mas"], fit["b_mas_per_day"])
    logger.info("  Drift Dec: %.1f + %.4f × t mas/day", fit["c_mas"], fit["d_mas_per_day"])
    logger.info("  Pre-enc obs: %d, Post-enc obs: %d", n_pre, n_post)

    # Save
    tag = f"{args.perturber:06d}_{args.target:06d}_linear"
    out_path = args.output or Path("data/output") / f"fit_{tag}.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    summary = {
        "method": "linear_horizons_detrend",
        "perturber": args.perturber,
        "target": args.target,
        "encounter_date": args.date,
        "m_est_kg": m_est,
        "fitted_mass_kg": mass_kg,
        "fitted_mass_sigma_kg": mass_sig,
        "alpha": alpha,
        "alpha_sigma": alpha_sig,
        "n_obs_total": obs.height,
        "n_pre": n_pre,
        "n_post": n_post,
        "max_signal_mas": max_signal,
        **fit,
    }
    with out_path.open("w") as fh:
        json.dump(summary, fh, indent=2)
    logger.info("Wrote → %s", out_path)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
