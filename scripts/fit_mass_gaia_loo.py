"""Leave-One-Out mass fit: Gaia-calibrated orbit + encounter-window residuals.

The fundamental improvement over fit_perturber_mass.py and fit_mass_linear.py:

  Problem with Horizons approach (fit_mass_linear.py):
    Horizons orbits for small MBAs have ~100-600 mas systematic drift relative
    to Gaia observations. The perturbation signal (5-40 mas) is buried 10-100×.

  Solution — LOO orbit fit:
    1. Fetch ALL Gaia DR3 observations of the target.
    2. Phase A (LOO orbit): fit 6 orbital elements using observations OUTSIDE
       the ±LOO_WINDOW_DAYS window around the encounter. With hundreds of
       observations and an N-body model including all 8 planets, the orbit is
       determined to Gaia's own precision (~few mas).
    3. Phase B (mass fit): using the LOO orbit, compute PREDICTED positions
       inside the encounter window. Residuals = Gaia - prediction. Fit perturber
       mass to explain the post-encounter excess above the pre-encounter baseline.

  Why this works:
    - Gaia calibrates its own orbit with the same data it will test. No Horizons
      drift. Residuals in the encounter window are at the Gaia noise level.
    - The LOO split ensures the orbit is NOT biased by any perturbation signal:
      pre-encounter data (before the velocity kick) fits a "clean" orbit; the
      post-encounter deviation is then purely the perturbation.
    - With a mid-mission encounter (2015-2016), hundreds of pre+post observations
      outside the LOO window give a tightly constrained baseline.

  Best candidates for this approach:
    - Encounters in the MIDDLE of the Gaia mission (2015-2016) so that both pre
      and post-encounter non-window observations are plentiful.
    - Targets with many Gaia transits (bright, well-tracked MBAs).
    - Loreley (2014-12-08) is the HARDEST case: only 5 months pre-encounter in
      the Gaia mission, giving few non-window observations. Still useful as a
      validation that the machinery is correct.

Usage
-----
    docker compose run --rm pipeline python -m scripts.fit_mass_gaia_loo \\
        --perturber 165 --target 31067 --date 2014-12-08
    docker compose run --rm pipeline python -m scripts.fit_mass_gaia_loo \\
        --perturber 57 --target 216875 --date 2016-08-26
"""

from __future__ import annotations

import argparse
import json
import logging
import math
from pathlib import Path

import numpy as np
import polars as pl
from astropy.time import Time
from astroquery.utils.tap.core import TapPlus
from scipy.optimize import least_squares

from src.astrometry.forward_model import forward_model, residuals_mas
from src.ingest.mpcorb import parse_mpcorb
from src.utils.config import load_config

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)

_J2010_TCB_JD = 2455197.5
_GAIA_START_JD_TCB = 2456863.5  # 2014-07-25
_GAIA_END_JD_TCB   = 2457910.5  # 2017-05-28


def fetch_gaia_full(archive_url: str, target: int) -> pl.DataFrame:
    d_min = _GAIA_START_JD_TCB - _J2010_TCB_JD
    d_max = _GAIA_END_JD_TCB   - _J2010_TCB_JD
    adql = (
        "SELECT number_mp, epoch, ra, dec, "
        "ra_error_systematic, dec_error_systematic, ra_dec_correlation_systematic, "
        "ra_error_random, dec_error_random, ra_dec_correlation_random, "
        "position_angle_scan, x_gaia, y_gaia, z_gaia "
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


def al_residuals_and_weights(
    dra_mas: np.ndarray,
    ddec_mas: np.ndarray,
    pa_scan_deg: np.ndarray,
    ra_err_sys: np.ndarray,
    dec_err_sys: np.ndarray,
    corr_sys: np.ndarray,
    ra_err_rand: np.ndarray,
    dec_err_rand: np.ndarray,
    corr_rand: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Project residuals onto the Gaia along-scan direction and compute σ_AL.

    Gaia SSO measures positions primarily along the scan direction (AL) with
    precision ~0.2-2 mas. The across-scan (AC) direction has ~300-600 mas
    uncertainty.  Fitting with unweighted (dRA² + dDec²) means each observation
    is dominated by AC noise; fitting with AL-projected, σ_AL-weighted residuals
    recovers the actual Gaia precision.

    Parameters
    ----------
    dra_mas, ddec_mas:
        (N,) residuals in mas (obs − pred).
    pa_scan_deg:
        (N,) scan position angle in deg, measured North toward East.
    ra_err_sys, dec_err_sys, corr_sys:
        (N,) systematic uncertainty components in mas and their correlation.
    ra_err_rand, dec_err_rand, corr_rand:
        (N,) random (photon noise) uncertainty components and their correlation.

    Returns
    -------
    r_al : (N,) along-scan residuals in mas
    sigma_al : (N,) per-transit AL uncertainties in mas
    """
    pa = np.radians(pa_scan_deg)
    # AL unit vector: (RA component, Dec component)
    e_al_ra  =  np.sin(pa)   # East-positive
    e_al_dec =  np.cos(pa)   # North-positive

    r_al = dra_mas * e_al_ra + ddec_mas * e_al_dec

    # σ_AL from the covariance ellipse projected onto scan direction.
    # Cov = Cov_sys + Cov_rand; σ²_AL = ê_AL^T Cov ê_AL.
    def _cov_al(s_ra, s_dec, rho):
        return (e_al_ra ** 2 * s_ra ** 2
                + 2.0 * e_al_ra * e_al_dec * rho * s_ra * s_dec
                + e_al_dec ** 2 * s_dec ** 2)

    var_al = (_cov_al(ra_err_sys, dec_err_sys, corr_sys)
              + _cov_al(ra_err_rand, dec_err_rand, corr_rand))
    sigma_al = np.sqrt(np.maximum(var_al, 1e-6))  # floor at 0.001 mas
    return r_al, sigma_al


def load_element_rows(snapshot: Path, numbers: list[int]) -> dict[int, dict]:
    """Parse MPCORB once and return element dicts for all requested numbers."""
    df = parse_mpcorb(str(snapshot), semimajor_min_au=0.0, semimajor_max_au=50.0)
    result: dict[int, dict] = {}
    for n in numbers:
        sub = df.filter(pl.col("number") == n)
        if sub.height == 0:
            raise ValueError(f"asteroid {n} not in {snapshot}")
        result[n] = sub.row(0, named=True)
    return result


def _mass_from_h(h: float | None) -> float:
    if h is None:
        return 1.0e18
    d_km = (1329.0 / math.sqrt(0.14)) * 10.0 ** (-h / 5.0)
    r = 0.5 * d_km * 1000.0
    return 1500.0 * (4.0 / 3.0) * math.pi * r ** 3


def fit_orbit_loo(
    target_elements: dict,
    perturber_elements: dict,
    obs_jd_tdb: np.ndarray,
    gaia_xyz: np.ndarray,
    ra_obs: np.ndarray,
    dec_obs: np.ndarray,
    pa_scan: np.ndarray,
    ra_err_sys: np.ndarray,
    dec_err_sys: np.ndarray,
    corr_sys: np.ndarray,
    ra_err_rand: np.ndarray,
    dec_err_rand: np.ndarray,
    corr_rand: np.ndarray,
    big4_elements: dict[str, dict] | None = None,
    reg_sigma: np.ndarray | None = None,
    max_nfev: int = 800,
    dt_days: float = 1.0,
    integrator: str = "whfast",
) -> tuple[dict, dict]:
    """Fit target orbit using observations OUTSIDE the encounter window.

    Uses AL-projected, σ_AL-weighted residuals (Gaia's actual precision is
    along-scan; across-scan is ~600 mas noise that must be ignored).

    Returns (fitted_elements_dict, diagnostics_dict).
    """
    a0 = target_elements["a_au"]
    e0 = target_elements["e"]
    i0 = target_elements["i_deg"]
    O0 = target_elements["Omega_deg"]
    o0 = target_elements["omega_deg"]
    M0 = target_elements["M_deg"]

    x0 = np.array([a0, e0, i0, O0, o0, M0])
    lo = np.array([a0 - 0.05, max(0.0, e0 - 0.05), max(0.0, i0 - 2.0),
                   O0 - 10.0, o0 - 10.0, M0 - 15.0])
    hi = np.array([a0 + 0.05, min(0.999, e0 + 0.05), i0 + 2.0,
                   O0 + 10.0, o0 + 10.0, M0 + 15.0])

    if reg_sigma is None:
        reg_sigma = np.array([5e-4, 5e-4, 0.1, 0.3, 0.3, 0.5])

    use_big4 = bool(big4_elements)

    def residual_func(params: np.ndarray) -> np.ndarray:
        tgt = dict(target_elements)
        tgt["a_au"]      = float(params[0])
        tgt["e"]         = float(params[1])
        tgt["i_deg"]     = float(params[2])
        tgt["Omega_deg"] = float(params[3])
        tgt["omega_deg"] = float(params[4])
        tgt["M_deg"]     = float(params[5])
        ra_pred, dec_pred = forward_model(
            tgt, perturber_elements, 0.0, obs_jd_tdb, gaia_xyz,
            include_big4=use_big4, big4_elements=big4_elements,
            dt_days=dt_days, integrator=integrator,
        )
        dra, ddec = residuals_mas(ra_obs, dec_obs, ra_pred, dec_pred)
        r_al, sigma_al = al_residuals_and_weights(
            dra, ddec, pa_scan,
            ra_err_sys, dec_err_sys, corr_sys,
            ra_err_rand, dec_err_rand, corr_rand,
        )
        reg = (params - x0) / reg_sigma
        return np.concatenate([r_al / sigma_al, reg])

    n_obs = len(obs_jd_tdb)
    logger.info(
        "[LOO orbit] fitting 6 params using %d out-of-window obs (AL-weighted)",
        n_obs,
    )
    res = least_squares(
        residual_func, x0, method="trf", bounds=(lo, hi),
        max_nfev=max_nfev, ftol=1e-10, xtol=1e-10, gtol=1e-10,
    )
    logger.info("[LOO orbit] done: nfev=%d  status=%d  message=%s",
                res.nfev, res.status, res.message)

    fitted = dict(target_elements)
    for key, idx in [("a_au", 0), ("e", 1), ("i_deg", 2),
                     ("Omega_deg", 3), ("omega_deg", 4), ("M_deg", 5)]:
        fitted[key] = float(res.x[idx])

    # Compute final AL residuals and σ_AL for diagnostics
    ra_pred_f, dec_pred_f = forward_model(
        fitted, perturber_elements, 0.0, obs_jd_tdb, gaia_xyz,
        include_big4=use_big4, big4_elements=big4_elements,
        dt_days=dt_days, integrator=integrator,
    )
    dra_f, ddec_f = residuals_mas(ra_obs, dec_obs, ra_pred_f, dec_pred_f)
    r_al_f, sigma_al_f = al_residuals_and_weights(
        dra_f, ddec_f, pa_scan,
        ra_err_sys, dec_err_sys, corr_sys,
        ra_err_rand, dec_err_rand, corr_rand,
    )
    # Normalised AL chi²_red
    chi2_norm = float(np.sum((r_al_f / sigma_al_f) ** 2))
    chi2_red  = chi2_norm / max(1, n_obs - 6)
    rms_al    = float(np.sqrt(np.mean(r_al_f ** 2)))
    rms_2d    = float(np.sqrt(np.mean(dra_f ** 2 + ddec_f ** 2)))
    return fitted, {
        "chi2_red": chi2_red,
        "rms_al_mas": rms_al,
        "rms_2d_mas": rms_2d,
        "nfev": res.nfev,
        "status": res.status,
        "n_obs": n_obs,
    }


def fit_mass_from_window(
    fitted_elements: dict,
    perturber_elements: dict,
    obs_jd_tdb: np.ndarray,
    gaia_xyz: np.ndarray,
    ra_obs: np.ndarray,
    dec_obs: np.ndarray,
    enc_jd_tdb: float,
    initial_mass_kg: float,
    pa_scan: np.ndarray,
    ra_err_sys: np.ndarray,
    dec_err_sys: np.ndarray,
    corr_sys: np.ndarray,
    ra_err_rand: np.ndarray,
    dec_err_rand: np.ndarray,
    corr_rand: np.ndarray,
    blackout_days: float = 7.0,
    big4_elements: dict[str, dict] | None = None,
    dt_days: float = 1.0,
    integrator: str = "whfast",
) -> dict:
    """Fit perturber mass to encounter-window residuals.

    The orbit (fitted_elements) comes from LOO Phase A. Residuals in the
    encounter window represent the perturbation signal, which is fit here.

    Splits into pre/post windows:
    - Pre-encounter residuals calibrate any remaining orbit drift (as a
      per-axis constant offset absorbed into a mean shift).
    - Post-encounter residuals contain the growing perturbation signal.
    """
    days = obs_jd_tdb - enc_jd_tdb
    pre_mask  = days < -blackout_days
    post_mask = days > blackout_days
    n_pre  = int(pre_mask.sum())
    n_post = int(post_mask.sum())

    logger.info("[mass fit] encounter window: %d pre + %d post obs", n_pre, n_post)

    use_big4 = bool(big4_elements)

    # Compute residuals under zero-mass prediction (orbit only)
    ra_pred0, dec_pred0 = forward_model(
        fitted_elements, perturber_elements, 0.0, obs_jd_tdb, gaia_xyz,
        include_big4=use_big4, big4_elements=big4_elements,
        dt_days=dt_days, integrator=integrator,
    )
    dra0, ddec0 = residuals_mas(ra_obs, dec_obs, ra_pred0, dec_pred0)
    r_al0, sig_al0 = al_residuals_and_weights(
        dra0, ddec0, pa_scan,
        ra_err_sys, dec_err_sys, corr_sys,
        ra_err_rand, dec_err_rand, corr_rand,
    )

    rms_pre_al  = float(np.sqrt(np.mean(r_al0[pre_mask] ** 2)))  if n_pre > 0 else float("nan")
    rms_post_al = float(np.sqrt(np.mean(r_al0[post_mask] ** 2))) if n_post > 0 else float("nan")
    shift_al = float(np.mean(r_al0[post_mask]) - np.mean(r_al0[pre_mask])) if (n_pre > 0 and n_post > 0) else float("nan")
    # Also keep 2D RMS for reference
    rms_pre_2d  = float(np.sqrt(np.mean(dra0[pre_mask]**2 + ddec0[pre_mask]**2))) if n_pre > 0 else float("nan")
    rms_post_2d = float(np.sqrt(np.mean(dra0[post_mask]**2 + ddec0[post_mask]**2))) if n_post > 0 else float("nan")
    logger.info(
        "[mass fit] zero-mass AL residuals:  pre=%.2f mas  post=%.2f mas  shift=%.2f mas",
        rms_pre_al, rms_post_al, shift_al,
    )
    logger.info(
        "[mass fit] zero-mass 2D RMS:  pre=%.1f mas  post=%.1f mas",
        rms_pre_2d, rms_post_2d,
    )

    if n_post < 3:
        return {
            "mass_kg": float("nan"), "mass_sigma_kg": float("nan"),
            "mass_sigma_inflated_kg": float("nan"),
            "log10_mass": float("nan"), "log10_mass_sigma": float("nan"),
            "log10_mass_sigma_inflated": float("nan"),
            "chi2_red": float("nan"), "n_pre": n_pre, "n_post": n_post,
            "rms_pre_al_mas": rms_pre_al, "rms_post_al_mas": rms_post_al,
            "shift_al_mas": shift_al,
            "rms_pre_2d_mas": rms_pre_2d, "rms_post_2d_mas": rms_post_2d,
            "note": "too few post-encounter observations",
        }

    log_m0 = math.log10(initial_mass_kg)

    def residual_func(params: np.ndarray) -> np.ndarray:
        mass = float(10.0 ** params[0])
        ra_p, dec_p = forward_model(
            fitted_elements, perturber_elements, mass, obs_jd_tdb[post_mask], gaia_xyz[post_mask],
            include_big4=use_big4, big4_elements=big4_elements,
            dt_days=dt_days, integrator=integrator,
        )
        dra, ddec = residuals_mas(
            ra_obs[post_mask], dec_obs[post_mask], ra_p, dec_p,
        )
        r_al, s_al = al_residuals_and_weights(
            dra, ddec, pa_scan[post_mask],
            ra_err_sys[post_mask], dec_err_sys[post_mask], corr_sys[post_mask],
            ra_err_rand[post_mask], dec_err_rand[post_mask], corr_rand[post_mask],
        )
        # Subtract pre-encounter mean AL drift to remove baseline offset
        if n_pre > 0:
            r_al = r_al - float(np.mean(r_al0[pre_mask]))
        return r_al / s_al

    res = least_squares(
        residual_func,
        np.array([log_m0]),
        method="trf",
        bounds=([14.0], [23.0]),
        max_nfev=500,
        ftol=1e-10, xtol=1e-10, gtol=1e-10,
        verbose=2,
    )

    log_m_fit = float(res.x[0])
    # res.fun contains normalised AL residuals (r_al / sigma_al)
    chi2_norm = float(np.sum(res.fun ** 2))
    n_data = len(res.fun)
    chi2_red = chi2_norm / max(1, n_data - 1)
    try:
        cov = np.linalg.inv(res.jac.T @ res.jac) * chi2_red
        log_m_sig = float(np.sqrt(cov[0, 0]))
    except np.linalg.LinAlgError:
        log_m_sig = float("nan")
    mass_fit = 10.0 ** log_m_fit
    mass_sig  = mass_fit * log_m_sig * math.log(10.0)
    # Inflate formal uncertainty by sqrt(chi2_red) to account for residual systematics
    sqrt_chi2 = math.sqrt(chi2_red) if math.isfinite(chi2_red) and chi2_red > 0 else float("nan")
    mass_sig_inflated   = mass_sig   * sqrt_chi2
    log_m_sig_inflated  = log_m_sig  * sqrt_chi2

    return {
        "mass_kg": mass_fit,
        "mass_sigma_kg": mass_sig,
        "mass_sigma_inflated_kg": mass_sig_inflated,
        "log10_mass": log_m_fit,
        "log10_mass_sigma": log_m_sig,
        "log10_mass_sigma_inflated": log_m_sig_inflated,
        "chi2_red": chi2_red,
        "n_pre": n_pre,
        "n_post": n_post,
        "rms_pre_al_mas": rms_pre_al,
        "rms_post_al_mas": rms_post_al,
        "shift_al_mas": shift_al,
        "rms_pre_2d_mas": rms_pre_2d,
        "rms_post_2d_mas": rms_post_2d,
        "note": "",
    }


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--config",    default="config.yaml")
    p.add_argument("--perturber", type=int, required=True)
    p.add_argument("--target",    type=int, required=True)
    p.add_argument("--date",      required=True, help="Encounter date, ISO UTC")
    p.add_argument(
        "--loo-window-days", type=float, default=180.0,
        help="Half-width of encounter window excluded from orbit fit (days).",
    )
    p.add_argument(
        "--blackout-days", type=float, default=7.0,
        help="Inner blackout around encounter (days).",
    )
    p.add_argument(
        "--mpcorb", type=Path,
        default=Path("data/raw/mpcorb_archive/MPCORB_20120918.DAT"),
    )
    p.add_argument(
        "--dt-days", type=float, default=1.0,
        help="WHFast timestep (days). For encounters with b < 0.001 AU use 0.05; "
             "or use --integrator ias15.",
    )
    p.add_argument(
        "--integrator", default="whfast",
        choices=["whfast", "ias15"],
        help="REBOUND integrator. Use ias15 for sub-day close encounters.",
    )
    p.add_argument("--output", type=Path, default=None)
    args = p.parse_args()

    cfg = load_config(args.config)
    archive_url = cfg.sources.gaia_sso.archive_url

    enc_jd_tdb = float(Time(args.date, scale="utc").tdb.jd)

    # ── Fetch all Gaia observations ──────────────────────────────────────────
    logger.info("Fetching all Gaia DR3 observations of target %d…", args.target)
    obs = fetch_gaia_full(archive_url, args.target)
    logger.info("  → %d transits over full DR3 window", obs.height)
    if obs.height < 15:
        logger.error("Too few transits (%d) for LOO fit.", obs.height)
        return 1

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
    pa_scan    = obs["position_angle_scan"].to_numpy().astype(float)
    ra_err_sys = obs["ra_error_systematic"].to_numpy().astype(float)
    dec_err_sys = obs["dec_error_systematic"].to_numpy().astype(float)
    corr_sys    = obs["ra_dec_correlation_systematic"].to_numpy().astype(float)
    ra_err_rand = obs["ra_error_random"].to_numpy().astype(float)
    dec_err_rand = obs["dec_error_random"].to_numpy().astype(float)
    corr_rand    = obs["ra_dec_correlation_random"].to_numpy().astype(float)

    days_from_enc = jd_tdb - enc_jd_tdb

    # ── Split: pre-encounter orbit baseline vs mass-fit window ──────────────
    # The orbit fit uses ONLY pre-encounter observations that are outside the
    # encounter window. Post-encounter data is perturbed and must NOT be used
    # for the orbit baseline — that would absorb the perturbation signal.
    # The mass fit window covers everything from (encounter - loo_window_days)
    # onwards (encounter window + all post-encounter data).
    out_mask = days_from_enc < -args.loo_window_days   # pre-enc, outside window
    win_mask = days_from_enc > -args.loo_window_days   # encounter + post-enc
    n_out = int(out_mask.sum())
    n_win = int(win_mask.sum())
    logger.info(
        "LOO split: %d pre-encounter baseline obs (orbit fit) | %d obs in mass window",
        n_out, n_win,
    )

    if n_out < 8:
        logger.error(
            "Only %d out-of-window observations — need ≥ 8 to fit 6 orbital params.", n_out,
        )
        return 1

    # ── Load orbital elements ────────────────────────────────────────────────
    logger.info("Loading MPCORB elements…")
    _BIG4_NUMBERS = {1: "ceres", 2: "pallas", 4: "vesta", 10: "hygiea"}
    big4_nums_to_load = [n for n in _BIG4_NUMBERS if n != args.perturber]
    all_nums = [args.target, args.perturber] + big4_nums_to_load
    elements_map = load_element_rows(args.mpcorb, all_nums)
    target_el    = elements_map[args.target]
    perturber_el = elements_map[args.perturber]
    big4_elements = {
        _BIG4_NUMBERS[n]: elements_map[n]
        for n in big4_nums_to_load
    }
    initial_mass = _mass_from_h(perturber_el.get("H"))
    logger.info(
        "Target %d: a=%.4f e=%.4f  Perturber %d: H=%.2f → M_est=%.2e kg",
        args.target, target_el["a_au"], target_el["e"],
        args.perturber, perturber_el.get("H", float("nan")), initial_mass,
    )
    logger.info("Big-4 loaded: %s", list(big4_elements.keys()))

    # ── Phase A: LOO orbit fit ───────────────────────────────────────────────
    fitted_elements, loo_diag = fit_orbit_loo(
        target_elements=target_el,
        perturber_elements=perturber_el,
        obs_jd_tdb=jd_tdb[out_mask],
        gaia_xyz=gaia_xyz[out_mask],
        ra_obs=ra_obs[out_mask],
        dec_obs=dec_obs[out_mask],
        pa_scan=pa_scan[out_mask],
        ra_err_sys=ra_err_sys[out_mask],
        dec_err_sys=dec_err_sys[out_mask],
        corr_sys=corr_sys[out_mask],
        ra_err_rand=ra_err_rand[out_mask],
        dec_err_rand=dec_err_rand[out_mask],
        corr_rand=corr_rand[out_mask],
        big4_elements=big4_elements,
        dt_days=args.dt_days,
        integrator=args.integrator,
    )
    logger.info(
        "[LOO orbit] AL-RMS = %.2f mas  2D-RMS = %.1f mas  chi²_red = %.2f  (n_obs=%d, nfev=%d)",
        loo_diag["rms_al_mas"], loo_diag["rms_2d_mas"],
        loo_diag["chi2_red"], loo_diag["n_obs"], loo_diag["nfev"],
    )

    # ── Phase B: mass fit on in-window residuals ─────────────────────────────
    if n_win < 6:
        logger.warning("Only %d in-window observations — mass fit will be weak.", n_win)

    mass_result = fit_mass_from_window(
        fitted_elements=fitted_elements,
        perturber_elements=perturber_el,
        obs_jd_tdb=jd_tdb[win_mask],
        gaia_xyz=gaia_xyz[win_mask],
        ra_obs=ra_obs[win_mask],
        dec_obs=dec_obs[win_mask],
        enc_jd_tdb=enc_jd_tdb,
        initial_mass_kg=initial_mass,
        pa_scan=pa_scan[win_mask],
        ra_err_sys=ra_err_sys[win_mask],
        dec_err_sys=dec_err_sys[win_mask],
        corr_sys=corr_sys[win_mask],
        ra_err_rand=ra_err_rand[win_mask],
        dec_err_rand=dec_err_rand[win_mask],
        corr_rand=corr_rand[win_mask],
        blackout_days=args.blackout_days,
        big4_elements=big4_elements,
        dt_days=args.dt_days,
        integrator=args.integrator,
    )

    logger.info("")
    logger.info("=== LOO FIT RESULT ===")
    logger.info(
        "  LOO orbit: AL-RMS=%.2f mas  2D-RMS=%.1f mas  chi²_red=%.2f  n_out=%d",
        loo_diag["rms_al_mas"], loo_diag["rms_2d_mas"],
        loo_diag["chi2_red"], n_out,
    )
    logger.info(
        "  Window AL residuals:  pre=%.2f mas  post=%.2f mas",
        mass_result["rms_pre_al_mas"], mass_result["rms_post_al_mas"],
    )
    logger.info(
        "  AL shift (post−pre): %.2f mas  (2D pre=%.1f mas  post=%.1f mas)",
        mass_result["shift_al_mas"],
        mass_result["rms_pre_2d_mas"], mass_result["rms_post_2d_mas"],
    )
    if math.isfinite(mass_result["mass_kg"]):
        logger.info(
            "  Fitted mass = %.3e ± %.2e kg  (log10=%.2f ± %.2f)",
            mass_result["mass_kg"], mass_result["mass_sigma_kg"],
            mass_result["log10_mass"], mass_result["log10_mass_sigma"],
        )
    else:
        logger.info("  Mass fit: %s", mass_result.get("note", "failed"))

    # ── Save ─────────────────────────────────────────────────────────────────
    tag = f"{args.perturber:06d}_{args.target:06d}_loo"
    out_path = args.output or Path("data/output") / f"fit_{tag}.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    summary = {
        "method": "gaia_loo_al_weighted",
        "perturber": args.perturber,
        "target": args.target,
        "encounter_date": args.date,
        "loo_window_days": args.loo_window_days,
        "include_big4": True,
        "big4_used": list(big4_elements.keys()),
        "n_obs_total": obs.height,
        "n_loo_orbit": n_out,
        "n_window": n_win,
        "loo_orbit_al_rms_mas": loo_diag["rms_al_mas"],
        "loo_orbit_2d_rms_mas": loo_diag["rms_2d_mas"],
        "loo_orbit_chi2_red": loo_diag["chi2_red"],
        "loo_orbit_nfev": loo_diag["nfev"],
        "fitted_elements": fitted_elements,
        **mass_result,
    }
    with out_path.open("w") as fh:
        json.dump(summary, fh, indent=2)
    logger.info("Wrote → %s", out_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
