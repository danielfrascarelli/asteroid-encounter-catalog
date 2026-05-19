"""Fit a perturber's mass from Gaia DR3 astrometry of its target.

This is the heart of the mini-roadmap. For a given (perturber, target,
encounter date) tuple:

  1. Pull Gaia DR3 observations of the target in a window around the
     encounter (default ±180 d).
  2. Load both bodies' osculating elements from MPCORB.
  3. Run scipy.optimize.least_squares over the 7 parameters:
         (a, e, i, Omega, omega, M₀, log10_mass)
     - First 6 = target's osculating elements (small deviations from MPCORB).
     - log10_mass = perturber mass (log-space for non-negativity).
  4. Residuals: (obs RA − pred RA) · cos(dec), (obs Dec − pred Dec) in mas.
  5. Report fitted mass with uncertainty from the Jacobian.

Output: ``data/output/fit_<perturber>_<target>.csv`` plus a summary JSON.

Usage
-----
    docker compose run --rm pipeline python -m scripts.fit_perturber_mass \
        --perturber 165 --target 31067 --date 2014-12-08
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


def fetch_gaia_observations(
    archive_url: str, target: int, d_min: float, d_max: float
) -> pl.DataFrame:
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


def load_element_row(snapshot: Path, number: int) -> dict:
    df = parse_mpcorb(str(snapshot), semimajor_min_au=0.0, semimajor_max_au=50.0)
    sub = df.filter(pl.col("number") == number)
    if sub.height == 0:
        raise ValueError(f"asteroid {number} not in {snapshot}")
    return sub.row(0, named=True)


def _mass_from_diameter(d_km: float, rho_kg_m3: float = 1500.0) -> float:
    """Sphere mass estimate."""
    if d_km is None or (isinstance(d_km, float) and math.isnan(d_km)) or d_km <= 0.0:
        return 1.0e18  # generic ~100 km MBA fallback
    radius_m = 0.5 * d_km * 1000.0
    return rho_kg_m3 * (4.0 / 3.0) * math.pi * radius_m ** 3


def fit_orbit_only(
    target_elements: dict,
    perturber_elements: dict,
    obs_jd_tdb: np.ndarray,
    gaia_xyz: np.ndarray,
    ra_obs: np.ndarray,
    dec_obs: np.ndarray,
) -> tuple[dict, dict]:
    """Phase A of the two-phase fit: fit only the target's orbital elements.

    Uses *pre-encounter* observations only. The perturber is treated as
    massless (mass = 0) so it doesn't affect the trajectory, just provides
    the structure the forward model expects.

    Returns
    -------
    (fitted_elements, diag)
        ``fitted_elements`` is a dict with the same MPCORB keys as the input
        but updated to the best-fit values. ``diag`` has chi²_red, nfev, etc.
    """
    a0 = target_elements["a_au"]
    e0 = target_elements["e"]
    i0 = target_elements["i_deg"]
    Omega0 = target_elements["Omega_deg"]  # noqa: N806
    omega0 = target_elements["omega_deg"]
    M0 = target_elements["M_deg"]  # noqa: N806

    x0 = np.array([a0, e0, i0, Omega0, omega0, M0])
    lo = np.array([a0 - 0.01, max(0.0, e0 - 0.01), max(0.0, i0 - 0.5),
                   Omega0 - 5.0, omega0 - 5.0, M0 - 5.0])
    hi = np.array([a0 + 0.01, min(0.999, e0 + 0.01), i0 + 0.5,
                   Omega0 + 5.0, omega0 + 5.0, M0 + 5.0])

    def residual_func(params: np.ndarray) -> np.ndarray:
        tgt = dict(target_elements)
        tgt["a_au"] = float(params[0])
        tgt["e"] = float(params[1])
        tgt["i_deg"] = float(params[2])
        tgt["Omega_deg"] = float(params[3])
        tgt["omega_deg"] = float(params[4])
        tgt["M_deg"] = float(params[5])
        ra_pred, dec_pred = forward_model(
            tgt, perturber_elements, 0.0, obs_jd_tdb, gaia_xyz,  # mass=0
        )
        dra, ddec = residuals_mas(ra_obs, dec_obs, ra_pred, dec_pred)
        return np.concatenate([dra, ddec])

    logger.info(
        "[Phase A] orbit-only fit: 6 params, %d pre-encounter obs (%d residuals)",
        len(obs_jd_tdb), 2 * len(obs_jd_tdb),
    )
    res = least_squares(
        residual_func, x0, method="trf", bounds=(lo, hi), max_nfev=100, verbose=2,
    )
    logger.info("[Phase A] done: nfev=%d cost=%.3e", res.nfev, res.cost)

    fitted = dict(target_elements)
    fitted["a_au"] = float(res.x[0])
    fitted["e"] = float(res.x[1])
    fitted["i_deg"] = float(res.x[2])
    fitted["Omega_deg"] = float(res.x[3])
    fitted["omega_deg"] = float(res.x[4])
    fitted["M_deg"] = float(res.x[5])

    chi2 = float(np.sum(res.fun ** 2))
    n_data = len(res.fun)
    chi2_red = chi2 / max(1, n_data - 6)
    return fitted, {
        "chi2": chi2,
        "chi2_red": chi2_red,
        "nfev": res.nfev,
        "fitted_elements": fitted,
        "residual_rms_mas": float(np.sqrt(chi2 / max(1, n_data))),
    }


def fit_mass_only(
    target_elements: dict,
    perturber_elements: dict,
    obs_jd_tdb: np.ndarray,
    gaia_xyz: np.ndarray,
    ra_obs: np.ndarray,
    dec_obs: np.ndarray,
    initial_mass_kg: float,
) -> dict:
    """Phase B of the two-phase fit: fix the orbit, fit only log10_mass.

    Uses *post-encounter* observations only. Orbit is fixed (presumably from
    Phase A) so the only thing that can explain post-encounter residuals is
    the perturber's gravity.
    """
    log_m0 = math.log10(initial_mass_kg)
    x0 = np.array([log_m0])
    lo = np.array([15.0])
    hi = np.array([22.0])

    def residual_func(params: np.ndarray) -> np.ndarray:
        mass = float(10.0 ** params[0])
        ra_pred, dec_pred = forward_model(
            target_elements, perturber_elements, mass, obs_jd_tdb, gaia_xyz,
        )
        dra, ddec = residuals_mas(ra_obs, dec_obs, ra_pred, dec_pred)
        return np.concatenate([dra, ddec])

    logger.info(
        "[Phase B] mass-only fit: 1 param, %d post-encounter obs (%d residuals)",
        len(obs_jd_tdb), 2 * len(obs_jd_tdb),
    )
    res = least_squares(
        residual_func, x0, method="trf", bounds=(lo, hi), max_nfev=100, verbose=2,
    )
    logger.info("[Phase B] done: nfev=%d cost=%.3e", res.nfev, res.cost)

    log_m_fit = float(res.x[0])
    chi2 = float(np.sum(res.fun ** 2))
    n_data = len(res.fun)
    chi2_red = chi2 / max(1, n_data - 1)
    try:
        cov = np.linalg.inv(res.jac.T @ res.jac) * chi2_red
        log_m_sig = float(np.sqrt(cov[0, 0]))
    except np.linalg.LinAlgError:
        log_m_sig = float("nan")
    mass_fit = 10.0 ** log_m_fit
    mass_sig = mass_fit * log_m_sig * math.log(10.0)
    return {
        "mass_kg": mass_fit,
        "mass_sigma_kg": mass_sig,
        "log10_mass": log_m_fit,
        "log10_mass_sigma": log_m_sig,
        "chi2": chi2,
        "chi2_red": chi2_red,
        "n_data": n_data,
        "nfev": res.nfev,
        "residual_rms_mas": float(np.sqrt(chi2 / max(1, n_data))),
        "post_residuals_mas": res.fun.tolist(),
    }


def fit_mass_two_phase(
    target_elements: dict,
    perturber_elements: dict,
    obs_jd_tdb: np.ndarray,
    gaia_xyz: np.ndarray,
    ra_obs: np.ndarray,
    dec_obs: np.ndarray,
    encounter_jd_tdb: float,
    initial_mass_kg: float,
    blackout_days: float = 7.0,
    min_obs_per_side: int = 6,
) -> dict:
    """Two-phase mass fit: orbit from pre-encounter, mass from post-encounter.

    This separates the perturbation signal from any orbit-fit ambiguity: the
    orbit is determined by pre-encounter data alone (where the perturber has
    not yet acted), and the mass is then the *only* parameter free to
    explain post-encounter residuals.
    """
    days_from_enc = obs_jd_tdb - encounter_jd_tdb
    pre_mask = days_from_enc < -blackout_days
    post_mask = days_from_enc > blackout_days
    n_pre = int(pre_mask.sum())
    n_post = int(post_mask.sum())
    if n_pre < min_obs_per_side or n_post < min_obs_per_side:
        raise ValueError(
            f"Need ≥ {min_obs_per_side} obs each side; got pre={n_pre}, post={n_post}"
        )
    logger.info(
        "Two-phase fit: %d pre + %d post observations (blackout ±%.1f d)",
        n_pre, n_post, blackout_days,
    )

    # Phase A: orbit from pre-encounter
    fitted_elements, phase_a = fit_orbit_only(
        target_elements=target_elements,
        perturber_elements=perturber_elements,
        obs_jd_tdb=obs_jd_tdb[pre_mask],
        gaia_xyz=gaia_xyz[pre_mask],
        ra_obs=ra_obs[pre_mask],
        dec_obs=dec_obs[pre_mask],
    )
    logger.info(
        "[Phase A] residual RMS = %.1f mas (chi²_red=%.1f)",
        phase_a["residual_rms_mas"], phase_a["chi2_red"],
    )

    # Phase B: mass from post-encounter, with orbit fixed
    phase_b = fit_mass_only(
        target_elements=fitted_elements,
        perturber_elements=perturber_elements,
        obs_jd_tdb=obs_jd_tdb[post_mask],
        gaia_xyz=gaia_xyz[post_mask],
        ra_obs=ra_obs[post_mask],
        dec_obs=dec_obs[post_mask],
        initial_mass_kg=initial_mass_kg,
    )
    logger.info(
        "[Phase B] residual RMS = %.1f mas (chi²_red=%.1f)",
        phase_b["residual_rms_mas"], phase_b["chi2_red"],
    )

    return {
        "mass_kg": phase_b["mass_kg"],
        "mass_sigma_kg": phase_b["mass_sigma_kg"],
        "log10_mass": phase_b["log10_mass"],
        "log10_mass_sigma": phase_b["log10_mass_sigma"],
        "phase_a_chi2_red": phase_a["chi2_red"],
        "phase_a_rms_mas": phase_a["residual_rms_mas"],
        "phase_a_nfev": phase_a["nfev"],
        "phase_b_chi2_red": phase_b["chi2_red"],
        "phase_b_rms_mas": phase_b["residual_rms_mas"],
        "phase_b_nfev": phase_b["nfev"],
        "n_pre": n_pre,
        "n_post": n_post,
        "fitted_elements_phase_a": fitted_elements,
        "post_residuals_mas": phase_b["post_residuals_mas"],
    }


def fit_mass(
    target_elements: dict,
    perturber_elements: dict,
    obs_jd_tdb: np.ndarray,
    gaia_xyz: np.ndarray,
    ra_obs: np.ndarray,
    dec_obs: np.ndarray,
    initial_mass_kg: float,
    fit_orbit: bool = True,
) -> dict:
    """Run scipy.optimize.least_squares for (orbit, mass).

    Parameters
    ----------
    target_elements / perturber_elements:
        MPCORB-style dicts. Target's elements are the initial guess for the
        6 orbital parameters; perturber's are fixed.
    obs_jd_tdb, gaia_xyz, ra_obs, dec_obs:
        Gaia observations.
    initial_mass_kg:
        Initial guess for the perturber's mass.
    fit_orbit:
        If True, fit (a, e, i, Omega, omega, M₀, log10_mass) together (7 dof).
        If False, fix the orbit and only fit log10_mass (1 dof).
    """
    a0 = target_elements["a_au"]
    e0 = target_elements["e"]
    i0 = target_elements["i_deg"]
    Omega0 = target_elements["Omega_deg"]
    omega0 = target_elements["omega_deg"]
    M0 = target_elements["M_deg"]
    log_m0 = math.log10(initial_mass_kg)

    if fit_orbit:
        x0 = np.array([a0, e0, i0, Omega0, omega0, M0, log_m0])
        # Reasonable bounds: small deviations on orbital elements + mass in [1e15, 1e22]
        lo = np.array([a0 - 0.01, max(0.0, e0 - 0.01), max(0.0, i0 - 0.5),
                       Omega0 - 5.0, omega0 - 5.0, M0 - 5.0, 15.0])
        hi = np.array([a0 + 0.01, min(0.999, e0 + 0.01), i0 + 0.5,
                       Omega0 + 5.0, omega0 + 5.0, M0 + 5.0, 22.0])
    else:
        x0 = np.array([log_m0])
        lo = np.array([15.0])
        hi = np.array([22.0])

    def _build_target_elements(params: np.ndarray) -> tuple[dict, float]:
        if fit_orbit:
            a, e, i, Omega, omega, M, log_m = params  # noqa: N806
        else:
            (log_m,) = params
            a, e, i, Omega, omega, M = a0, e0, i0, Omega0, omega0, M0  # noqa: N806
        tgt = dict(target_elements)
        tgt["a_au"] = float(a)
        tgt["e"] = float(e)
        tgt["i_deg"] = float(i)
        tgt["Omega_deg"] = float(Omega)
        tgt["omega_deg"] = float(omega)
        tgt["M_deg"] = float(M)
        return tgt, float(10.0 ** log_m)

    def residual_func(params: np.ndarray) -> np.ndarray:
        tgt, mass = _build_target_elements(params)
        ra_pred, dec_pred = forward_model(
            tgt, perturber_elements, mass, obs_jd_tdb, gaia_xyz,
        )
        dra, ddec = residuals_mas(ra_obs, dec_obs, ra_pred, dec_pred)
        return np.concatenate([dra, ddec])

    logger.info(
        "Starting least_squares: %d free params, %d observations, %d residuals",
        len(x0), len(obs_jd_tdb), 2 * len(obs_jd_tdb),
    )
    res = least_squares(
        residual_func,
        x0,
        method="trf",
        bounds=(lo, hi),
        max_nfev=100,
        verbose=2,
    )
    logger.info("Done: nfev=%d cost=%.3e", res.nfev, res.cost)

    # Uncertainty: (Jᵀ J)⁻¹ × χ²_red, diagonal = variance per param
    jac = res.jac
    chi2 = float(np.sum(res.fun ** 2))
    n_data = len(res.fun)
    n_params = len(res.x)
    chi2_red = chi2 / max(1, n_data - n_params)
    try:
        cov = np.linalg.inv(jac.T @ jac) * chi2_red
        sigma = np.sqrt(np.diag(cov))
    except np.linalg.LinAlgError:
        sigma = np.full(n_params, float("nan"))

    if fit_orbit:
        log_m_fit = res.x[-1]
        log_m_sig = sigma[-1]
    else:
        log_m_fit = res.x[0]
        log_m_sig = sigma[0]

    mass_fit = 10.0 ** log_m_fit
    # Propagate log → linear uncertainty
    mass_sig = mass_fit * log_m_sig * math.log(10.0)

    return {
        "mass_kg": mass_fit,
        "mass_sigma_kg": mass_sig,
        "log10_mass": log_m_fit,
        "log10_mass_sigma": log_m_sig,
        "chi2": chi2,
        "chi2_red": chi2_red,
        "n_data": n_data,
        "n_params": n_params,
        "nfev": res.nfev,
        "residuals_mas": res.fun.tolist(),
        "fitted_params": res.x.tolist(),
        "param_sigmas": sigma.tolist(),
    }


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--config", default="config.yaml")
    p.add_argument("--perturber", type=int, required=True)
    p.add_argument("--target", type=int, required=True)
    p.add_argument("--date", required=True, help="Encounter date, ISO UTC")
    p.add_argument("--half-window-days", type=float, default=180.0)
    p.add_argument(
        "--initial-mass-kg",
        type=float,
        default=None,
        help="Initial mass guess (kg). If unset, derived from perturber diameter.",
    )
    p.add_argument(
        "--fit-orbit",
        action="store_true",
        help="Fit orbital elements jointly with mass (default: mass-only).",
    )
    p.add_argument(
        "--two-phase",
        action="store_true",
        help="Two-phase fit: orbit from pre-encounter, mass from post-encounter. "
             "Mutually exclusive with --fit-orbit.",
    )
    p.add_argument(
        "--blackout-days",
        type=float,
        default=7.0,
        help="Half-width of blackout window around encounter (days).",
    )
    p.add_argument(
        "--mpcorb",
        type=Path,
        default=Path("data/raw/mpcorb_archive/MPCORB_20120918.DAT"),
        help="MPCORB snapshot to load elements from.",
    )
    args = p.parse_args()

    cfg = load_config(args.config)
    archive_url = cfg.sources.gaia_sso.archive_url

    enc_jd_tdb = float(Time(args.date, scale="utc").tdb.jd)
    enc_days = float(Time(args.date, scale="utc").tcb.jd) - _J2010_TCB_JD

    logger.info("Fetching Gaia observations of target %d…", args.target)
    obs = fetch_gaia_observations(
        archive_url,
        args.target,
        enc_days - args.half_window_days,
        enc_days + args.half_window_days,
    )
    if obs.height < 6:
        logger.error("Only %d Gaia transits — too few to fit", obs.height)
        return 1
    logger.info("  → %d transits", obs.height)

    epochs_days = obs["epoch"].to_numpy()
    jd_tcb = epochs_days + _J2010_TCB_JD
    jd_tdb = Time(jd_tcb, format="jd", scale="tcb").tdb.jd.astype(float)
    gaia_xyz = np.column_stack([
        obs["x_gaia"].to_numpy(),
        obs["y_gaia"].to_numpy(),
        obs["z_gaia"].to_numpy(),
    ]).astype(float)
    ra_obs = obs["ra"].to_numpy().astype(float)
    dec_obs = obs["dec"].to_numpy().astype(float)

    logger.info("Loading orbital elements from %s", args.mpcorb)
    target_el = load_element_row(args.mpcorb, args.target)
    perturber_el = load_element_row(args.mpcorb, args.perturber)
    logger.info(
        "Target %d: a=%.4f e=%.4f i=%.2f° epoch=%.1f",
        args.target, target_el["a_au"], target_el["e"], target_el["i_deg"],
        target_el["epoch_jd"],
    )
    logger.info(
        "Perturber %d: a=%.4f e=%.4f i=%.2f° epoch=%.1f H=%.2f",
        args.perturber, perturber_el["a_au"], perturber_el["e"],
        perturber_el["i_deg"], perturber_el["epoch_jd"], perturber_el.get("H", 0.0),
    )

    # Initial mass: from perturber diameter (H → D via ρ=0.14 albedo, ρ=1.5 g/cm³)
    if args.initial_mass_kg is not None:
        initial_mass = args.initial_mass_kg
    else:
        h = perturber_el.get("H", None)
        if h is not None:
            d_km = (1329.0 / math.sqrt(0.14)) * 10.0 ** (-h / 5.0)
        else:
            d_km = 100.0
        initial_mass = _mass_from_diameter(d_km)
    logger.info("Initial mass guess: %.2e kg", initial_mass)

    if args.two_phase and args.fit_orbit:
        logger.error("--two-phase and --fit-orbit are mutually exclusive")
        return 2

    if args.two_phase:
        fit = fit_mass_two_phase(
            target_el, perturber_el,
            obs_jd_tdb=jd_tdb, gaia_xyz=gaia_xyz,
            ra_obs=ra_obs, dec_obs=dec_obs,
            encounter_jd_tdb=enc_jd_tdb,
            initial_mass_kg=initial_mass,
            blackout_days=args.blackout_days,
        )
        logger.info("")
        logger.info("=== TWO-PHASE FIT RESULT ===")
        logger.info(
            "  fitted mass = %.3e ± %.2e kg  (%.2f log10)",
            fit["mass_kg"], fit["mass_sigma_kg"], fit["log10_mass"],
        )
        logger.info(
            "  Phase A (orbit, pre-encounter): %d obs, RMS=%.1f mas, chi²_red=%.1f",
            fit["n_pre"], fit["phase_a_rms_mas"], fit["phase_a_chi2_red"],
        )
        logger.info(
            "  Phase B (mass, post-encounter): %d obs, RMS=%.1f mas, chi²_red=%.1f",
            fit["n_post"], fit["phase_b_rms_mas"], fit["phase_b_chi2_red"],
        )
        n_obs_for_save = fit["n_pre"] + fit["n_post"]
        chi2_save = fit["phase_b_chi2_red"]
        n_params_save = 1
        nfev_save = fit["phase_a_nfev"] + fit["phase_b_nfev"]
        fitted_params_save = fit["fitted_elements_phase_a"]
        param_sigmas_save = None
    else:
        fit = fit_mass(
            target_el, perturber_el,
            obs_jd_tdb=jd_tdb, gaia_xyz=gaia_xyz,
            ra_obs=ra_obs, dec_obs=dec_obs,
            initial_mass_kg=initial_mass,
            fit_orbit=args.fit_orbit,
        )
        logger.info("")
        logger.info("=== FIT RESULT ===")
        logger.info(
            "  fitted mass = %.3e ± %.2e kg  (%.2f log10)",
            fit["mass_kg"], fit["mass_sigma_kg"], fit["log10_mass"],
        )
        logger.info("  chi2 = %.2e   chi2_red = %.2f", fit["chi2"], fit["chi2_red"])
        logger.info("  n_obs = %d   n_params = %d   nfev = %d",
                    fit["n_data"] // 2, fit["n_params"], fit["nfev"])
        n_obs_for_save = fit["n_data"] // 2
        chi2_save = fit["chi2_red"]
        n_params_save = fit["n_params"]
        nfev_save = fit["nfev"]
        fitted_params_save = fit["fitted_params"]
        param_sigmas_save = fit["param_sigmas"]

    # Save
    out_dir = Path("data/output")
    out_dir.mkdir(parents=True, exist_ok=True)
    tag = f"{args.perturber:06d}_{args.target:06d}"
    if args.two_phase:
        tag = tag + "_2phase"
    summary_path = out_dir / f"fit_{tag}.json"
    summary = {
        "perturber": args.perturber,
        "target": args.target,
        "encounter_date": args.date,
        "fitted_mass_kg": fit["mass_kg"],
        "fitted_mass_sigma_kg": fit["mass_sigma_kg"],
        "log10_mass": fit["log10_mass"],
        "log10_mass_sigma": fit["log10_mass_sigma"],
        "n_obs": n_obs_for_save,
        "n_params": n_params_save,
        "chi2_red": chi2_save,
        "nfev": nfev_save,
        "method": "two_phase" if args.two_phase else ("joint_fit" if args.fit_orbit else "mass_only"),
    }
    if args.two_phase:
        summary["phase_a_chi2_red"] = fit["phase_a_chi2_red"]
        summary["phase_a_rms_mas"] = fit["phase_a_rms_mas"]
        summary["phase_b_chi2_red"] = fit["phase_b_chi2_red"]
        summary["phase_b_rms_mas"] = fit["phase_b_rms_mas"]
        summary["n_pre"] = fit["n_pre"]
        summary["n_post"] = fit["n_post"]
        summary["fitted_elements_phase_a"] = fit["fitted_elements_phase_a"]
    else:
        summary["fitted_params"] = fitted_params_save
        summary["param_sigmas"] = param_sigmas_save

    with summary_path.open("w") as fh:
        json.dump(summary, fh, indent=2)
    logger.info("Wrote summary → %s", summary_path)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
