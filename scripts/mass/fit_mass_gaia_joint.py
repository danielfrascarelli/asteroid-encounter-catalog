"""Joint orbit-drift + perturber-mass fit for one Gaia close encounter.

This wrapper keeps the existing LOO orbit fit as an initialisation step, then
fits perturber mass and six residual target-orbit deltas in one least-squares
problem over the encounter window.

Usage
-----
    docker compose run --rm pipeline python -m scripts.mass.fit_mass_gaia_joint \\
        --perturber 57 --target 216875 --date 2016-08-26
"""

from __future__ import annotations

import argparse
import json
import logging
import math
from pathlib import Path

import numpy as np
from astropy.time import Time

from scripts.mass.fit_mass_gaia_loo import (
    _J2010_TCB_JD,
    _MPCORB_ARCHIVE_DIR,
    _best_mpcorb_snapshot,
    _mass_from_h,
    fetch_gaia_full,
    fit_orbit_loo,
    load_element_rows,
)
from src.mass.forward_model_joint import (
    DEFAULT_PRIORS,
    GaiaObservationBundle,
    JointFitPriors,
    fit_joint,
)
from src.propagate.nbody import _MAJOR_ASTEROIDS
from src.utils.config import load_config

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
logger = logging.getLogger(__name__)


_PARAM_NAMES = [
    "log10_mass",
    "da_rel",
    "de",
    "di_deg",
    "dOmega_deg",
    "domega_deg",
    "dM_deg",
]


def _build_obs_bundle(
    *,
    jd_tdb: np.ndarray,
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
    mask: np.ndarray,
) -> GaiaObservationBundle:
    return GaiaObservationBundle(
        jd_tdb=jd_tdb[mask],
        gaia_xyz_bary=gaia_xyz[mask],
        ra_deg=ra_obs[mask],
        dec_deg=dec_obs[mask],
        position_angle_scan_deg=pa_scan[mask],
        ra_error_systematic_mas=ra_err_sys[mask],
        dec_error_systematic_mas=dec_err_sys[mask],
        ra_dec_correlation_systematic=corr_sys[mask],
        ra_error_random_mas=ra_err_rand[mask],
        dec_error_random_mas=dec_err_rand[mask],
        ra_dec_correlation_random=corr_rand[mask],
    )


def _active_bounds(params: np.ndarray, priors: JointFitPriors) -> list[str]:
    active: list[str] = []
    lo = priors.lower_bounds
    hi = priors.upper_bounds
    span = hi - lo
    for name, value, low, high, width in zip(_PARAM_NAMES, params, lo, hi, span, strict=True):
        tolerance = max(width * 0.01, 1e-12)
        if abs(value - low) <= tolerance or abs(high - value) <= tolerance:
            active.append(name)
    return active


def _mass_uncertainty(result, chi2_red: float) -> tuple[float, float]:
    try:
        cov = np.linalg.inv(result.jac.T @ result.jac) * chi2_red
        log10_sigma = float(np.sqrt(cov[0, 0]))
    except np.linalg.LinAlgError:
        log10_sigma = float("nan")
    mass = 10.0 ** float(result.x[0])
    mass_sigma = mass * log10_sigma * math.log(10.0)
    return log10_sigma, mass_sigma


def _jtj_condition(result) -> float:
    try:
        return float(np.linalg.cond(result.jac.T @ result.jac))
    except np.linalg.LinAlgError:
        return float("nan")


def _json_default(value):
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--perturber", type=int, required=True)
    parser.add_argument("--target", type=int, required=True)
    parser.add_argument("--date", required=True, help="Encounter date, ISO UTC")
    parser.add_argument("--loo-window-days", type=float, default=180.0)
    parser.add_argument("--blackout-days", type=float, default=7.0)
    parser.add_argument("--mpcorb", type=Path, default=None)
    parser.add_argument("--dt-days", type=float, default=1.0)
    parser.add_argument("--integrator", default="whfast", choices=["whfast", "ias15"])
    parser.add_argument("--background-n", type=int, default=20)
    parser.add_argument("--loo-max-nfev", type=int, default=800)
    parser.add_argument("--max-nfev", type=int, default=800)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()

    cfg = load_config(args.config)
    archive_url = cfg.sources.gaia_sso.archive_url
    enc_jd_tdb = float(Time(args.date, scale="utc").tdb.jd)

    if args.mpcorb is None:
        args.mpcorb = _best_mpcorb_snapshot(_MPCORB_ARCHIVE_DIR, enc_jd_tdb)
        logger.info("Auto-selected MPCORB snapshot: %s", args.mpcorb.name)

    logger.info("Fetching all Gaia DR3 observations of target %d", args.target)
    obs = fetch_gaia_full(archive_url, args.target)
    if obs.height < 15:
        logger.error("Too few transits (%d) for joint fit.", obs.height)
        return 1

    jd_tcb = obs["epoch"].to_numpy() + _J2010_TCB_JD
    jd_tdb = Time(jd_tcb, format="jd", scale="tcb").tdb.jd.astype(float)
    gaia_xyz = np.column_stack(
        [obs["x_gaia"].to_numpy(), obs["y_gaia"].to_numpy(), obs["z_gaia"].to_numpy()]
    ).astype(float)
    ra_obs = obs["ra"].to_numpy().astype(float)
    dec_obs = obs["dec"].to_numpy().astype(float)
    pa_scan = obs["position_angle_scan"].to_numpy().astype(float)
    ra_err_sys = obs["ra_error_systematic"].to_numpy().astype(float)
    dec_err_sys = obs["dec_error_systematic"].to_numpy().astype(float)
    corr_sys = obs["ra_dec_correlation_systematic"].to_numpy().astype(float)
    ra_err_rand = obs["ra_error_random"].to_numpy().astype(float)
    dec_err_rand = obs["dec_error_random"].to_numpy().astype(float)
    corr_rand = obs["ra_dec_correlation_random"].to_numpy().astype(float)

    days_from_enc = jd_tdb - enc_jd_tdb
    loo_mask = days_from_enc < -args.loo_window_days
    joint_mask = (days_from_enc > -args.loo_window_days) & (
        np.abs(days_from_enc) >= args.blackout_days
    )
    n_loo = int(loo_mask.sum())
    n_joint = int(joint_mask.sum())
    n_pre = int(((days_from_enc < -args.blackout_days) & joint_mask).sum())
    n_post = int(((days_from_enc > args.blackout_days) & joint_mask).sum())
    logger.info("LOO obs=%d | joint obs=%d (%d pre, %d post)", n_loo, n_joint, n_pre, n_post)
    if n_loo < 8 or n_joint < 8:
        logger.error("Insufficient observations for joint fit.")
        return 1

    background_n = max(0, min(args.background_n, len(_MAJOR_ASTEROIDS)))
    background_registry = dict(
        sorted(_MAJOR_ASTEROIDS.items(), key=lambda kv: kv[1][1], reverse=True)
    )
    background_names: dict[str, tuple[int, float]] = {}
    for name, (num, gm) in background_registry.items():
        if len(background_names) >= background_n:
            break
        if num != args.perturber:
            background_names[name] = (num, gm)

    all_nums = [args.target, args.perturber] + [num for num, _ in background_names.values()]
    elements_map = load_element_rows(args.mpcorb, all_nums)
    target_el = elements_map[args.target]
    perturber_el = elements_map[args.perturber]
    background_elements = {
        name: elements_map[num]
        for name, (num, _gm) in background_names.items()
        if num in elements_map
    }

    fitted_elements, loo_diag = fit_orbit_loo(
        target_elements=target_el,
        perturber_elements=perturber_el,
        obs_jd_tdb=jd_tdb[loo_mask],
        gaia_xyz=gaia_xyz[loo_mask],
        ra_obs=ra_obs[loo_mask],
        dec_obs=dec_obs[loo_mask],
        pa_scan=pa_scan[loo_mask],
        ra_err_sys=ra_err_sys[loo_mask],
        dec_err_sys=dec_err_sys[loo_mask],
        corr_sys=corr_sys[loo_mask],
        ra_err_rand=ra_err_rand[loo_mask],
        dec_err_rand=dec_err_rand[loo_mask],
        corr_rand=corr_rand[loo_mask],
        background_elements=background_elements,
        max_nfev=args.loo_max_nfev,
        dt_days=args.dt_days,
        integrator=args.integrator,
    )

    initial_mass = _mass_from_h(perturber_el.get("H"))
    initial_log10_mass = math.log10(initial_mass if initial_mass > 0 else 1e18)
    joint_obs = _build_obs_bundle(
        jd_tdb=jd_tdb,
        gaia_xyz=gaia_xyz,
        ra_obs=ra_obs,
        dec_obs=dec_obs,
        pa_scan=pa_scan,
        ra_err_sys=ra_err_sys,
        dec_err_sys=dec_err_sys,
        corr_sys=corr_sys,
        ra_err_rand=ra_err_rand,
        dec_err_rand=dec_err_rand,
        corr_rand=corr_rand,
        mask=joint_mask,
    )
    result = fit_joint(
        fitted_elements,
        perturber_el,
        joint_obs,
        initial_log10_mass=initial_log10_mass,
        priors=DEFAULT_PRIORS,
        background_elements=background_elements,
        dt_days=args.dt_days,
        integrator=args.integrator,
        max_nfev=args.max_nfev,
    )

    chi2_data = float(np.sum(result.fun[:n_joint] ** 2))
    chi2_red = chi2_data / max(1, n_joint - len(_PARAM_NAMES))
    log10_sigma, mass_sigma = _mass_uncertainty(result, chi2_red)
    mass_kg = 10.0 ** float(result.x[0])
    active_bounds = _active_bounds(result.x, DEFAULT_PRIORS)
    condition = _jtj_condition(result)

    logger.info(
        "Joint fit done: success=%s nfev=%d chi2_red=%.2f log10M=%.3f active_bounds=%s",
        result.success,
        result.nfev,
        chi2_red,
        result.x[0],
        active_bounds,
    )

    tag = f"{args.perturber:06d}_{args.target:06d}_joint"
    out_path = args.output or Path("data/output") / f"fit_{tag}.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    summary = {
        "method": "gaia_joint_orbit_mass_al_weighted",
        "perturber": args.perturber,
        "target": args.target,
        "encounter_date": args.date,
        "mpcorb": str(args.mpcorb),
        "loo_window_days": args.loo_window_days,
        "blackout_days": args.blackout_days,
        "loo_max_nfev": args.loo_max_nfev,
        "joint_max_nfev": args.max_nfev,
        "background_n": len(background_elements),
        "background_used": list(background_elements.keys()),
        "n_obs_total": obs.height,
        "n_loo_orbit": n_loo,
        "n_joint": n_joint,
        "n_pre": n_pre,
        "n_post": n_post,
        "loo_orbit_al_rms_mas": loo_diag["rms_al_mas"],
        "loo_orbit_chi2_red": loo_diag["chi2_red"],
        "joint_success": bool(result.success),
        "joint_status": int(result.status),
        "joint_message": result.message,
        "joint_nfev": int(result.nfev),
        "chi2_red_joint": chi2_red,
        "jtj_condition": condition,
        "active_bounds": active_bounds,
        "mass_kg": mass_kg,
        "mass_sigma_kg": mass_sigma,
        "log10_mass": float(result.x[0]),
        "log10_mass_sigma": log10_sigma,
        "da_rel": float(result.x[1]),
        "de": float(result.x[2]),
        "di_deg": float(result.x[3]),
        "dOmega_deg": float(result.x[4]),
        "domega_deg": float(result.x[5]),
        "dM_deg": float(result.x[6]),
        "fitted_elements_initial": fitted_elements,
    }
    out_path.write_text(json.dumps(summary, indent=2, default=_json_default))
    logger.info("Wrote %s", out_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
