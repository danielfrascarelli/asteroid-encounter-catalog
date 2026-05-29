"""Multi-target joint mass fit for a single perturber across N Gaia encounters.

The fit shares one log10_M_perturber across all targets and keeps six orbital
deltas free per target, breaking the M ↔ deltas degeneracy that biased the
single-target Stage 4 validation (see ``docs/mass_layer_validation.md``).

Usage
-----
    docker compose run --rm pipeline python -m scripts.mass.fit_mass_gaia_multitarget \\
        --perturber 2 \\
        --targets-csv data/output/stage4_validation_summary.csv \\
        --likelihood mahalanobis2d \\
        --priors default \\
        --output data/output/multitarget/fit_000002_pallas.json
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import math
from dataclasses import dataclass
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
from src.mass.forward_model_joint import PRIOR_PRESETS, GaiaObservationBundle, resolve_priors
from src.mass.forward_model_joint_multitarget import (
    TargetBundle,
    fit_joint_multitarget,
)
from src.propagate.nbody import _MAJOR_ASTEROIDS
from src.utils.config import load_config

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class TargetSpec:
    target: int
    date_utc: str

    @property
    def jd_tdb(self) -> float:
        return float(Time(self.date_utc, scale="utc").tdb.jd)


def _read_targets_from_csv(path: Path, perturber: int) -> list[TargetSpec]:
    out: list[TargetSpec] = []
    seen: set[int] = set()
    with path.open() as fh:
        reader = csv.DictReader(fh)
        if "perturber" not in reader.fieldnames or "target" not in reader.fieldnames:
            raise ValueError(
                f"{path}: expected columns 'perturber', 'target', 'encounter_date'; "
                f"got {reader.fieldnames}"
            )
        date_col = "encounter_date" if "encounter_date" in reader.fieldnames else "date"
        for row in reader:
            try:
                if int(row["perturber"]) != perturber:
                    continue
                t = int(row["target"])
                if t in seen:
                    continue
                seen.add(t)
                out.append(TargetSpec(target=t, date_utc=row[date_col]))
            except (KeyError, ValueError):
                continue
    return out


def _read_targets_from_json(path: Path) -> list[TargetSpec]:
    data = json.loads(path.read_text())
    if not isinstance(data, list):
        raise ValueError(f"{path}: expected JSON list of {{target, date}} entries")
    return [TargetSpec(target=int(d["target"]), date_utc=str(d["date"])) for d in data]


def _build_bundle(
    *,
    target: int,
    encounter_jd_tdb: float,
    archive_url: str,
    target_elements: dict,
    perturber_elements: dict,
    background_elements: dict[str, dict],
    loo_window_days: float,
    blackout_days: float,
    loo_max_nfev: int,
    dt_days: float,
    integrator: str,
) -> tuple[TargetBundle, dict] | None:
    """Fetch Gaia obs, run LOO orbit fit, and assemble a TargetBundle.

    Returns ``None`` if there are not enough observations for the joint window.
    """
    obs = fetch_gaia_full(archive_url, target)
    if obs.height < 15:
        logger.warning("  target %d: too few transits (%d), skipping", target, obs.height)
        return None

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

    days_from_enc = jd_tdb - encounter_jd_tdb
    loo_mask = days_from_enc < -loo_window_days
    joint_mask = (days_from_enc > -loo_window_days) & (np.abs(days_from_enc) >= blackout_days)
    n_loo = int(loo_mask.sum())
    n_joint = int(joint_mask.sum())
    if n_loo < 8 or n_joint < 8:
        logger.warning(
            "  target %d: insufficient obs (loo=%d joint=%d), skipping", target, n_loo, n_joint
        )
        return None

    fitted_elements, loo_diag = fit_orbit_loo(
        target_elements=target_elements,
        perturber_elements=perturber_elements,
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
        max_nfev=loo_max_nfev,
        dt_days=dt_days,
        integrator=integrator,
    )

    bundle_obs = GaiaObservationBundle(
        jd_tdb=jd_tdb[joint_mask],
        gaia_xyz_bary=gaia_xyz[joint_mask],
        ra_deg=ra_obs[joint_mask],
        dec_deg=dec_obs[joint_mask],
        position_angle_scan_deg=pa_scan[joint_mask],
        ra_error_systematic_mas=ra_err_sys[joint_mask],
        dec_error_systematic_mas=dec_err_sys[joint_mask],
        ra_dec_correlation_systematic=corr_sys[joint_mask],
        ra_error_random_mas=ra_err_rand[joint_mask],
        dec_error_random_mas=dec_err_rand[joint_mask],
        ra_dec_correlation_random=corr_rand[joint_mask],
    )
    bundle = TargetBundle(target_number=target, elements=fitted_elements, obs=bundle_obs)
    diag = {
        "target": target,
        "n_obs_total": obs.height,
        "n_loo_orbit": n_loo,
        "n_joint": n_joint,
        "loo_orbit_al_rms_mas": loo_diag["rms_al_mas"],
        "loo_orbit_chi2_red": loo_diag["chi2_red"],
    }
    return bundle, diag


def _json_default(value):
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


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


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--perturber", type=int, required=True)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--targets-csv", type=Path, help="CSV with perturber,target,encounter_date")
    group.add_argument("--targets-json", type=Path, help="JSON list of {target,date} objects")
    parser.add_argument("--loo-window-days", type=float, default=180.0)
    parser.add_argument("--blackout-days", type=float, default=7.0)
    parser.add_argument("--mpcorb", type=Path, default=None)
    parser.add_argument("--dt-days", type=float, default=1.0)
    parser.add_argument("--integrator", default="whfast", choices=["whfast", "ias15"])
    parser.add_argument("--background-n", type=int, default=20)
    parser.add_argument("--loo-max-nfev", type=int, default=800)
    parser.add_argument("--max-nfev", type=int, default=1500)
    parser.add_argument("--likelihood", choices=["al", "mahalanobis2d"], default="mahalanobis2d")
    parser.add_argument("--priors", choices=sorted(PRIOR_PRESETS), default="default")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    priors = resolve_priors(args.priors)
    cfg = load_config(args.config)
    archive_url = cfg.sources.gaia_sso.archive_url

    if args.targets_csv is not None:
        target_specs = _read_targets_from_csv(args.targets_csv, args.perturber)
    else:
        target_specs = _read_targets_from_json(args.targets_json)
    if not target_specs:
        logger.error("No targets found for perturber %d", args.perturber)
        return 1
    logger.info(
        "Perturber %d: %d targets — %s",
        args.perturber,
        len(target_specs),
        [s.target for s in target_specs],
    )

    if args.mpcorb is None:
        mid_jd = float(np.mean([s.jd_tdb for s in target_specs]))
        args.mpcorb = _best_mpcorb_snapshot(_MPCORB_ARCHIVE_DIR, mid_jd)
        logger.info("Auto-selected MPCORB snapshot: %s", args.mpcorb.name)

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

    target_numbers = [s.target for s in target_specs]
    all_nums = [args.perturber] + target_numbers + [num for num, _ in background_names.values()]
    elements_map = load_element_rows(args.mpcorb, all_nums)
    perturber_el = elements_map[args.perturber]
    background_elements = {
        name: elements_map[num]
        for name, (num, _gm) in background_names.items()
        if num in elements_map
    }

    bundles: list[TargetBundle] = []
    diagnostics_per_target: list[dict] = []
    for spec in target_specs:
        logger.info(
            "Building bundle: perturber=%d target=%d date=%s",
            args.perturber,
            spec.target,
            spec.date_utc,
        )
        if spec.target not in elements_map:
            logger.warning("  target %d: no MPCORB elements found, skipping", spec.target)
            continue
        result = _build_bundle(
            target=spec.target,
            encounter_jd_tdb=spec.jd_tdb,
            archive_url=archive_url,
            target_elements=elements_map[spec.target],
            perturber_elements=perturber_el,
            background_elements=background_elements,
            loo_window_days=args.loo_window_days,
            blackout_days=args.blackout_days,
            loo_max_nfev=args.loo_max_nfev,
            dt_days=args.dt_days,
            integrator=args.integrator,
        )
        if result is None:
            continue
        bundle, diag = result
        diag["encounter_date"] = spec.date_utc
        bundles.append(bundle)
        diagnostics_per_target.append(diag)

    if len(bundles) < 2:
        logger.error(
            "Multi-target fit needs ≥2 valid bundles; got %d. Aborting.", len(bundles)
        )
        return 1

    initial_mass = _mass_from_h(perturber_el.get("H"))
    initial_log10_mass = math.log10(initial_mass if initial_mass > 0 else 1e18)

    logger.info(
        "Running multi-target joint fit: N=%d targets, params=%d, likelihood=%s, priors=%s",
        len(bundles),
        1 + 6 * len(bundles),
        args.likelihood,
        args.priors,
    )
    result = fit_joint_multitarget(
        bundles,
        perturber_el,
        initial_log10_mass=initial_log10_mass,
        priors=priors,
        background_elements=background_elements,
        dt_days=args.dt_days,
        integrator=args.integrator,
        max_nfev=args.max_nfev,
        likelihood=args.likelihood,
    )

    n_obs_per_target = [len(b.obs.jd_tdb) for b in bundles]
    n_obs_total = sum(n_obs_per_target)
    n_astrometric = n_obs_total if args.likelihood == "al" else 2 * n_obs_total
    n_params = 1 + 6 * len(bundles)
    chi2_data = float(np.sum(result.fun[:n_astrometric] ** 2))
    chi2_red = chi2_data / max(1, n_astrometric - n_params)
    log10_sigma, mass_sigma = _mass_uncertainty(result, chi2_red)
    mass_kg = 10.0 ** float(result.x[0])
    condition = _jtj_condition(result)

    per_target_deltas = []
    for i, bundle in enumerate(bundles):
        s = 1 + 6 * i
        per_target_deltas.append(
            {
                "target": bundle.target_number,
                "da_rel": float(result.x[s + 0]),
                "de": float(result.x[s + 1]),
                "di_deg": float(result.x[s + 2]),
                "dOmega_deg": float(result.x[s + 3]),
                "domega_deg": float(result.x[s + 4]),
                "dM_deg": float(result.x[s + 5]),
            }
        )

    summary = {
        "method": "gaia_joint_orbit_mass_multitarget",
        "likelihood": args.likelihood,
        "priors_preset": args.priors,
        "perturber": args.perturber,
        "n_targets": len(bundles),
        "targets": [b.target_number for b in bundles],
        "n_params": n_params,
        "n_obs_total": n_obs_total,
        "n_astrometric_residuals": n_astrometric,
        "mpcorb": str(args.mpcorb),
        "loo_window_days": args.loo_window_days,
        "blackout_days": args.blackout_days,
        "background_n": len(background_elements),
        "background_used": list(background_elements.keys()),
        "joint_success": bool(result.success),
        "joint_status": int(result.status),
        "joint_message": result.message,
        "joint_nfev": int(result.nfev),
        "chi2_red_joint": chi2_red,
        "jtj_condition": condition,
        "mass_kg": mass_kg,
        "mass_sigma_kg": mass_sigma,
        "log10_mass": float(result.x[0]),
        "log10_mass_sigma": log10_sigma,
        "per_target_deltas": per_target_deltas,
        "per_target_diagnostics": diagnostics_per_target,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(summary, indent=2, default=_json_default))
    logger.info(
        "Wrote %s | log10_M=%.3f σ=%.3f | χ²_red=%.2f | success=%s",
        args.output,
        result.x[0],
        log10_sigma,
        chi2_red,
        result.success,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
