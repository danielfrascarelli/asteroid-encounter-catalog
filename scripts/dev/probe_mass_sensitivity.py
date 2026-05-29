"""Probe the astrometric sensitivity of the joint forward model to perturber mass.

Decisive check for the closing-the-loop result: the noiseless fit returned the
H-based initial mass unchanged (χ²≈0 at the *wrong* mass). Either (a) these
encounters carry no measurable deflection signal, or (b) the mass parameter is
not propagating into the predicted RA/Dec (a bug).

This builds one real target bundle and evaluates ``forward_model`` over a grid
of perturber masses, reporting the max on-sky displacement (mas) relative to a
near-zero mass baseline. If injecting the perturber's literature mass moves
positions by ≫ the ~mas Gaia precision, signal exists and the optimiser is at
fault; if ≪ mas, the encounters are genuinely weak.

Usage
-----
    docker compose run --rm pipeline python -m scripts.dev.probe_mass_sensitivity \\
        --perturber 2 --target 28036 --encounter-date 2017-03-01 \\
        --masses 1e10 1e19 5e19 1.17e20 2.05e20 4e20
"""

from __future__ import annotations

import argparse
import logging

import numpy as np

from scripts.mass.fit_mass_gaia_loo import (
    _MPCORB_ARCHIVE_DIR,
    _best_mpcorb_snapshot,
    load_element_rows,
)
from scripts.mass.fit_mass_gaia_multitarget import TargetSpec, _build_bundle
from src.astrometry.forward_model import forward_model, residuals_mas
from src.mass.forward_model_joint import al_residuals_and_weights
from src.propagate.nbody import _MAJOR_ASTEROIDS
from src.utils.config import load_config

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
logger = logging.getLogger(__name__)

_PLANETS = (
    "sun",
    "mercury",
    "venus",
    "earth",
    "mars",
    "jupiter",
    "saturn",
    "uranus",
    "neptune",
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--perturber", type=int, required=True)
    parser.add_argument("--target", type=int, required=True)
    parser.add_argument("--encounter-date", required=True, help="UTC ISO date")
    parser.add_argument(
        "--masses",
        type=float,
        nargs="+",
        default=[1e10, 1e19, 5e19, 1.17e20, 2.05e20, 4e20],
    )
    parser.add_argument("--loo-window-days", type=float, default=180.0)
    parser.add_argument("--blackout-days", type=float, default=7.0)
    parser.add_argument("--dt-days", type=float, default=1.0)
    parser.add_argument("--integrator", default="whfast", choices=["whfast", "ias15"])
    parser.add_argument("--background-n", type=int, default=20)
    args = parser.parse_args()

    cfg = load_config(args.config)
    archive_url = cfg.sources.gaia_sso.archive_url
    spec = TargetSpec(target=args.target, date_utc=args.encounter_date)

    mpcorb = _best_mpcorb_snapshot(_MPCORB_ARCHIVE_DIR, spec.jd_tdb)
    logger.info("MPCORB snapshot: %s", mpcorb.name)

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

    all_nums = [args.perturber, args.target] + [num for num, _ in background_names.values()]
    elements_map = load_element_rows(mpcorb, all_nums)
    perturber_el = elements_map[args.perturber]
    background_elements = {
        name: elements_map[num]
        for name, (num, _gm) in background_names.items()
        if num in elements_map
    }

    built = _build_bundle(
        target=args.target,
        encounter_jd_tdb=spec.jd_tdb,
        archive_url=archive_url,
        target_elements=elements_map[args.target],
        perturber_elements=perturber_el,
        background_elements=background_elements,
        loo_window_days=args.loo_window_days,
        blackout_days=args.blackout_days,
        loo_max_nfev=800,
        dt_days=args.dt_days,
        integrator=args.integrator,
    )
    if built is None:
        logger.error("Could not build bundle for target %d", args.target)
        return 1
    bundle, _diag = built
    obs = bundle.obs

    def predict(mass_kg: float) -> tuple[np.ndarray, np.ndarray]:
        return forward_model(
            bundle.elements,
            perturber_el,
            float(mass_kg),
            obs.jd_tdb,
            obs.gaia_xyz_bary,
            include_planets=_PLANETS,
            include_background=bool(background_elements),
            background_elements=background_elements,
            dt_days=args.dt_days,
            integrator=args.integrator,
        )

    masses = sorted(args.masses)
    base_mass = masses[0]
    ra0, dec0 = predict(base_mass)
    # Gaia SSO astrometry is precise only along-scan (AL); across-scan (AC) is
    # ~100s of mas. The meaningful noise floor is the AL-projected 1-sigma the
    # fit actually uses, not the AC-dominated RA/Dec quadrature sum.
    zeros = np.zeros_like(obs.ra_deg)
    _r_al, sigma_al = al_residuals_and_weights(
        zeros,
        zeros,
        obs.position_angle_scan_deg,
        obs.ra_error_systematic_mas,
        obs.dec_error_systematic_mas,
        obs.ra_dec_correlation_systematic,
        obs.ra_error_random_mas,
        obs.dec_error_random_mas,
        obs.ra_dec_correlation_random,
    )
    noise_floor = float(np.median(sigma_al))
    logger.info(
        "Bundle target=%d: %d joint obs. Baseline mass=%.2e kg. "
        "Median along-scan (AL) 1-sigma: %.3f mas.",
        args.target,
        len(obs.jd_tdb),
        base_mass,
        noise_floor,
    )
    print(f"\n{'mass_kg':>12}  {'rms_2d mas':>12}  {'rms_AL mas':>11}  {'chi_AL (frozen δ)':>17}")
    print("-" * 60)
    for m in masses:
        ra, dec = predict(m)
        dra, ddec = residuals_mas(ra, dec, ra0, dec0)  # (m) - (baseline), in mas
        rms = float(np.sqrt(np.mean(dra**2 + ddec**2)))
        r_al, _s = al_residuals_and_weights(
            dra,
            ddec,
            obs.position_angle_scan_deg,
            obs.ra_error_systematic_mas,
            obs.dec_error_systematic_mas,
            obs.ra_dec_correlation_systematic,
            obs.ra_error_random_mas,
            obs.dec_error_random_mas,
            obs.ra_dec_correlation_random,
        )
        rms_al = float(np.sqrt(np.mean(r_al**2)))
        # Whitened along-scan signal norm with deltas frozen = intrinsic mass
        # detectability (√Δχ²); detection needs ≳3.
        chi_al = float(np.sqrt(np.sum((r_al / sigma_al) ** 2)))
        print(f"{m:12.3e}  {rms:12.4f}  {rms_al:11.4f}  {chi_al:17.3f}")
    print(f"\nMedian along-scan 1-sigma: {noise_floor:.3f} mas; N_obs={len(obs.jd_tdb)}")
    print("chi_AL = whitened AL signal vs baseline with orbital deltas frozen.")
    print("It is the intrinsic mass detectability (√Δχ²); a 3σ detection needs ≳3.\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
