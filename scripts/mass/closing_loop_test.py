"""Closing-the-loop test for the joint perturber-mass forward model.

Motivation
----------
Track A Stage 1 (tight priors) and Stage 2 (multi-target joint fit) both
failed the calibrator gate, and A2 showed the fitted orbital deltas are ~1e-7
— they absorb no mass signal. That rules out the M↔deltas degeneracy as the
cause of the structural fit/lit bias and points at the deflection forward
model (or the DR3 dataset) itself.

This script discriminates between those two:

1. Build real bundles for a calibrator (real Gaia geometry + LOO-fitted
   elements), exactly as the multi-target CLI does.
2. **Inject** a known perturber mass: regenerate each target's RA/Dec with
   ``forward_model`` at ``M_inject`` and zero deltas. This is synthetic
   "truth" that lives on the real observation epochs and Gaia positions.
3. Optionally add astrometric noise.
4. Re-fit ``log10_M`` + deltas starting from the H-based initial mass (NOT
   the truth), and check whether the fit recovers ``M_inject``.

Interpretation
--------------
- **Noiseless run recovers M_inject (ratio≈1):** the forward model + optimiser
  are self-consistent. The real-data bias is physical / dataset-limited
  (orbit drift, unmodeled perturbers, Gaia systematics) → A3/DR4 territory.
- **Noiseless run does NOT recover M_inject:** the bias is a bug in the
  forward model / fit machinery; neither A2 nor A3 would fix it. Fix the model.
- The realistic-noise run quantifies how much DR3-level noise degrades the
  recovery (leverage of the dataset).

Usage
-----
    docker compose run --rm pipeline python -m scripts.mass.closing_loop_test \\
        --perturber 2 --targets-csv data/output/stage4_validation_summary.csv \\
        --inject-mass-kg 2.05e20 --noise none \\
        --output data/output/multitarget/closing_loop_pallas_noiseless.json
"""

from __future__ import annotations

import argparse
import json
import logging
import math
from pathlib import Path

import numpy as np

from scripts.mass.fit_mass_gaia_loo import (
    _MPCORB_ARCHIVE_DIR,
    _best_mpcorb_snapshot,
    load_element_rows,
)
from scripts.mass.fit_mass_gaia_multitarget import (
    _build_bundle,
    _json_default,
    _mass_from_h,
    _read_targets_from_csv,
    _read_targets_from_json,
)
from src.astrometry.forward_model import forward_model
from src.mass.forward_model_joint import (
    PRIOR_PRESETS,
    GaiaObservationBundle,
    resolve_priors,
)
from src.mass.forward_model_joint_multitarget import TargetBundle, fit_joint_multitarget
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


def _inject_truth_obs(
    bundle: TargetBundle,
    perturber_elements: dict,
    mass_inject_kg: float,
    background_elements: dict[str, dict],
    dt_days: float,
    integrator: str,
    rng: np.random.Generator,
    noise: str,
) -> GaiaObservationBundle:
    """Replace a bundle's RA/Dec with forward_model truth at *mass_inject_kg*.

    Deltas are zero (truth elements == bundle.elements). With ``noise='none'``
    the synthetic obs are exact; with ``noise='realistic'`` per-axis Gaussian
    noise is drawn from the quadrature sum of the systematic and random error
    columns (correlation ignored — conservative, slightly over-weights).
    """
    obs = bundle.obs
    ra_pred, dec_pred = forward_model(
        bundle.elements,
        perturber_elements,
        float(mass_inject_kg),
        obs.jd_tdb,
        obs.gaia_xyz_bary,
        include_planets=_PLANETS,
        include_background=bool(background_elements),
        background_elements=background_elements,
        dt_days=dt_days,
        integrator=integrator,
    )
    ra = np.asarray(ra_pred, dtype=float).copy()
    dec = np.asarray(dec_pred, dtype=float).copy()

    if noise == "realistic":
        sigma_ra_mas = np.hypot(obs.ra_error_systematic_mas, obs.ra_error_random_mas)
        sigma_dec_mas = np.hypot(obs.dec_error_systematic_mas, obs.dec_error_random_mas)
        deg = np.pi / 180.0
        # RA residual convention carries a cos(dec) factor (see residuals_mas):
        # invert it so the injected noise has the intended on-sky magnitude.
        dra_mas = rng.normal(0.0, sigma_ra_mas)
        ddec_mas = rng.normal(0.0, sigma_dec_mas)
        ra = ra + dra_mas / (np.cos(dec * deg) * 3_600_000.0)
        dec = dec + ddec_mas / 3_600_000.0

    return GaiaObservationBundle(
        jd_tdb=obs.jd_tdb,
        gaia_xyz_bary=obs.gaia_xyz_bary,
        ra_deg=ra,
        dec_deg=dec,
        position_angle_scan_deg=obs.position_angle_scan_deg,
        ra_error_systematic_mas=obs.ra_error_systematic_mas,
        dec_error_systematic_mas=obs.dec_error_systematic_mas,
        ra_dec_correlation_systematic=obs.ra_dec_correlation_systematic,
        ra_error_random_mas=obs.ra_error_random_mas,
        dec_error_random_mas=obs.dec_error_random_mas,
        ra_dec_correlation_random=obs.ra_dec_correlation_random,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--perturber", type=int, required=True)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--targets-csv", type=Path)
    group.add_argument("--targets-json", type=Path)
    parser.add_argument("--inject-mass-kg", type=float, required=True)
    parser.add_argument("--noise", choices=["none", "realistic"], default="none")
    parser.add_argument("--seed", type=int, default=20260529)
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
    parser.add_argument(
        "--scan-mass",
        action="store_true",
        help="Instead of fitting, scan χ²(log10_M) with deltas frozen at 0 over the "
        "synthetic-truth data. Shows whether a mass minimum exists (optimiser failure) "
        "or the likelihood is flat (no leverage).",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    priors = resolve_priors(args.priors)
    cfg = load_config(args.config)
    archive_url = cfg.sources.gaia_sso.archive_url
    rng = np.random.default_rng(args.seed)

    if args.targets_csv is not None:
        target_specs = _read_targets_from_csv(args.targets_csv, args.perturber)
    else:
        target_specs = _read_targets_from_json(args.targets_json)
    if not target_specs:
        logger.error("No targets found for perturber %d", args.perturber)
        return 1

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

    real_bundles: list[TargetBundle] = []
    for spec in target_specs:
        if spec.target not in elements_map:
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
        bundle, _diag = result
        real_bundles.append(bundle)

    if len(real_bundles) < 2:
        logger.error("Need ≥2 bundles; got %d", len(real_bundles))
        return 1

    logger.info(
        "Injecting M=%.3e kg into %d bundles (noise=%s)",
        args.inject_mass_kg,
        len(real_bundles),
        args.noise,
    )
    synth_bundles = [
        TargetBundle(
            target_number=b.target_number,
            elements=b.elements,
            obs=_inject_truth_obs(
                b,
                perturber_el,
                args.inject_mass_kg,
                background_elements,
                args.dt_days,
                args.integrator,
                rng,
                args.noise,
            ),
        )
        for b in real_bundles
    ]

    initial_mass = _mass_from_h(perturber_el.get("H"))
    initial_log10_mass = math.log10(initial_mass if initial_mass > 0 else 1e18)

    if args.scan_mass:
        from src.mass.forward_model_joint_multitarget import residuals_joint_multitarget

        truth_log10 = math.log10(args.inject_mass_kg)
        n_targets = len(synth_bundles)
        grid = np.linspace(truth_log10 - 1.2, truth_log10 + 0.6, 19)
        logger.info(
            "Scanning χ²(log10_M) with deltas=0 | truth=%.3f H-init=%.3f",
            truth_log10,
            initial_log10_mass,
        )
        print(f"\n{'log10_M':>10}  {'M_kg':>12}  {'chi2 (δ=0)':>14}  note")
        print("-" * 54)
        chi2_at = {}
        for g in grid:
            params = np.concatenate([[g], np.zeros(6 * n_targets)])
            r = residuals_joint_multitarget(
                params,
                synth_bundles,
                perturber_el,
                priors=priors,
                background_elements=background_elements,
                dt_days=args.dt_days,
                integrator=args.integrator,
                likelihood=args.likelihood,
            )
            n_astro = len(r) - 6 * n_targets
            chi2 = float(np.sum(r[:n_astro] ** 2))
            chi2_at[round(float(g), 4)] = chi2
            note = ""
            if abs(g - truth_log10) < 1e-9:
                note = "<- truth"
            elif abs(g - initial_log10_mass) < 0.06:
                note = "<- ~H-init"
            print(f"{g:10.4f}  {10.0**g:12.3e}  {chi2:14.4f}  {note}")
        gmin = min(chi2_at, key=chi2_at.get)
        print(f"\nχ² minimum on grid at log10_M={gmin:.4f} (M={10.0**gmin:.3e} kg)")
        print(f"truth log10_M={truth_log10:.4f}; H-init log10_M={initial_log10_mass:.4f}\n")
        return 0

    logger.info(
        "Fitting from initial log10_M=%.3f (truth log10_M=%.3f)",
        initial_log10_mass,
        math.log10(args.inject_mass_kg),
    )

    result = fit_joint_multitarget(
        synth_bundles,
        perturber_el,
        initial_log10_mass=initial_log10_mass,
        priors=priors,
        background_elements=background_elements,
        dt_days=args.dt_days,
        integrator=args.integrator,
        max_nfev=args.max_nfev,
        likelihood=args.likelihood,
    )
    logger.info("least_squares status=%s message=%s", result.status, result.message)

    mass_fit = 10.0 ** float(result.x[0])
    ratio = mass_fit / args.inject_mass_kg
    n_obs_total = sum(len(b.obs.jd_tdb) for b in synth_bundles)
    n_astrometric = n_obs_total if args.likelihood == "al" else 2 * n_obs_total
    n_params = 1 + 6 * len(synth_bundles)
    chi2_red = float(np.sum(result.fun[:n_astrometric] ** 2)) / max(1, n_astrometric - n_params)

    per_target_deltas = []
    for i, bundle in enumerate(synth_bundles):
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
        "test": "closing_the_loop",
        "perturber": args.perturber,
        "n_targets": len(synth_bundles),
        "targets": [b.target_number for b in synth_bundles],
        "noise": args.noise,
        "seed": args.seed,
        "likelihood": args.likelihood,
        "priors_preset": args.priors,
        "inject_mass_kg": args.inject_mass_kg,
        "inject_log10_mass": math.log10(args.inject_mass_kg),
        "initial_log10_mass": initial_log10_mass,
        "mass_fit_kg": mass_fit,
        "log10_mass_fit": float(result.x[0]),
        "ratio_fit_over_inject": ratio,
        "chi2_red": chi2_red,
        "joint_success": bool(result.success),
        "joint_nfev": int(result.nfev),
        "per_target_deltas": per_target_deltas,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(summary, indent=2, default=_json_default))
    logger.info(
        "Wrote %s | inject=%.3e fit=%.3e ratio=%.4f χ²_red=%.2f",
        args.output,
        args.inject_mass_kg,
        mass_fit,
        ratio,
        chi2_red,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
