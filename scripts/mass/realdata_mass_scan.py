"""Real-data profiled χ²(mass) scan — Track A Stage 2.6 diagnostic.

Background
----------
Stage 2.5 fixed the optimiser bug (closing-loop recovers the injected mass
exactly) but the *real* Hygiea fit lands at ~9× the literature mass with a good
χ²_red≈1.2 (see ``docs/mass_layer_stage_a2_5_profiled.md``). With the optimiser
and the M↔deltas degeneracy both ruled out, the remaining suspect is a
data/model systematic in the signal regime. This script probes it directly.

For each target of a perturber *individually* (and for the joint bundle), it
scans the profiled astrometric χ² over ``log10_M`` on the **real** Gaia
observations — at every trial mass the 6 orbital deltas are re-optimised
(``_profiled_delta_fit``), so the curve is the true profiled likelihood, not a
deltas-frozen slice. It answers three questions the followup plan poses:

1. **Is the real-data χ²(mass) unimodal?** A clean single minimum means the
   optimiser is trustworthy and the high mass is what the data prefer; a flat
   or multi-modal curve means non-identifiability.
2. **Where is each target's minimum vs literature?** Per-target minima cluster
   near 9× → a coherent systematic (orbit drift / Gaia along-scan bias);
   minima scattered with one or two pulling high → a few bad encounters.
3. **What do the per-observation residuals look like at the optimum?** A
   post-encounter temporal ramp in the along-scan residual is the signature of
   unmodelled orbital drift the mass is absorbing.

Usage
-----
    docker compose run --rm pipeline python -m scripts.mass.realdata_mass_scan \\
        --perturber 10 --targets-csv data/output/stage4_validation_summary.csv \\
        --lit-mass-kg 8.3e19 \\
        --out-prefix data/output/stage2_6/hygiea
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
    _mass_from_h,
    load_element_rows,
)
from scripts.mass.fit_mass_gaia_multitarget import (
    _build_bundle,
    _json_default,
    _read_targets_from_csv,
    _read_targets_from_json,
)
from src.mass.forward_model_joint import PRIOR_PRESETS, resolve_priors
from src.mass.forward_model_joint_multitarget import (
    TargetBundle,
    _astrometric_residuals_one_target,
    _profiled_delta_fit,
)
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


def _profiled_chi2(
    log10_mass: float,
    bundles: list[TargetBundle],
    perturber_el: dict,
    *,
    priors,
    background_elements,
    dt_days: float,
    integrator: str,
    likelihood: str,
    inner_max_nfev: int,
) -> tuple[float, np.ndarray]:
    """Profiled astrometric χ² at *log10_mass* + the optimised delta vector."""
    n_obs = sum(len(b.obs.jd_tdb) for b in bundles)
    n_astro = n_obs if likelihood == "al" else 2 * n_obs
    res = _profiled_delta_fit(
        log10_mass,
        bundles,
        perturber_el,
        priors=priors,
        background_elements=background_elements,
        include_planets=_PLANETS,
        dt_days=dt_days,
        integrator=integrator,
        likelihood=likelihood,
        max_nfev=inner_max_nfev,
    )
    chi2 = float(np.sum(res.fun[:n_astro] ** 2))
    return chi2, np.asarray(res.x, dtype=float)


def _per_obs_chi2(
    bundle: TargetBundle,
    perturber_el: dict,
    log10_mass: float,
    deltas6: np.ndarray,
    *,
    background_elements,
    dt_days: float,
    integrator: str,
    likelihood: str,
) -> np.ndarray:
    """χ² contribution of each observation at a given mass+deltas."""
    params7 = np.concatenate([[log10_mass], deltas6])
    whitened = _astrometric_residuals_one_target(
        params7,
        bundle.elements,
        perturber_el,
        bundle.obs,
        background_elements,
        _PLANETS,
        dt_days,
        integrator,
        likelihood,
    )
    if likelihood == "mahalanobis2d":
        w = whitened.reshape(-1, 2)
        return np.sum(w**2, axis=1)
    return whitened**2


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--perturber", type=int, required=True)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--targets-csv", type=Path)
    group.add_argument("--targets-json", type=Path)
    parser.add_argument("--lit-mass-kg", type=float, required=True)
    parser.add_argument("--grid-lo-dex", type=float, default=-1.6, help="grid start, dex below lit")
    parser.add_argument("--grid-hi-dex", type=float, default=1.6, help="grid end, dex above lit")
    parser.add_argument("--grid-n", type=int, default=21)
    parser.add_argument("--loo-window-days", type=float, default=180.0)
    parser.add_argument("--blackout-days", type=float, default=7.0)
    parser.add_argument(
        "--joint-window-days",
        type=float,
        default=None,
        help="If set, bound each bundle's joint obs to a symmetric ±W-day band around "
        "the encounter (the default _build_bundle window is one-sided and lets the "
        "post-encounter arc run to the end of DR3, accumulating orbit drift).",
    )
    parser.add_argument("--mpcorb", type=Path, default=None)
    parser.add_argument("--dt-days", type=float, default=1.0)
    parser.add_argument("--integrator", default="whfast", choices=["whfast", "ias15"])
    parser.add_argument("--background-n", type=int, default=20)
    parser.add_argument("--loo-max-nfev", type=int, default=800)
    parser.add_argument("--inner-max-nfev", type=int, default=1500)
    parser.add_argument("--likelihood", choices=["al", "mahalanobis2d"], default="mahalanobis2d")
    parser.add_argument("--priors", choices=sorted(PRIOR_PRESETS), default="default")
    parser.add_argument("--out-prefix", type=Path, required=True)
    parser.add_argument(
        "--release",
        default=None,
        help="Gaia release ('dr3' | 'fpr'). Defaults to config's gaia_sso.release.",
    )
    args = parser.parse_args()

    priors = resolve_priors(args.priors)
    cfg = load_config(args.config)
    gaia = cfg.sources.gaia_sso
    if args.release is not None:
        gaia.release = args.release
    release_cfg = gaia.active()
    archive_url = gaia.archive_url
    logger.info("Gaia release: %s (table %s)", gaia.release, release_cfg.table)

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

    bundles: list[TargetBundle] = []
    for spec in target_specs:
        if spec.target not in elements_map:
            continue
        logger.info("Building bundle: perturber=%d target=%d", args.perturber, spec.target)
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
            joint_window_days=args.joint_window_days,
            release_cfg=release_cfg,
        )
        if result is None:
            continue
        bundle, _diag = result
        if args.joint_window_days is not None:
            logger.info(
                "  target %d: ±%.0fd joint window keeps %d obs",
                spec.target,
                args.joint_window_days,
                len(bundle.obs.jd_tdb),
            )
        bundles.append(bundle)

    if not bundles:
        logger.error("No valid bundles built")
        return 1

    lit_log10 = math.log10(args.lit_mass_kg)
    h_mass = _mass_from_h(perturber_el.get("H"))
    h_log10 = math.log10(h_mass if h_mass > 0 else 1e18)
    grid = np.linspace(lit_log10 + args.grid_lo_dex, lit_log10 + args.grid_hi_dex, args.grid_n)

    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)

    # ---- per-target + joint scans -------------------------------------------
    scan_groups: dict[str, list[TargetBundle]] = {f"target_{b.target_number}": [b] for b in bundles}
    if len(bundles) > 1:
        scan_groups["joint"] = bundles

    common = dict(
        priors=priors,
        background_elements=background_elements,
        dt_days=args.dt_days,
        integrator=args.integrator,
        likelihood=args.likelihood,
        inner_max_nfev=args.inner_max_nfev,
    )

    scan_rows: list[dict] = []
    summary_groups: list[dict] = []
    for gname, gbundles in scan_groups.items():
        logger.info("Scanning %s (%d obs)", gname, sum(len(b.obs.jd_tdb) for b in gbundles))
        chi2_curve = np.empty(len(grid))
        for k, g in enumerate(grid):
            chi2, _x = _profiled_chi2(float(g), gbundles, perturber_el, **common)
            chi2_curve[k] = chi2
            scan_rows.append(
                {"group": gname, "log10_M": float(g), "M_kg": float(10.0**g), "chi2_astro": chi2}
            )

        kmin = int(np.argmin(chi2_curve))
        # χ² at the literature mass (interpolated onto the grid via nearest scan)
        chi2_lit, _x_lit = _profiled_chi2(lit_log10, gbundles, perturber_el, **common)
        n_obs = sum(len(b.obs.jd_tdb) for b in gbundles)
        n_astro = n_obs if args.likelihood == "al" else 2 * n_obs
        n_params = 1 + 6 * len(gbundles)
        dof = max(1, n_astro - n_params)
        # unimodality: count interior local minima on the grid
        c = chi2_curve
        n_local_min = int(
            sum(1 for i in range(1, len(c) - 1) if c[i] < c[i - 1] and c[i] < c[i + 1])
        )
        summary_groups.append(
            {
                "group": gname,
                "n_targets": len(gbundles),
                "n_obs": n_obs,
                "log10_M_min_grid": float(grid[kmin]),
                "M_min_grid_kg": float(10.0 ** grid[kmin]),
                "ratio_min_over_lit": float(10.0 ** (grid[kmin] - lit_log10)),
                "chi2_at_min": float(chi2_curve[kmin]),
                "chi2_red_at_min": float(chi2_curve[kmin] / dof),
                "chi2_at_lit": float(chi2_lit),
                "chi2_red_at_lit": float(chi2_lit / dof),
                "delta_chi2_lit_minus_min": float(chi2_lit - chi2_curve[kmin]),
                "sigma_pref_for_min": float(math.sqrt(max(0.0, chi2_lit - chi2_curve[kmin]))),
                "n_interior_local_minima": n_local_min,
                "unimodal": n_local_min <= 1,
            }
        )

    # ---- write scan + summary first (so a per-obs error can't lose the scan) -
    import csv as _csv

    scan_csv = args.out_prefix.with_name(args.out_prefix.name + "_scan.csv")
    with scan_csv.open("w", newline="") as fh:
        w = _csv.DictWriter(fh, fieldnames=["group", "log10_M", "M_kg", "chi2_astro"])
        w.writeheader()
        w.writerows(scan_rows)

    summary = {
        "perturber": args.perturber,
        "lit_mass_kg": args.lit_mass_kg,
        "lit_log10_M": lit_log10,
        "h_prior_log10_M": h_log10,
        "likelihood": args.likelihood,
        "priors_preset": args.priors,
        "grid_log10_M": grid.tolist(),
        "groups": summary_groups,
    }
    summary_json = args.out_prefix.with_name(args.out_prefix.name + "_summary.json")
    summary_json.write_text(json.dumps(summary, indent=2, default=_json_default))

    # ---- per-observation residuals (single-target groups only) --------------
    # _profiled_delta_fit returns the 6 optimised deltas (mass is held fixed),
    # so x_min / x_lit are length-6 delta vectors, not the full (1+6) vector.
    per_obs_rows: list[dict] = []
    for b in bundles:
        gb = [b]
        # mass at this target's own profiled minimum
        chi2_curve = np.array(
            [r["chi2_astro"] for r in scan_rows if r["group"] == f"target_{b.target_number}"]
        )
        kmin = int(np.argmin(chi2_curve))
        log10_min = float(grid[kmin])
        _, deltas_min = _profiled_chi2(log10_min, gb, perturber_el, **common)
        _, deltas_lit = _profiled_chi2(lit_log10, gb, perturber_el, **common)
        chi2_obs_min = _per_obs_chi2(
            b,
            perturber_el,
            log10_min,
            deltas_min[:6],
            background_elements=background_elements,
            dt_days=args.dt_days,
            integrator=args.integrator,
            likelihood=args.likelihood,
        )
        chi2_obs_lit = _per_obs_chi2(
            b,
            perturber_el,
            lit_log10,
            deltas_lit[:6],
            background_elements=background_elements,
            dt_days=args.dt_days,
            integrator=args.integrator,
            likelihood=args.likelihood,
        )
        jd = b.obs.jd_tdb
        for j in range(len(jd)):
            per_obs_rows.append(
                {
                    "target": b.target_number,
                    "jd_tdb": float(jd[j]),
                    "days_rel": float(jd[j] - np.median(jd)),
                    "chi2_at_target_min": float(chi2_obs_min[j]),
                    "chi2_at_lit": float(chi2_obs_lit[j]),
                }
            )

    perobs_csv = args.out_prefix.with_name(args.out_prefix.name + "_perobs.csv")
    with perobs_csv.open("w", newline="") as fh:
        w = _csv.DictWriter(
            fh,
            fieldnames=["target", "jd_tdb", "days_rel", "chi2_at_target_min", "chi2_at_lit"],
        )
        w.writeheader()
        w.writerows(per_obs_rows)

    # ---- console report ------------------------------------------------------
    print(f"\n=== Real-data profiled χ²(mass) scan — perturber {args.perturber} ===")
    print(
        f"lit mass = {args.lit_mass_kg:.3e} kg (log10={lit_log10:.3f}); "
        f"H-prior log10={h_log10:.3f}\n"
    )
    hdr = f"{'group':>14}  {'M_min/lit':>10}  {'χ²r@min':>9}  {'χ²r@lit':>9}  {'Δχ²(lit-min)':>13}  {'σ':>6}  uni"
    print(hdr)
    print("-" * len(hdr))
    for s in summary_groups:
        print(
            f"{s['group']:>14}  {s['ratio_min_over_lit']:10.2f}  "
            f"{s['chi2_red_at_min']:9.2f}  {s['chi2_red_at_lit']:9.2f}  "
            f"{s['delta_chi2_lit_minus_min']:13.1f}  {s['sigma_pref_for_min']:6.1f}  "
            f"{'Y' if s['unimodal'] else 'N'}"
        )
    print(f"\nWrote:\n  {scan_csv}\n  {perobs_csv}\n  {summary_json}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
