"""Per-observation diagnosis of the Stage 2 mass-fit outliers — Track B Stage 1.

Stage 2 of the mass layer (2D Mahalanobis joint fit) left two pairs with an
anomalously high reduced chi-squared:

* (124) Alkeste -> 57942  : chi2_red_joint = 84.45
* (389) Industria -> 176865 : chi2_red jumped 4.86 (AL-only) -> 31.87 (2D)

See ``docs/mass_layer_stage2_diagnostic.md``. Three hypotheses were posed:
unmodelled secondary perturber, a Gaia across-scan (AC) catalogue systematic
that the 2D likelihood exposes but the AL-only fit hid, or a single bad transit.

This script discriminates between them. For each pair it rebuilds the exact
Stage 2 joint observation set, re-runs the single-target joint fit, and at the
optimum decomposes every observation into its along-scan (AL) and across-scan
(AC) residuals with their own sigmas, plus the full 2D Mahalanobis chi-squared.
It then asks:

1. **Is the excess in AC, not AL?** If the AL residuals are well-behaved
   (chi2_red_AL ~ 1) but the AC residuals carry the excess, the misfit is a
   2D/AC systematic the AL-only Stage 1 ignored, not a mass/orbit signal.
2. **Is it a temporal ramp?** A post-encounter slope in the AL residual is the
   signature of unmodelled orbital drift the fit cannot absorb.
3. **Is it one bad transit?** A few observations with chi2 >> the rest point to
   an outlier transit (cosmic ray, edge of FOV).
4. **Is there a secondary perturber?** Query the hybrid catalogue for other
   encounters of the target (<0.3 AU, +/-90 d) and flag any massive body.

Usage
-----
    docker compose run --rm pipeline python -m scripts.mass.diagnose_stage2_outliers \\
        --out-prefix data/output/stage2_outliers/diag
"""

from __future__ import annotations

import argparse
import csv as _csv
import json
import logging
import math
from pathlib import Path

import numpy as np
import polars as pl
from astropy.time import Time

from scripts.mass.fit_mass_gaia_loo import (
    _MPCORB_ARCHIVE_DIR,
    _best_mpcorb_snapshot,
    _mass_from_h,
    load_element_rows,
)
from scripts.mass.fit_mass_gaia_multitarget import _build_bundle
from src.astrometry.forward_model import forward_model, residuals_mas
from src.mass.forward_model_joint import (
    DEFAULT_PRIORS,
    al_residuals_and_weights,
    apply_target_deltas,
    fit_joint,
)
from src.mass.likelihood_al import mahalanobis_residuals_2d
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

# The two Stage 2 outliers (perturber, target, encounter date UTC).
_DEFAULT_PAIRS: tuple[tuple[int, int, str], ...] = (
    (124, 57942, "2015-12-26"),
    (389, 176865, "2015-10-17"),
)


def _ac_residuals_and_weights(
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
    """Project tangential residuals onto the across-scan axis with sigma_AC.

    The along-scan unit vector is ``(sin PA, cos PA)`` in ``(RA*, Dec)``; the
    across-scan axis is orthogonal, ``(cos PA, -sin PA)``. Mirrors
    :func:`src.mass.forward_model_joint.al_residuals_and_weights`.
    """
    pa = np.radians(pa_scan_deg)
    e_ac_ra = np.cos(pa)
    e_ac_dec = -np.sin(pa)
    r_ac = dra_mas * e_ac_ra + ddec_mas * e_ac_dec

    def _projected_var(s_ra: np.ndarray, s_dec: np.ndarray, rho: np.ndarray) -> np.ndarray:
        out = (
            e_ac_ra**2 * s_ra**2
            + 2.0 * e_ac_ra * e_ac_dec * rho * s_ra * s_dec
            + e_ac_dec**2 * s_dec**2
        )
        return np.asarray(out, dtype=float)

    var_ac = _projected_var(ra_err_sys, dec_err_sys, corr_sys) + _projected_var(
        ra_err_rand, dec_err_rand, corr_rand
    )
    sigma_ac = np.sqrt(np.maximum(var_ac, 1e-6))
    return r_ac, sigma_ac


def _tangential_residuals(
    params7: np.ndarray,
    target_elements: dict,
    perturber_el: dict,
    obs,
    background_elements: dict,
    dt_days: float,
    integrator: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Raw (dra*, ddec) tangential residuals in mas at the given parameters."""
    adjusted = apply_target_deltas(target_elements, params7)
    mass_kg = float(10.0 ** params7[0])
    ra_pred, dec_pred = forward_model(
        adjusted,
        perturber_el,
        mass_kg,
        obs.jd_tdb,
        obs.gaia_xyz_bary,
        include_planets=_PLANETS,
        include_background=bool(background_elements),
        background_elements=background_elements,
        dt_days=dt_days,
        integrator=integrator,
    )
    return residuals_mas(obs.ra_deg, obs.dec_deg, ra_pred, dec_pred)


def _search_secondary_perturbers(
    catalog_path: Path,
    target: int,
    primary: int,
    encounter_jd_tdb: float,
    *,
    window_days: float,
    max_dist_au: float,
) -> dict:
    """Find other close encounters of *target* near its primary encounter.

    The hybrid catalogue is keyed by ``jd_tdb`` (no ``date_utc`` column), so the
    temporal window is applied directly on the Julian date.
    """
    if not catalog_path.exists():
        return {"available": False, "catalog": str(catalog_path)}
    lo, hi = encounter_jd_tdb - window_days, encounter_jd_tdb + window_days
    # Lazy scan: only the small filtered result is materialised (the full
    # hybrid catalogue is ~1 GB, and loading it eagerly per call can OOM).
    sub = (
        pl.scan_parquet(catalog_path)
        .filter(
            (pl.col("dist_au") <= max_dist_au)
            & (pl.col("jd_tdb") >= lo)
            & (pl.col("jd_tdb") <= hi)
            & ((pl.col("number_1") == target) | (pl.col("number_2") == target))
        )
        .collect()
    )
    massive_numbers = {num for _name, (num, _gm) in _MAJOR_ASTEROIDS.items()}
    partners: list[dict] = []
    for row in sub.iter_rows(named=True):
        n1, n2 = int(row["number_1"]), int(row["number_2"])
        partner = n2 if n1 == target else n1
        if partner == primary or partner == target:
            continue
        days_from_primary = float(row["jd_tdb"]) - encounter_jd_tdb
        partners.append(
            {
                "partner": partner,
                "dist_au": float(row["dist_au"]),
                "days_from_primary_encounter": days_from_primary,
                "is_massive_registry": partner in massive_numbers,
            }
        )
    partners.sort(key=lambda r: r["dist_au"])
    return {
        "available": True,
        "n_candidates": len(partners),
        "n_massive": sum(1 for p in partners if p["is_massive_registry"]),
        "closest": partners[:10],
    }


def _diagnose_pair(
    perturber: int,
    target: int,
    date_utc: str,
    *,
    archive_url: str,
    mpcorb: Path,
    catalog_path: Path,
    loo_window_days: float,
    blackout_days: float,
    dt_days: float,
    integrator: str,
    background_n: int,
    max_nfev: int,
) -> tuple[dict, list[dict]]:
    """Run the per-observation diagnosis for one (perturber, target) pair."""
    enc_jd_tdb = float(Time(date_utc, scale="utc").tdb.jd)

    background_n = max(0, min(background_n, len(_MAJOR_ASTEROIDS)))
    background_registry = dict(
        sorted(_MAJOR_ASTEROIDS.items(), key=lambda kv: kv[1][1], reverse=True)
    )
    background_names: dict[str, tuple[int, float]] = {}
    for name, (num, gm) in background_registry.items():
        if len(background_names) >= background_n:
            break
        if num != perturber:
            background_names[name] = (num, gm)

    all_nums = [perturber, target] + [num for num, _ in background_names.values()]
    elements_map = load_element_rows(mpcorb, all_nums)
    perturber_el = elements_map[perturber]
    background_elements = {
        name: elements_map[num]
        for name, (num, _gm) in background_names.items()
        if num in elements_map
    }

    built = _build_bundle(
        target=target,
        encounter_jd_tdb=enc_jd_tdb,
        archive_url=archive_url,
        target_elements=elements_map[target],
        perturber_elements=perturber_el,
        background_elements=background_elements,
        loo_window_days=loo_window_days,
        blackout_days=blackout_days,
        loo_max_nfev=800,
        dt_days=dt_days,
        integrator=integrator,
        joint_window_days=None,  # reproduce the one-sided Stage 2 window
    )
    if built is None:
        return {"perturber": perturber, "target": target, "status": "no_bundle"}, []
    bundle, diag = built

    h_mass = _mass_from_h(perturber_el.get("H"))
    h_log10 = math.log10(h_mass if h_mass > 0 else 1e18)
    result = fit_joint(
        bundle.elements,
        perturber_el,
        bundle.obs,
        initial_log10_mass=h_log10,
        priors=DEFAULT_PRIORS,
        background_elements=background_elements,
        dt_days=dt_days,
        integrator=integrator,
        max_nfev=max_nfev,
        likelihood="mahalanobis2d",
    )
    params7 = np.asarray(result.x, dtype=float)

    obs = bundle.obs
    dra, ddec = _tangential_residuals(
        params7, bundle.elements, perturber_el, obs, background_elements, dt_days, integrator
    )
    r_al, sigma_al = al_residuals_and_weights(
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
    r_ac, sigma_ac = _ac_residuals_and_weights(
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
    _whitened, chi2_2d = mahalanobis_residuals_2d(
        dra,
        ddec,
        obs.ra_error_systematic_mas,
        obs.dec_error_systematic_mas,
        obs.ra_dec_correlation_systematic,
        obs.ra_error_random_mas,
        obs.dec_error_random_mas,
        obs.ra_dec_correlation_random,
    )
    # Per-observation 1-D chi2 in the along-scan and across-scan directions,
    # each weighted by its own projected sigma. These are *not* the additive
    # components of the 2D Mahalanobis chi2 (which folds in the AL-AC
    # correlation); they are independent goodness-of-fit checks per axis. The
    # 2D Mahalanobis chi2 is reported separately as the Stage 2 headline.
    chi2_al = (r_al / sigma_al) ** 2
    chi2_ac = (r_ac / sigma_ac) ** 2
    days_rel = obs.jd_tdb - enc_jd_tdb

    n_obs = len(obs.jd_tdb)
    n_params = 7
    dof = max(1, 2 * n_obs - n_params)
    dof_al = max(1, n_obs - n_params)

    # temporal ramp: slope of AL residual (mas) vs days, post-encounter arc
    post = days_rel > 0
    slope_post = float("nan")
    if int(post.sum()) >= 3:
        slope_post = float(np.polyfit(days_rel[post], r_al[post], 1)[0])

    per_obs_rows = [
        {
            "perturber": perturber,
            "target": target,
            "jd_tdb": float(obs.jd_tdb[j]),
            "days_rel": float(days_rel[j]),
            "r_al_mas": float(r_al[j]),
            "sigma_al_mas": float(sigma_al[j]),
            "r_ac_mas": float(r_ac[j]),
            "sigma_ac_mas": float(sigma_ac[j]),
            "chi2_al": float(chi2_al[j]),
            "chi2_ac": float(chi2_ac[j]),
            "chi2_2d": float(chi2_2d[j]),
        }
        for j in range(n_obs)
    ]

    secondary = _search_secondary_perturbers(
        catalog_path,
        target,
        perturber,
        enc_jd_tdb,
        window_days=90.0,
        max_dist_au=0.3,
    )

    summary = {
        "perturber": perturber,
        "target": target,
        "encounter_date": date_utc,
        "status": "ok",
        "n_obs_joint": n_obs,
        "n_obs_post": int(post.sum()),
        "mass_fit_kg": float(10.0 ** params7[0]),
        "log10_mass_fit": float(params7[0]),
        "deltas6": params7[1:].tolist(),
        "loo_orbit_chi2_red": diag["loo_orbit_chi2_red"],
        "chi2_2d_total": float(np.sum(chi2_2d)),
        "chi2_red_2d": float(np.sum(chi2_2d) / dof),
        "chi2_red_al": float(np.sum(chi2_al) / dof_al),
        "chi2_red_ac": float(np.sum(chi2_ac) / dof_al),
        "median_abs_al_pull": float(np.median(np.abs(r_al / sigma_al))),
        "median_abs_ac_pull": float(np.median(np.abs(r_ac / sigma_ac))),
        "sigma_al_median_mas": float(np.median(sigma_al)),
        "sigma_ac_median_mas": float(np.median(sigma_ac)),
        "frac_obs_chi2_2d_gt25": float(np.mean(chi2_2d > 25.0)),
        "n_obs_chi2_2d_gt25": int(np.sum(chi2_2d > 25.0)),
        "n_obs_chi2_2d_gt100": int(np.sum(chi2_2d > 100.0)),
        "max_chi2_2d": float(np.max(chi2_2d)),
        "top1_frac_of_total": float(np.max(chi2_2d) / max(np.sum(chi2_2d), 1e-12)),
        "top5_frac_of_total": float(
            np.sum(np.sort(chi2_2d)[::-1][:5]) / max(np.sum(chi2_2d), 1e-12)
        ),
        "al_residual_slope_post_mas_per_day": slope_post,
        "rms_al_mas": float(np.sqrt(np.mean(r_al**2))),
        "rms_ac_mas": float(np.sqrt(np.mean(r_ac**2))),
        "secondary_perturber_search": secondary,
    }
    return summary, per_obs_rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument(
        "--catalog",
        type=Path,
        default=Path("data/output/encounters_catalog_hybrid_stageb.parquet"),
    )
    parser.add_argument("--mpcorb", type=Path, default=None)
    parser.add_argument("--loo-window-days", type=float, default=180.0)
    parser.add_argument("--blackout-days", type=float, default=7.0)
    parser.add_argument("--dt-days", type=float, default=1.0)
    parser.add_argument("--integrator", default="whfast", choices=["whfast", "ias15"])
    parser.add_argument("--background-n", type=int, default=20)
    parser.add_argument("--max-nfev", type=int, default=800)
    parser.add_argument("--out-prefix", type=Path, default=Path("data/output/stage2_outliers/diag"))
    args = parser.parse_args()

    cfg = load_config(args.config)
    archive_url = cfg.sources.gaia_sso.archive_url

    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)

    all_summaries: list[dict] = []
    all_perobs: list[dict] = []
    for perturber, target, date_utc in _DEFAULT_PAIRS:
        mpcorb = args.mpcorb
        if mpcorb is None:
            enc_jd_tdb = float(Time(date_utc, scale="utc").tdb.jd)
            mpcorb = _best_mpcorb_snapshot(_MPCORB_ARCHIVE_DIR, enc_jd_tdb)
            logger.info("target %d: MPCORB snapshot %s", target, mpcorb.name)
        logger.info("Diagnosing perturber=%d target=%d", perturber, target)
        summary, per_obs = _diagnose_pair(
            perturber,
            target,
            date_utc,
            archive_url=archive_url,
            mpcorb=mpcorb,
            catalog_path=args.catalog,
            loo_window_days=args.loo_window_days,
            blackout_days=args.blackout_days,
            dt_days=args.dt_days,
            integrator=args.integrator,
            background_n=args.background_n,
            max_nfev=args.max_nfev,
        )
        all_summaries.append(summary)
        all_perobs.extend(per_obs)

    perobs_csv = args.out_prefix.with_name(args.out_prefix.name + "_perobs.csv")
    if all_perobs:
        with perobs_csv.open("w", newline="") as fh:
            w = _csv.DictWriter(fh, fieldnames=list(all_perobs[0].keys()))
            w.writeheader()
            w.writerows(all_perobs)

    summary_json = args.out_prefix.with_name(args.out_prefix.name + "_summary.json")
    summary_json.write_text(json.dumps({"pairs": all_summaries}, indent=2))

    print("\n=== Stage 2 outlier diagnosis ===")
    hdr = (
        f"{'pair':>16}  {'n':>4}  {'chi2r_2d':>9}  {'|AL_pull|':>9}  {'|AC_pull|':>9}  "
        f"{'sigAL':>6}  {'sigAC':>7}  {'top1%':>6}  {'sec_mass':>8}"
    )
    print(hdr)
    print("-" * len(hdr))
    for s in all_summaries:
        if s.get("status") != "ok":
            print(f"{s['perturber']}->{s['target']:>10}  {s.get('status')}")
            continue
        sec = s["secondary_perturber_search"]
        n_massive = sec.get("n_massive", "?") if sec.get("available") else "n/a"
        print(
            f"{str(s['perturber']) + '->' + str(s['target']):>16}  {s['n_obs_joint']:>4}  "
            f"{s['chi2_red_2d']:9.2f}  {s['median_abs_al_pull']:9.2f}  {s['median_abs_ac_pull']:9.2f}  "
            f"{s['sigma_al_median_mas']:6.1f}  {s['sigma_ac_median_mas']:7.1f}  "
            f"{100 * s['top1_frac_of_total']:6.1f}  {n_massive:>8}"
        )
    print(f"\nWrote:\n  {perobs_csv}\n  {summary_json}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
