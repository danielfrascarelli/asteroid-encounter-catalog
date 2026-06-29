"""Closing-loop sobre la geometría REAL: ¿el motor es insesgado en estos datos?

Toma los objetivos reales de un perturbador (tiempos de tránsito, posiciones de
Gaia, ángulos de barrido y σ_AL reales), pero **reemplaza la astrometría observada
por la predicción del propio motor** a la masa verdadera (la de la efeméride) con
la órbita semilla MPCORB, más ruido realista: N(0, σ_AL) por CCD + un offset común
por cruce FOV de amplitud ``s_c`` (el piso calibrado). Luego corre el mismo ajuste
(calibración de piso + stacking) y mide el ratio recuperado.

Interpretación:
  * ratio ≈ 1.0  → el motor es insesgado y tiene leverage suficiente en esta
    geometría: el sobre-tiro ~20% de los datos reales es un **sistemático de la
    astrometría real**, no un bug del motor.
  * ratio ≠ 1.0  → hay un problema de leverage/condicionamiento o un sesgo del motor.
"""

from __future__ import annotations

import argparse
import dataclasses
from pathlib import Path

import numpy as np

from scripts.mass.fit_mass_gaia_loo import (
    _MPCORB_ARCHIVE_DIR,
    _best_mpcorb_snapshot,
    load_element_rows,
)
from scripts.mass.fit_mass_gaia_multitarget import _read_targets_from_csv
from scripts.mass.orbdet_fit_realdata import (
    _BIG4_NAME_BY_NUMBER,
    _build_target_obs,
    _epoch_ref,
    _fetch_target,
)
from src.orbdet.constants import M_SUN_KG
from src.orbdet.dynamics_assist import big_asteroid_perturbers
from src.orbdet.frames import ecliptic_to_equatorial
from src.orbdet.mass_determination import (
    _assist_positions,
    _ModelConfig,
    calibrate_sys_floor,
    determine_shared_mass,
)
from src.orbdet.observation import light_time_correct, radec_from_positions
from src.utils.config import load_config


def _synthesize(tobs, mass_msun, cfg, s_c, rng):
    """Devuelve un TargetObservations con RA/Dec sintéticos a *mass_msun* + ruido."""
    gaia = np.asarray(tobs.gaia_bary_icrs, dtype=float)
    obs_jd = np.asarray(tobs.obs_jd_tdb, dtype=float)
    el = tobs.initial_elements

    def bary_ecl_at(jd):
        return _assist_positions(el, mass_msun, np.atleast_1d(jd), cfg)

    jd_ret, _ = light_time_correct(bary_ecl_at, obs_jd, gaia, n_iter=cfg.n_lighttime_iter)
    pos = _assist_positions(el, mass_msun, jd_ret, cfg)
    ra, dec = radec_from_positions(ecliptic_to_equatorial(pos), gaia)

    pa = np.radians(np.asarray(tobs.pa_scan_deg, dtype=float))
    sig = np.asarray(tobs.sigma_al_mas, dtype=float)
    fov = np.asarray(tobs.fov_group)
    # Ruido along-scan: aleatorio por CCD + común por FOV (amplitud s_c).
    al_noise = rng.normal(0.0, sig)
    for g in np.unique(fov):
        idx = fov == g
        al_noise[idx] += rng.normal(0.0, s_c)
    cos_dec = np.cos(np.radians(dec))
    ra_obs = ra + (al_noise * np.sin(pa)) / 3.6e6 / np.maximum(cos_dec, 1e-6)
    dec_obs = dec + (al_noise * np.cos(pa)) / 3.6e6
    return dataclasses.replace(tobs, ra_obs_deg=ra_obs, dec_obs_deg=dec_obs)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--perturber", type=int, default=1)
    ap.add_argument("--release", default="fpr")
    ap.add_argument(
        "--s-c", type=float, default=2.16, help="amplitud del offset común por FOV (mas)"
    )
    ap.add_argument("--seed", type=int, default=12345)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    cfg_top = load_config("config.yaml")
    gaia = cfg_top.sources.gaia_sso
    gaia.release = args.release
    release_cfg = gaia.active()

    pname = _BIG4_NAME_BY_NUMBER[args.perturber]
    specs = _read_targets_from_csv(
        Path("data/output/stage4_validation_summary.csv"), args.perturber
    )
    mid_jd = float(np.mean([s.jd_tdb for s in specs]))
    snapshot = _best_mpcorb_snapshot(_MPCORB_ARCHIVE_DIR, mid_jd)
    elements_map = load_element_rows(snapshot, [s.target for s in specs])
    common_epoch = float(next(iter(elements_map.values()))["epoch_jd"])

    studied = big_asteroid_perturbers(common_epoch, names=(pname,))[0]
    background = big_asteroid_perturbers(common_epoch, exclude=(pname,))
    true_mass = studied.mass_msun
    cfg = _ModelConfig(
        epoch_jd_tdb=common_epoch,
        perturber_elements=studied.elements,
        perturber_name=pname,
        background_perturbers=background,
        perturbers=(),
        integrator="ias15",
        dt_days=1.0,
        n_lighttime_iter=3,
        gm_rel_delta=1e-3,
        backend="assist",
        gr=True,
    )
    print(f"{pname}: masa verdadera inyectada {true_mass * M_SUN_KG:.4e} kg, s_c={args.s_c} mas")

    targets = []
    for spec in specs:
        if spec.target not in elements_map:
            continue
        raw = _fetch_target(gaia.archive_url, spec.target, release_cfg)
        if raw is None:
            continue
        tobs = _build_target_obs(
            raw, elements_map[spec.target], common_epoch, _epoch_ref(release_cfg), background
        )
        targets.append(_synthesize(tobs, true_mass, cfg, args.s_c, rng))

    floor, _ = calibrate_sys_floor(
        targets,
        true_mass,
        [t.initial_elements for t in targets],
        common_epoch,
        perturber_elements=studied.elements,
        perturber_name=pname,
        background_perturbers=background,
        backend="assist",
        gr=True,
    )
    mass, _fit, result = determine_shared_mass(
        targets,
        true_mass,
        studied.elements,
        common_epoch,
        perturber_name=pname,
        background_perturbers=background,
        backend="assist",
        gr=True,
        sys_floor_mas=floor,
        max_iter=40,
    )
    ratio = mass / true_mass
    z = (mass - true_mass) / float(np.sqrt(result.covariance[0, 0]))
    print(
        f"\nRECUPERADO: mass={mass * M_SUN_KG:.4e} kg  ratio={ratio:.4f}  "
        f"z(vs verdad)={z:.2f}  χ²_red={result.chi2_reduced:.3f}  s_c_cal={floor:.2f} mas"
    )
    print("→ ratio≈1 ⇒ motor insesgado; el sobre-tiro real es sistemático de los datos.")


if __name__ == "__main__":
    main()
