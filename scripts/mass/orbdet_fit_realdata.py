"""Ajuste conjunto órbita+masa sobre datos Gaia **reales** vía el motor ``orbdet``.

Este es el script IO que cierra los gates literales de T8 (χ²_red≈1 + rechazo de
outliers sobre datos reales) y T9 (un Big-4 end-to-end) del
``planning/MASS_DETERMINATION_PLAN.md``: carga la astrometría por-tránsito de Gaia
(DR3 o FPR) de cada test-asteroid de un perturbador, la pasa por
:func:`orbdet.gaia_adapter.build_target_observations` y resuelve la masa del
perturbador **conjuntamente** con las órbitas de todos los objetivos vía
:func:`orbdet.mass_determination.determine_shared_mass` (backend ASSIST: DE440 +
GR + 16 perturbadores asteroidales).

El IO vive **fuera** de ``orbdet`` por el contrato de aislamiento: este módulo lee
config/MPCORB y consulta el TAP de Gaia (polars/astroquery), y solo entrega arrays
numpy al motor.

Diseño del ajuste (estrategia Fuentes-Muñoz / OrbFit / JPL)
----------------------------------------------------------
- **Perturbador** (uno de los 4 calibradores Big-4, todos en la efeméride DE441):
  su órbita se toma de la **misma efeméride** que el fondo
  (:func:`orbdet.dynamics_assist.big_asteroid_perturbers`) y se mantiene **fija**;
  solo su masa es parámetro libre. La masa semilla es la de la efeméride (≈ la de
  literatura para el Big-4), de modo que σ(masa) y la información de Fisher leen
  directamente si los datos **constrainen** la masa o no (la pregunta abierta del
  cierre Track A).
- **Objetivos** (test-asteroids pequeños, no en la efeméride): órbita semilla de
  MPCORB, 6 elementos libres por objetivo, todos compartiendo la masa.
- **Época común** del ajuste = época del snapshot MPCORB más cercano al arco de
  datos (buen condicionamiento). El perturbador y el fondo se evalúan en esa época
  desde la efeméride; los objetivos se llevan ahí con el propio N-cuerpos.

Rechazo de outliers
-------------------
Tras converger, se hace sigma-clipping iterativo sobre los residuos blanqueados
(|r| > ``--reject-sigma``) y se re-ajusta, hasta ``--reject-passes`` pasadas o sin
descartes. Es el rechazo estándar de OD que lleva χ²_red hacia 1 cuando hay unos
pocos tránsitos contaminados.

Uso
---
    docker compose run --rm pipeline python -m scripts.mass.orbdet_fit_realdata \\
        --perturber 1 --release fpr \\
        --targets-csv data/output/stage4_validation_summary.csv \\
        --out data/output/orbdet/ceres_fpr.json

    # Los 4 calibradores de una:
    docker compose run --rm pipeline python -m scripts.mass.orbdet_fit_realdata \\
        --perturber big4 --release fpr --out-dir data/output/orbdet
"""

from __future__ import annotations

import argparse
import csv as _csv
import dataclasses
import json
import logging
import math
import os
import time
from pathlib import Path

import numpy as np

from scripts.mass.fit_mass_gaia_loo import (
    _MPCORB_ARCHIVE_DIR,
    _best_mpcorb_snapshot,
    fetch_gaia_full,
    load_element_rows,
)
from scripts.mass.fit_mass_gaia_multitarget import (
    TargetSpec,
    _read_targets_from_csv,
    _read_targets_from_json,
)
from src.orbdet.constants import M_SUN_KG
from src.orbdet.dynamics_assist import BIG_ASTEROIDS, big_asteroid_perturbers
from src.orbdet.gaia_adapter import (
    build_target_observations,
    elements_from_mpcorb,
    propagate_elements,
)
from src.orbdet.mass_determination import (
    TargetObservations,
    calibrate_sys_floor,
    determine_shared_mass,
)
from src.utils.config import GaiaReleaseConfig, load_config

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
logger = logging.getLogger(__name__)

# Nombre de la efeméride (DE441) por número MPC para los 16 perturbadores grandes
# (los únicos cuya órbita expone la efeméride sb441-n16; ver BIG_ASTEROIDS).
_EPHEM_NAME_BY_NUMBER: dict[int, str] = {
    1: "Ceres",
    2: "Pallas",
    3: "Juno",
    4: "Vesta",
    7: "Iris",
    10: "Hygiea",
    15: "Eunomia",
    16: "Psyche",
    31: "Euphrosyne",
    52: "Europa",
    65: "Cybele",
    87: "Sylvia",
    88: "Thisbe",
    107: "Camilla",
    511: "Davida",
    704: "Interamnia",
}
# Calibradores con masa de literatura de alta precisión (subconjunto de validación).
_BIG4_NAME_BY_NUMBER: dict[int, str] = {1: "Ceres", 2: "Pallas", 4: "Vesta", 10: "Hygiea"}

# Mínimo de tránsitos por objetivo para que aporte al ajuste.
_MIN_OBS_PER_TARGET: int = 8


def _ephem_name_for_perturber(number: int, csv_name: str | None) -> str:
    """Nombre tal como lo expone la efeméride DE441 (ASSIST) para *number*.

    Prefiere el mapa explícito de los 16 grandes; si no, intenta el nombre del CSV
    (capitalizado) y verifica que ASSIST lo conozca.
    """
    if number in _EPHEM_NAME_BY_NUMBER:
        return _EPHEM_NAME_BY_NUMBER[number]
    if csv_name:
        cand = csv_name.strip().capitalize()
        if cand in BIG_ASTEROIDS:
            return cand
    raise ValueError(
        f"Perturbador {number} ({csv_name!r}) no está en la efeméride DE441 "
        f"(BIG_ASTEROIDS={BIG_ASTEROIDS}). Este motor requiere la órbita del "
        "perturbador de la efeméride; solo los 16 grandes están soportados."
    )


def _read_targets_from_catalog(
    catalog_path: Path,
    perturber: int,
    *,
    top_n: int,
    max_target_number: int = 100_000,
    max_dist_au: float = 0.05,
) -> list[TargetSpec]:
    """Selecciona los *top_n* objetivos más cercanos al perturbador desde el catálogo.

    Misma lógica que `run_stage4_validation._select_top_targets` pero devuelve
    ``TargetSpec`` directamente desde el catálogo de encuentros (no del CSV de
    validación, que sólo trae 8 por cuerpo y arrastra columnas del método LOO viejo).
    Usar MUCHOS objetivos por perturbador es lo que aprieta σ(masa) y promedia los
    sistemáticos por-encuentro. ``date_utc`` se deriva del ``jd_tdb`` del encuentro.
    """
    import polars as pl
    from astropy.time import Time

    cat = pl.read_parquet(catalog_path, columns=["number_1", "number_2", "dist_au", "jd_tdb"])
    sub = cat.filter((pl.col("number_1") == perturber) | (pl.col("number_2") == perturber))
    sub = sub.with_columns(
        pl.when(pl.col("number_1") == perturber)
        .then(pl.col("number_2"))
        .otherwise(pl.col("number_1"))
        .alias("target")
    )
    sub = sub.filter((pl.col("target") < max_target_number) & (pl.col("dist_au") < max_dist_au))
    sub = sub.sort("dist_au").unique(subset=["target"], keep="first").sort("dist_au").head(top_n)
    out: list[TargetSpec] = []
    for row in sub.iter_rows(named=True):
        date_utc = str(Time(row["jd_tdb"], format="jd", scale="tdb").utc.isot)
        out.append(TargetSpec(target=int(row["target"]), date_utc=date_utc))
    return out


def _read_perturber_meta(csv_path: Path, perturber: int) -> dict:
    """Lee nombre y masa de literatura del perturbador desde el CSV de validación.

    Devuelve ``{"name": str|None, "mass_lit_kg": float|None,
    "mass_lit_sigma_kg": float|None, "source": str|None}`` (campos ausentes → None).
    """
    meta: dict = {"name": None, "mass_lit_kg": None, "mass_lit_sigma_kg": None, "source": None}
    if not csv_path.exists():
        return meta
    with csv_path.open() as fh:
        for row in _csv.DictReader(fh):
            try:
                if int(row.get("perturber", "")) != perturber:
                    continue
            except ValueError:
                continue
            meta["name"] = row.get("perturber_name") or meta["name"]
            for key, col in (
                ("mass_lit_kg", "mass_lit_kg"),
                ("mass_lit_sigma_kg", "mass_lit_sigma_kg"),
            ):
                val = row.get(col)
                if val and meta[key] is None:
                    try:
                        meta[key] = float(val)
                    except ValueError:
                        pass
            meta["source"] = row.get("literature_source") or meta["source"]
            if meta["name"] and meta["mass_lit_kg"] is not None:
                break
    return meta


def _fetch_target(
    archive_url: str, target: int, release_cfg: GaiaReleaseConfig | None
) -> dict[str, np.ndarray] | None:
    """Trae los tránsitos Gaia de *target* y los devuelve como arrays numpy.

    Devuelve ``None`` si el objetivo tiene menos de :data:`_MIN_OBS_PER_TARGET`
    tránsitos utilizables.
    """
    df = fetch_gaia_full(archive_url, target, release_cfg)
    if df.height < _MIN_OBS_PER_TARGET:
        logger.warning(
            "target %d: solo %d tránsitos (<%d) — descartado",
            target,
            df.height,
            _MIN_OBS_PER_TARGET,
        )
        return None
    cols = (
        "epoch",
        "ra",
        "dec",
        "position_angle_scan",
        "ra_error_systematic",
        "dec_error_systematic",
        "ra_dec_correlation_systematic",
        "ra_error_random",
        "dec_error_random",
        "ra_dec_correlation_random",
        "x_gaia",
        "y_gaia",
        "z_gaia",
    )
    return {c: df[c].to_numpy().astype(float) for c in cols}


def _build_target_obs(
    raw: dict[str, np.ndarray],
    mpcorb_row: dict,
    common_epoch_jd_tdb: float,
    epoch_ref_jd_tcb: float,
    background_perturbers: tuple,
) -> TargetObservations:
    """Convierte un objetivo (obs crudas + fila MPCORB) en :class:`TargetObservations`.

    La semilla orbital del objetivo se lleva de la época del snapshot MPCORB a la
    época común del ajuste con el propio N-cuerpos (no-op si coinciden).
    """
    el_mpc = elements_from_mpcorb(
        mpcorb_row["a_au"],
        mpcorb_row["e"],
        mpcorb_row["i_deg"],
        mpcorb_row["Omega_deg"],
        mpcorb_row["omega_deg"],
        mpcorb_row["M_deg"],
    )
    el0 = propagate_elements(
        el_mpc,
        float(mpcorb_row["epoch_jd"]),
        common_epoch_jd_tdb,
        backend="assist",
        asteroid_perturbers=background_perturbers,
    )
    return build_target_observations(
        initial_elements=el0,
        epoch_jd_tdb=common_epoch_jd_tdb,
        epoch_days_tcb=raw["epoch"],
        ra_deg=raw["ra"],
        dec_deg=raw["dec"],
        pa_scan_deg=raw["position_angle_scan"],
        ra_err_sys=raw["ra_error_systematic"],
        dec_err_sys=raw["dec_error_systematic"],
        corr_sys=raw["ra_dec_correlation_systematic"],
        ra_err_rand=raw["ra_error_random"],
        dec_err_rand=raw["dec_error_random"],
        corr_rand=raw["ra_dec_correlation_random"],
        x_gaia=raw["x_gaia"],
        y_gaia=raw["y_gaia"],
        z_gaia=raw["z_gaia"],
        epoch_ref_jd_tcb=epoch_ref_jd_tcb,
    )


def _split_by_targets(values: np.ndarray, targets: list[TargetObservations]) -> list[np.ndarray]:
    """Parte un vector concatenado en orden de objetivos (mismas longitudes)."""
    out: list[np.ndarray] = []
    start = 0
    for t in targets:
        n = int(np.asarray(t.obs_jd_tdb).size)
        out.append(np.asarray(values[start : start + n]))
        start += n
    return out


def _mask_target(tobs: TargetObservations, keep: np.ndarray) -> TargetObservations:
    """Devuelve un :class:`TargetObservations` con solo los tránsitos en *keep*."""
    fov = tobs.fov_group
    return dataclasses.replace(
        tobs,
        obs_jd_tdb=np.asarray(tobs.obs_jd_tdb)[keep],
        ra_obs_deg=np.asarray(tobs.ra_obs_deg)[keep],
        dec_obs_deg=np.asarray(tobs.dec_obs_deg)[keep],
        pa_scan_deg=np.asarray(tobs.pa_scan_deg)[keep],
        sigma_al_mas=np.asarray(tobs.sigma_al_mas)[keep],
        gaia_bary_icrs=np.asarray(tobs.gaia_bary_icrs)[keep],
        fov_group=(None if fov is None else np.asarray(fov)[keep]),
    )


def _fit_with_rejection(
    targets: list[TargetObservations],
    seed_mass_msun: float,
    perturber_elements,
    common_epoch_jd_tdb: float,
    *,
    perturber_name: str,
    background_perturbers: tuple,
    max_iter: int,
    reject_sigma: float,
    reject_passes: int,
    sys_floor_mas: float | None,
    n_workers: int = 1,
):
    """Ajuste conjunto con sigma-clipping iterativo y piso sistemático por FOV.

    ``sys_floor_mas`` fija el piso de error correlacionado intra-tránsito (mas) que
    modela la covarianza en bloques; ``None`` → se **autocalibra** tras la primera
    pasada para que χ²_red ≈ 1 (cuenta la correlación de los ~7 CCDs por cruce, que
    de otro modo subestima σ(masa) ~√7). Devuelve
    ``(mass_msun, fitted_elements, result, targets_final, floor, passes)``.
    """
    cur = list(targets)
    seed = float(seed_mass_msun)
    floor = 0.0 if sys_floor_mas is None else float(sys_floor_mas)
    passes_log: list[dict] = []

    # Pre-paso: si el piso es automático, un ajuste a s_c=0 y se calibra el piso
    # para χ²_red≈1 (cuenta la correlación intra-FOV). Sirve de warm-start del bucle.
    if sys_floor_mas is None:
        t0 = time.time()
        m0, el0, r0 = determine_shared_mass(
            cur,
            seed,
            perturber_elements,
            common_epoch_jd_tdb,
            perturber_name=perturber_name,
            background_perturbers=background_perturbers,
            backend="assist",
            gr=True,
            sys_floor_mas=0.0,
            max_iter=max_iter,
            n_workers=n_workers,
        )
        floor, chi2_at_floor = calibrate_sys_floor(
            cur,
            m0,
            el0,
            common_epoch_jd_tdb,
            perturber_elements=perturber_elements,
            perturber_name=perturber_name,
            background_perturbers=background_perturbers,
            backend="assist",
            gr=True,
            n_workers=n_workers,
        )
        seed = float(m0)
        passes_log.append(
            {
                "pass": -1,
                "n_obs": int(r0.residuals.size),
                "chi2_red": float(r0.chi2_reduced),
                "mass_kg": float(m0 * M_SUN_KG),
                "sys_floor_mas": 0.0,
                "converged": bool(r0.converged),
                "seconds": round(time.time() - t0, 1),
            }
        )
        logger.info(
            "calibración (s_c=0): χ²_red=%.3f mass=%.4e kg → piso s_c=%.3f mas (χ²_red→%.3f)",
            r0.chi2_reduced,
            m0 * M_SUN_KG,
            floor,
            chi2_at_floor,
        )

    last = None
    for p in range(reject_passes + 1):
        t0 = time.time()
        mass_msun, fitted, result = determine_shared_mass(
            cur,
            seed,
            perturber_elements,
            common_epoch_jd_tdb,
            perturber_name=perturber_name,
            background_perturbers=background_perturbers,
            backend="assist",
            gr=True,
            sys_floor_mas=floor,
            max_iter=max_iter,
            n_workers=n_workers,
        )
        n_obs = int(result.residuals.size)
        passes_log.append(
            {
                "pass": p,
                "n_obs": n_obs,
                "chi2_red": float(result.chi2_reduced),
                "mass_kg": float(mass_msun * M_SUN_KG),
                "sys_floor_mas": float(floor),
                "converged": bool(result.converged),
                "seconds": round(time.time() - t0, 1),
            }
        )
        logger.info(
            "pass %d: n_obs=%d χ²_red=%.3f mass=%.4e kg s_c=%.2f conv=%s (%.0fs)",
            p,
            n_obs,
            result.chi2_reduced,
            mass_msun * M_SUN_KG,
            floor,
            result.converged,
            passes_log[-1]["seconds"],
        )
        last = (mass_msun, fitted, result, list(cur))
        if p == reject_passes or reject_sigma <= 0:
            break
        # Sigma-clip: descarta |r| > reject_sigma y re-ajusta con la órbita ya hallada.
        per_target = _split_by_targets(np.abs(result.residuals), cur)
        new_targets: list[TargetObservations] = []
        n_dropped = 0
        for tobs, el_fit, r_abs in zip(cur, fitted, per_target):
            keep = r_abs <= reject_sigma
            n_dropped += int((~keep).sum())
            if int(keep.sum()) < _MIN_OBS_PER_TARGET:
                # No descartar por debajo del mínimo: conservar el objetivo completo.
                new_targets.append(dataclasses.replace(tobs, initial_elements=el_fit))
                continue
            new_targets.append(
                _mask_target(dataclasses.replace(tobs, initial_elements=el_fit), keep)
            )
        if n_dropped == 0:
            logger.info("pass %d: sin outliers >%.1fσ — fin", p, reject_sigma)
            break
        logger.info("pass %d: descartados %d tránsitos >%.1fσ", p, n_dropped, reject_sigma)
        cur = new_targets
        seed = float(mass_msun)
    assert last is not None
    return (*last, floor, passes_log)


def _run_perturber(
    perturber: int,
    *,
    args: argparse.Namespace,
    archive_url: str,
    release_cfg: GaiaReleaseConfig | None,
    release_name: str,
) -> dict:
    """Corre el ajuste conjunto para un perturbador y devuelve el dict de resultado."""
    meta = _read_perturber_meta(args.targets_csv, perturber) if args.targets_csv else {}
    pname = _ephem_name_for_perturber(perturber, meta.get("name"))
    lit_kg = args.lit_mass_kg if args.lit_mass_kg is not None else meta.get("mass_lit_kg")
    lit_sigma_kg = meta.get("mass_lit_sigma_kg")

    if args.from_catalog is not None:
        specs: list[TargetSpec] = _read_targets_from_catalog(
            args.from_catalog,
            perturber,
            top_n=args.top_per_perturber,
            max_dist_au=args.max_dist_au,
        )
    elif args.targets_json is not None:
        specs = _read_targets_from_json(args.targets_json)
    else:
        specs = _read_targets_from_csv(args.targets_csv, perturber)
    if args.max_targets:
        specs = specs[: args.max_targets]
    if not specs:
        raise ValueError(f"Sin objetivos para el perturbador {perturber}")
    logger.info("Perturbador %d (%s): %d objetivos candidatos", perturber, pname, len(specs))

    # Época común = snapshot MPCORB más cercano al arco de encuentros.
    mid_jd = float(np.mean([s.jd_tdb for s in specs]))
    snapshot = args.mpcorb or _best_mpcorb_snapshot(_MPCORB_ARCHIVE_DIR, mid_jd)
    logger.info("Snapshot MPCORB: %s", Path(snapshot).name)
    target_numbers = [s.target for s in specs]
    elements_map = load_element_rows(snapshot, target_numbers)
    common_epoch = float(next(iter(elements_map.values()))["epoch_jd"])
    logger.info("Época común del ajuste (JD TDB): %.5f", common_epoch)

    # Perturbador (órbita + masa-semilla) y fondo, desde la efeméride DE441.
    studied = big_asteroid_perturbers(common_epoch, names=(pname,))[0]
    perturber_elements = studied.elements
    seed_mass_msun = (
        (args.seed_mass_kg / M_SUN_KG) if args.seed_mass_kg is not None else studied.mass_msun
    )
    background = big_asteroid_perturbers(common_epoch, exclude=(pname,))
    logger.info(
        "Fondo: %d perturbadores asteroidales; masa-semilla %.4e kg",
        len(background),
        seed_mass_msun * M_SUN_KG,
    )

    targets: list[TargetObservations] = []
    used_numbers: list[int] = []
    for spec in specs:
        if spec.target not in elements_map:
            logger.warning("target %d no está en el snapshot MPCORB — saltado", spec.target)
            continue
        raw = _fetch_target(archive_url, spec.target, release_cfg)
        if raw is None:
            continue
        tobs = _build_target_obs(
            raw, elements_map[spec.target], common_epoch, _epoch_ref(release_cfg), background
        )
        targets.append(tobs)
        used_numbers.append(spec.target)
        logger.info("target %d: %d tránsitos", spec.target, int(tobs.obs_jd_tdb.size))

    if not targets:
        raise ValueError(f"Sin objetivos con datos para el perturbador {perturber}")

    mass_msun, fitted, result, targets_final, sys_floor, passes_log = _fit_with_rejection(
        targets,
        seed_mass_msun,
        perturber_elements,
        common_epoch,
        perturber_name=pname,
        background_perturbers=background,
        max_iter=args.max_iter,
        reject_sigma=args.reject_sigma,
        reject_passes=args.reject_passes,
        sys_floor_mas=args.sys_floor,
        n_workers=args.workers,
    )

    mass_kg = float(mass_msun * M_SUN_KG)
    var_mass = (
        float(result.covariance[0, 0]) if np.all(np.isfinite(result.covariance)) else math.nan
    )
    sigma_mass_msun = math.sqrt(var_mass) if var_mass > 0 else math.nan
    sigma_kg = float(sigma_mass_msun * M_SUN_KG)
    ratio = (mass_kg / lit_kg) if lit_kg else None
    z = None
    if lit_kg and math.isfinite(sigma_kg):
        denom = math.sqrt(sigma_kg**2 + (lit_sigma_kg or 0.0) ** 2)
        z = (mass_kg - lit_kg) / denom if denom > 0 else None

    out = {
        "perturber": perturber,
        "perturber_name": pname,
        "release": release_name,
        "common_epoch_jd_tdb": common_epoch,
        "mpcorb_snapshot": Path(snapshot).name,
        "n_targets": len(targets_final),
        "target_numbers": used_numbers,
        "n_obs_final": int(result.residuals.size),
        "seed_mass_kg": float(seed_mass_msun * M_SUN_KG),
        "mass_fit_kg": mass_kg,
        "mass_fit_sigma_kg": sigma_kg,
        "mass_fit_msun": float(mass_msun),
        "chi2": float(result.chi2),
        "dof": int(result.dof),
        "chi2_red": float(result.chi2_reduced),
        "converged": bool(result.converged),
        "n_iter": int(result.n_iter),
        "sys_floor_mas": float(sys_floor),
        "mass_lit_kg": lit_kg,
        "mass_lit_sigma_kg": lit_sigma_kg,
        "literature_source": meta.get("source"),
        "ratio_fit_over_lit": ratio,
        "z_score": z,
        "reject_sigma": args.reject_sigma,
        "passes": passes_log,
    }
    return out


def _epoch_ref(release_cfg: GaiaReleaseConfig | None) -> float:
    """Época de referencia TCB del release (J2010 por defecto / DR3 legacy)."""
    from src.orbdet.time_scales import J2010_TCB_JD

    return release_cfg.epoch_ref_jd_tcb if release_cfg is not None else J2010_TCB_JD


def _report(out: dict) -> None:
    print(f"\n=== orbdet joint mass fit — {out['perturber_name']} ({out['release']}) ===")
    print(f"  objetivos:        {out['n_targets']}  ({out['n_obs_final']} tránsitos)")
    print(f"  masa ajustada:    {out['mass_fit_kg']:.4e} ± {out['mass_fit_sigma_kg']:.2e} kg")
    if out["mass_lit_kg"]:
        print(
            f"  literatura:       {out['mass_lit_kg']:.4e} kg "
            f"({out.get('literature_source') or '?'})"
        )
        if out["ratio_fit_over_lit"] is not None:
            print(f"  ratio fit/lit:    {out['ratio_fit_over_lit']:.3f}")
        if out["z_score"] is not None:
            print(f"  z-score:          {out['z_score']:.2f}   (gate |z|<3)")
    print(f"  χ²_red:           {out['chi2_red']:.3f}   (gate ≈1)")
    print(f"  piso sist. s_c:   {out['sys_floor_mas']:.3f} mas (covarianza en bloques por FOV)")
    print(f"  convergió:        {out['converged']} en {out['n_iter']} iter")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument(
        "--perturber",
        required=True,
        help="número MPC del perturbador, lista separada por comas, o 'big4' "
        "(=1,2,4,10: Ceres/Pallas/Vesta/Hygiea)",
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--targets-csv",
        type=Path,
        default=Path("data/output/stage4_validation_summary.csv"),
    )
    group.add_argument("--targets-json", type=Path)
    parser.add_argument(
        "--from-catalog",
        type=Path,
        default=None,
        help="selecciona los objetivos más cercanos directamente del catálogo de "
        "encuentros (parquet) en vez del CSV — permite usar muchos más por perturbador",
    )
    parser.add_argument(
        "--top-per-perturber",
        type=int,
        default=40,
        help="nº de objetivos más cercanos a usar con --from-catalog",
    )
    parser.add_argument(
        "--max-dist-au",
        type=float,
        default=0.05,
        help="distancia máxima de encuentro (--from-catalog)",
    )
    parser.add_argument("--release", default=None, help="'dr3' | 'fpr' (default: el del config)")
    parser.add_argument(
        "--lit-mass-kg", type=float, default=None, help="override de masa de literatura"
    )
    parser.add_argument("--seed-mass-kg", type=float, default=None, help="override de masa-semilla")
    parser.add_argument("--mpcorb", type=Path, default=None, help="snapshot MPCORB explícito")
    parser.add_argument(
        "--max-targets", type=int, default=None, help="limita nº de objetivos (debug)"
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=max(1, (os.cpu_count() or 2) - 2),
        help="procesos paralelos para evaluar objetivos (default: nº de cores − 2)",
    )
    parser.add_argument("--max-iter", type=int, default=40)
    parser.add_argument("--reject-sigma", type=float, default=4.0, help="umbral sigma-clip (0=off)")
    parser.add_argument("--reject-passes", type=int, default=2)
    parser.add_argument(
        "--sys-floor",
        type=float,
        default=None,
        help="piso de error sistemático correlacionado intra-FOV (mas); "
        "omitir → autocalibrar para χ²_red≈1 (recomendado)",
    )
    parser.add_argument("--out", type=Path, default=None, help="JSON de salida (perturbador único)")
    parser.add_argument(
        "--out-dir", type=Path, default=None, help="directorio de salida (varios perturbadores)"
    )
    args = parser.parse_args()

    if args.targets_json is not None:
        args.targets_csv = None

    cfg = load_config(args.config)
    gaia = cfg.sources.gaia_sso
    if args.release is not None:
        gaia.release = args.release
    release_cfg = gaia.active()
    archive_url = gaia.archive_url
    logger.info("Gaia release: %s (tabla %s)", gaia.release, release_cfg.table)

    if args.perturber.lower() == "big4":
        perturbers = [1, 2, 4, 10]
    else:
        perturbers = [int(p) for p in args.perturber.split(",") if p.strip()]

    results: list[dict] = []
    for pert in perturbers:
        try:
            out = _run_perturber(
                pert,
                args=args,
                archive_url=archive_url,
                release_cfg=release_cfg,
                release_name=gaia.release,
            )
        except Exception as exc:  # noqa: BLE001 — un perturbador no debe tumbar el resto
            logger.exception("Perturbador %d falló: %s", pert, exc)
            continue
        results.append(out)
        _report(out)

        if args.out_dir is not None:
            args.out_dir.mkdir(parents=True, exist_ok=True)
            dest = args.out_dir / f"{out['perturber_name'].lower()}_{gaia.release}.json"
            dest.write_text(json.dumps(out, indent=2))
            logger.info("Escrito %s", dest)
        elif args.out is not None and len(perturbers) == 1:
            args.out.parent.mkdir(parents=True, exist_ok=True)
            args.out.write_text(json.dumps(out, indent=2))
            logger.info("Escrito %s", args.out)

    if not results:
        logger.error("Ningún perturbador produjo resultado")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
