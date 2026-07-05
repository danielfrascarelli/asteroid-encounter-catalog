"""Injection-recovery end-to-end de la capa de detección (tribunal 2026-07-04, M10/B4).

Genera pares sintéticos de asteroides con encuentros de geometría controlada y
verifica que ``detect_encounters`` (scan KD-tree grueso + refinamiento Kepler, con
los parámetros de producción) los recupera con la distancia y la época correctas.
Es el test que habría detectado B1 (ventana de refinamiento menor que el semipaso
grueso): los mínimos se colocan con **fase uniforme respecto de la grilla gruesa**,
así que cualquier recorte de ventana aparece como sesgo en ``d`` y ``t``.

Construcción de cada par
------------------------
1. Punto de encuentro ``P`` a radio r ∈ [1.6, 3.5] AU y época t* uniforme en la
   ventana (⇒ fase uniforme vs la grilla de 12 h).
2. Cuerpo 1 en ``P`` con velocidad ~circular perturbada (excentricidades hasta
   ~0.4). Cuerpo 2 en ``P + b`` con ``v₂ = v₁ + Δv``; |Δv| log-uniforme en
   [0.5, v_max] km/s y ``b ⊥ Δv`` con |b| log-uniforme en [10⁻⁴, 0.045] AU
   (⇒ el mínimo real ≈ |b| cerca de t*).
3. Estados → elementos con ``state_to_elements`` (mismo GM que el propagador).
4. Verdad de terreno: scan denso a 5 s (±0.5 d alrededor de t*), independiente
   del refinador.

Gates (los del plan de remediación, Tarea 5/18)
-----------------------------------------------
- Recuperación ≥ 99 % de los pares inyectados con d_min real < 0.045 AU.
- |d_pipeline − d_true| ≤ max(1 μAU, 10⁻³·d_true) por par (ratio d ≈ 1).
- |t_pipeline − t_true| ≤ paso fino.

Uso
---
    docker compose run --rm pipeline python -m scripts.validate.injection_recovery_detection \
        --n-pairs 200 --seed 42
"""

from __future__ import annotations

import argparse
import logging

import numpy as np
import polars as pl

from src.detect.pipeline import detect_encounters
from src.orbdet.constants import GM_SUN
from src.orbdet.kepler import state_to_elements
from src.propagate.grid import make_time_grid
from src.propagate.kepler import kepler_to_cartesian

logger = logging.getLogger(__name__)

_DEG = 180.0 / np.pi
_KMS_TO_AUD = 86400.0 / 149_597_870.7  # km/s → AU/día

# Parámetros de producción (config.yaml)
_COARSE_STEP_HOURS = 12.0
_FINE_STEP_SECONDS = 120.0
_WINDOW_HOURS = 6.0
_THRESHOLD_AU = 0.05
_MAX_REL_VEL_KM_S = 25.0  # cota usada para ensanchar el query radius

_T0 = 2457000.0  # JD TDB, inicio de la ventana sintética


def _random_unit(rng: np.random.Generator) -> np.ndarray:
    v = rng.normal(size=3)
    return v / np.linalg.norm(v)


def _make_pair(
    rng: np.random.Generator, t_star: float, v_rel_max_km_s: float
) -> tuple[dict, dict, float]:
    """Construye un par sintético con encuentro ~|b| cerca de t_star.

    Returns
    -------
    (elem_row_1, elem_row_2, b_au)
    """
    while True:
        r_mag = rng.uniform(1.6, 3.5)
        # Dirección de r con inclinación moderada (belt-like)
        u_r = _random_unit(rng)
        u_r[2] *= 0.3
        u_r /= np.linalg.norm(u_r)
        r_vec = r_mag * u_r

        # v1: ~circular con perturbación (e hasta ~0.4)
        v_circ = np.sqrt(GM_SUN / r_mag)
        u_t = np.cross(np.array([0.0, 0.0, 1.0]), u_r)
        u_t /= np.linalg.norm(u_t)
        v1 = v_circ * rng.uniform(0.9, 1.1) * u_t
        v1 += v_circ * 0.15 * rng.normal(size=3) * np.array([1.0, 1.0, 0.5])

        # Δv y miss vector b ⊥ Δv
        dv_kms = np.exp(rng.uniform(np.log(0.5), np.log(v_rel_max_km_s)))
        u_dv = _random_unit(rng)
        v2 = v1 + dv_kms * _KMS_TO_AUD * u_dv
        b_mag = np.exp(rng.uniform(np.log(1e-4), np.log(0.045)))
        u_b = _random_unit(rng)
        u_b -= np.dot(u_b, u_dv) * u_dv  # ⊥ Δv → mínimo ≈ |b| en t*
        u_b /= np.linalg.norm(u_b)
        r2_vec = r_vec + b_mag * u_b

        try:
            el1 = state_to_elements(r_vec, v1)
            el2 = state_to_elements(r2_vec, v2)
        except (ValueError, ZeroDivisionError):
            continue
        # Órbitas acotadas y belt-like (evitar hipérbolas y e extremas)
        if not (1.2 < el1.a < 5.0 and 1.2 < el2.a < 5.0):
            continue
        if el1.e >= 0.9 or el2.e >= 0.9:
            continue

        def _row(el, num: int) -> dict:
            return {
                "number": num,
                "designation": f"SYN{num}",
                "a_au": el.a,
                "e": el.e,
                "i_deg": el.i * _DEG,
                "Omega_deg": el.Omega * _DEG,
                "omega_deg": el.omega * _DEG,
                "M_deg": el.M * _DEG,
                "epoch_jd": t_star,
            }

        return _row(el1, 0), _row(el2, 0), float(b_mag)


def _dense_minimum(elems: pl.DataFrame, i: int, j: int, t_center: float) -> tuple[float, float]:
    """Mínimo por scan denso a 5 s en ±0.5 d — verdad de terreno sin refinador."""
    t = np.arange(t_center - 0.5, t_center + 0.5, 5.0 / 86400.0)

    def _pos(k: int) -> np.ndarray:
        row = elems.row(k, named=True)
        return kepler_to_cartesian(
            a_au=row["a_au"],
            e=row["e"],
            i_rad=np.radians(row["i_deg"]),
            Omega_rad=np.radians(row["Omega_deg"]),
            omega_rad=np.radians(row["omega_deg"]),
            M0_rad=np.radians(row["M_deg"]),
            epoch_jd=row["epoch_jd"],
            t_jd=t,
        )

    d = np.linalg.norm(_pos(i) - _pos(j), axis=1)
    k = int(np.argmin(d))
    return float(t[k]), float(d[k])


def run_injection_recovery(
    n_pairs: int = 200,
    seed: int = 42,
    window_days: float = 30.0,
    v_rel_max_km_s: float = _MAX_REL_VEL_KM_S,
) -> dict:
    """Corre la inyección-recuperación; devuelve el resumen con los gates."""
    rng = np.random.default_rng(seed)

    rows: list[dict] = []
    truths: list[tuple[int, int, float, float]] = []  # (idx1, idx2, t_true, d_true)
    for _ in range(n_pairs):
        t_star = rng.uniform(_T0 + 2.0, _T0 + window_days - 2.0)
        r1, r2, _b = _make_pair(rng, t_star, v_rel_max_km_s)
        idx1 = len(rows)
        r1["number"], r2["number"] = idx1 + 1, idx1 + 2  # numbers 1-based únicos
        rows.extend([r1, r2])
        truths.append((idx1, idx1 + 1, t_star, np.nan))

    elems = pl.DataFrame(rows).with_columns(pl.col("number").cast(pl.Int32))

    # Verdad de terreno (independiente del refinador)
    truths = [(i, j, *_dense_minimum(elems, i, j, t_star)) for (i, j, t_star, _) in truths]
    injected = [(i, j, t, d) for (i, j, t, d) in truths if d < 0.045]
    logger.info(
        "%d pares inyectados, %d con d_min real < 0.045 AU (evaluables)",
        n_pairs,
        len(injected),
    )

    grid = make_time_grid(_T0, _T0 + window_days, step_hours=_COARSE_STEP_HOURS)
    query_radius = _THRESHOLD_AU + (
        v_rel_max_km_s * _KMS_TO_AUD * (_COARSE_STEP_HOURS / 24.0) / 2.0
    )
    result = detect_encounters(
        elems,
        grid,
        threshold_au=_THRESHOLD_AU,
        semimajor_diff_max_au=0.5,
        inclination_diff_max_deg=30.0,
        leaf_size=30,
        fine_step_seconds=_FINE_STEP_SECONDS,
        window_hours=_WINDOW_HOURS,
        prefilter_enabled=False,  # pares todos-contra-todos, como la corrida grande
        refinement_enabled=True,
        n_workers=1,
        chunk_size_days=30.0,
        query_radius_au=query_radius,
    )
    by_pair = {
        (min(a, b), max(a, b)): (t, d)
        for a, b, t, d in zip(
            result["number_1"].to_list(),
            result["number_2"].to_list(),
            result["jd_tdb"].to_list(),
            result["dist_au"].to_list(),
        )
    }

    n_found = 0
    n_d_fail = 0
    d_errors: list[float] = []
    t_errors_s: list[float] = []
    misses: list[tuple[int, int, float]] = []
    for i, j, t_true, d_true in injected:
        key = (i + 1, j + 1)
        if key not in by_pair:
            misses.append((i + 1, j + 1, d_true))
            continue
        n_found += 1
        t_pipe, d_pipe = by_pair[key]
        err = abs(d_pipe - d_true)
        d_errors.append(err)
        if err > max(1e-6, 1e-3 * d_true):
            n_d_fail += 1
        t_errors_s.append(abs(t_pipe - t_true) * 86400.0)

    recovery = n_found / len(injected) if injected else float("nan")
    d_ok = n_found > 0 and n_d_fail == 0
    t_ok = (max(t_errors_s) <= _FINE_STEP_SECONDS) if t_errors_s else False

    summary = {
        "n_pairs": n_pairs,
        "n_evaluable": len(injected),
        "n_recovered": n_found,
        "recovery_frac": recovery,
        "max_abs_d_error_au": max(d_errors) if d_errors else None,
        "median_abs_d_error_au": float(np.median(d_errors)) if d_errors else None,
        "max_abs_t_error_s": max(t_errors_s) if t_errors_s else None,
        "misses": misses,
        "gate_recovery_ge_99pct": recovery >= 0.99,
        "gate_distance_ratio_1": d_ok,
        "gate_epoch_within_fine_step": t_ok,
    }
    return summary


def main() -> int:
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s"
    )
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--n-pairs", type=int, default=200)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--window-days", type=float, default=30.0)
    ap.add_argument("--v-rel-max-km-s", type=float, default=_MAX_REL_VEL_KM_S)
    args = ap.parse_args()

    s = run_injection_recovery(
        n_pairs=args.n_pairs,
        seed=args.seed,
        window_days=args.window_days,
        v_rel_max_km_s=args.v_rel_max_km_s,
    )
    print()
    print(f"pares inyectados evaluables : {s['n_evaluable']} / {s['n_pairs']}")
    print(f"recuperados                 : {s['n_recovered']} ({s['recovery_frac'] * 100:.1f} %)")
    print(
        f"|Δd| máx / mediana          : {s['max_abs_d_error_au']:.2e} / {s['median_abs_d_error_au']:.2e} AU"
    )
    print(
        f"|Δt| máx                    : {s['max_abs_t_error_s']:.1f} s (paso fino {_FINE_STEP_SECONDS:.0f} s)"
    )
    if s["misses"]:
        print(f"perdidos: {s['misses']}")
    ok = (
        s["gate_recovery_ge_99pct"]
        and s["gate_distance_ratio_1"]
        and s["gate_epoch_within_fine_step"]
    )
    print(
        f"gates: recovery≥99% {'OK' if s['gate_recovery_ge_99pct'] else 'FAIL'} | "
        f"d-ratio≈1 {'OK' if s['gate_distance_ratio_1'] else 'FAIL'} | "
        f"t≤fino {'OK' if s['gate_epoch_within_fine_step'] else 'FAIL'} → "
        f"{'PASS' if ok else 'FAIL'}"
    )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
