"""Construye el catálogo de masas a partir de los ajustes por-perturbador.

Toma los JSON que produce ``orbdet_fit_realdata`` (uno por perturbador) y arma un
catálogo con el **modelo de error correcto para una medición limitada por
sistemáticos**.

Por qué hace falta un piso sistemático
--------------------------------------
Con muchos objetivos la σ formal (Fisher, de la covarianza) baja como 1/√N y se
vuelve diminuta (<1%), pero la **exactitud real está limitada por sistemáticos
por-encuentro** (imperfección de la órbita del objetivo, sistemáticos astrométricos
locales, perturbadores menores) que NO se reducen apilando. Se ve en que los
calibradores Big-4 —cuya masa verdadera se conoce (DAWN/Goffin/Vernazza)— quedan a
unos pocos % de la verdad aunque σ_formal sea ~0.2%.

El piso ``f_sys`` se **calibra con los Big-4**: es la dispersión RMS de
``masa_fit/masa_lit − 1`` sobre los calibradores. Cada masa se reporta entonces con
``σ_total = √(σ_stat² + (f_sys·|masa|)²)`` y un z basado en σ_total. Es el tratamiento
estándar de "incertidumbre externa" cuando la formal subestima por sistemáticos no
modelados.

Diagnósticos de validez de σ_jack (tribunal 2026-07-04, B6)
-----------------------------------------------------------
Con N réplicas jackknife, σ_jack puede estar dominada por una sola réplica (leverage
top-1 de hasta el 92 % de la varianza en este dataset → ~1 grado de libertad efectivo).
El catálogo reporta por fila:

- ``jack_leverage_top1``: fracción de la varianza jackknife aportada por la réplica
  más influyente (``max((m_i − m̄)²) / Σ(m_i − m̄)²``).
- ``sigma_jack_excl_top1_kg`` / ``snr_jack_excl_top1``: σ y SNR recalculados sin esa
  réplica — diagnóstico de estabilidad de la clasificación.
- ``sigma_jack_defensible``: True sólo si ``n_targets ≥ --min-n-jack`` (default 10) y
  ``jack_leverage_top1 ≤ 0.5``. Una masa sólo se declara ``measured`` si su σ_jack es
  defendible; masas ajustadas ≤ 0 se marcan ``non_physical``.

Uso
---
    python -m scripts.mass.build_mass_catalog \\
        --in-dir data/output/orbdet/expanded \\
        --out data/output/orbdet/mass_catalog.csv
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

from src.orbdet.identifiability import (
    delta_chi2_quadratic,
    false_alarm_probability,
    threshold_for_nsigma,
)

_BIG4 = {1, 2, 4, 10}


def _load(in_dir: Path) -> list[dict]:
    rows = []
    for p in sorted(in_dir.glob("*.json")):
        try:
            d = json.loads(p.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        if "mass_fit_kg" in d and "perturber" in d:
            rows.append(d)
    return rows


def _jackknife_diagnostics(
    jack_masses: list[float] | None,
) -> tuple[float | None, float | None]:
    """Leverage top-1 y σ_jack sin la réplica dominante (B6).

    ``σ_jack = √((n−1)/n · Σ(m_i − m̄)²)`` — la misma fórmula del fitter (verificada
    contra ``mass_fit_sigma_jack_kg`` de los JSON). El leverage top-1 es la fracción
    de ``Σ(m_i − m̄)²`` que aporta la réplica más desviada; con leverage > 0.5 la
    σ_jack tiene ~1 grado de libertad efectivo y no es defendible como σ.

    Returns
    -------
    (leverage_top1, sigma_jack_excl_top1_kg) — ``(None, None)`` si hay < 3 réplicas.
    """
    if not jack_masses or len(jack_masses) < 3:
        return None, None
    n = len(jack_masses)
    mean = sum(jack_masses) / n
    devs2 = [(x - mean) ** 2 for x in jack_masses]
    ss = sum(devs2)
    if ss <= 0.0:
        return 0.0, 0.0
    k_top = devs2.index(max(devs2))
    leverage = devs2[k_top] / ss
    rest = [x for i, x in enumerate(jack_masses) if i != k_top]
    m = len(rest)
    mean_r = sum(rest) / m
    ss_r = sum((x - mean_r) ** 2 for x in rest)
    sigma_excl = math.sqrt((m - 1) / m * ss_r)
    return leverage, sigma_excl


def _calibrate_f_sys(rows: list[dict], min_n: int) -> tuple[float, float, int]:
    """Piso sistemático fraccional desde los calibradores Big-4 **bien muestreados**.

    Sólo usa calibradores con ``n_targets >= min_n``: con pocos objetivos la
    dispersión por-encuentro domina (p. ej. Pallas con N=6 se desvía +24% por
    estadística de muestra chica, no por un sistemático) y contaminaría el piso. Con
    N≳20 la estimación ya promedió esa dispersión y mide la **exactitud real**.

    Devuelve ``(f_sys, mean_offset, n_cal)``: dispersión RMS y sesgo medio de
    ``ratio − 1`` sobre los calibradores bien muestreados.
    """
    devs = [
        d["ratio_fit_over_lit"] - 1.0
        for d in rows
        if d.get("perturber") in _BIG4
        and d.get("ratio_fit_over_lit")
        and (d.get("n_targets") or 0) >= min_n
    ]
    if not devs:
        return 0.0, 0.0, 0
    mean = sum(devs) / len(devs)
    rms = math.sqrt(sum(x * x for x in devs) / len(devs))
    return rms, mean, len(devs)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--in-dir", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument(
        "--f-sys",
        type=float,
        default=None,
        help="piso sistemático fraccional (default: calibrado de los Big-4 bien muestreados)",
    )
    ap.add_argument(
        "--min-n-cal",
        type=int,
        default=15,
        help="nº mínimo de objetivos para que un calibrador defina el piso f_sys",
    )
    ap.add_argument(
        "--min-n-reliable",
        type=int,
        default=15,
        help="nº mínimo de objetivos para marcar una masa como fiable (no target-limited)",
    )
    ap.add_argument(
        "--min-snr-jack",
        type=float,
        default=3.0,
        help="SNR de identificabilidad mínimo (masa/σ_jack) para declarar la masa "
        "como medida en vez de cota (F2). Requiere --jackknife en el ajuste",
    )
    ap.add_argument(
        "--identif-nsigma",
        type=float,
        default=3.0,
        help="significancia (en σ) del criterio de identificabilidad por verosimilitud "
        "perfilada (M13/T21): umbral Δχ²(M=0) = nσ² (default 3σ → Δχ²=9). Se computa "
        "en paralelo al snr_jack como criterio alternativo",
    )
    ap.add_argument(
        "--min-n-jack",
        type=int,
        default=10,
        help="nº mínimo de réplicas jackknife para que σ_jack sea defendible como σ "
        "(B6: con N chico y leverage alto, σ_jack tiene ~1 gdl efectivo)",
    )
    ap.add_argument(
        "--max-leverage",
        type=float,
        default=0.5,
        help="leverage top-1 máximo para que σ_jack sea defendible (B6)",
    )
    args = ap.parse_args()

    rows = _load(args.in_dir)
    if not rows:
        print(f"Sin JSON de masas en {args.in_dir}")
        return 1

    f_sys_cal, mean_off, n_cal = _calibrate_f_sys(rows, args.min_n_cal)
    f_sys = args.f_sys if args.f_sys is not None else f_sys_cal
    print(
        f"Piso sistemático f_sys = {f_sys:.3f} ({f_sys * 100:.1f}%)"
        + (
            f" — calibrado de {n_cal} Big-4 (sesgo medio {mean_off * 100:+.1f}%)"
            if args.f_sys is None
            else " (manual)"
        )
    )

    delta_chi2_thr = threshold_for_nsigma(args.identif_nsigma)

    out_rows = []
    for d in sorted(rows, key=lambda x: x["perturber"]):
        m = d["mass_fit_kg"]
        s_stat = d.get("mass_fit_sigma_kg") or 0.0
        # B6 — σ leverage-robusta. Cuando el jackknife está dominado por una réplica
        # (leverage top-1 > max_leverage) y hay σ bootstrap disponible, se usa la mayor
        # de ambas como σ estadística: el bootstrap no sufre el ~1 gdl efectivo del
        # jackknife (un objetivo de alto leverage no entra en toda remuestra).
        lev_top1, s_jack_excl = _jackknife_diagnostics(d.get("jackknife_masses_kg"))
        s_boot = d.get("mass_fit_sigma_boot_kg")
        use_boot = (
            s_boot is not None
            and s_boot > 0
            and lev_top1 is not None
            and lev_top1 > args.max_leverage
        )
        if use_boot:
            s_stat = max(s_stat, s_boot)
        # |m|: con masa ajustada negativa (no física, p. ej. Davida N=3) el piso
        # proporcional debe seguir siendo una σ positiva (bug de signo, B6).
        s_sys = f_sys * abs(m)
        s_tot = math.sqrt(s_stat**2 + s_sys**2)
        # Referencia: masa de literatura si la hay (Big-4); si no, la masa semilla de
        # la efeméride DE441 (también un valor publicado) para los demás perturbadores.
        lit = d.get("mass_lit_kg")
        ref = lit if lit else d.get("seed_mass_kg")
        ref_src = d.get("literature_source") if lit else "DE441 ephemeris (seed)"
        ratio = (m / ref) if ref else None
        s_lit = d.get("mass_lit_sigma_kg") or 0.0
        z_tot = None
        if ref:
            denom = math.sqrt(s_tot**2 + s_lit**2)
            z_tot = (m - ref) / denom if denom > 0 else None

        # F2 — identificabilidad por jackknife. La masa es una determinación genuina
        # sólo si supera su σ externa (que captura la regresión masa↔órbita); si no,
        # la deflexión queda bajo el ruido por-encuentro y se reporta como cota.
        # B6 — σ_jack sólo es defendible con suficientes réplicas y sin una réplica
        # dominante; masas ≤ 0 son no físicas y nunca "measured".
        s_formal = d.get("mass_fit_sigma_formal_kg")
        s_jack = d.get("mass_fit_sigma_jack_kg")
        snr_jack = (m / s_jack) if (s_jack and s_jack > 0) else None
        snr_boot = (m / s_boot) if (use_boot and s_boot and s_boot > 0) else None

        # M13/T21 — criterio alternativo por verosimilitud perfilada. Δχ²(M=0) es la
        # curvatura del χ² perfilado sobre la órbita; bajo la aproximación cuadrática
        # (exacta en el límite lineal-Gaussiano) Δχ² = (M̂/σ_formal)², computable
        # directo de los JSON sin re-fits. Umbral Δχ² = nσ² (9 ≈ 3σ). Se reporta en
        # paralelo al snr_jack (no lo reemplaza). σ_formal es el denominador correcto
        # de la curvatura (a diferencia de σ_jack, que sufre el ~1 gdl efectivo).
        delta_chi2_m0 = (
            delta_chi2_quadratic(m, s_formal) if (m > 0 and s_formal and s_formal > 0) else None
        )
        p_false_alarm = (
            false_alarm_probability(delta_chi2_m0) if delta_chi2_m0 is not None else None
        )
        identifiable_profile = delta_chi2_m0 is not None and delta_chi2_m0 >= delta_chi2_thr
        snr_jack_excl = (m / s_jack_excl) if (s_jack_excl and s_jack_excl > 0) else None
        n_targets = d.get("n_targets") or 0
        sigma_jack_defensible = (
            s_jack is not None
            and n_targets >= args.min_n_jack
            and lev_top1 is not None
            and lev_top1 <= args.max_leverage
        )
        # "measured" exige estabilidad de la clasificación (gate B6): el SNR debe
        # superar el umbral tanto con σ_jack completa como sin la réplica dominante,
        # con N mínimo de réplicas y masa física. La defensibilidad de σ_jack como
        # barra publicable es ortogonal y va en `sigma_jack_defensible` (leverage).
        if m <= 0:
            mass_status = "non_physical"
        elif s_jack is None:
            mass_status = "unknown"  # ajuste sin --jackknife: identificabilidad no evaluada
        elif use_boot:
            # Caso con leverage alto: el bootstrap es la σ defendible; la decisión de
            # identificabilidad usa snr_boot (no snr_jack, que subestima bajo leverage).
            mass_status = (
                "measured"
                if (
                    snr_boot is not None
                    and snr_boot >= args.min_snr_jack
                    and n_targets >= args.min_n_jack
                )
                else "not_identifiable"
            )
        elif (
            snr_jack is not None
            and snr_jack >= args.min_snr_jack
            and (snr_jack_excl is None or snr_jack_excl >= args.min_snr_jack)
            and n_targets >= args.min_n_jack
        ):
            mass_status = "measured"
        else:
            mass_status = "not_identifiable"
        out_rows.append(
            {
                "perturber": d["perturber"],
                "name": d.get("perturber_name", ""),
                "n_targets": d.get("n_targets"),
                "n_obs": d.get("n_obs_final"),
                "mass_fit_kg": m,
                "sigma_stat_kg": s_stat,
                "sigma_formal_kg": s_formal,
                "sigma_jack_kg": s_jack,
                "jack_leverage_top1": lev_top1,
                "sigma_jack_excl_top1_kg": s_jack_excl,
                "snr_jack_excl_top1": snr_jack_excl,
                "sigma_jack_defensible": sigma_jack_defensible,
                "sigma_boot_kg": s_boot,
                "snr_boot": snr_boot,
                "used_bootstrap_sigma": use_boot,
                "sigma_sys_kg": s_sys,
                "sigma_total_kg": s_tot,
                "sigma_total_frac": s_tot / abs(m) if m else None,
                "snr_jack": snr_jack,
                "delta_chi2_M0": delta_chi2_m0,
                "identif_p_false_alarm": p_false_alarm,
                "identifiable_profile": identifiable_profile,
                "mass_status": mass_status,
                "chi2_red": d.get("chi2_red"),
                "sys_floor_mas": d.get("sys_floor_mas"),
                "ref_mass_kg": ref,
                "ref_source": ref_src,
                "ratio_fit_over_ref": ratio,
                "z_total": z_tot,
                "is_calibrator": d["perturber"] in _BIG4,
                "reliable": (d.get("n_targets") or 0) >= args.min_n_reliable,
            }
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()))
        w.writeheader()
        w.writerows(out_rows)
    print(f"Escrito {args.out} ({len(out_rows)} masas)")

    print(
        f"\n{'name':10} {'N':>3} {'mass (kg)':>12} {'σ_tot':>10} {'σ%':>5} "
        f"{'ratio':>6} {'z_tot':>6} {'snrJ':>5} {'Δχ²0':>8} {'idP':>4} {'lev1':>5} "
        f"{'status':>16} {'cal':>4} {'rel':>4}"
    )
    for r in out_rows:
        ztxt = f"{r['z_total']:+.2f}" if r["z_total"] is not None else "  -  "
        rtxt = f"{r['ratio_fit_over_ref']:.3f}" if r["ratio_fit_over_ref"] else "  -  "
        sfrac = f"{r['sigma_total_frac'] * 100:.1f}" if r["sigma_total_frac"] else "-"
        snrtxt = f"{r['snr_jack']:.1f}" if r["snr_jack"] is not None else "  -  "
        dchtxt = f"{r['delta_chi2_M0']:.1f}" if r["delta_chi2_M0"] is not None else "  -  "
        idptxt = "Y" if r["identifiable_profile"] else "no"
        levtxt = (
            f"{r['jack_leverage_top1']:.2f}" if r["jack_leverage_top1"] is not None else "  -  "
        )
        print(
            f"{r['name']:10} {r['n_targets'] or 0:>3} {r['mass_fit_kg']:>12.3e} "
            f"{r['sigma_total_kg']:>10.2e} {sfrac:>5} {rtxt:>6} {ztxt:>6} {snrtxt:>5} "
            f"{dchtxt:>8} {idptxt:>4} {levtxt:>5} {r['mass_status']:>16} "
            f"{'Y' if r['is_calibrator'] else '':>4} {'Y' if r['reliable'] else 'no':>4}"
        )
    n_meas = sum(1 for r in out_rows if r["mass_status"] == "measured")
    n_bound = sum(1 for r in out_rows if r["mass_status"] == "not_identifiable")
    n_unk = sum(1 for r in out_rows if r["mass_status"] == "unknown")
    n_nonphys = sum(1 for r in out_rows if r["mass_status"] == "non_physical")
    n_lever = sum(
        1
        for r in out_rows
        if r["jack_leverage_top1"] is not None and r["jack_leverage_top1"] > args.max_leverage
    )
    if n_unk == len(out_rows):
        print("\n[F2] identificabilidad no evaluada (ajustes sin --jackknife)")
    else:
        print(
            f"\n[F2] {n_meas} medidas (snr_jack≥{args.min_snr_jack:g} con y sin réplica "
            f"dominante, N≥{args.min_n_jack}), "
            f"{n_bound} no identificables (cota), {n_nonphys} no físicas (M≤0), "
            f"{n_unk} sin jackknife"
        )
        if n_lever:
            print(
                f"[B6] {n_lever} perturbadores con leverage top-1 > {args.max_leverage:g}: "
                "σ_jack dominada por una réplica (~1 gdl efectivo) — no defendible como σ; "
                "requiere bootstrap/delete-d o más encuentros"
            )

    # M13/T21 — resumen del criterio de verosimilitud perfilada (paralelo al snr_jack).
    n_id_prof = sum(1 for r in out_rows if r["identifiable_profile"])
    n_prof_eval = sum(1 for r in out_rows if r["delta_chi2_M0"] is not None)
    n_disagree = sum(
        1
        for r in out_rows
        if r["delta_chi2_M0"] is not None
        and r["identifiable_profile"] != (r["mass_status"] == "measured")
    )
    print(
        f"\n[M13/T21] identificabilidad por verosimilitud perfilada "
        f"(Δχ²(M=0) = (M̂/σ_formal)² > {delta_chi2_thr:g} ≈ {args.identif_nsigma:g}σ): "
        f"{n_id_prof}/{n_prof_eval} identificables; "
        f"{n_disagree} discrepan del veredicto snr_jack (measured)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
