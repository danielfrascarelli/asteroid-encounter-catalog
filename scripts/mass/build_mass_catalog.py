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
``σ_total = √(σ_stat² + (f_sys·masa)²)`` y un z basado en σ_total. Es el tratamiento
estándar de "incertidumbre externa" cuando la formal subestima por sistemáticos no
modelados.

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

    out_rows = []
    for d in sorted(rows, key=lambda x: x["perturber"]):
        m = d["mass_fit_kg"]
        s_stat = d.get("mass_fit_sigma_kg") or 0.0
        s_sys = f_sys * m
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
        s_formal = d.get("mass_fit_sigma_formal_kg")
        s_jack = d.get("mass_fit_sigma_jack_kg")
        snr_jack = (m / s_jack) if (s_jack and s_jack > 0) else None
        if s_jack is None:
            mass_status = "unknown"  # ajuste sin --jackknife: identificabilidad no evaluada
        elif snr_jack is not None and snr_jack >= args.min_snr_jack:
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
                "sigma_sys_kg": s_sys,
                "sigma_total_kg": s_tot,
                "sigma_total_frac": s_tot / m if m else None,
                "snr_jack": snr_jack,
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
        f"{'ratio':>6} {'z_tot':>6} {'snrJ':>5} {'status':>16} {'cal':>4} {'rel':>4}"
    )
    for r in out_rows:
        ztxt = f"{r['z_total']:+.2f}" if r["z_total"] is not None else "  -  "
        rtxt = f"{r['ratio_fit_over_ref']:.3f}" if r["ratio_fit_over_ref"] else "  -  "
        sfrac = f"{r['sigma_total_frac'] * 100:.1f}" if r["sigma_total_frac"] else "-"
        snrtxt = f"{r['snr_jack']:.1f}" if r["snr_jack"] is not None else "  -  "
        print(
            f"{r['name']:10} {r['n_targets'] or 0:>3} {r['mass_fit_kg']:>12.3e} "
            f"{r['sigma_total_kg']:>10.2e} {sfrac:>5} {rtxt:>6} {ztxt:>6} {snrtxt:>5} "
            f"{r['mass_status']:>16} {'Y' if r['is_calibrator'] else '':>4} "
            f"{'Y' if r['reliable'] else 'no':>4}"
        )
    n_meas = sum(1 for r in out_rows if r["mass_status"] == "measured")
    n_bound = sum(1 for r in out_rows if r["mass_status"] == "not_identifiable")
    n_unk = sum(1 for r in out_rows if r["mass_status"] == "unknown")
    if n_unk == len(out_rows):
        print("\n[F2] identificabilidad no evaluada (ajustes sin --jackknife)")
    else:
        print(
            f"\n[F2] {n_meas} medidas (snr_jack≥{args.min_snr_jack:g}), "
            f"{n_bound} no identificables (cota), {n_unk} sin jackknife"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
