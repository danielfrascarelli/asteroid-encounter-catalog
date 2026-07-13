"""F3 — cómputo de f_sys (RMS de ratio−1) sobre calibradores con N≥20.

Compara el fondo estándar de 16 (leído del catálogo jackknife) contra el fondo
extendido (leído de los JSON de ``orbdet_fit_realdata`` con ``--extra-background``).
No hace cómputo pesado; solo agrega ratios ya calculados.

Uso:
    docker compose run --rm pipeline python -m scripts.mass.f3_fsys \
        --jack-catalog data/output/orbdet/mass_catalog_jack.csv \
        --ext-json data/output/orbdet/f3_extbg/ceres_extbg20_fpr.json \
        --ext-json data/output/orbdet/f3_extbg/vesta_extbg20_fpr.json \
        --ext-json data/output/orbdet/f3_extbg/hygiea_extbg20_fpr.json
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

# Calibradores con N≥20 que entran a f_sys (Pallas queda fuera, N=6).
_CALIBRATORS_N20 = {1, 4, 10}


def _rms(vals: list[float]) -> float:
    return math.sqrt(sum(v * v for v in vals) / len(vals)) if vals else float("nan")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--jack-catalog", type=Path, required=True)
    ap.add_argument("--ext-json", type=Path, action="append", default=[])
    args = ap.parse_args()

    base: dict[int, dict] = {}
    with args.jack_catalog.open() as fh:
        for r in csv.DictReader(fh):
            num = int(r["perturber"])
            if num in _CALIBRATORS_N20:
                base[num] = {
                    "name": r["name"],
                    "n": int(r["n_targets"]),
                    "ratio": float(r["ratio_fit_over_ref"]),
                    "chi2": float(r["chi2_red"]),
                    "ref_kg": float(r["ref_mass_kg"]),
                }

    ext: dict[int, dict] = {}
    for jp in args.ext_json:
        d = json.loads(jp.read_text())
        num = int(d["perturber"])
        ref = d.get("mass_lit_kg") or (base.get(num, {}).get("ref_kg"))
        ratio = (d["mass_fit_kg"] / ref) if ref else None
        ext[num] = {
            "name": d["perturber_name"],
            "n": int(d["n_targets"]),
            "ratio": ratio,
            "chi2": float(d["chi2_red"]),
            "n_background": d.get("n_background"),
        }

    print(
        f"{'#':>4} {'name':<10} {'N':>4}  {'ratio16':>8} {'chi16':>6}  "
        f"{'ratioEXT':>8} {'chiEXT':>6}  {'nbg':>4}"
    )
    base_ratios: list[float] = []
    ext_ratios: list[float] = []
    for num in sorted(_CALIBRATORS_N20):
        b = base.get(num, {})
        e = ext.get(num, {})
        if b:
            base_ratios.append(b["ratio"] - 1.0)
        er = e.get("ratio")
        if er is not None:
            ext_ratios.append(er - 1.0)
        print(
            f"{num:>4} {b.get('name', e.get('name','?')):<10} "
            f"{b.get('n', e.get('n','')):>4}  "
            f"{b.get('ratio', float('nan')):>8.4f} {b.get('chi2', float('nan')):>6.3f}  "
            f"{(er if er is not None else float('nan')):>8.4f} "
            f"{e.get('chi2', float('nan')):>6.3f}  {str(e.get('n_background','')):>4}"
        )
    print("-" * 60)
    print(
        f"f_sys (fondo 16, N calib={len(base_ratios)}):  {_rms(base_ratios)*100:.3f} % (informativo)"
    )
    if ext_ratios:
        print(
            f"f_sys (fondo EXT, N calib={len(ext_ratios)}): {_rms(ext_ratios)*100:.3f} % (informativo)"
        )

    # Gate (reformulado 2026-07-04, tribunal menor 14): la comparación de f_sys
    # de 3 puntos (±41 % de incertidumbre relativa) no puede resolver el efecto
    # del fondo extendido — 4.158 % vs 4.257 % es ruido de redondeo. El test con
    # potencia es la **Δmasa pareada por calibrador**: mismo calibrador, mismos
    # datos, solo cambia el fondo. Si extender el fondo moviera la masa, se ve
    # acá directamente.
    print()
    deltas: list[tuple[str, float]] = []
    for num in sorted(_CALIBRATORS_N20):
        b, e = base.get(num), ext.get(num)
        if not b or not e or e.get("ratio") is None:
            continue
        m16 = b["ratio"] * b["ref_kg"]
        mext = e["ratio"] * b["ref_kg"]
        delta = abs(mext - m16) / m16
        deltas.append((b["name"], delta))
        print(f"Δmasa pareada {b['name']:<10}: {delta*100:.4f} %")
    if deltas:
        worst = max(d for _, d in deltas)
        ok = worst < 0.0025
        print(f"gate: max Δmasa pareada = {worst*100:.4f} % < 0.25 % → {'PASS' if ok else 'FAIL'}")
        return 0 if ok else 1
    print("gate: sin pares base/EXT comparables — no evaluado")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
