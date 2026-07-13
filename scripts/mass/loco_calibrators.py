"""LOCO — validación leave-one-calibrator-out del piso sistemático f_sys (B8).

El tribunal (2026-07-04, B8) señaló que validar los calibradores con un σ_total
que incluye ``f_sys·M``, cuando f_sys es el RMS de esos mismos 3 calibradores, es
circular: con 3 puntos, ``|desvío|/RMS ≤ √3`` por construcción. Este script rompe
la circularidad: para cada calibrador ``c``, calcula ``f_sys^(−c)`` con los
*otros* calibradores y reporta el z de ``c`` contra ese piso independiente.

También reporta la incertidumbre relativa de f_sys (≈ 1/√(2(n−1)) para un RMS de
n puntos) para que nadie lea el piso como una constante exacta.

Uso:
    docker compose run --rm pipeline python -m scripts.mass.loco_calibrators \
        --catalog data/output/orbdet/mass_catalog_jack.csv
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

# Calibradores que definen el piso (N≥15; Pallas N=6 queda fuera, como en el pipeline).
_FLOOR_CALIBRATORS = {1, 4, 10}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--catalog", type=Path, required=True, help="mass_catalog CSV (jackknife)")
    args = ap.parse_args()

    cal: dict[int, dict] = {}
    with args.catalog.open() as fh:
        for r in csv.DictReader(fh):
            num = int(r["perturber"])
            if num in _FLOOR_CALIBRATORS and r.get("ratio_fit_over_ref"):
                cal[num] = {
                    "name": r["name"],
                    "ratio": float(r["ratio_fit_over_ref"]),
                    "m": float(r["mass_fit_kg"]),
                    "s_stat": float(r["sigma_stat_kg"] or 0.0),
                    "ref": float(r["ref_mass_kg"]),
                }

    if len(cal) < 3:
        print(f"Solo {len(cal)} calibradores con ratio en {args.catalog} — LOCO necesita ≥3")
        return 1

    devs_all = [c["ratio"] - 1.0 for c in cal.values()]
    f_all = math.sqrt(sum(d * d for d in devs_all) / len(devs_all))
    rel_unc = 1.0 / math.sqrt(2.0 * (len(devs_all) - 1))
    print(
        f"f_sys (todos, n={len(devs_all)}): {f_all * 100:.2f} % "
        f"± {f_all * rel_unc * 100:.2f} % (incertidumbre relativa ≈ {rel_unc * 100:.0f} %)"
    )
    print()
    print(f"{'cal':<8} {'ratio':>7} {'f_sys(-c)':>10} {'z_circular':>10} {'z_LOCO':>8}")

    worst_z = 0.0
    for num, c in sorted(cal.items()):
        others = [cal[k]["ratio"] - 1.0 for k in cal if k != num]
        f_loco = math.sqrt(sum(d * d for d in others) / len(others))

        def _z(f_sys: float) -> float:
            s_tot = math.sqrt(c["s_stat"] ** 2 + (f_sys * abs(c["m"])) ** 2)
            return (c["m"] - c["ref"]) / s_tot if s_tot > 0 else float("nan")

        z_circ = _z(f_all)
        z_loco = _z(f_loco)
        worst_z = max(worst_z, abs(z_loco))
        print(
            f"{c['name']:<8} {c['ratio']:>7.4f} {f_loco * 100:>9.2f}% "
            f"{z_circ:>+10.2f} {z_loco:>+8.2f}"
        )

    print()
    ok = worst_z < 3.0
    print(f"gate LOCO: max |z_LOCO| = {worst_z:.2f} < 3 → {'PASS' if ok else 'FAIL'}")
    print(
        "Nota: el z 'circular' (piso que incluye al propio calibrador) está acotado a "
        "√3 por construcción y NO es evidencia de exactitud; el z_LOCO sí lo es."
    )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
