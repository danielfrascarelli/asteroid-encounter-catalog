"""Diagnóstico: ¿cuánto deriva un perturbador integrado libremente vs su efeméride?

El ajuste de masa siembra los 16 asteroides grandes (incluido el perturbador bajo
estudio) desde DE441 en la época común y luego los **integra libremente** en
rebound (para poder variar la masa). Este script mide el error posicional entre
esa órbita integrada y la efeméride DE441 a lo largo del arco FPR, para cada Big-4.
Si la deriva es grande, la geometría del encuentro queda mal modelada y la masa
ajustada lo absorbe (sesgo).
"""

from __future__ import annotations

import numpy as np

from src.orbdet.dynamics_assist import (
    BIG_ASTEROIDS,
    _build_sim,
    big_asteroid_perturbers,
    load_ephem,
)
from src.orbdet.kepler import KeplerElements

AU_KM = 149_597_870.7
COMMON_EPOCH = 2457200.5  # snapshot MPCORB 20150524, época común del Big-4 FPR
ARC_DAYS = 500.0  # ~ media ventana FPR a cada lado


def _ephem_pos(ephem, name, jd):
    p = ephem.get_particle(name, jd - ephem.jd_ref)
    return np.array([p.x, p.y, p.z], dtype=float)


def main() -> None:
    ephem = load_ephem()
    # Dummy test element (no se usa su salida, sólo arrastra la integración).
    dummy = KeplerElements(a=2.5, e=0.1, i=0.1, Omega=0.2, omega=0.3, M=0.4)
    # Fondo = los 16 asteroides desde la efeméride en la época común.
    all16 = big_asteroid_perturbers(COMMON_EPOCH, names=BIG_ASTEROIDS)
    name_to_idx = {ap.name: i for i, ap in enumerate(all16)}

    sample_jd = COMMON_EPOCH + np.linspace(-ARC_DAYS, ARC_DAYS, 11)

    print(f"Época común JD TDB = {COMMON_EPOCH}, arco ±{ARC_DAYS:.0f} d")
    print(f"{'body':10} {'|dr| max (km)':>14} {'|dr| end (km)':>14} {'|dr| max (mas@2AU)':>20}")
    for name in ("Ceres", "Pallas", "Vesta", "Hygiea"):
        # Integra el conjunto completo de 16 (igual que el ajuste real).
        sim, _ = _build_sim(COMMON_EPOCH, dummy, all16, ephem, gr=True)
        idx = name_to_idx[name]
        errs_km = []
        for jd in sample_jd:
            sim.integrate(float(jd - ephem.jd_ref))
            p = sim.particles[idx]
            r_int = np.array([p.x, p.y, p.z], dtype=float)  # baricéntrico ecuatorial
            r_eph = _ephem_pos(ephem, name, jd)
            errs_km.append(np.linalg.norm(r_int - r_eph) * AU_KM)
        errs_km = np.array(errs_km)
        max_km = errs_km.max()
        end_km = errs_km[-1]
        # Ángulo subtendido a ~2 AU (escala típica Gaia-asteroide), en mas.
        max_mas = np.degrees((max_km / AU_KM) / 2.0) * 3.6e6
        print(f"{name:10} {max_km:>14.1f} {end_km:>14.1f} {max_mas:>20.3f}")


if __name__ == "__main__":
    main()
