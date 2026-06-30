"""Diagnóstico: estructura temporal de los tránsitos Gaia de un objetivo.

Gaia observa cada objeto en "tránsitos de CCD": al cruzar el plano focal, hasta
~9 CCDs lo miden en ~40 s, todos compartiendo el mismo error de actitud/calibración
→ altamente correlacionados. Si el ajuste los trata como independientes,
sobre-cuenta la información (σ subestimada, χ²_red sesgado). Este script mide cuántos
tránsitos caen en cada "visita" (cluster temporal) para cuantificar la correlación.
"""

from __future__ import annotations

import sys

import numpy as np

from scripts.mass.fit_mass_gaia_loo import fetch_gaia_full
from src.utils.config import load_config


def main() -> None:
    target = int(sys.argv[1]) if len(sys.argv) > 1 else 18937
    release = sys.argv[2] if len(sys.argv) > 2 else "fpr"
    cfg = load_config("config.yaml")
    gaia = cfg.sources.gaia_sso
    gaia.release = release
    df = fetch_gaia_full(gaia.archive_url, target, gaia.active())
    ep = np.sort(df["epoch"].to_numpy().astype(float))  # días desde época ref TCB
    n = ep.size
    dt = np.diff(ep)
    dt_days = dt  # ya en días

    # Cluster: gap < 0.05 d (~72 min, una vuelta de FOV es ~106 min; CCDs en ~40 s).
    # Usamos varios umbrales para ver la jerarquía CCD→FOV→visita.
    print(f"target {target} ({release}): {n} tránsitos")
    print(f"  arco: {ep.max() - ep.min():.1f} d")
    for thr in (1e-3, 0.01, 0.05, 0.5, 5.0):
        # nº de clusters separados por gaps > thr
        n_clusters = 1 + int((dt_days > thr).sum())
        sizes = []
        cur = 1
        for g in dt_days:
            if g > thr:
                sizes.append(cur)
                cur = 1
            else:
                cur += 1
        sizes.append(cur)
        sizes = np.array(sizes)
        print(
            f"  gap>{thr:>6.3g} d: {n_clusters:>4} clusters, "
            f"tamaño medio {sizes.mean():.2f}, máx {sizes.max()}, "
            f"(N_indep/N = {n_clusters / n:.3f})"
        )
    # Histograma de gaps pequeños (estructura CCD).
    small = dt_days[dt_days < 0.01] * 86400.0  # s
    if small.size:
        print(
            f"  gaps intra-cluster (<864 s): n={small.size}, "
            f"mediana {np.median(small):.1f} s, rango [{small.min():.1f}, {small.max():.1f}] s"
        )


if __name__ == "__main__":
    main()
