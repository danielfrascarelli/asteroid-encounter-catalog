"""¿Queda correlación de residuos MÁS ALLÁ del cruce FOV (a escala de visita)?

El fix de covarianza en bloques decorrelaciona los CCDs DENTRO de cada cruce FOV.
Pero una "visita" de Gaia contiene ~2 cruces FOV (los dos campos de visión, ~106 min
aparte) que comparten el mismo período de actitud y podrían seguir correlacionados.
Si así fuera, σ(masa) sigue subestimada y conviene blanquear en bloques a escala de
visita (gap mayor).

Este script ajusta con el piso ya calibrado y mide el ICC de los residuos blanqueados
a varias escalas temporales de agrupación. ICC≈0 a la escala de visita ⇒ el bloqueo
por FOV basta; ICC>0 ⇒ hay correlación residual a esa escala que justifica agrupar más.
"""

from __future__ import annotations

import argparse
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
from src.orbdet.mass_determination import calibrate_sys_floor, determine_shared_mass
from src.utils.config import load_config


def _icc_at_gap(resid: np.ndarray, epochs: np.ndarray, gap: float) -> tuple[float, float, int]:
    """ICC (ANOVA one-way) de *resid* agrupado por gaps > *gap*. (icc, m_bar, n_clusters)."""
    order = np.argsort(epochs)
    r = resid[order]
    ep = epochs[order]
    cid = np.zeros(r.size, dtype=int)
    c = 0
    for i in range(1, r.size):
        if ep[i] - ep[i - 1] > gap:
            c += 1
        cid[i] = c
    groups = [r[cid == k] for k in range(c + 1) if (cid == k).sum() >= 2]
    if len(groups) < 2:
        return 0.0, float(r.size), c + 1
    grand = r.mean()
    n, k = r.size, len(groups)
    ssb = sum(g.size * (g.mean() - grand) ** 2 for g in groups)
    ssw = sum(((g - g.mean()) ** 2).sum() for g in groups)
    msb, msw = ssb / (k - 1), ssw / max(n - k, 1)
    sizes = np.array([g.size for g in groups])
    m0 = (n - (sizes**2).sum() / n) / (k - 1)
    denom = msb + (m0 - 1) * msw
    icc = (msb - msw) / denom if denom > 0 else 0.0
    return float(icc), float(n / (c + 1)), c + 1


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--perturber", type=int, default=2)
    ap.add_argument("--release", default="fpr")
    args = ap.parse_args()

    cfg = load_config("config.yaml")
    gaia = cfg.sources.gaia_sso
    gaia.release = args.release
    release_cfg = gaia.active()
    pname = _BIG4_NAME_BY_NUMBER[args.perturber]
    specs = _read_targets_from_csv(
        Path("data/output/stage4_validation_summary.csv"), args.perturber
    )
    mid = float(np.mean([s.jd_tdb for s in specs]))
    snap = _best_mpcorb_snapshot(_MPCORB_ARCHIVE_DIR, mid)
    emap = load_element_rows(snap, [s.target for s in specs])
    epoch = float(next(iter(emap.values()))["epoch_jd"])
    studied = big_asteroid_perturbers(epoch, names=(pname,))[0]
    bg = big_asteroid_perturbers(epoch, exclude=(pname,))

    targets, eps = [], []
    for s in specs:
        if s.target not in emap:
            continue
        raw = _fetch_target(gaia.archive_url, s.target, release_cfg)
        if raw is None:
            continue
        t = _build_target_obs(raw, emap[s.target], epoch, _epoch_ref(release_cfg), bg)
        targets.append(t)
        eps.append(np.asarray(t.obs_jd_tdb, dtype=float))

    # Ajuste inicial a piso 0, luego calibrar el piso en las órbitas YA ajustadas
    # (no en la semilla MPCORB) y re-ajustar — igual que la producción.
    m0, el0, _r0 = determine_shared_mass(
        targets,
        studied.mass_msun,
        studied.elements,
        epoch,
        perturber_name=pname,
        background_perturbers=bg,
        backend="assist",
        gr=True,
        sys_floor_mas=0.0,
        max_iter=40,
    )
    floor, _ = calibrate_sys_floor(
        targets,
        m0,
        el0,
        epoch,
        perturber_elements=studied.elements,
        perturber_name=pname,
        background_perturbers=bg,
        backend="assist",
        gr=True,
    )
    mass, _f, res = determine_shared_mass(
        targets,
        m0,
        studied.elements,
        epoch,
        perturber_name=pname,
        background_perturbers=bg,
        backend="assist",
        gr=True,
        sys_floor_mas=floor,
        max_iter=40,
    )
    print(
        f"{pname}: mass={mass * M_SUN_KG:.3e} ratio={mass / studied.mass_msun:.3f} "
        f"χ²_red={res.chi2_reduced:.3f} s_c={floor:.2f} mas"
    )

    r = np.asarray(res.residuals)
    print("\nICC de residuos blanqueados (post-fix) por escala de agrupación:")
    print(f"  {'gap (d)':>10} {'escala':>8} {'ICC medio':>10} {'m̄':>6}")
    start = 0
    for gap, label in [(0.01, "FOV"), (0.1, "~visita"), (0.5, "visita"), (2.0, "multidía")]:
        iccs, mbars = [], []
        s2 = start
        for ep in eps:
            n = ep.size
            icc, mbar, _ = _icc_at_gap(r[s2 : s2 + n], ep, gap)
            iccs.append(icc)
            mbars.append(mbar)
            s2 += n
        print(f"  {gap:>10.2f} {label:>8} {np.mean(iccs):>+10.3f} {np.mean(mbars):>6.1f}")
    print(
        "\n→ ICC≈0 a escala visita ⇒ el bloqueo por FOV es suficiente; el residual "
        "no es correlación a esa escala."
    )


if __name__ == "__main__":
    main()
