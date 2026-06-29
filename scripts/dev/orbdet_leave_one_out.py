"""Leave-one-out por objetivo: ¿el sobre-tiro de masa lo arrastra un encuentro?

Ajusta la masa con los N objetivos y luego N veces quitando uno. Si quitar un
objetivo concreto baja mucho la masa hacia la verdad, ese encuentro domina (mala
geometría / sistemático local) y conviene vetarlo; si el sobre-tiro es uniforme, es
un sistemático común de los datos, no un outlier.
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


def _fit(targets, studied, bg, epoch):
    # fit a piso 0 → calibrar piso en órbitas ajustadas → re-ajustar (como producción).
    m0, el0, _r0 = determine_shared_mass(
        targets,
        studied.mass_msun,
        studied.elements,
        epoch,
        perturber_name=studied.name,
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
        perturber_name=studied.name,
        background_perturbers=bg,
        backend="assist",
        gr=True,
    )
    mass, _f, res = determine_shared_mass(
        targets,
        m0,
        studied.elements,
        epoch,
        perturber_name=studied.name,
        background_perturbers=bg,
        backend="assist",
        gr=True,
        sys_floor_mas=floor,
        max_iter=40,
    )
    sig = float(np.sqrt(res.covariance[0, 0]))
    return mass, sig, res.chi2_reduced, floor


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--perturber", type=int, default=2)
    ap.add_argument("--release", default="fpr")
    args = ap.parse_args()

    cfg = load_config("config.yaml")
    gaia = cfg.sources.gaia_sso
    gaia.release = args.release
    rc = gaia.active()
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
    true = studied.mass_msun * M_SUN_KG

    targets, nums = [], []
    for s in specs:
        if s.target not in emap:
            continue
        raw = _fetch_target(gaia.archive_url, s.target, rc)
        if raw is None:
            continue
        targets.append(_build_target_obs(raw, emap[s.target], epoch, _epoch_ref(rc), bg))
        nums.append(s.target)

    m, sig, chi, floor = _fit(targets, studied, bg, epoch)
    print(
        f"{pname}: verdad={true:.3e} kg | ALL N={len(targets)}: "
        f"mass={m * M_SUN_KG:.3e} ratio={m * M_SUN_KG / true:.3f} z={(m * M_SUN_KG - true) / (sig * M_SUN_KG):.2f} χ²={chi:.3f}"
    )
    print("\nLeave-one-out (quita el objetivo indicado):")
    print(f"  {'drop':>8} {'ratio':>7} {'z':>6} {'Δratio':>7}")
    base_ratio = m * M_SUN_KG / true
    for i, num in enumerate(nums):
        sub = [t for j, t in enumerate(targets) if j != i]
        mi, sigi, chii, _ = _fit(sub, studied, bg, epoch)
        ri = mi * M_SUN_KG / true
        zi = (mi * M_SUN_KG - true) / (sigi * M_SUN_KG)
        print(f"  {num:>8} {ri:>7.3f} {zi:>6.2f} {ri - base_ratio:>+7.3f}")


if __name__ == "__main__":
    main()
