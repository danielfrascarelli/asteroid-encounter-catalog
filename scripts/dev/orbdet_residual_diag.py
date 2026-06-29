"""Diagnóstico profundo de residuos del ajuste conjunto sobre datos reales.

Corre el ajuste de un perturbador (default Ceres/FPR), extrae los residuos
blanqueados por tránsito y mide empíricamente:
  1. χ²_red global.
  2. Correlación intra-FOV de los residuos (ICC): ¿los ~7 CCDs de un mismo cruce
     comparten residuo (correlados → σ subestimada) o están dispersos (independientes)?
     Da el N efectivo de grados de libertad.
  3. Perfil χ²(masa) re-ajustando órbitas en cada masa fija (verosimilitud perfilada):
     ¿el mínimo es agudo y dónde cae respecto a la verdad?
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
from src.orbdet.mass_determination import determine_shared_mass
from src.utils.config import load_config


def _intra_fov_icc(resid: np.ndarray, epochs: np.ndarray, gap: float = 0.01) -> tuple[float, int]:
    """Correlación intra-cluster (ICC) de *resid* agrupado por cruces FOV.

    Devuelve ``(icc, n_clusters)``. ICC≈1 → residuos idénticos dentro del FOV
    (totalmente correlados); ICC≈0 → independientes. N_efectivo ≈ N/(1+(m̄-1)·ICC).
    """
    order = np.argsort(epochs)
    r = resid[order]
    ep = epochs[order]
    # Asigna cluster por gaps.
    cluster = np.zeros(r.size, dtype=int)
    c = 0
    for i in range(1, r.size):
        if ep[i] - ep[i - 1] > gap:
            c += 1
        cluster[i] = c
    groups = [r[cluster == k] for k in range(c + 1)]
    groups = [g for g in groups if g.size >= 2]
    if not groups:
        return 0.0, c + 1
    grand = r.mean()
    n = r.size
    k = len(groups)
    # ANOVA one-way: ICC = (MSB - MSW)/(MSB + (m0-1) MSW)
    ssb = sum(g.size * (g.mean() - grand) ** 2 for g in groups)
    ssw = sum(((g - g.mean()) ** 2).sum() for g in groups)
    msb = ssb / (k - 1) if k > 1 else 0.0
    msw = ssw / (n - k) if n > k else 0.0
    sizes = np.array([g.size for g in groups])
    m0 = (n - (sizes**2).sum() / n) / (k - 1) if k > 1 else sizes.mean()
    icc = (msb - msw) / (msb + (m0 - 1) * msw) if (msb + (m0 - 1) * msw) > 0 else 0.0
    return float(icc), c + 1


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--perturber", type=int, default=1)
    ap.add_argument("--release", default="fpr")
    ap.add_argument("--max-targets", type=int, default=None)
    ap.add_argument("--profile", action="store_true", help="perfil χ²(masa) (caro)")
    args = ap.parse_args()

    cfg = load_config("config.yaml")
    gaia = cfg.sources.gaia_sso
    gaia.release = args.release
    release_cfg = gaia.active()
    archive_url = gaia.archive_url

    pname = _BIG4_NAME_BY_NUMBER[args.perturber]
    csv = Path("data/output/stage4_validation_summary.csv")
    specs = _read_targets_from_csv(csv, args.perturber)
    if args.max_targets:
        specs = specs[: args.max_targets]

    mid_jd = float(np.mean([s.jd_tdb for s in specs]))
    snapshot = _best_mpcorb_snapshot(_MPCORB_ARCHIVE_DIR, mid_jd)
    elements_map = load_element_rows(snapshot, [s.target for s in specs])
    common_epoch = float(next(iter(elements_map.values()))["epoch_jd"])

    studied = big_asteroid_perturbers(common_epoch, names=(pname,))[0]
    background = big_asteroid_perturbers(common_epoch, exclude=(pname,))
    seed_mass = studied.mass_msun
    print(f"{pname}: seed mass {seed_mass * M_SUN_KG:.4e} kg (efeméride/verdad)")

    targets = []
    epochs_per_t = []
    for spec in specs:
        if spec.target not in elements_map:
            continue
        raw = _fetch_target(archive_url, spec.target, release_cfg)
        if raw is None:
            continue
        tobs = _build_target_obs(
            raw, elements_map[spec.target], common_epoch, _epoch_ref(release_cfg), background
        )
        targets.append(tobs)
        epochs_per_t.append(np.asarray(tobs.obs_jd_tdb, dtype=float))

    mass, fitted, result = determine_shared_mass(
        targets,
        seed_mass,
        studied.elements,
        common_epoch,
        perturber_name=pname,
        background_perturbers=background,
        backend="assist",
        gr=True,
        max_iter=40,
    )
    print(
        f"\nfit: mass={mass * M_SUN_KG:.4e} kg  ratio={mass / seed_mass:.3f}  "
        f"χ²_red={result.chi2_reduced:.3f}  conv={result.converged}"
    )

    # ICC por objetivo.
    res = np.asarray(result.residuals)
    start = 0
    iccs, neffs = [], []
    print("\nresiduos por objetivo (whitened):")
    for tobs, ep in zip(targets, epochs_per_t):
        n = ep.size
        r = res[start : start + n]
        start += n
        icc, ncl = _intra_fov_icc(r, ep)
        m_bar = n / ncl
        n_eff = n / (1 + (m_bar - 1) * max(icc, 0.0))
        iccs.append(icc)
        neffs.append(n_eff)
        print(
            f"  n={n:4d} FOV={ncl:3d} m̄={m_bar:4.1f}  RMS={np.sqrt((r**2).mean()):.2f}  "
            f"ICC={icc:+.3f}  N_eff={n_eff:5.1f}"
        )
    icc_mean = float(np.mean(iccs))
    n_tot = sum(e.size for e in epochs_per_t)
    n_eff_tot = sum(neffs)
    print(
        f"\nICC medio={icc_mean:+.3f}  N={n_tot}  N_eff≈{n_eff_tot:.0f}  "
        f"→ inflar σ por √(N/N_eff)={np.sqrt(n_tot / n_eff_tot):.2f}"
    )
    sig0 = float(np.sqrt(result.covariance[0, 0])) * M_SUN_KG
    print(f"σ(masa) formal={sig0:.3e} kg → corregida≈{sig0 * np.sqrt(n_tot / n_eff_tot):.3e} kg")

    if args.profile:
        print("\nperfil χ²(masa) [orbitas re-ajustadas a masa fija]:")
        # Re-ajusta órbitas con masa fija escaneando ratios.
        from src.orbdet.mass_determination import _ModelConfig, _target_resid_and_blocks

        for ratio in (0.6, 0.8, 1.0, 1.2, 1.4, 1.6):
            m_fix = seed_mass * ratio
            cfg_m = _ModelConfig(
                epoch_jd_tdb=common_epoch,
                perturber_elements=studied.elements,
                perturber_name=pname,
                background_perturbers=background,
                perturbers=(),
                integrator="ias15",
                dt_days=1.0,
                n_lighttime_iter=3,
                gm_rel_delta=1e-3,
                backend="assist",
                gr=True,
            )
            chi2 = 0.0
            ndof = 0
            for tobs, el in zip(targets, fitted):
                r, _jm, _je = _target_resid_and_blocks(m_fix, el, tobs, cfg_m)
                chi2 += float((r**2).sum())
                ndof += r.size
            print(
                f"  ratio={ratio:.2f} mass={m_fix * M_SUN_KG:.3e}  χ²={chi2:.1f} "
                f"χ²_red={chi2 / (ndof - 1):.3f}"
            )


if __name__ == "__main__":
    main()
