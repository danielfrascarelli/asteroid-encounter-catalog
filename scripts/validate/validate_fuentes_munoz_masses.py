"""Cross-check our orbdet mass catalog against Fuentes-Muñoz et al. (2025) masses.

This is the **mass** cross-check against an independent Gaia-FPR mass study (see
``docs/mass_determination_results.md``). The companion script
``validate_fuentes_munoz_2025.py`` validates *encounter pairs*; this one validates
the *determined masses*.

Source
------
Fuentes-Muñoz, Farnocchia, Giorgini & Park (2025), "Asteroid Mass Estimation by
Mutual Perturbations during Close Encounters after Gaia FPR", AJ 170, 353
(DOI 10.3847/1538-3881/ae0cc9). Machine-readable Table 5 carries, per perturber:

    GMfin   (km^3 s^-2)  final mass parameter (their FPR posterior)
    e_GMfin (km^3 s^-2)  1-sigma uncertainty

We convert their GM to mass via M = GM / G with G = 6.67430e-20 km^3 kg^-1 s^-2
and compare against our ``mass_fit_kg`` per perturber, reporting the ratio and a
z-score that folds both uncertainties:

    z = (M_ours - M_fm) / sqrt(sigma_ours^2 + sigma_fm^2)

Caveat
------
For the Big-4 calibrators Fuentes-Muñoz pin GMfin to the literature/SB441 seed
(tiny e_GMfin), so those rows are NOT an independent comparison — they recover the
same DAWN/Goffin/Vernazza values we already calibrate against. The meaningful,
independent comparisons are the **non-calibrator** perturbers (Psyche, Sylvia, …)
where FM ran their own FPR fit. The summary reports both populations separately.

Usage
-----
    docker compose run --rm pipeline python -m scripts.ingest.download_fuentes_munoz
    docker compose run --rm pipeline python -m scripts.validate.validate_fuentes_munoz_masses
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import re
from pathlib import Path

import polars as pl

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
logger = logging.getLogger(__name__)

_MRT_PATH = Path("data/raw/fuentes_munoz_2025/ajae0cc9t5_mrt.txt")
_CATALOG_PATH = Path("data/output/orbdet/mass_catalog.csv")
_OUT_DIR = Path("data/output/literature_validation")

# Newtonian constant in the MRT's GM units: G = 6.67430e-11 m^3 kg^-1 s^-2, and
# 1 m^3 = 1e-9 km^3, so G = 6.67430e-20 km^3 kg^-1 s^-2 (CODATA 2018).
_G_KM3_KG_S2 = 6.67430e-20

# Fixed-width byte ranges from the MRT byte-by-byte description (1-indexed,
# converted to 0-indexed half-open slices).
_COL_ASTEROID = (0, 22)
_COL_GMFIN = (65, 76)
_COL_E_GMFIN = (77, 88)

_PROVISIONAL = re.compile(r"^\d{4}\s[A-Z]{1,2}\d*$")
_LEADING_NUMBER = re.compile(r"^(\d+)(?:\s|$)")


def _to_float(text: str) -> float | None:
    """Parse a possibly-blank MRT numeric field to float (None if blank)."""
    text = text.strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def parse_table5_masses(path: Path) -> pl.DataFrame:
    """Parse the MRT into per-perturber masses.

    Returns
    -------
    pl.DataFrame
        Columns ``perturber`` (Int64, MPC number), ``fm_gm_fin`` (km^3/s^2),
        ``fm_e_gm_fin`` (km^3/s^2), ``fm_mass_kg``, ``fm_sigma_kg`` — one row per
        numbered perturber with a parseable final mass.
    """
    lines = path.read_text().splitlines()
    div = [i for i, ln in enumerate(lines) if ln.startswith("---")]
    data = lines[div[-1] + 1 :]

    rows: list[dict[str, float | int]] = []
    for ln in data:
        if not ln.strip():
            continue
        asteroid = ln[_COL_ASTEROID[0] : _COL_ASTEROID[1]].strip()
        m = _LEADING_NUMBER.match(asteroid)
        if _PROVISIONAL.match(asteroid) or not m:
            continue  # numbered perturbers only
        gm_fin = _to_float(ln[_COL_GMFIN[0] : _COL_GMFIN[1]])
        if gm_fin is None or gm_fin <= 0:
            continue
        e_gm_fin = _to_float(ln[_COL_E_GMFIN[0] : _COL_E_GMFIN[1]])
        rows.append(
            {
                "perturber": int(m.group(1)),
                "fm_gm_fin": gm_fin,
                "fm_e_gm_fin": e_gm_fin if e_gm_fin is not None else float("nan"),
                "fm_mass_kg": gm_fin / _G_KM3_KG_S2,
                "fm_sigma_kg": (e_gm_fin / _G_KM3_KG_S2) if e_gm_fin is not None else float("nan"),
            }
        )

    return pl.DataFrame(rows).unique(subset=["perturber"], keep="first")


def main() -> int:
    parser = argparse.ArgumentParser(description="Cross-check orbdet masses vs Fuentes-Muñoz 2025")
    parser.add_argument("--mrt", default=str(_MRT_PATH))
    parser.add_argument("--catalog", default=str(_CATALOG_PATH))
    parser.add_argument("--out-dir", default=str(_OUT_DIR))
    args = parser.parse_args()

    mrt = Path(args.mrt)
    cat_path = Path(args.catalog)
    if not mrt.exists():
        logger.error("MRT not found at %s — run scripts.ingest.download_fuentes_munoz first.", mrt)
        return 1
    if not cat_path.exists():
        logger.error(
            "Mass catalog not found at %s — run scripts.mass.build_mass_catalog first.", cat_path
        )
        return 1

    fm = parse_table5_masses(mrt)
    logger.info("Parsed %d numbered FM perturber masses from Table 5.", fm.height)

    ours = pl.read_csv(cat_path)
    joined = (
        ours.join(fm, on="perturber", how="inner")
        .with_columns(
            (pl.col("mass_fit_kg") / pl.col("fm_mass_kg")).alias("ratio_ours_over_fm"),
            (
                (pl.col("mass_fit_kg") - pl.col("fm_mass_kg"))
                / (pl.col("sigma_total_kg") ** 2 + pl.col("fm_sigma_kg") ** 2).sqrt()
            ).alias("z_vs_fm"),
        )
        .sort("perturber")
    )

    if joined.is_empty():
        logger.error("No perturbers overlap between our catalog and FM Table 5.")
        return 1

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_csv = out_dir / "fuentes_munoz_2025_mass_comparison.csv"
    joined.select(
        "perturber",
        "name",
        "is_calibrator",
        "reliable",
        "n_targets",
        "mass_fit_kg",
        "sigma_total_kg",
        "fm_mass_kg",
        "fm_sigma_kg",
        "ratio_ours_over_fm",
        "z_vs_fm",
    ).write_csv(out_csv)

    # Independent comparison = non-calibrators (FM pins calibrators to the seed).
    indep = joined.filter(~pl.col("is_calibrator"))
    reliable_indep = indep.filter(pl.col("reliable"))

    def _stats(df: pl.DataFrame) -> dict:
        if df.is_empty():
            return {"n": 0}
        ratios = df["ratio_ours_over_fm"].to_list()
        zs = [z for z in df["z_vs_fm"].to_list() if z is not None and not math.isnan(z)]
        return {
            "n": df.height,
            "median_ratio": float(df["ratio_ours_over_fm"].median()),
            "mean_ratio": float(df["ratio_ours_over_fm"].mean()),
            "min_ratio": min(ratios),
            "max_ratio": max(ratios),
            "n_within_z3": sum(1 for z in zs if abs(z) < 3),
            "n_with_z": len(zs),
        }

    summary = {
        "source": "Fuentes-Muñoz et al. 2025, AJ 170, 353, Table 5 (GMfin)",
        "catalog": str(cat_path),
        "n_overlap": joined.height,
        "all_overlap": _stats(joined),
        "independent_noncalibrators": _stats(indep),
        "independent_reliable": _stats(reliable_indep),
        "note": (
            "Calibrators (Ceres/Vesta/Pallas/Hygiea) are pinned by FM to the "
            "literature/SB441 seed, so they are not an independent comparison. The "
            "independent_* blocks (non-calibrators) carry the real cross-check."
        ),
        "outputs": {"comparison_csv": str(out_csv)},
    }
    (out_dir / "fuentes_munoz_2025_mass_summary.json").write_text(json.dumps(summary, indent=2))

    print("\n" + "=" * 88)
    print("FUENTES-MUÑOZ et al. 2025 (AJ 170, 353) — MASS cross-check vs orbdet catalog")
    print("=" * 88)
    print(
        f"{'#':>5} {'name':<12} {'cal':>3} {'rel':>3} {'M_ours':>10} {'M_FM':>10} {'ratio':>7} {'z':>7}"
    )
    print("-" * 88)
    for r in joined.iter_rows(named=True):
        z = r["z_vs_fm"]
        z_str = f"{z:+.2f}" if z is not None and not math.isnan(z) else "  n/a"
        print(
            f"{r['perturber']:>5} {r['name'][:12]:<12} "
            f"{'Y' if r['is_calibrator'] else '.':>3} {'Y' if r['reliable'] else '.':>3} "
            f"{r['mass_fit_kg']:>10.3e} {r['fm_mass_kg']:>10.3e} "
            f"{r['ratio_ours_over_fm']:>7.3f} {z_str:>7}"
        )
    print("-" * 88)
    si = summary["independent_noncalibrators"]
    sr = summary["independent_reliable"]
    if si["n"]:
        print(
            f"Independent (non-calibrator) perturbers: {si['n']}  "
            f"median ratio={si['median_ratio']:.3f}  "
            f"({si['n_within_z3']}/{si['n_with_z']} within |z|<3)"
        )
    if sr["n"]:
        print(
            f"  of which 'reliable': {sr['n']}  median ratio={sr['median_ratio']:.3f}  "
            f"({sr['n_within_z3']}/{sr['n_with_z']} within |z|<3)"
        )
    print("Note: calibrators are pinned by FM to the seed → not an independent check.")
    print("=" * 88)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
