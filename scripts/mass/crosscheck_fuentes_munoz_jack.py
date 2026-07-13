"""Cross-check jackknife mass catalog vs Fuentes-Muñoz (2025), σ_jack vs σ_formal.

This is frente-P3 item **F5**: re-run the Fuentes-Muñoz mass cross-check using the
**external jackknife σ (F1)** instead of the formal Fisher σ, and report how many
perturbers stay consistent (|z| < 3) with the independent reference under the
larger, more honest error bar.

Why redo the cross-check
------------------------
``scripts.validate.validate_fuentes_munoz_masses`` compares our masses to
Fuentes-Muñoz using ``sigma_total_kg`` from the *non-jack* catalog, whose
statistical part is the Fisher/formal σ (covariance diagonal). That σ shrinks as
1/√N and **underestimates** the true error because it ignores the mass↔orbit
regression captured only by leaving encounters out. The jackknife σ (F1) folds
that regression back in, so a z built on it is the defendable one. This script
puts both z side by side so the impact is explicit.

Two z-scores per perturber (reference = Fuentes-Muñoz GMfin, Table 5):

    z_formal = (M_ours - M_fm) / sqrt(sigma_formal^2 + sigma_fm^2)
    z_total  = (M_ours - M_fm) / sqrt(sigma_total^2  + sigma_fm^2)

where ``sigma_total`` already combines σ_jack with the systematic floor
(``sigma_total_kg`` in the jack catalog). ``sigma_formal`` is the raw Fisher σ.

Identifiability (F2)
--------------------
``mass_status`` distinguishes ``measured`` (snr_jack ≥ 3, a genuine
determination), ``not_identifiable`` (deflection under the per-encounter noise —
the value is an upper bound, **not** a measurement) and ``unknown`` (fit run
without jackknife). Only ``measured`` rows count as consistency tests; the
``not_identifiable`` rows are reported as bounds and excluded from the |z| < 3
tally.

Caveat on the reference
-----------------------
For the Big-4 calibrators (Ceres/Vesta/Pallas/Hygiea) Fuentes-Muñoz pin GMfin to
the literature/SB441 seed, so those rows are **not** independent — they recover
the same DAWN/Goffin/Vernazza values we calibrate against. The independent
cross-check is the non-calibrators. Both populations are reported separately.

Usage
-----
    docker compose run --rm pipeline python -m scripts.mass.crosscheck_fuentes_munoz_jack \\
        --catalog data/output/orbdet/mass_catalog_jack.csv \\
        --out-dir data/output/literature_validation
"""

from __future__ import annotations

import argparse
import json
import logging
import math
from pathlib import Path

import polars as pl

from scripts.validate.validate_fuentes_munoz_masses import parse_table5_masses

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
logger = logging.getLogger(__name__)

_MRT_PATH = Path("data/raw/fuentes_munoz_2025/ajae0cc9t5_mrt.txt")
_CATALOG_PATH = Path("data/output/orbdet/mass_catalog_jack.csv")
_OUT_DIR = Path("data/output/literature_validation")


def _z_score(m: float, ref: float, sigma_ours: float | None, sigma_ref: float) -> float | None:
    """z folding our σ and the reference σ; ``None`` if σ_ours is missing/zero.

    Parameters
    ----------
    m : float
        Our fitted mass (kg).
    ref : float
        Reference mass (kg).
    sigma_ours : float or None
        Our 1-σ uncertainty (kg). ``None`` or non-finite yields ``None``.
    sigma_ref : float
        Reference 1-σ uncertainty (kg); treated as 0 if non-finite.

    Returns
    -------
    float or None
        Standardized residual ``(m - ref) / sqrt(σ_ours² + σ_ref²)``.
    """
    if sigma_ours is None or not math.isfinite(sigma_ours) or sigma_ours <= 0:
        return None
    s_ref = sigma_ref if (sigma_ref is not None and math.isfinite(sigma_ref)) else 0.0
    denom = math.sqrt(sigma_ours**2 + s_ref**2)
    return (m - ref) / denom if denom > 0 else None


def crosscheck(catalog: pl.DataFrame, fm: pl.DataFrame) -> pl.DataFrame:
    """Join the jackknife catalog to FM Table 5 and add both z-scores.

    Parameters
    ----------
    catalog : pl.DataFrame
        The jackknife mass catalog (``mass_catalog_jack.csv``).
    fm : pl.DataFrame
        Fuentes-Muñoz Table 5 masses (from :func:`parse_table5_masses`).

    Returns
    -------
    pl.DataFrame
        One row per overlapping perturber, sorted by MPC number, with columns
        ``ratio_ours_over_fm``, ``z_formal`` and ``z_total`` added.
    """
    joined = catalog.join(fm, on="perturber", how="inner").sort("perturber")

    ratios: list[float] = []
    z_formal: list[float | None] = []
    z_total: list[float | None] = []
    min_det: list[float] = []
    for r in joined.iter_rows(named=True):
        m = r["mass_fit_kg"]
        ref = r["fm_mass_kg"]
        s_ref = r["fm_sigma_kg"]
        ratios.append(m / ref if ref else float("nan"))
        z_formal.append(_z_score(m, ref, r.get("sigma_formal_kg"), s_ref))
        z_total.append(_z_score(m, ref, r.get("sigma_total_kg"), s_ref))
        # Potencia (B7): desviación fraccional mínima detectable a 3σ. Un |z|<3
        # solo excluye desviaciones mayores que esto; si la desviación observada
        # es menor, el test no tiene potencia para verla.
        s_tot = r.get("sigma_total_kg") or 0.0
        min_det.append(
            3.0 * math.sqrt(s_tot**2 + (s_ref or 0.0) ** 2) / ref if ref else float("nan")
        )

    return joined.with_columns(
        pl.Series("ratio_ours_over_fm", ratios),
        pl.Series("z_formal", z_formal, dtype=pl.Float64),
        pl.Series("z_total", z_total, dtype=pl.Float64),
        pl.Series("min_detectable_frac_3sigma", min_det, dtype=pl.Float64),
    )


def _count_within(df: pl.DataFrame, col: str) -> tuple[int, int]:
    """Count ``measured`` rows with finite ``|col| < 3`` over those with finite z."""
    meas = df.filter(pl.col("mass_status") == "measured")
    zs = [z for z in meas[col].to_list() if z is not None and math.isfinite(z)]
    return sum(1 for z in zs if abs(z) < 3), len(zs)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
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
        logger.error("Jackknife catalog not found at %s.", cat_path)
        return 1

    fm = parse_table5_masses(mrt)
    logger.info("Parsed %d numbered FM perturber masses from Table 5.", fm.height)

    ours = pl.read_csv(cat_path)
    joined = crosscheck(ours, fm)
    if joined.is_empty():
        logger.error("No perturbers overlap between the jack catalog and FM Table 5.")
        return 1
    logger.info("Overlap: %d perturbers.", joined.height)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_csv = out_dir / "fuentes_munoz_jack_mass_comparison.csv"
    joined.select(
        "perturber",
        "name",
        "is_calibrator",
        "reliable",
        "mass_status",
        "n_targets",
        "mass_fit_kg",
        "sigma_formal_kg",
        "sigma_jack_kg",
        "sigma_total_kg",
        "fm_mass_kg",
        "fm_sigma_kg",
        "ratio_ours_over_fm",
        "z_formal",
        "z_total",
        "min_detectable_frac_3sigma",
    ).write_csv(out_csv)

    # Tallies. Independent = non-calibrators (FM pins calibrators to the seed).
    indep = joined.filter(~pl.col("is_calibrator"))
    all_f = _count_within(joined, "z_formal")
    all_j = _count_within(joined, "z_total")
    ind_f = _count_within(indep, "z_formal")
    ind_j = _count_within(indep, "z_total")

    summary = {
        "reference": "Fuentes-Muñoz et al. 2025, AJ 170, 353, Table 5 (GMfin)",
        "catalog": str(cat_path),
        "n_overlap": joined.height,
        "note": (
            "z_formal uses the Fisher σ (sigma_formal_kg); z_total uses sigma_total_kg "
            "which folds the external jackknife σ (F1) with the systematic floor. Only "
            "mass_status=='measured' rows count toward the |z|<3 tally; "
            "'not_identifiable' rows are bounds, not measurements."
        ),
        "measured_within_z3": {
            "all_overlap": {"formal": all_f[0], "jack": all_j[0], "n": all_f[1]},
            "independent_noncalibrators": {"formal": ind_f[0], "jack": ind_j[0], "n": ind_f[1]},
        },
        "outputs": {"comparison_csv": str(out_csv)},
    }
    (out_dir / "fuentes_munoz_jack_mass_summary.json").write_text(json.dumps(summary, indent=2))

    print("\n" + "=" * 104)
    print("FUENTES-MUÑOZ 2025 (AJ 170, 353) — jackknife cross-check: σ_formal vs σ_jack")
    print("=" * 104)
    print(
        f"{'#':>5} {'name':<12} {'cal':>3} {'st':>4} {'M_ours':>10} {'M_FM':>10} "
        f"{'ratio':>7} {'z_form':>7} {'z_tot':>7} {'pow3σ':>6} {'status':>16}"
    )
    print("-" * 104)

    def _zt(z: float | None) -> str:
        return f"{z:+.2f}" if z is not None and math.isfinite(z) else "  n/a"

    for r in joined.iter_rows(named=True):
        st = {"measured": "M", "not_identifiable": "bnd", "unknown": "?"}.get(r["mass_status"], "?")
        print(
            f"{r['perturber']:>5} {r['name'][:12]:<12} "
            f"{'Y' if r['is_calibrator'] else '.':>3} {st:>4} "
            f"{r['mass_fit_kg']:>10.3e} {r['fm_mass_kg']:>10.3e} "
            f"{r['ratio_ours_over_fm']:>7.3f} {_zt(r['z_formal']):>7} {_zt(r['z_total']):>7} "
            f"{r['min_detectable_frac_3sigma'] * 100:>5.0f}% "
            f"{r['mass_status']:>16}"
        )
    print("-" * 104)
    print(
        f"measured within |z|<3  —  ALL overlap: formal {all_f[0]}/{all_f[1]}, "
        f"total {all_j[0]}/{all_j[1]}"
    )
    print(
        f"                          non-calibrators: formal {ind_f[0]}/{ind_f[1]}, "
        f"total {ind_j[0]}/{ind_j[1]}"
    )

    # Test de signo (B7, resultado principal para no-calibradores): con barras del
    # 20-70 %, |z|<3 no tiene potencia frente a los sesgos observados (14-30 %).
    # El signo de ratio-1 sí la tiene: bajo H0 (sin sesgo) es Binomial(n, 1/2).
    meas_ind = indep.filter(pl.col("mass_status") == "measured")
    rr = [x for x in meas_ind["ratio_ours_over_fm"].to_list() if x and math.isfinite(x)]
    if rr:
        n_below = sum(1 for x in rr if x < 1.0)
        n = len(rr)
        # p-value exacto de dos colas del test de signo
        k = min(n_below, n - n_below)
        p_two = sum(math.comb(n, i) for i in range(0, k + 1)) * 2 / 2**n
        p_two = min(1.0, p_two)
        geo = math.exp(sum(math.log(x) for x in rr) / n)
        print(
            f"sign test (measured non-calibrators): {n_below}/{n} below FM, "
            f"geometric-mean ratio {geo:.3f}, two-sided p = {p_two:.3f}"
        )
        summary["sign_test_noncalibrators"] = {
            "n_below_fm": n_below,
            "n": n,
            "geometric_mean_ratio": geo,
            "p_two_sided": p_two,
        }
        (out_dir / "fuentes_munoz_jack_mass_summary.json").write_text(json.dumps(summary, indent=2))
    print("Note: calibrators are pinned by FM to the seed → not an independent check.")
    print("      'not_identifiable' rows are bounds, excluded from the tally.")
    print("      FM 2025 fits the same Gaia FPR astrometry → errors are NOT independent;")
    print("      pow3σ is the minimum fractional deviation this test could detect at 3σ.")
    print("=" * 104)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
