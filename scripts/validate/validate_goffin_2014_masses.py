"""Cross-check our orbdet mass catalog against Goffin (2014) ground-based masses.

This is item F7 of the mass layer: an **independent ground-based** cross-check.
Goffin (2014) derived asteroid masses from close encounters observed from the
ground (roughly 1900–2012), entirely separate from Gaia. Comparing our Gaia-FPR
orbdet masses against Goffin's tests our absolute mass scale against a fully
independent data set and method.

Source
------
Goffin E. (2014), "New determination of asteroid masses from close encounters",
A&A 565, A56 (DOI 10.1051/0004-6361/201322766). VizieR catalog J/A+A/565/A56,
``table5`` carries the derived masses. Run
``scripts.ingest.download_goffin_2014_masses`` first to produce
``data/raw/goffin_2014_masses.parquet`` (M and e_M already converted from
1e-10 solar masses to kg with ``src.orbdet.constants.M_SUN_KG``).

We join Goffin's masses with our ``mass_catalog.csv`` on perturber number and
report, per perturber, the ratio and a z-score folding both uncertainties:

    ratio = M_ours / M_goffin
    z     = (M_ours - M_goffin) / sqrt(sigma_ours^2 + sigma_goffin^2)

with ``sigma_ours = sigma_total_kg`` (our external error incl. systematic floor)
and ``sigma_goffin`` Goffin's own published 1-sigma.

Usage
-----
    docker compose run --rm pipeline python -m scripts.ingest.download_goffin_2014_masses
    docker compose run --rm pipeline python -m scripts.validate.validate_goffin_2014_masses
"""

from __future__ import annotations

import argparse
import json
import logging
import math
from pathlib import Path

import polars as pl

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
logger = logging.getLogger(__name__)

_GOFFIN_PATH = Path("data/raw/goffin_2014_masses.parquet")
_CATALOG_PATH = Path("data/output/orbdet/mass_catalog.csv")
_OUT_DIR = Path("data/output/literature_validation")

# Our 16 orbdet perturbers (MPC numbers), for reporting overlap coverage.
_OUR_PERTURBERS = {1, 2, 3, 4, 7, 10, 15, 16, 31, 52, 65, 87, 88, 107, 511, 704}


def load_goffin_masses(path: Path) -> pl.DataFrame:
    """Load the parsed Goffin (2014) mass table.

    Returns
    -------
    pl.DataFrame
        Columns ``perturber`` (Int64), ``goffin_mass_kg``, ``goffin_sigma_kg``
        and the auxiliary Goffin columns, one row per numbered asteroid.
    """
    return pl.read_parquet(path).select(
        "perturber",
        pl.col("name").alias("goffin_name"),
        "goffin_mass_kg",
        "goffin_sigma_kg",
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Cross-check orbdet masses vs Goffin 2014")
    parser.add_argument("--goffin", default=str(_GOFFIN_PATH))
    parser.add_argument("--catalog", default=str(_CATALOG_PATH))
    parser.add_argument("--out-dir", default=str(_OUT_DIR))
    args = parser.parse_args()

    goffin_path = Path(args.goffin)
    cat_path = Path(args.catalog)
    if not goffin_path.exists():
        logger.error(
            "Goffin masses not found at %s — run scripts.ingest.download_goffin_2014_masses first.",
            goffin_path,
        )
        return 1
    if not cat_path.exists():
        logger.error(
            "Mass catalog not found at %s — run scripts.mass.build_mass_catalog first.", cat_path
        )
        return 1

    goffin = load_goffin_masses(goffin_path)
    logger.info("Loaded %d Goffin masses.", goffin.height)

    ours = pl.read_csv(cat_path)
    joined = (
        ours.join(goffin, on="perturber", how="inner")
        .with_columns(
            (pl.col("mass_fit_kg") / pl.col("goffin_mass_kg")).alias("ratio_ours_over_goffin"),
            (
                (pl.col("mass_fit_kg") - pl.col("goffin_mass_kg"))
                / (pl.col("sigma_total_kg") ** 2 + pl.col("goffin_sigma_kg") ** 2).sqrt()
            ).alias("z_vs_goffin"),
        )
        .sort("perturber")
    )

    if joined.is_empty():
        logger.error("No perturbers overlap between our catalog and Goffin table5.")
        return 1

    overlap_nums = set(joined["perturber"].to_list())
    missing = sorted(_OUR_PERTURBERS - overlap_nums)
    logger.info(
        "%d/%d of our perturbers overlap with Goffin; missing from Goffin: %s",
        len(overlap_nums & _OUR_PERTURBERS),
        len(_OUR_PERTURBERS),
        missing,
    )

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_csv = out_dir / "goffin_2014_mass_comparison.csv"
    joined.select(
        "perturber",
        "name",
        "is_calibrator",
        "reliable",
        "n_targets",
        "mass_fit_kg",
        "sigma_total_kg",
        "goffin_mass_kg",
        "goffin_sigma_kg",
        "ratio_ours_over_goffin",
        "z_vs_goffin",
    ).write_csv(out_csv)

    # Independent comparison = non-calibrators (calibrators are tied to a seed/
    # literature mass and so are not a clean independent check of our pipeline).
    indep = joined.filter(~pl.col("is_calibrator"))
    reliable_indep = indep.filter(pl.col("reliable"))

    def _stats(df: pl.DataFrame) -> dict:
        if df.is_empty():
            return {"n": 0}
        ratios = df["ratio_ours_over_goffin"].to_list()
        zs = [z for z in df["z_vs_goffin"].to_list() if z is not None and not math.isnan(z)]
        return {
            "n": df.height,
            "median_ratio": float(df["ratio_ours_over_goffin"].median()),
            "mean_ratio": float(df["ratio_ours_over_goffin"].mean()),
            "min_ratio": min(ratios),
            "max_ratio": max(ratios),
            "n_within_z3": sum(1 for z in zs if abs(z) < 3),
            "n_with_z": len(zs),
        }

    summary = {
        "source": "Goffin E. 2014, A&A 565, A56, table5 (ground-based masses)",
        "catalog": str(cat_path),
        "n_overlap": joined.height,
        "our_perturbers_overlapping_goffin": sorted(overlap_nums & _OUR_PERTURBERS),
        "our_perturbers_missing_from_goffin": missing,
        "all_overlap": _stats(joined),
        "independent_noncalibrators": _stats(indep),
        "independent_reliable": _stats(reliable_indep),
        "unit_note": (
            "Goffin masses converted from 1e-10 solar masses to kg with "
            "M_SUN_KG from src.orbdet.constants; z uses Goffin's own e_M as sigma_ref."
        ),
        "note": (
            "Goffin is fully ground-based and independent of Gaia. Calibrators "
            "(Ceres/Pallas/Vesta/Hygiea) are tied to seed/literature values in our "
            "fit, so the independent_* blocks (non-calibrators) carry the real check."
        ),
        "outputs": {"comparison_csv": str(out_csv)},
    }
    (out_dir / "goffin_2014_mass_summary.json").write_text(json.dumps(summary, indent=2))

    print("\n" + "=" * 88)
    print("GOFFIN (2014) A&A 565, A56 — ground-based MASS cross-check vs orbdet catalog")
    print("=" * 88)
    print(
        f"{'#':>5} {'name':<12} {'cal':>3} {'rel':>3} "
        f"{'M_ours':>10} {'M_Goffin':>10} {'ratio':>7} {'z':>7}"
    )
    print("-" * 88)
    for r in joined.iter_rows(named=True):
        z = r["z_vs_goffin"]
        z_str = f"{z:+.2f}" if z is not None and not math.isnan(z) else "  n/a"
        print(
            f"{r['perturber']:>5} {str(r['name'])[:12]:<12} "
            f"{'Y' if r['is_calibrator'] else '.':>3} {'Y' if r['reliable'] else '.':>3} "
            f"{r['mass_fit_kg']:>10.3e} {r['goffin_mass_kg']:>10.3e} "
            f"{r['ratio_ours_over_goffin']:>7.3f} {z_str:>7}"
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
    print(f"Overlap with our 16 perturbers: {sorted(overlap_nums & _OUR_PERTURBERS)}")
    print("=" * 88)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
