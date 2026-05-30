"""Cross-match the frozen encounter catalog against Fuentes-Muñoz et al. (2025).

Source
------
Fuentes-Muñoz, Farnocchia, Giorgini & Park (2025), "Asteroid Mass Estimation by
Mutual Perturbations during Close Encounters after Gaia FPR", AJ 170, 353
(DOI 10.3847/1538-3881/ae0cc9; first presented as LPSC 2024 #2388). Machine-
readable Table 5 ("Asteroid initial masses and uncertainties"):

    https://content.cld.iop.org/journals/1538-3881/170/6/353/revision2/ajae0cc9t5_mrt.txt

Each row is a *perturber* asteroid; the ``List`` column is a pipe-delimited list
of the (first 100, highest-signal) *test* asteroids whose astrometry showed a
measurable mass signal — i.e. that had a dynamically significant close encounter
with the perturber. This is the closest thing to a machine-readable
perturber→target pair list the mass-determination literature provides (Goffin
2014's VizieR catalog J/A+A/565/A56 carries only mass tables, no encounter list).

What this validates — and what it does NOT
-------------------------------------------
Fuentes-Muñoz fit orbits over the **full Gaia FPR baseline plus all archival
astrometry**, so the encounter that produces each pair's signal can fall at *any*
epoch — frequently outside our Gaia DR3 observation window (2014-07-25 →
2017-05-28). Our catalog only contains encounters < 0.05 AU *inside* that window.
So the raw overlap is a **lower bound**, not a recall: a Fuentes-Muñoz pair absent
from our catalog usually means its closest approach happened outside the DR3
window (or beyond 0.05 AU, or outside our a∈[1.5,4.0] numbered-MBA scope), NOT
that we missed an in-window encounter. The meaningful, defensible statement is the
*positive* one: every pair present in both is an independent confirmation that our
DR3-window geometry agrees with a pair an independent Gaia-FPR mass study flagged.

Usage
-----
    # 1. fetch the table (official source)
    docker compose run --rm pipeline python -m scripts.ingest.download_fuentes_munoz
    # 2. cross-match
    docker compose run --rm pipeline python -m scripts.validate.validate_fuentes_munoz_2025
"""

from __future__ import annotations

import argparse
import json
import logging
import re
from pathlib import Path

import polars as pl

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
logger = logging.getLogger(__name__)

_MRT_PATH = Path("data/raw/fuentes_munoz_2025/ajae0cc9t5_mrt.txt")
_CATALOG_PATH = Path("data/output/encounters_catalog_rebound_005au.parquet")
_OUT_DIR = Path("data/output/literature_validation")

# Fixed-width byte ranges from the MRT byte-by-byte description (1-indexed).
# NB: trailing columns (Bibcode/Shortbib/SPKID) are NOT reliably aligned when a
# perturber's List is short, so the perturber number is taken from the Asteroid
# field, not the SPK-ID column.
_COL_ASTEROID = (0, 22)
_COL_LIST = (100, 776)

# A provisional designation: 4-digit year, space, 1-2 letters, optional digits
# (e.g. "2013 KY18", "2007 VQ345"). Numbered asteroids appear as "<n> <Name>"
# ("216 Kleopatra") or bare "<n>" when unnamed ("29943"). The provisional form
# also starts with digits, so it MUST be tested before the leading-number rule.
_PROVISIONAL = re.compile(r"^\d{4}\s[A-Z]{1,2}\d*$")
_LEADING_NUMBER = re.compile(r"^(\d+)(?:\s|$)")


def parse_table5(path: Path) -> tuple[pl.DataFrame, dict]:
    """Parse the MRT into (perturber, target) numbered pairs.

    Returns
    -------
    pairs_df:
        DataFrame with ``perturber`` (Int64), ``target`` (Int64), ``lo``, ``hi``
        (ordered key) — one row per unique numbered pair.
    stats:
        Counts for provenance (perturbers, dropped provisional designations, …).
    """
    lines = path.read_text().splitlines()
    # Data begins after the LAST '---' divider (which closes the Notes block).
    div = [i for i, ln in enumerate(lines) if ln.startswith("---")]
    data = lines[div[-1] + 1 :]

    pairs: set[tuple[int, int]] = set()
    n_perturbers_numbered = 0
    n_perturbers_provisional = 0
    n_target_tokens = 0
    n_target_provisional = 0

    for ln in data:
        if not ln.strip():
            continue
        asteroid = ln[_COL_ASTEROID[0] : _COL_ASTEROID[1]].strip()
        m = _LEADING_NUMBER.match(asteroid)
        if _PROVISIONAL.match(asteroid) or not m:
            n_perturbers_provisional += 1
            continue  # our catalog is numbered-only
        pnum = int(m.group(1))
        n_perturbers_numbered += 1
        list_field = ln[_COL_LIST[0] : _COL_LIST[1]].strip()
        for tok in (t.strip() for t in list_field.split("|")):
            if not tok or tok == "...":
                continue
            n_target_tokens += 1
            if re.fullmatch(r"\d+", tok):
                tnum = int(tok)
                pairs.add((min(pnum, tnum), max(pnum, tnum)))
            else:
                n_target_provisional += 1

    pairs_df = pl.DataFrame(
        {
            "lo": [p[0] for p in pairs],
            "hi": [p[1] for p in pairs],
        },
        schema={"lo": pl.Int64, "hi": pl.Int64},
    )
    stats = {
        "n_perturbers_numbered": n_perturbers_numbered,
        "n_perturbers_provisional_dropped": n_perturbers_provisional,
        "n_target_tokens": n_target_tokens,
        "n_target_provisional_dropped": n_target_provisional,
        "n_unique_numbered_pairs": len(pairs),
    }
    return pairs_df, stats


def main() -> int:
    parser = argparse.ArgumentParser(description="Cross-match catalog vs Fuentes-Muñoz 2025")
    parser.add_argument("--mrt", default=str(_MRT_PATH))
    parser.add_argument("--catalog", default=str(_CATALOG_PATH))
    parser.add_argument("--out-dir", default=str(_OUT_DIR))
    args = parser.parse_args()

    mrt = Path(args.mrt)
    if not mrt.exists():
        logger.error("MRT not found at %s — run scripts.ingest.download_fuentes_munoz first.", mrt)
        return 1

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    pairs_df, stats = parse_table5(mrt)
    logger.info(
        "Parsed Fuentes-Muñoz Table 5: %d numbered perturbers, %d unique numbered pairs "
        "(%d provisional targets dropped)",
        stats["n_perturbers_numbered"],
        stats["n_unique_numbered_pairs"],
        stats["n_target_provisional_dropped"],
    )

    cat = (
        pl.scan_parquet(args.catalog)
        .select(["number_1", "number_2", "dist_au", "jd_tdb"])
        .with_columns(
            pl.min_horizontal("number_1", "number_2").cast(pl.Int64).alias("lo"),
            pl.max_horizontal("number_1", "number_2").cast(pl.Int64).alias("hi"),
        )
    )
    matched = pairs_df.lazy().join(cat, on=["lo", "hi"], how="inner").collect()

    n_fm = len(pairs_df)
    n_match = matched.height
    overlap = n_match / n_fm if n_fm else 0.0
    logger.info(
        "Overlap: %d / %d Fuentes-Muñoz numbered pairs present in the DR3 catalog (%.2f%%)",
        n_match,
        n_fm,
        100 * overlap,
    )

    matched_out = matched.select(
        pl.col("lo").alias("body_1"),
        pl.col("hi").alias("body_2"),
        "dist_au",
        "jd_tdb",
    ).sort("dist_au")
    matched_path = out_dir / "fuentes_munoz_2025_matches.parquet"
    matched_out.write_parquet(matched_path)

    dist = matched.select(
        pl.col("dist_au").min().alias("min"),
        pl.col("dist_au").median().alias("median"),
        pl.col("dist_au").max().alias("max"),
    ).to_dicts()[0]

    summary = {
        "source": "Fuentes-Muñoz et al. 2025, AJ 170, 353 (DOI 10.3847/1538-3881/ae0cc9), Table 5",
        "catalog": args.catalog,
        "parse_stats": stats,
        "overlap": {
            "n_fm_numbered_pairs": n_fm,
            "n_present_in_dr3_catalog": n_match,
            "fraction": overlap,
        },
        "matched_distance_au": dist,
        "caveat": (
            "Overlap is a LOWER BOUND, not recall: Fuentes-Muñoz fit over the full "
            "Gaia FPR + archival baseline, so most pairs' close approaches fall "
            "outside the DR3 window (2014-07-25..2017-05-28), beyond 0.05 AU, or "
            "outside our a∈[1.5,4.0] numbered-MBA scope. Present-in-both pairs are "
            "positive independent confirmations."
        ),
        "outputs": {"matches": str(matched_path)},
    }
    (out_dir / "fuentes_munoz_2025_summary.json").write_text(json.dumps(summary, indent=2))

    print("\n" + "=" * 70)
    print("FUENTES-MUÑOZ et al. 2025 (AJ 170, 353) — cross-match vs DR3 catalog")
    print("=" * 70)
    print(f"FM numbered (perturber,target) pairs: {n_fm}")
    print(f"Present in our DR3-window catalog:    {n_match}  ({100 * overlap:.2f}%)")
    print(
        f"  matched dist: min={dist['min']:.6f}  median={dist['median']:.6f}  max={dist['max']:.6f} AU"
    )
    print("Note: overlap is a lower bound (FM spans the full Gaia FPR baseline; many")
    print("      encounters fall outside the DR3 window). Matches are confirmations.")
    print("=" * 70)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
