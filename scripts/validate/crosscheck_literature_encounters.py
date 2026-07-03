"""Completeness cross-check: literature mass-determination pairs vs our catalog.

Scientific question
-------------------
The determination of asteroid masses from mutual close encounters relies on a
*hand-selected* set of (perturber, test-asteroid) pairs. For a close-encounter
dataset paper the natural referee question is one of **completeness**: of the
pairs that the literature actually used to measure masses, how many reappear as
sub-0.05 AU encounters in *our* catalog?

This script answers that for the primary, most complete literature source,
Fuentes-Muñoz et al. (2025), and (optionally, if extractable) for Goffin (2014).

Sources
-------
Fuentes-Muñoz, Farnocchia, Giorgini & Park (2025), "Asteroid Mass Estimation by
Mutual Perturbations during Close Encounters after Gaia FPR", AJ 170, 353. The
machine-readable Table 5 lists, per numbered perturber, a pipe-delimited ``List``
column (bytes 101-776) of the test asteroids *for which there was signal*. Per
the table's Note (5) this list is **truncated to the first 100 objects** (sorted
by decreasing signal) with a trailing ``...`` when longer. Our recovery
denominator is therefore "FM-listed pairs with a resolvable numbered target",
not FM's full ``Ntest`` population — this is a conservative, honest denominator
and is documented as such in the report.

Our catalog
-----------
``data/output/encounters_catalog_hybrid_stageb.parquet`` — the catalog consumed
by the mass engine. Columns ``number_1``, ``number_2`` (MPC numbers, either
order) and ``dist_au`` (minimum encounter distance). It already contains only
encounters below ~0.05 AU; we additionally enforce ``dist_au < THRESHOLD``.

Matching
--------
A FM pair (perturber P, target T), both numbered, is *recovered* iff the catalog
contains at least one row with {number_1, number_2} == {P, T} and
``dist_au < THRESHOLD``. Matching is order-insensitive.

Causes of non-recovery are classified as:

``provisional_target``
    The FM target is a provisional designation (e.g. ``2007 VQ345``), not a
    numbered object → out of scope for a numbered-asteroid catalog.
``encounter_above_threshold``
    The pair exists in the catalog but its minimum distance is >= THRESHOLD
    (i.e. the closest approach we found is farther than the encounter cut).
``absent``
    The pair does not appear in the catalog at all (prefilter/censoring, or the
    perturber/target was not in our propagated numbered subset).

Usage
-----
    docker compose run --rm pipeline python -m scripts.validate.crosscheck_literature_encounters
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
_CATALOG_PATH = Path("data/output/encounters_catalog_hybrid_stageb.parquet")
_OUT_DIR = Path("data/output/literature_validation")
_DEFAULT_THRESHOLD = 0.05

# Fixed-width byte ranges from the MRT byte-by-byte description (1-indexed;
# converted to 0-indexed half-open slices).
_COL_ASTEROID = (0, 22)
_COL_LIST = (100, 776)

_PROVISIONAL = re.compile(r"^\d{4}\s[A-Z]{1,2}\d*$")
_LEADING_NUMBER = re.compile(r"^(\d+)(?:\s|$)")


def _resolve_number(token: str) -> int | None:
    """Return the MPC number of a target token, or ``None`` if not numbered.

    Parameters
    ----------
    token : str
        One entry from the pipe-delimited FM ``List`` column, e.g. ``"59117"``
        or ``"2007 VQ345"``.

    Returns
    -------
    int or None
        The integer MPC number if the token is a bare number; ``None`` for
        provisional designations or the truncation marker ``...``.
    """
    token = token.strip()
    if not token or token == "...":
        return None
    if _PROVISIONAL.match(token):
        return None
    return int(token) if token.isdigit() else None


def parse_table5_targets(path: Path) -> pl.DataFrame:
    """Parse FM Table 5 into one row per (perturber, listed target).

    Parameters
    ----------
    path : Path
        Path to the FM machine-readable Table 5 (``ajae0cc9t5_mrt.txt``).

    Returns
    -------
    pl.DataFrame
        Columns:

        ``perturber`` : Int64
            Numbered perturber MPC number.
        ``perturber_name`` : String
            Perturber name (as printed after the number).
        ``target_token`` : String
            Raw target token from the pipe-delimited list.
        ``target_number`` : Int64 (nullable)
            Resolved numbered target, ``None`` for provisional designations.
        ``truncated`` : Boolean
            Whether this perturber's list was truncated (ended with ``...``).
    """
    lines = path.read_text().splitlines()
    div = [i for i, ln in enumerate(lines) if ln.startswith("---")]
    data = lines[div[-1] + 1 :]

    rows: list[dict[str, object]] = []
    for ln in data:
        if not ln.strip():
            continue
        asteroid = ln[_COL_ASTEROID[0] : _COL_ASTEROID[1]].strip()
        m = _LEADING_NUMBER.match(asteroid)
        if _PROVISIONAL.match(asteroid) or not m:
            continue  # numbered perturbers only
        perturber = int(m.group(1))
        # Name = asteroid string minus the leading number.
        name = asteroid[m.end() :].strip() or str(perturber)

        raw_list = ln[_COL_LIST[0] : _COL_LIST[1]]
        tokens = [t for t in (tok.strip() for tok in raw_list.split("|")) if t]
        truncated = bool(tokens) and tokens[-1] == "..."

        for tok in tokens:
            if tok == "...":
                continue
            rows.append(
                {
                    "perturber": perturber,
                    "perturber_name": name,
                    "target_token": tok,
                    "target_number": _resolve_number(tok),
                    "truncated": truncated,
                }
            )

    return pl.DataFrame(
        rows,
        schema={
            "perturber": pl.Int64,
            "perturber_name": pl.String,
            "target_token": pl.String,
            "target_number": pl.Int64,
            "truncated": pl.Boolean,
        },
    )


def matched_pairs_from_catalog(
    catalog_path: Path,
    pairs: pl.DataFrame,
    threshold: float,
) -> tuple[pl.DataFrame, set[int]]:
    """Find which FM numbered pairs appear in the catalog below ``threshold``.

    The catalog is scanned lazily and filtered *early* to rows whose both
    numbers fall in the small set of FM perturbers/targets, so the 72M-row
    catalog is never materialised in full. A second lazy pass computes the
    catalog *universe* — every asteroid number that appears in any catalog row
    — so that non-recovered pairs can be attributed to a missing object vs a
    genuinely absent close encounter.

    Parameters
    ----------
    catalog_path : Path
        Parquet catalog with ``number_1``, ``number_2``, ``dist_au``.
    pairs : pl.DataFrame
        FM numbered pairs; must have ``perturber`` and ``target_number``
        (non-null) columns.
    threshold : float
        Encounter distance cut in AU.

    Returns
    -------
    joined : pl.DataFrame
        Columns ``lo``, ``hi`` (the order-insensitive number pair), and
        ``min_dist_au`` = minimum ``dist_au`` found for that pair over *all*
        catalog rows (regardless of threshold), plus ``present`` (bool: pair
        exists at any distance in the scanned/relevant subset) and
        ``recovered`` (bool: exists with ``min_dist_au < threshold``).
    universe : set of int
        Every distinct asteroid number that appears anywhere in the catalog
        (as ``number_1`` or ``number_2``).
    """
    perts = pairs["perturber"].unique().to_list()
    tgts = pairs.filter(pl.col("target_number").is_not_null())["target_number"].unique().to_list()
    relevant = set(perts) | set(tgts)

    # Canonical order-insensitive FM pair keys.
    fm_keys = (
        pairs.filter(pl.col("target_number").is_not_null())
        .select(
            pl.min_horizontal("perturber", "target_number").alias("lo"),
            pl.max_horizontal("perturber", "target_number").alias("hi"),
        )
        .unique()
    )

    lf = pl.scan_parquet(catalog_path)
    # Early filter: keep only rows where BOTH numbers are in the FM-relevant set.
    # This collapses 72M rows to a few thousand before any pairwise work.
    sub = (
        lf.filter(pl.col("number_1").is_in(relevant) & pl.col("number_2").is_in(relevant))
        .select(
            pl.min_horizontal("number_1", "number_2").cast(pl.Int64).alias("lo"),
            pl.max_horizontal("number_1", "number_2").cast(pl.Int64).alias("hi"),
            pl.col("dist_au"),
        )
        .group_by("lo", "hi")
        .agg(pl.col("dist_au").min().alias("min_dist_au"))
        .collect(engine="streaming")
    )

    # Catalog universe (streamed union of both number columns).
    universe_df = (
        pl.concat(
            [
                lf.select(pl.col("number_1").cast(pl.Int64).alias("n")),
                lf.select(pl.col("number_2").cast(pl.Int64).alias("n")),
            ]
        )
        .unique()
        .collect(engine="streaming")
    )
    universe = set(universe_df["n"].to_list())

    joined = fm_keys.join(sub, on=["lo", "hi"], how="left").with_columns(
        pl.col("min_dist_au").is_not_null().alias("present"),
        (pl.col("min_dist_au") < threshold).fill_null(False).alias("recovered"),
    )
    return joined, universe


def build_report(
    pairs: pl.DataFrame,
    matched: pl.DataFrame,
    universe: set[int],
    threshold: float,
) -> tuple[pl.DataFrame, pl.DataFrame, dict]:
    """Assemble the per-pair classification, per-perturber summary, and totals.

    Parameters
    ----------
    pairs : pl.DataFrame
        Output of :func:`parse_table5_targets`.
    matched : pl.DataFrame
        Output of :func:`matched_pairs_from_catalog` (the ``joined`` element).
    universe : set of int
        Catalog universe (every number appearing in any catalog row).
    threshold : float
        Encounter distance cut in AU (for reporting).

    Returns
    -------
    pair_class : pl.DataFrame
        One row per FM (perturber, target) with its outcome classification.
    per_pert : pl.DataFrame
        Per-perturber recovery summary.
    totals : dict
        Global recovery statistics and cause breakdown.

    Notes
    -----
    Outcome categories (numbered targets):

    ``recovered``
        Pair present with ``min_dist_au < threshold``.
    ``encounter_above_threshold``
        Pair present in the catalog but its minimum distance is >= threshold.
    ``absent_perturber_missing``
        Not present; the perturber never appears in the catalog universe.
    ``absent_target_missing``
        Not present; the target never appears in the catalog universe (but the
        perturber does).
    ``absent_no_close_encounter``
        Not present; both objects are in the catalog universe, so the pair
        simply never came within the encounter cut in our propagation. This is
        the *expected* residual: FM's "signal" list is not distance-limited to
        0.05 AU, so many FM pairs are genuinely wider encounters.
    """
    universe_list = list(universe)

    # Numbered pairs only carry a canonical key we can join on.
    numbered = pairs.filter(pl.col("target_number").is_not_null()).with_columns(
        pl.min_horizontal("perturber", "target_number").alias("lo"),
        pl.max_horizontal("perturber", "target_number").alias("hi"),
    )
    numbered = numbered.join(
        matched.select("lo", "hi", "min_dist_au", "present", "recovered"),
        on=["lo", "hi"],
        how="left",
    ).with_columns(
        pl.col("perturber").is_in(universe_list).alias("_pert_in_universe"),
        pl.col("target_number").is_in(universe_list).alias("_tgt_in_universe"),
    )

    numbered = numbered.with_columns(
        pl.when(pl.col("recovered"))
        .then(pl.lit("recovered"))
        .when(pl.col("present"))
        .then(pl.lit("encounter_above_threshold"))
        .when(~pl.col("_pert_in_universe"))
        .then(pl.lit("absent_perturber_missing"))
        .when(~pl.col("_tgt_in_universe"))
        .then(pl.lit("absent_target_missing"))
        .otherwise(pl.lit("absent_no_close_encounter"))
        .alias("outcome")
    )

    provisional = pairs.filter(pl.col("target_number").is_null()).with_columns(
        pl.lit(None, dtype=pl.Int64).alias("lo"),
        pl.lit(None, dtype=pl.Int64).alias("hi"),
        pl.lit(None, dtype=pl.Float64).alias("min_dist_au"),
        pl.lit(False).alias("present"),
        pl.lit(False).alias("recovered"),
        pl.lit("provisional_target").alias("outcome"),
    )

    pair_class = pl.concat(
        [
            numbered.select(
                "perturber",
                "perturber_name",
                "target_token",
                "target_number",
                "truncated",
                "min_dist_au",
                "outcome",
            ),
            provisional.select(
                "perturber",
                "perturber_name",
                "target_token",
                "target_number",
                "truncated",
                "min_dist_au",
                "outcome",
            ),
        ]
    )

    # Per-perturber summary over NUMBERED targets (the resolvable denominator).
    per_pert = (
        numbered.group_by("perturber", "perturber_name")
        .agg(
            pl.len().alias("n_numbered_targets"),
            (pl.col("outcome") == "recovered").sum().alias("n_recovered"),
            (pl.col("outcome") == "encounter_above_threshold").sum().alias("n_above_thresh"),
            pl.col("outcome").str.starts_with("absent").sum().alias("n_absent"),
            pl.col("truncated").max().alias("list_truncated"),
        )
        .with_columns(
            (pl.col("n_recovered") / pl.col("n_numbered_targets") * 100)
            .round(1)
            .alias("pct_recovered")
        )
        .sort("n_numbered_targets", descending=True)
    )

    n_total_listed = pairs.height
    n_provisional = pairs.filter(pl.col("target_number").is_null()).height
    n_numbered = numbered.height
    counts = {
        row["outcome"]: row["len"]
        for row in numbered.group_by("outcome").len().iter_rows(named=True)
    }
    n_recovered = counts.get("recovered", 0)

    totals = {
        "threshold_au": threshold,
        "n_perturbers": int(pairs["perturber"].n_unique()),
        "n_listed_pairs_total": n_total_listed,
        "n_provisional_targets": n_provisional,
        "n_numbered_pairs": n_numbered,
        "n_recovered": n_recovered,
        "n_encounter_above_threshold": counts.get("encounter_above_threshold", 0),
        "n_absent_perturber_missing": counts.get("absent_perturber_missing", 0),
        "n_absent_target_missing": counts.get("absent_target_missing", 0),
        "n_absent_no_close_encounter": counts.get("absent_no_close_encounter", 0),
        "recovery_fraction_of_numbered": (n_recovered / n_numbered) if n_numbered else 0.0,
        "recovery_fraction_of_all_listed": (
            (n_recovered / n_total_listed) if n_total_listed else 0.0
        ),
    }
    return pair_class, per_pert, totals


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Completeness cross-check of literature mass pairs vs our catalog"
    )
    parser.add_argument("--mrt", default=str(_MRT_PATH))
    parser.add_argument("--catalog", default=str(_CATALOG_PATH))
    parser.add_argument("--out-dir", default=str(_OUT_DIR))
    parser.add_argument("--threshold", type=float, default=_DEFAULT_THRESHOLD)
    args = parser.parse_args()

    mrt = Path(args.mrt)
    catalog = Path(args.catalog)
    if not mrt.exists():
        logger.error("FM MRT not found at %s", mrt)
        return 1
    if not catalog.exists():
        logger.error("Catalog not found at %s", catalog)
        return 1

    pairs = parse_table5_targets(mrt)
    logger.info(
        "Parsed %d listed (perturber, target) pairs across %d numbered perturbers.",
        pairs.height,
        pairs["perturber"].n_unique(),
    )

    matched, universe = matched_pairs_from_catalog(catalog, pairs, args.threshold)
    logger.info(
        "Scanned catalog; %d FM numbered pairs present; %d distinct numbers in universe.",
        int(matched["present"].sum()),
        len(universe),
    )

    pair_class, per_pert, totals = build_report(pairs, matched, universe, args.threshold)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    pair_class.write_csv(out_dir / "completeness_pairs.csv")
    per_pert.write_csv(out_dir / "completeness_per_perturber.csv")
    (out_dir / "completeness_summary.json").write_text(json.dumps(totals, indent=2))

    print("\n" + "=" * 78)
    print(
        f"COMPLETENESS vs Fuentes-Muñoz et al. 2025 (numbered targets, dist < {args.threshold:.3f} AU)"
    )
    print("=" * 78)
    print(
        f"Numbered FM pairs: {totals['n_numbered_pairs']}   "
        f"recovered: {totals['n_recovered']}   "
        f"({totals['recovery_fraction_of_numbered'] * 100:.1f}%)"
    )
    print("  non-recovery causes (numbered pairs):")
    print(f"    encounter >= threshold in catalog : {totals['n_encounter_above_threshold']}")
    print(f"    perturber absent from catalog     : {totals['n_absent_perturber_missing']}")
    print(f"    target absent from catalog        : {totals['n_absent_target_missing']}")
    print(f"    both present, no <thr encounter   : {totals['n_absent_no_close_encounter']}")
    print(f"  provisional targets (out of scope)  : {totals['n_provisional_targets']}")
    print("-" * 78)
    print(
        f"{'#':>5} {'name':<14} {'Nnum':>5} {'rec':>5} {'>thr':>5} {'abs':>5} {'%rec':>6} {'trunc':>5}"
    )
    print("-" * 78)
    for r in per_pert.iter_rows(named=True):
        print(
            f"{r['perturber']:>5} {r['perturber_name'][:14]:<14} "
            f"{r['n_numbered_targets']:>5} {r['n_recovered']:>5} "
            f"{r['n_above_thresh']:>5} {r['n_absent']:>5} {r['pct_recovered']:>6.1f} "
            f"{'Y' if r['list_truncated'] else '.':>5}"
        )
    print("=" * 78)
    logger.info("Wrote outputs to %s", out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
