"""Filter novel asteroid encounters to identify mass-determination candidates.

Reads ``data/output/novel_encounters_not_in_literature.csv`` (output of the
cross-match against published literature) and applies astrophysical criteria
to retain only encounters with meaningful gravitational deflection potential.

Output: ``data/output/relevant_novel_encounters.csv``

Usage
-----
    docker compose run --rm pipeline python encounter_analysis/filter_candidates.py
    docker compose run --rm pipeline python encounter_analysis/filter_candidates.py --input path/to/other.csv
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import polars as pl

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Thresholds — edit here or override via CLI arguments
# ---------------------------------------------------------------------------

_GAIA_OBSERVABLE_ONLY: bool = True
_MAX_DIST_AU: float = 0.02
_MIN_DIAMETER_KM: float = 30.0
_MAX_VEL_KM_S: float = 8.0

_INPUT_DEFAULT = Path("data/output/novel_encounters_not_in_literature.csv")
_OUTPUT_DEFAULT = Path("data/output/relevant_novel_encounters.csv")


# ---------------------------------------------------------------------------
# Score
# ---------------------------------------------------------------------------


def deflection_score(df: pl.DataFrame) -> pl.Series:
    """Proxy for gravitational deflection strength.

    δ ∝ M / (v² · b)  where  M ∝ D³,  b = dist_au,  v = rel_vel_km_s
    → score = diameter_1_km³ / (dist_au · rel_vel_km_s²)

    Higher score = stronger expected deflection = better mass-determination candidate.
    """
    return pl.col("diameter_1_km") ** 3 / (pl.col("dist_au") * pl.col("rel_vel_km_s") ** 2)


# ---------------------------------------------------------------------------
# Main filter
# ---------------------------------------------------------------------------


def filter_candidates(
    df: pl.DataFrame,
    *,
    gaia_observable_only: bool = _GAIA_OBSERVABLE_ONLY,
    max_dist_au: float = _MAX_DIST_AU,
    min_diameter_km: float = _MIN_DIAMETER_KM,
    max_vel_km_s: float = _MAX_VEL_KM_S,
) -> pl.DataFrame:
    """Apply astrophysical filters and add derived columns.

    Parameters
    ----------
    df:
        Raw novel-encounters DataFrame.
    gaia_observable_only:
        Keep only encounters where ``gaia_observable == true``.
    max_dist_au:
        Maximum closest-approach distance (AU).
    min_diameter_km:
        Minimum diameter of the perturber (body_1) in km.
    max_vel_km_s:
        Maximum relative velocity (km/s).

    Returns
    -------
    pl.DataFrame
        Filtered and enriched DataFrame, sorted by ``deflection_score`` desc.
    """
    n_raw = len(df)

    # Normalise gaia_observable: may be bool or string "true"/"false"
    if df["gaia_observable"].dtype == pl.Utf8:
        df = df.with_columns(
            pl.col("gaia_observable").str.to_lowercase().eq("true").alias("gaia_observable")
        )

    mask = (
        (pl.col("dist_au") < max_dist_au)
        & (pl.col("diameter_1_km") > min_diameter_km)
        & (pl.col("rel_vel_km_s") < max_vel_km_s)
    )
    if gaia_observable_only:
        mask = mask & pl.col("gaia_observable")

    filtered = df.filter(mask)

    # Add deflection score
    filtered = filtered.with_columns(deflection_score(filtered).alias("deflection_score"))

    # Flag encounters where the perturber's mass is not in the literature
    # (diameter_1_km is present but no published mass — approximated by
    # checking against well-known massive bodies)
    _KNOWN_MASS_NUMBERS = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 16, 52, 65, 87, 88, 107, 121, 128}
    filtered = filtered.with_columns(
        pl.col("number_1").is_in(list(_KNOWN_MASS_NUMBERS)).not_().alias("mass_unknown")
    )

    filtered = filtered.sort("deflection_score", descending=True)

    logger.info(
        "Filtered %d → %d encounters (dist < %.3f AU, D₁ > %.0f km, v < %.1f km/s, gaia_observable=%s)",
        n_raw,
        len(filtered),
        max_dist_au,
        min_diameter_km,
        max_vel_km_s,
        gaia_observable_only,
    )
    return filtered


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--input", type=Path, default=_INPUT_DEFAULT)
    p.add_argument("--output", type=Path, default=_OUTPUT_DEFAULT)
    p.add_argument("--max-dist-au", type=float, default=_MAX_DIST_AU)
    p.add_argument("--min-diameter-km", type=float, default=_MIN_DIAMETER_KM)
    p.add_argument("--max-vel-km-s", type=float, default=_MAX_VEL_KM_S)
    p.add_argument(
        "--all-observability", action="store_true", help="Include encounters not observable by Gaia"
    )
    args = p.parse_args()

    if not args.input.exists():
        p.error(f"Input file not found: {args.input}")

    logger.info("Reading %s …", args.input)
    df = pl.read_csv(args.input)
    logger.info("Loaded %d encounters", len(df))

    result = filter_candidates(
        df,
        gaia_observable_only=not args.all_observability,
        max_dist_au=args.max_dist_au,
        min_diameter_km=args.min_diameter_km,
        max_vel_km_s=args.max_vel_km_s,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.write_csv(args.output)
    logger.info("Saved %d candidates → %s", len(result), args.output)

    # Quick summary
    cat_a = result.filter(~pl.col("mass_unknown"))
    cat_b = result.filter(pl.col("mass_unknown"))
    logger.info("  Category A (known mass): %d encounters", len(cat_a))
    logger.info("  Category B (unknown mass): %d encounters", len(cat_b))
    if len(result) > 0:
        top = result.row(0, named=True)
        logger.info(
            "  Top candidate: (%s) + (%s)  dist=%.5f AU  score=%.2e",
            top["designation_1"],
            top["designation_2"],
            top["dist_au"],
            top["deflection_score"],
        )


if __name__ == "__main__":
    main()
