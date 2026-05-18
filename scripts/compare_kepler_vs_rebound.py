"""Side-by-side comparison of Kepler vs N-body encounter catalogs.

Joins the two catalogs (``encounters_catalog_005au.parquet`` from the Kepler
2-body run and ``encounters_catalog_rebound_005au.parquet`` from the N-body
run) on (number_1, number_2, sign(jd_tdb)) and produces:

1. Pairs only in Kepler (false positives in 2-body model).
2. Pairs only in rebound (true encounters missed by 2-body).
3. Pairs in both — with Δdist and Δjd for each, so we can characterise the
   typical correction introduced by Jupiter+Saturn perturbations.

Also re-runs the Fienga and Galád spot-checks and reports, for each matched
literature event, the Kepler distance vs the rebound distance, JPL ground
truth, and the residual reduction (or otherwise).

Output: ``data/output/kepler_vs_rebound_comparison.csv`` + summary log.

Usage
-----
    docker compose run --rm pipeline python -m scripts.compare_kepler_vs_rebound
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import polars as pl

from src.utils.config import load_config

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)


def _ordered_pair(df: pl.DataFrame) -> pl.DataFrame:
    """Ensure (number_1, number_2) with number_1 < number_2 for join compatibility."""
    return df.with_columns(
        [
            pl.min_horizontal("number_1", "number_2").alias("n_lo"),
            pl.max_horizontal("number_1", "number_2").alias("n_hi"),
        ]
    )


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--config", default="config.yaml")
    p.add_argument("--kepler", default="encounters_catalog_005au.parquet")
    p.add_argument("--rebound", default="encounters_catalog_rebound_005au.parquet")
    args = p.parse_args()

    cfg = load_config(args.config)
    out_dir = Path(cfg.paths.output)
    kepler_path = out_dir / args.kepler
    rebound_path = out_dir / args.rebound

    if not kepler_path.exists():
        logger.error("Missing Kepler catalog: %s", kepler_path)
        return 1
    if not rebound_path.exists():
        logger.error("Missing rebound catalog: %s — run the rebound pipeline first.", rebound_path)
        return 1

    kep = _ordered_pair(pl.read_parquet(kepler_path)).rename(
        {"jd_tdb": "jd_kep", "dist_au": "dist_kep"}
    )
    reb = _ordered_pair(pl.read_parquet(rebound_path)).rename(
        {"jd_tdb": "jd_reb", "dist_au": "dist_reb"}
    )

    logger.info("Kepler catalog : %d encounters", len(kep))
    logger.info("rebound catalog: %d encounters", len(reb))

    # Outer join on the ordered pair: lets us count only-in-A, only-in-B, both.
    joined = kep.join(reb, on=["n_lo", "n_hi"], how="full", coalesce=True)
    only_kep = joined.filter(pl.col("dist_reb").is_null())
    only_reb = joined.filter(pl.col("dist_kep").is_null())
    both = joined.filter(pl.col("dist_kep").is_not_null() & pl.col("dist_reb").is_not_null())

    logger.info("")
    logger.info("=== Set intersection ===")
    logger.info("Only in Kepler  (potential 2-body false positives): %d", len(only_kep))
    logger.info("Only in rebound (true encounters Kepler missed)   : %d", len(only_reb))
    logger.info("In both                                            : %d", len(both))

    if len(both) > 0:
        delta_dist = (both["dist_reb"] - both["dist_kep"]).abs()
        delta_jd = (both["jd_reb"] - both["jd_kep"]).abs()

        def _f(v: object) -> float:
            return float(v) if isinstance(v, (int, float)) else 0.0

        logger.info("")
        logger.info("=== Distance deltas on shared pairs ===")
        logger.info(
            "  |Δdist| (AU): mean=%.5f  med=%.5f  p95=%.5f  max=%.5f",
            _f(delta_dist.mean()),
            _f(delta_dist.median()),
            _f(delta_dist.quantile(0.95)),
            _f(delta_dist.max()),
        )
        logger.info(
            "  |Δt|   (d) : mean=%.2f  med=%.2f  p95=%.2f  max=%.2f",
            _f(delta_jd.mean()),
            _f(delta_jd.median()),
            _f(delta_jd.quantile(0.95)),
            _f(delta_jd.max()),
        )

    # Save shared-pair comparison
    both_out = both.select(
        [
            pl.col("n_lo").alias("number_1"),
            pl.col("n_hi").alias("number_2"),
            pl.col("dist_kep").alias("dist_kepler_au"),
            pl.col("dist_reb").alias("dist_rebound_au"),
            (pl.col("dist_reb") - pl.col("dist_kep")).alias("delta_au"),
            pl.col("jd_kep").alias("jd_kepler_tdb"),
            pl.col("jd_reb").alias("jd_rebound_tdb"),
        ]
    ).sort(pl.col("delta_au").abs(), descending=True)

    cmp_path = out_dir / "kepler_vs_rebound_comparison.csv"
    both_out.write_csv(cmp_path)
    logger.info("Wrote %d shared-pair comparisons to %s", len(both_out), cmp_path)

    # Show the top 10 largest absolute deltas — these are the pairs most affected
    # by the addition of Jupiter+Saturn perturbations.
    if len(both_out) > 0:
        logger.info("")
        logger.info("=== Top 10 |Δdist| pairs (where rebound moves the answer most) ===")
        for row in both_out.head(10).iter_rows(named=True):
            logger.info(
                "  (%d, %d)  Kepler=%.5f AU  rebound=%.5f AU  Δ=%+.5f AU",
                row["number_1"],
                row["number_2"],
                row["dist_kepler_au"],
                row["dist_rebound_au"],
                row["delta_au"],
            )

    # ------------------------------------------------------------------------
    # Re-validate the literature-matched events against the rebound catalog
    # ------------------------------------------------------------------------
    fienga_csv = out_dir / "fienga_2003_matches.csv"
    galad_csv = out_dir / "galad_2002_matches.csv"

    def _lookup(df: pl.DataFrame, a: int, b: int) -> tuple[float, float] | None:
        lo, hi = min(a, b), max(a, b)
        hit = df.filter((pl.col("n_lo") == lo) & (pl.col("n_hi") == hi))
        if len(hit) == 0:
            return None
        r = hit.row(0, named=True)
        return float(r["dist_reb"]), float(r["jd_reb"])

    lit_rows: list[dict] = []
    for csv_path, label, dist_col in (
        (fienga_csv, "Fienga 2003", "fienga_impact_au"),
        (galad_csv, "Galád 2002", "galad_r_au"),
    ):
        if not csv_path.exists():
            continue
        for r in pl.read_csv(csv_path).iter_rows(named=True):
            a, b = int(r["perturber"]), int(r["target"])
            kep_dist = float(r["our_dist_au"])
            lit_dist = float(r[dist_col])
            reb_hit = _lookup(reb, a, b)
            reb_dist = reb_hit[0] if reb_hit else float("nan")
            lit_rows.append(
                {
                    "source": label,
                    "perturber": a,
                    "target": b,
                    "literature_au": lit_dist,
                    "kepler_au": kep_dist,
                    "rebound_au": reb_dist,
                    "delta_kepler_minus_lit": kep_dist - lit_dist,
                    "delta_rebound_minus_lit": (reb_dist - lit_dist) if reb_hit else float("nan"),
                }
            )

    if lit_rows:
        lit_df = pl.DataFrame(lit_rows)
        lit_path = out_dir / "literature_kepler_vs_rebound.csv"
        lit_df.write_csv(lit_path)
        logger.info("")
        logger.info("=== Literature events: Kepler vs rebound (vs literature) ===")
        for r in lit_df.iter_rows(named=True):
            logger.info(
                "  %-12s (%d, %d)  lit=%.5f  Kepler=%.5f (Δ=%+.5f)  rebound=%.5f (Δ=%+.5f)",
                r["source"],
                r["perturber"],
                r["target"],
                r["literature_au"],
                r["kepler_au"],
                r["delta_kepler_minus_lit"],
                r["rebound_au"],
                r["delta_rebound_minus_lit"],
            )
        logger.info("Saved → %s", lit_path)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
