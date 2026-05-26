"""Stratified sampling of the frozen catalog for Kepler-vs-N-body validation.

Reads ``data/output/encounters_catalog_rebound_005au.parquet`` and the MPCORB
snapshot that produced it (``MPCORB_20160217.DAT``).  Joins each catalog row
with the orbital elements of *both* bodies, then samples ~1000 pairs covering
the orbital-parameter space as uniformly as possible — not the "typical" pair,
but the corners (high e, high i, small Δa) where Kepler-2-body refinement is
most likely to diverge from a full N-body solution.

Stratification axes (5 bins each unless stated otherwise):
    a_1  : 1.5 – 4.0 AU                       (linear)
    e_1  : 0.0 – 0.7                          (linear; covers tail)
    i_1  : 0 – 25°                            (linear; covers tail)
    dist_au   : log10 0.0001 – 0.05           (4 bins)
    |Δa| : 0.0 – 0.5 AU                       (4 bins, prefilter cap)

Total combinatorial cells = 5⁵ × ... we cap to 200 *occupied* bins × 5 each
= ~1000 candidates.  Bin labels are stored as ``bin_id`` strings so the report
can drill in on the worst bins.

Output
------
``data/cache/nbody_validation/sample_1000.parquet`` with the catalog columns
plus the joined orbital elements for both bodies and a ``bin_id`` column.

Usage
-----
    docker compose run --rm pipeline python -m scripts.validate.sample_for_nbody_check \\
        --n-per-bin 5 --max-bins 200 --seed 42
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import polars as pl

from src.ingest.mpcorb import parse_mpcorb

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
logger = logging.getLogger(__name__)

CATALOG = Path("data/output/encounters_catalog_rebound_005au.parquet")
MPCORB = Path("data/raw/mpcorb_archive/MPCORB_20160217.DAT")
OUT = Path("data/cache/nbody_validation/sample_1000.parquet")

A_EDGES = np.array([1.5, 2.0, 2.5, 3.0, 3.5, 4.0])
E_EDGES = np.array([0.0, 0.10, 0.20, 0.30, 0.45, 0.70])
I_EDGES = np.array([0.0, 3.0, 6.0, 10.0, 15.0, 25.0])
DIST_LOG10_EDGES = np.array([-4.0, -3.0, -2.0, -1.5, np.log10(0.05) + 1e-9])
DA_EDGES = np.array([0.0, 0.05, 0.15, 0.30, 0.50 + 1e-9])


def _bin_index(values: np.ndarray, edges: np.ndarray) -> np.ndarray:
    """Return integer bin index in [0, len(edges)-2] (right-exclusive); -1 if outside."""
    idx = np.searchsorted(edges, values, side="right") - 1
    idx = np.where((values < edges[0]) | (values >= edges[-1]), -1, idx)
    return idx


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-per-bin", type=int, default=5)
    parser.add_argument("--max-bins", type=int, default=200)
    parser.add_argument("--pool-size", type=int, default=500_000,
                        help="Catalog rows to sample before binning (default 500k).")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--catalog", type=Path, default=CATALOG)
    parser.add_argument("--mpcorb", type=Path, default=MPCORB)
    parser.add_argument("--out", type=Path, default=OUT)
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)

    logger.info("Reading MPCORB elements: %s", args.mpcorb)
    elements = parse_mpcorb(args.mpcorb, only_numbered=True)
    logger.info("MPCORB elements: %d numbered asteroids", len(elements))

    logger.info("Reading frozen catalog: %s", args.catalog)
    cat = pl.read_parquet(args.catalog)
    logger.info("Catalog rows: %d", len(cat))

    # The catalog has 72M rows; joining the full table against 455k MPCORB
    # elements blows memory.  Pre-subsample to a much smaller pool that is
    # still huge relative to the ~1000-pair sample we ultimately need.  The
    # pool size needs to be large enough that high-e/high-i tail bins are
    # well-populated: 500k random rows at the empirical e>0.3 fraction
    # (~3%) gives ~15k high-e candidates — plenty.
    pool_n = min(args.pool_size, len(cat))
    cat = cat.sample(n=pool_n, seed=args.seed)
    logger.info("Subsampled catalog pool: %d rows", len(cat))

    elem_keep = elements.select(
        ["number", "a_au", "e", "i_deg", "Omega_deg", "omega_deg", "M_deg", "epoch_jd", "H"]
    )

    joined = (
        cat.join(
            elem_keep.rename(
                {
                    "number": "number_1",
                    "a_au": "a_1",
                    "e": "e_1",
                    "i_deg": "i_1",
                    "Omega_deg": "Omega_1",
                    "omega_deg": "omega_1",
                    "M_deg": "M_1",
                    "epoch_jd": "epoch_1",
                    "H": "H_1",
                }
            ),
            on="number_1",
            how="inner",
        )
        .join(
            elem_keep.rename(
                {
                    "number": "number_2",
                    "a_au": "a_2",
                    "e": "e_2",
                    "i_deg": "i_2",
                    "Omega_deg": "Omega_2",
                    "omega_deg": "omega_2",
                    "M_deg": "M_2",
                    "epoch_jd": "epoch_2",
                    "H": "H_2",
                }
            ),
            on="number_2",
            how="inner",
        )
        .with_columns(
            (pl.col("a_1") - pl.col("a_2")).abs().alias("delta_a_au"),
            (pl.col("a_1") * (1.0 - pl.col("e_1"))).alias("q_1"),
            (pl.col("a_2") * (1.0 - pl.col("e_2"))).alias("q_2"),
        )
    )
    logger.info("After joining MPCORB elements: %d rows", len(joined))

    a1 = joined["a_1"].to_numpy()
    e1 = joined["e_1"].to_numpy()
    i1 = joined["i_1"].to_numpy()
    dist = joined["dist_au"].to_numpy()
    da = joined["delta_a_au"].to_numpy()

    log10_dist = np.log10(np.clip(dist, 1e-6, None))

    bi_a = _bin_index(a1, A_EDGES)
    bi_e = _bin_index(e1, E_EDGES)
    bi_i = _bin_index(i1, I_EDGES)
    bi_d = _bin_index(log10_dist, DIST_LOG10_EDGES)
    bi_da = _bin_index(da, DA_EDGES)

    valid = (bi_a >= 0) & (bi_e >= 0) & (bi_i >= 0) & (bi_d >= 0) & (bi_da >= 0)
    logger.info("Rows in valid bin range: %d / %d", int(valid.sum()), len(joined))

    # Bin id = packed string (compact, debug-friendly).
    bin_codes = np.where(
        valid,
        np.char.add(
            np.char.add(
                np.char.add(
                    np.char.add(
                        np.char.add(
                            np.char.add(bi_a.astype(str), "_"),
                            bi_e.astype(str),
                        ),
                        "_",
                    ),
                    bi_i.astype(str),
                ),
                "_",
            ),
            np.char.add(np.char.add(bi_d.astype(str), "_"), bi_da.astype(str)),
        ),
        "INVALID",
    )

    joined = joined.with_columns(pl.Series("bin_id", bin_codes))
    joined = joined.filter(pl.col("bin_id") != "INVALID")

    counts = joined.group_by("bin_id").len().rename({"len": "n"}).sort("n", descending=True)
    logger.info("Occupied bins: %d", len(counts))
    logger.info(
        "Top-5 most-populated bins:\n%s",
        counts.head(5),
    )

    # Sample per bin (up to args.n_per_bin), then trim to args.max_bins.
    occupied_bins = counts["bin_id"].to_list()
    rng.shuffle(occupied_bins)
    selected_bins = occupied_bins[: args.max_bins]
    logger.info("Sampling from %d bins", len(selected_bins))

    parts: list[pl.DataFrame] = []
    rejected_too_small = 0
    for bid in selected_bins:
        sub = joined.filter(pl.col("bin_id") == bid)
        if len(sub) < 1:
            rejected_too_small += 1
            continue
        # Sample without replacement (or take all if bin is smaller than ask).
        n_take = min(args.n_per_bin, len(sub))
        sub = sub.sample(n=n_take, seed=int(rng.integers(0, 2**31 - 1)))
        parts.append(sub)

    if not parts:
        logger.error("No rows sampled — aborting.")
        return 1

    sample = pl.concat(parts)
    logger.info(
        "Sampled %d candidates across %d bins (rejected %d empty bins)",
        len(sample),
        len(parts),
        rejected_too_small,
    )

    # Sanity coverage stats.
    logger.info("Coverage:")
    logger.info("  e_1   p50=%.3f  p95=%.3f  max=%.3f", *np.quantile(sample["e_1"].to_numpy(), [0.5, 0.95, 1.0]))
    logger.info("  i_1   p50=%.3f  p95=%.3f  max=%.3f", *np.quantile(sample["i_1"].to_numpy(), [0.5, 0.95, 1.0]))
    logger.info("  dist  p50=%.5f  p95=%.5f  max=%.5f", *np.quantile(sample["dist_au"].to_numpy(), [0.5, 0.95, 1.0]))
    logger.info("  |Δa|  p50=%.3f  p95=%.3f  max=%.3f", *np.quantile(sample["delta_a_au"].to_numpy(), [0.5, 0.95, 1.0]))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    sample.write_parquet(args.out)
    logger.info("Wrote %s (%d rows, %d columns)", args.out, len(sample), len(sample.columns))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
