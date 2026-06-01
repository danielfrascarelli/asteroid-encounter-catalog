"""Measure the systematic error of the frozen run's N-body perturber set.

The frozen catalog's coarse scan used a 3-body N-body model (Sun + Jupiter +
Saturn) and the Kepler-refined minimum distances were validated against *that
same* model (FROZEN_RUN.md Stage A/B).  That budget is therefore *internal*: it
says nothing about how much the truncated perturber set itself shifts the
closest-approach geometry.

This script quantifies that ceiling.  It draws a stratified sample of pairs from
the frozen catalog and re-refines each one under N-body twice:

* **baseline** — the frozen perturber set ``(Sun, Jupiter, Saturn)``;
* **full**     — all eight planets ``(Sun … Neptune)``.

The distribution of ``Δdist = dist_full − dist_baseline`` is the systematic
uncertainty attributable to the missing perturbers (Uranus, Neptune and the
terrestrials), reported overall and per orbital stratum.  Major asteroids are
excluded from both configs (matching the frozen scan) so the comparison
isolates the *planetary* perturber set.

Usage
-----
    docker compose run --rm pipeline python -m scripts.validate.measure_nbody_perturber_ceiling \\
        --sample-per-stratum 20 --output data/output/nbody_perturber_ceiling
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import numpy as np
import polars as pl

from scripts.validate.refine_pair_nbody import refine_pair_nbody
from src.ingest.mpcorb import parse_mpcorb

logger = logging.getLogger(__name__)

CATALOG = Path("data/output/encounters_catalog_rebound_005au.parquet")
MPCORB = Path("data/raw/mpcorb_archive/MPCORB_20160217.DAT")

_BASELINE_PLANETS = ("sun", "jupiter", "saturn")
_FULL_PLANETS = (
    "sun",
    "mercury",
    "venus",
    "earth",
    "mars",
    "jupiter",
    "saturn",
    "uranus",
    "neptune",
)


def _classify(q_min: float, e_max: float) -> str:
    """Assign a pair to an orbital stratum (mirrors the FROZEN_RUN high-error cut)."""
    if q_min < 1.8:
        return "low_q"
    if e_max > 0.3:
        return "high_e"
    return "cold"


def _build_sample(
    n_per_stratum: int, *, subsample_every: int, seed: int
) -> tuple[pl.DataFrame, dict[int, dict]]:
    """Return (sampled pairs with elements joined, elements-by-number lookup)."""
    logger.info("Parsing MPCORB snapshot %s …", MPCORB.name)
    elements = parse_mpcorb(MPCORB)
    el = elements.select(["number", "a_au", "e"]).with_columns(
        (pl.col("a_au") * (1.0 - pl.col("e"))).alias("q_au")
    )

    logger.info("Subsampling catalog (every %d rows) …", subsample_every)
    cat = (
        pl.scan_parquet(CATALOG)
        .gather_every(subsample_every)
        .select(["number_1", "number_2", "jd_tdb", "dist_au"])
        .collect()
    )

    # Join elements for both bodies to derive q_min / e_max per pair.
    cat = (
        cat.join(
            el.rename({"number": "number_1", "a_au": "a1", "e": "e1", "q_au": "q1"}),
            on="number_1",
            how="inner",
        )
        .join(
            el.rename({"number": "number_2", "a_au": "a2", "e": "e2", "q_au": "q2"}),
            on="number_2",
            how="inner",
        )
        .with_columns(
            pl.min_horizontal("q1", "q2").alias("q_min"),
            pl.max_horizontal("e1", "e2").alias("e_max"),
        )
    )
    cat = cat.with_columns(
        pl.struct(["q_min", "e_max"])
        .map_elements(lambda s: _classify(s["q_min"], s["e_max"]), return_dtype=pl.Utf8)
        .alias("stratum")
    )

    rng = np.random.default_rng(seed)
    parts: list[pl.DataFrame] = []
    for stratum in ("cold", "high_e", "low_q"):
        pool = cat.filter(pl.col("stratum") == stratum)
        if len(pool) == 0:
            logger.warning("Stratum %s has no pairs in the subsample", stratum)
            continue
        take = min(n_per_stratum, len(pool))
        idx = rng.choice(len(pool), size=take, replace=False)
        parts.append(pool[idx.tolist()])
        logger.info("Stratum %-7s: %d available, sampled %d", stratum, len(pool), take)

    sample = pl.concat(parts)
    lookup = {int(r["number"]): r for r in elements.to_dicts()}
    return sample, lookup


def _refine(row: dict, lookup: dict[int, dict], planets: tuple[str, ...], window_hours: float):
    """Refine one pair under the given perturber set; None on failure."""
    e1 = lookup.get(int(row["number_1"]))
    e2 = lookup.get(int(row["number_2"]))
    if e1 is None or e2 is None:
        return None
    return refine_pair_nbody(
        elements_1=e1,
        elements_2=e2,
        t_center_jd=float(row["jd_tdb"]),
        window_hours=window_hours,
        include_planets=planets,
        include_major_asteroids=False,  # matches the frozen scan (majors off)
    )


def _percentiles(values: np.ndarray) -> dict[str, float]:
    if len(values) == 0:
        return {}
    a = np.abs(values)
    return {
        "n": int(len(values)),
        "median_abs": float(np.median(a)),
        "p95_abs": float(np.percentile(a, 95)),
        "p99_abs": float(np.percentile(a, 99)),
        "max_abs": float(a.max()),
        "signed_median": float(np.median(values)),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-per-stratum", type=int, default=20)
    parser.add_argument("--subsample-every", type=int, default=2000)
    parser.add_argument("--window-hours", type=float, default=6.0)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output", type=Path, default=Path("data/output/nbody_perturber_ceiling"))
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s"
    )

    sample, lookup = _build_sample(
        args.sample_per_stratum, subsample_every=args.subsample_every, seed=args.seed
    )
    logger.info("Refining %d pairs × 2 perturber sets …", len(sample))

    records = []
    for i, row in enumerate(sample.iter_rows(named=True)):
        base = _refine(row, lookup, _BASELINE_PLANETS, args.window_hours)
        full = _refine(row, lookup, _FULL_PLANETS, args.window_hours)
        if base is None or full is None:
            logger.warning("Skipping (%s, %s): missing elements", row["number_1"], row["number_2"])
            continue
        d_base = base.dist_au_nbody
        d_full = full.dist_au_nbody
        records.append(
            {
                "number_1": row["number_1"],
                "number_2": row["number_2"],
                "stratum": row["stratum"],
                "q_min": row["q_min"],
                "e_max": row["e_max"],
                "dist_kepler_catalog": row["dist_au"],
                "dist_baseline_sjs": d_base,
                "dist_full_planets": d_full,
                "delta_dist_au": d_full - d_base,
                "delta_t_min_days": full.t_min_nbody_jd - base.t_min_nbody_jd,
                "base_converged": base.converged,
                "full_converged": full.converged,
            }
        )
        if (i + 1) % 10 == 0:
            logger.info("  %d/%d done", i + 1, len(sample))

    if not records:
        logger.error("No pairs refined successfully.")
        return 1

    df = pl.DataFrame(records)
    args.output.mkdir(parents=True, exist_ok=True)
    csv_path = args.output / "perturber_ceiling_pairs.csv"
    df.write_csv(csv_path)

    deltas = df["delta_dist_au"].to_numpy()
    summary = {
        "n_pairs": len(df),
        "baseline_perturbers": list(_BASELINE_PLANETS),
        "full_perturbers": list(_FULL_PLANETS),
        "window_hours": args.window_hours,
        "seed": args.seed,
        "delta_dist_au": _percentiles(deltas),
        "by_stratum": {
            s: _percentiles(df.filter(pl.col("stratum") == s)["delta_dist_au"].to_numpy())
            for s in df["stratum"].unique().to_list()
        },
    }
    json_path = args.output / "perturber_ceiling_summary.json"
    json_path.write_text(json.dumps(summary, indent=2))

    logger.info("Wrote %s and %s", csv_path, json_path)
    logger.info(
        "Δdist (full − baseline) |median|=%.3g AU  |p95|=%.3g AU  |max|=%.3g AU",
        summary["delta_dist_au"]["median_abs"],
        summary["delta_dist_au"]["p95_abs"],
        summary["delta_dist_au"]["max_abs"],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
