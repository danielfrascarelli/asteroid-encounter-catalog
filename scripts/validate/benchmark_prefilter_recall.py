"""Quantify the recall of the orbital prefilter on the adverse high-e/high-i tail.

The orbital prefilter (:func:`src.detect.prefilter.compatible_pairs`) drops any
pair whose semimajor axes differ by more than ``semimajor_diff_max_au`` (0.5 AU)
OR whose inclinations differ by more than ``inclination_diff_max_deg`` (30°)
*before* the KD-tree spatial scan, for cost.  It is a heuristic: high-eccentricity
or high-inclination crossing orbits can come within the encounter threshold even
when Δa or Δi is large, so a real encounter can be dropped.  Its recall on that
adverse tail has never been measured — which is why the frozen catalog cannot
claim completeness (FROZEN_RUN.md caveat #2 / audit blocker #2).

Method
------
The prefilter is a **pure, deterministic mask on orbital elements**: a pair is
kept iff ``|Δa| ≤ 0.5 AU AND |Δi| ≤ 30°``.  Running detection *with* the
prefilter is therefore *identical* to running it *without* the prefilter and
then intersecting the resulting encounter set with the mask — the KD-tree scan
and the refinement are byte-for-byte the same on the surviving pairs.  We
exploit this:

1. Build the adverse subset: numbered asteroids in the frozen-catalog scope
   (``a ∈ [1.5, 4.0] AU``, ``only_numbered``) restricted to ``e > 0.3`` OR
   ``i > 15°`` — the regime where the Δa/Δi heuristic is most likely to fail.
2. Run detection **without** the prefilter over the full Gaia DR3 window
   (``pairs=None`` → KD-tree spatial query only, O(N log N) per step, so the
   full ~52 k-body adverse population is tractable without the O(N²) pair
   materialisation).  This is the ground-truth encounter set ``F`` (all pairs
   whose Kepler-refined minimum distance ≤ 0.05 AU).
3. Apply the prefilter mask analytically to ``F`` → ``P`` (what the prefilter
   keeps).  ``recall = |P| / |F|``.
4. ``--cross-check N``: additionally run the *real* pipeline with the prefilter
   enabled on a ≤5000-body sub-sample and assert it reproduces the analytic
   mask, proving the equivalence holds in the production code path.
5. Characterise the missed pairs (``F \\ P``) by (Δa, Δi, e_max, i_max) and
   report recall per band with a Wilson binomial confidence interval.

Propagation uses the **Kepler 2-body** model throughout (``method=kepler``).
This matches the model that produces the *final* reported distances in the
frozen catalog (the rebound trajectory is used only for the coarse scan there;
the refinement — hence every reported minimum distance — is Kepler).  The
recall *ratio* is robust to the propagation model regardless, since both the
with- and without-prefilter sets are computed under the same model.

Usage
-----
    docker compose run --rm pipeline python -m scripts.validate.benchmark_prefilter_recall
    docker compose run --rm pipeline python -m scripts.validate.benchmark_prefilter_recall \
        --max-bodies 8000 --cross-check 4000
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import time
from pathlib import Path

import polars as pl
from astropy.time import Time

from src.detect.pipeline import detect_encounters
from src.ingest.mpcorb import parse_mpcorb
from src.ingest.mpcorb_archive import discover_snapshots, select_for_window
from src.propagate.grid import make_time_grid
from src.utils.config import load_config

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)

_AU_KM = 1.495_978_707e8


def _wilson_ci(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    """Wilson score interval for a binomial proportion k/n (default 95 %)."""
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    denom = 1.0 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = (z / denom) * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (max(0.0, centre - half), min(1.0, centre + half))


def _select_adverse(
    snapshot_path: Path,
    *,
    a_min: float,
    a_max: float,
    e_min: float,
    i_min_deg: float,
    max_bodies: int | None,
    seed: int,
) -> pl.DataFrame:
    """Numbered MBAs in [a_min, a_max] with e > e_min OR i > i_min_deg."""
    df = parse_mpcorb(
        snapshot_path,
        only_numbered=True,
        semimajor_min_au=a_min,
        semimajor_max_au=a_max,
    )
    adverse = df.filter((pl.col("e") > e_min) | (pl.col("i_deg") > i_min_deg))
    logger.info(
        "Adverse subset: %d / %d numbered bodies in a∈[%.1f,%.1f] (e>%.2f OR i>%.0f°)",
        len(adverse),
        len(df),
        a_min,
        a_max,
        e_min,
        i_min_deg,
    )
    if max_bodies is not None and len(adverse) > max_bodies:
        adverse = adverse.sample(n=max_bodies, seed=seed)
        logger.info("Sub-sampled adverse subset to %d bodies (seed=%d)", max_bodies, seed)
    # Stable order so prefilter row indices are deterministic.
    return adverse.sort("number")


def _pair_key_frame(df: pl.DataFrame) -> pl.DataFrame:
    """Add an ordered (lo, hi) MPC-number key column for set membership."""
    return df.with_columns(
        pl.min_horizontal("number_1", "number_2").alias("_lo"),
        pl.max_horizontal("number_1", "number_2").alias("_hi"),
    )


def _annotate_orbital_deltas(encounters: pl.DataFrame, elements: pl.DataFrame) -> pl.DataFrame:
    """Join per-body (a, e, i) onto each encounter and compute Δa, Δi, e_max, i_max."""
    elem = elements.select(
        pl.col("number"),
        pl.col("a_au"),
        pl.col("e"),
        pl.col("i_deg"),
    )
    out = (
        encounters.join(
            elem.rename({"number": "number_1", "a_au": "a1", "e": "e1", "i_deg": "i1"}),
            on="number_1",
            how="left",
        )
        .join(
            elem.rename({"number": "number_2", "a_au": "a2", "e": "e2", "i_deg": "i2"}),
            on="number_2",
            how="left",
        )
        .with_columns(
            (pl.col("a1") - pl.col("a2")).abs().alias("delta_a_au"),
            (pl.col("i1") - pl.col("i2")).abs().alias("delta_i_deg"),
            pl.max_horizontal("e1", "e2").alias("e_max"),
            pl.max_horizontal("i1", "i2").alias("i_max_deg"),
        )
    )
    return out


def _passes_prefilter(df: pl.DataFrame, da_max: float, di_max: float) -> pl.Series:
    return (pl.col("delta_a_au") <= da_max) & (pl.col("delta_i_deg") <= di_max)


def _recall_by_band(
    annotated: pl.DataFrame,
    kept_mask: pl.Series,
    band_col: str,
    edges: list[float],
) -> list[dict]:
    """Recall within bands of *band_col* defined by *edges* (right-open)."""
    rows = []
    vals = annotated[band_col].to_numpy()
    kept = kept_mask.to_numpy()
    for lo, hi in zip(edges[:-1], edges[1:]):
        in_band = (vals >= lo) & (vals < hi)
        n = int(in_band.sum())
        k = int((in_band & kept).sum())
        ci_lo, ci_hi = _wilson_ci(k, n)
        rows.append(
            {
                "band": f"[{lo:g}, {hi:g})",
                "n_encounters": n,
                "n_kept": k,
                "n_missed": n - k,
                "recall": (k / n) if n else float("nan"),
                "ci95_lo": ci_lo,
                "ci95_hi": ci_hi,
            }
        )
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description="Measure orbital prefilter recall")
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--e-min", type=float, default=0.3, help="adverse threshold on e")
    parser.add_argument("--i-min", type=float, default=15.0, help="adverse threshold on i (deg)")
    parser.add_argument(
        "--max-bodies",
        type=int,
        default=None,
        help="cap the adverse subset (default: all). Sampled with the config seed.",
    )
    parser.add_argument(
        "--cross-check",
        type=int,
        default=0,
        help="also run the real pipeline WITH the prefilter on an N≤5000 sub-sample "
        "and assert it equals the analytic mask (0 = skip).",
    )
    parser.add_argument("--out-dir", default="data/output/prefilter_recall")
    args = parser.parse_args()

    cfg = load_config(args.config)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    da_max = cfg.detection.prefilter.semimajor_diff_max_au
    di_max = cfg.detection.prefilter.inclination_diff_max_deg

    # --- Time window + frozen snapshot ---
    tw = cfg.time_window
    t_start = Time(tw.start, scale=tw.scale).tdb.jd
    t_end = Time(tw.end, scale=tw.scale).tdb.jd
    snapshots = discover_snapshots(Path(cfg.paths.raw))
    snap = select_for_window(snapshots, t_start, t_end)
    logger.info("Using MPCORB snapshot %s", snap.path.name)

    sub = cfg.subset
    adverse = _select_adverse(
        snap.path,
        a_min=sub.semimajor_axis_au.min,
        a_max=sub.semimajor_axis_au.max,
        e_min=args.e_min,
        i_min_deg=args.i_min,
        max_bodies=args.max_bodies,
        seed=cfg.run.seed,
    )
    n_bodies = len(adverse)

    # --- Time grid (tiered: coarse bulk + Kepler fine refinement) ---
    fine_step_hours = cfg.propagation.time_step_hours
    coarse_step_hours = cfg.propagation.coarse_step_hours or fine_step_hours
    grid = make_time_grid(t_start, t_end, step_hours=coarse_step_hours)

    det = cfg.detection
    v_max_au_per_day = det.max_relative_velocity_km_s * 86_400.0 / _AU_KM
    query_radius_au = det.threshold_au + v_max_au_per_day * (coarse_step_hours / 24.0)

    n_workers = cfg.parallel.n_workers if cfg.parallel.enabled else 1

    def _run(prefilter_enabled: bool, elements: pl.DataFrame) -> pl.DataFrame:
        return detect_encounters(
            elements,
            grid,
            threshold_au=det.threshold_au,
            semimajor_diff_max_au=da_max,
            inclination_diff_max_deg=di_max,
            leaf_size=det.kdtree.leaf_size,
            fine_step_seconds=det.refinement.fine_time_step_seconds,
            window_hours=det.refinement.window_hours,
            prefilter_enabled=prefilter_enabled,
            refinement_enabled=det.refinement.enabled,
            n_workers=n_workers,
            chunk_size_days=cfg.parallel.chunk_size_days,
            positions=None,  # Kepler propagation
            query_radius_au=query_radius_au,
        )

    # --- Ground truth: detection WITHOUT the prefilter over the full adverse set ---
    logger.info(
        "Running detection WITHOUT prefilter on %d adverse bodies "
        "(%d coarse steps, threshold=%.4f AU)…",
        n_bodies,
        len(grid),
        det.threshold_au,
    )
    t0 = time.monotonic()
    full = _run(prefilter_enabled=False, elements=adverse)
    logger.info("No-prefilter detection: %d encounters in %.1fs", len(full), time.monotonic() - t0)

    annotated = _annotate_orbital_deltas(full, adverse)
    kept_mask = annotated.select(_passes_prefilter(annotated, da_max, di_max)).to_series()

    n_total = len(annotated)
    n_kept = int(kept_mask.sum())
    n_missed = n_total - n_kept
    recall = n_kept / n_total if n_total else float("nan")
    ci_lo, ci_hi = _wilson_ci(n_kept, n_total)

    logger.info(
        "RECALL (analytic mask): %d/%d kept = %.4f%%  [95%% CI %.4f%%, %.4f%%]; %d missed",
        n_kept,
        n_total,
        100 * recall,
        100 * ci_lo,
        100 * ci_hi,
        n_missed,
    )

    missed = annotated.filter(~kept_mask)

    # --- Recall by band ---
    bands = {
        "delta_a_au": _recall_by_band(
            annotated, kept_mask, "delta_a_au", [0.0, 0.1, 0.25, 0.5, 0.75, 1.0, 2.0, 10.0]
        ),
        "delta_i_deg": _recall_by_band(
            annotated, kept_mask, "delta_i_deg", [0.0, 5.0, 15.0, 30.0, 45.0, 90.0, 180.0]
        ),
        "e_max": _recall_by_band(
            annotated, kept_mask, "e_max", [0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0]
        ),
        "i_max_deg": _recall_by_band(
            annotated, kept_mask, "i_max_deg", [0.0, 5.0, 15.0, 30.0, 45.0, 90.0, 180.0]
        ),
    }

    # --- Optional cross-check: real pipeline with prefilter == analytic mask ---
    cross_check: dict | None = None
    if args.cross_check and args.cross_check > 0:
        cc_n = min(args.cross_check, 5000)
        cc_elements = adverse.sample(n=min(cc_n, n_bodies), seed=cfg.run.seed).sort("number")
        logger.info("Cross-check: running BOTH paths on %d bodies…", len(cc_elements))
        cc_full = _run(prefilter_enabled=False, elements=cc_elements)
        cc_pref = _run(prefilter_enabled=True, elements=cc_elements)

        cc_full_keys = set(_pair_key_frame(cc_full).select("_lo", "_hi").iter_rows())
        cc_pref_keys = set(_pair_key_frame(cc_pref).select("_lo", "_hi").iter_rows())
        cc_annot = _annotate_orbital_deltas(cc_full, cc_elements)
        cc_mask = cc_annot.select(_passes_prefilter(cc_annot, da_max, di_max)).to_series()
        cc_analytic_keys = set(
            _pair_key_frame(cc_annot.filter(cc_mask)).select("_lo", "_hi").iter_rows()
        )

        only_pipeline = cc_pref_keys - cc_analytic_keys
        only_analytic = cc_analytic_keys - cc_pref_keys
        subset_ok = cc_pref_keys.issubset(cc_full_keys)
        equal = (len(only_pipeline) == 0) and (len(only_analytic) == 0)
        cross_check = {
            "n_bodies": len(cc_elements),
            "n_no_prefilter": len(cc_full),
            "n_with_prefilter": len(cc_pref),
            "with_prefilter_subset_of_no_prefilter": bool(subset_ok),
            "pipeline_equals_analytic_mask": bool(equal),
            "n_only_in_pipeline": len(only_pipeline),
            "n_only_in_analytic": len(only_analytic),
        }
        logger.info(
            "Cross-check: pipeline-with-prefilter == analytic mask? %s "
            "(only_pipeline=%d, only_analytic=%d, subset_ok=%s)",
            equal,
            len(only_pipeline),
            len(only_analytic),
            subset_ok,
        )

    # --- Persist ---
    full_path = out_dir / "adverse_no_prefilter_encounters.parquet"
    annotated.write_parquet(full_path)
    missed_path = out_dir / "missed_pairs.parquet"
    missed.write_parquet(missed_path)

    summary = {
        "snapshot": snap.path.name,
        "scope": {
            "a_min_au": sub.semimajor_axis_au.min,
            "a_max_au": sub.semimajor_axis_au.max,
            "e_min": args.e_min,
            "i_min_deg": args.i_min,
            "only_numbered": True,
            "n_adverse_bodies": n_bodies,
            "max_bodies_cap": args.max_bodies,
        },
        "config": {
            "threshold_au": det.threshold_au,
            "prefilter_semimajor_diff_max_au": da_max,
            "prefilter_inclination_diff_max_deg": di_max,
            "coarse_step_hours": coarse_step_hours,
            "fine_step_seconds": det.refinement.fine_time_step_seconds,
            "window_hours": det.refinement.window_hours,
            "query_radius_au": query_radius_au,
            "propagation_method": "kepler",
        },
        "result": {
            "n_encounters_total": n_total,
            "n_kept_by_prefilter": n_kept,
            "n_missed_by_prefilter": n_missed,
            "recall": recall,
            "recall_ci95_lo": ci_lo,
            "recall_ci95_hi": ci_hi,
        },
        "recall_by_band": bands,
        "cross_check": cross_check,
        "outputs": {
            "encounters": str(full_path),
            "missed_pairs": str(missed_path),
        },
    }
    summary_path = out_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))
    logger.info("Wrote %s, %s, %s", full_path, missed_path, summary_path)

    # --- Console summary ---
    print("\n" + "=" * 70)
    print(
        f"PREFILTER RECALL — adverse subset ({n_bodies} bodies, e>{args.e_min} OR i>{args.i_min}°)"
    )
    print("=" * 70)
    print(f"Encounters <{det.threshold_au} AU (no prefilter): {n_total}")
    print(f"Kept by prefilter:   {n_kept}")
    print(f"MISSED by prefilter: {n_missed}")
    print(f"RECALL = {100 * recall:.4f}%  [95% CI {100 * ci_lo:.4f}%, {100 * ci_hi:.4f}%]")
    if n_missed:
        print(f"\nMissed pairs — why (Δa, Δi vs limits {da_max:.2f} AU / {di_max:.0f}°):")
        diag = missed.select(
            (pl.col("delta_a_au") > da_max).sum().alias("over_da"),
            (pl.col("delta_i_deg") > di_max).sum().alias("over_di"),
            ((pl.col("delta_a_au") > da_max) & (pl.col("delta_i_deg") > di_max))
            .sum()
            .alias("both"),
        )
        print(diag)
        print("\nMissed pairs by e_max / i_max:")
        print(
            missed.select(
                pl.col("e_max").min().alias("e_min"),
                pl.col("e_max").median().alias("e_med"),
                pl.col("e_max").max().alias("e_max"),
                pl.col("i_max_deg").min().alias("i_min"),
                pl.col("i_max_deg").median().alias("i_med"),
                pl.col("i_max_deg").max().alias("i_max"),
            )
        )
    print("=" * 70)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
