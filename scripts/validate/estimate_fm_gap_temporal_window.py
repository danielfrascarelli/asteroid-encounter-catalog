"""Quantify how much of the FM-vs-catalog gap is explained by the DR3 window.

Context (Tarea 11 / tribunal M1)
---------------------------------
``crosscheck_literature_encounters.py`` classifies 25,962 Fuentes-Muñoz et al.
(2025) (perturber, target) pairs as ``absent_no_close_encounter``: both bodies
appear in our catalog's numbered universe, but the pair never approached
< ~0.05 AU *within* the DR3 window (2014-07-25 -> 2017-05-28). FM fits orbits
over the full archival baseline (decades, up to Nov 2024) using the SAME
0.05 AU threshold for MBAs, so the leading hypothesis (tribunal M1) is that
these pairs' close approach simply falls outside our 2.8-year window, not
that the threshold differs.

Does FM publish a per-pair encounter epoch?
--------------------------------------------
No. Table 5's byte-by-byte header (verified directly against
``data/raw/fuentes_munoz_2025/ajae0cc9t5_mrt.txt``) has exactly these fields:
Asteroid, Group, GMini, e_GMini, r_GMini, GMfin, e_GMfin, Ntest, Nsign, List,
Bibcode, Shortbib, SPKID. There is no epoch/date/MJD column anywhere, and the
pipe-delimited ``List`` is just target identifiers (truncated to top-100 by
signal), not epochs. So a per-pair FM close-approach date cannot be recovered
by parsing -- it requires propagating each pair ourselves.

Method (tractable estimate on a sample, not the full 25,962)
--------------------------------------------------------------
1. Load the ``absent_no_close_encounter`` pairs from
   ``data/output/literature_validation/completeness_pairs.csv``, dedupe to
   unordered (lo, hi) keys.
2. Draw a random sample (seed from config.yaml, default 42) of ``--n-sample``
   unique pairs.
3. For each pair, propagate BOTH bodies with the project's two-body Kepler
   propagator (``src/propagate/kepler.py``), using the same MPCORB elements
   snapshot that fed the frozen catalog
   (``data/cache/nbody_validation/mpcorb_stageb_elements.parquet``, single
   epoch JD ~2457400 ~= 2016-02-16), over a daily grid spanning
   ``--start``..``--end`` (default 1990-01-01 .. 2024-12-31, a 35-year
   archival-scale baseline), then refine the minimum with an hourly sub-grid
   around the coarse minimum.
4. Classify each pair by its GLOBAL minimum separation over the full baseline:
   - ``outside_window``: min < threshold AND its epoch falls outside the DR3
     window -> confirms the temporal-window explanation for this pair.
   - ``inside_window``: min < threshold AND its epoch falls INSIDE the DR3
     window -> inconsistent with the catalog's absence call; flagged as an
     anomaly (two-body-vs-n-body model mismatch, or a borderline case).
   - ``no_encounter_in_range``: no epoch in [start, end] drops below
     threshold -> the temporal-window explanation does NOT account for this
     pair (could be a genuinely wide encounter that never gets close in this
     baseline, an encounter beyond the propagated range, or Kepler drift
     error masking/creating a false negative).
5. Report the fraction ``outside_window / n_sample`` with a Wilson 95% CI as
   the quantified "gap explained by temporal window" estimate.

Caveat -- read before citing this number
-----------------------------------------
Two-body Kepler with elements fixed at a single 2016 epoch is NOT an
n-body integration: it ignores secular/resonant planetary perturbations,
whose effect on the true orbit grows with |t - epoch|. Over the ~10-year
excursions used here this is a coarse (but standard, and explicitly endorsed
by the task as tractable) approximation good enough to bucket a close
approach into "which multi-year window did the deepest approach happen in",
not to pin down its exact epoch or distance to sub-0.001 AU precision. Some
individual classifications (particularly ones near the 0.05 AU threshold or
near a window boundary) may be wrong; the sample-level fraction, not any
single pair, is the deliverable.

Usage
-----
    docker compose run --rm pipeline python -m \\
        scripts.validate.estimate_fm_gap_temporal_window --n-sample 500
"""

from __future__ import annotations

import argparse
import json
import logging
import random
from pathlib import Path

import numpy as np
import polars as pl
from astropy.time import Time

from src.propagate.kepler import kepler_to_cartesian

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
logger = logging.getLogger(__name__)

_PAIRS_CSV = Path("data/output/literature_validation/completeness_pairs.csv")
_ELEMENTS_PATH = Path("data/cache/nbody_validation/mpcorb_stageb_elements.parquet")
_OUT_DIR = Path("data/output/literature_validation")

_DR3_WINDOW_START = "2014-07-25T00:00:00"
_DR3_WINDOW_END = "2017-05-28T00:00:00"
_DEFAULT_THRESHOLD = 0.05
_DEFAULT_SEED = 42


def _jd_tdb(iso: str) -> float:
    return float(Time(iso, format="isot", scale="utc").tdb.jd)


def load_absent_pairs(path: Path) -> list[tuple[int, int]]:
    """Return unique (lo, hi) keys classified `absent_no_close_encounter`."""
    df = pl.read_csv(
        path,
        schema_overrides={"target_token": pl.String, "min_dist_au": pl.Float64},
        infer_schema_length=10000,
    )
    sub = df.filter(pl.col("outcome") == "absent_no_close_encounter")
    keys = sub.select(
        pl.min_horizontal("perturber", "target_number").alias("lo"),
        pl.max_horizontal("perturber", "target_number").alias("hi"),
    ).unique()
    return list(zip(keys["lo"].to_list(), keys["hi"].to_list(), strict=True))


def load_elements(path: Path) -> dict[int, dict[str, float]]:
    """Return {number: {a_au, e, i_deg, Omega_deg, omega_deg, M_deg, epoch_jd}}."""
    df = pl.read_parquet(path).rename(
        {"node_deg": "Omega_deg", "argperi_deg": "omega_deg", "mean_anomaly_deg": "M_deg"}
    )
    out: dict[int, dict[str, float]] = {}
    for row in df.select(
        "number", "a_au", "e", "i_deg", "Omega_deg", "omega_deg", "M_deg", "epoch_jd"
    ).iter_rows(named=True):
        out[int(row["number"])] = row
    return out


def closest_approach(
    el_p: dict[str, float],
    el_t: dict[str, float],
    t_grid_jd: np.ndarray,
) -> tuple[float, float]:
    """Return (min_dist_au, epoch_jd_of_min) of body P vs body T over t_grid_jd."""
    pos_p = kepler_to_cartesian(
        el_p["a_au"],
        el_p["e"],
        np.deg2rad(el_p["i_deg"]),
        np.deg2rad(el_p["Omega_deg"]),
        np.deg2rad(el_p["omega_deg"]),
        np.deg2rad(el_p["M_deg"]),
        el_p["epoch_jd"],
        t_grid_jd,
    )
    pos_t = kepler_to_cartesian(
        el_t["a_au"],
        el_t["e"],
        np.deg2rad(el_t["i_deg"]),
        np.deg2rad(el_t["Omega_deg"]),
        np.deg2rad(el_t["omega_deg"]),
        np.deg2rad(el_t["M_deg"]),
        el_t["epoch_jd"],
        t_grid_jd,
    )
    dist = np.linalg.norm(pos_p - pos_t, axis=-1)
    i = int(np.argmin(dist))
    return float(dist[i]), float(t_grid_jd[i])


def refine(
    el_p: dict[str, float],
    el_t: dict[str, float],
    t_center_jd: float,
    coarse_step_days: float,
) -> tuple[float, float]:
    """Refine around a coarse minimum with an hourly sub-grid spanning +/- one coarse step."""
    span = max(coarse_step_days, 2.0)
    fine = np.arange(t_center_jd - span, t_center_jd + span, 1.0 / 24.0)
    return closest_approach(el_p, el_t, fine)


def wilson_ci(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    """Wilson score interval for a binomial proportion."""
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    denom = 1 + z**2 / n
    centre = p + z**2 / (2 * n)
    half = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2))
    return ((centre - half) / denom, (centre + half) / denom)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pairs-csv", default=str(_PAIRS_CSV))
    ap.add_argument("--elements", default=str(_ELEMENTS_PATH))
    ap.add_argument("--out-dir", default=str(_OUT_DIR))
    ap.add_argument("--threshold", type=float, default=_DEFAULT_THRESHOLD)
    ap.add_argument("--n-sample", type=int, default=500)
    ap.add_argument("--seed", type=int, default=_DEFAULT_SEED)
    ap.add_argument("--start", default="1990-01-01T00:00:00")
    ap.add_argument("--end", default="2024-12-31T00:00:00")
    ap.add_argument("--coarse-step-days", type=float, default=1.0)
    args = ap.parse_args()

    pairs = load_absent_pairs(Path(args.pairs_csv))
    logger.info("Loaded %d unique absent_no_close_encounter pairs.", len(pairs))

    elements = load_elements(Path(args.elements))
    logger.info("Loaded elements for %d numbered bodies.", len(elements))

    rng = random.Random(args.seed)
    pairs_shuffled = pairs[:]
    rng.shuffle(pairs_shuffled)

    t_start = _jd_tdb(args.start)
    t_end = _jd_tdb(args.end)
    win_start = _jd_tdb(_DR3_WINDOW_START)
    win_end = _jd_tdb(_DR3_WINDOW_END)
    t_grid = np.arange(t_start, t_end, args.coarse_step_days)
    logger.info(
        "Propagation baseline: JD [%.1f, %.1f] (%s .. %s), %d coarse points; "
        "DR3 window JD [%.1f, %.1f].",
        t_start,
        t_end,
        args.start,
        args.end,
        len(t_grid),
        win_start,
        win_end,
    )

    rows: list[dict] = []
    n_skipped_missing_elements = 0
    for lo, hi in pairs_shuffled:
        if len(rows) >= args.n_sample:
            break
        el_lo = elements.get(int(lo))
        el_hi = elements.get(int(hi))
        if el_lo is None or el_hi is None:
            n_skipped_missing_elements += 1
            continue

        min_dist, t_min = closest_approach(el_lo, el_hi, t_grid)
        min_dist, t_min = refine(el_lo, el_hi, t_min, args.coarse_step_days)

        has_encounter = min_dist < args.threshold
        in_window = win_start <= t_min <= win_end
        if has_encounter and not in_window:
            outcome = "outside_window"
        elif has_encounter and in_window:
            outcome = "inside_window"
        else:
            outcome = "no_encounter_in_range"

        iso_min = Time(t_min, format="jd", scale="tdb").isot
        rows.append(
            {
                "lo": int(lo),
                "hi": int(hi),
                "min_dist_au": min_dist,
                "t_min_jd_tdb": t_min,
                "t_min_iso": iso_min,
                "in_dr3_window": in_window,
                "outcome": outcome,
            }
        )

    result = pl.DataFrame(rows)
    n = result.height
    counts = {
        k: int((result["outcome"] == k).sum())
        for k in ("outside_window", "inside_window", "no_encounter_in_range")
    }
    k_outside = counts["outside_window"]
    frac = k_outside / n if n else 0.0
    lo_ci, hi_ci = wilson_ci(k_outside, n)

    summary = {
        "method": "two-body Kepler propagation (single 2016-02-16 MPCORB elements "
        "epoch), daily coarse grid + hourly local refinement, over "
        f"[{args.start}, {args.end}]",
        "threshold_au": args.threshold,
        "dr3_window": [_DR3_WINDOW_START, _DR3_WINDOW_END],
        "propagation_baseline": [args.start, args.end],
        "n_absent_no_close_encounter_total": len(pairs),
        "n_sample_requested": args.n_sample,
        "n_sample_evaluated": n,
        "n_skipped_missing_elements": n_skipped_missing_elements,
        "seed": args.seed,
        "counts": counts,
        "fraction_outside_window": frac,
        "fraction_outside_window_wilson95_ci": [lo_ci, hi_ci],
        "fraction_inside_window_anomaly": counts["inside_window"] / n if n else 0.0,
        "fraction_no_encounter_in_range": counts["no_encounter_in_range"] / n if n else 0.0,
    }

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    result.write_csv(out_dir / "fm_gap_temporal_window_sample.csv")
    (out_dir / "fm_gap_temporal_window_summary.json").write_text(json.dumps(summary, indent=2))

    print("\n" + "=" * 78)
    print("FM gap x temporal window -- sample-based estimate (Tarea 11 / tribunal M1)")
    print("=" * 78)
    print(f"Sample: {n} / {len(pairs)} absent_no_close_encounter pairs (seed={args.seed})")
    print(f"Baseline propagated: {args.start} .. {args.end} (two-body Kepler)")
    print(
        f"  outside_window (encounter <{args.threshold} AU, outside DR3 window): "
        f"{counts['outside_window']} ({100*frac:.1f}%, 95% CI "
        f"[{100*lo_ci:.1f}, {100*hi_ci:.1f}]%)"
    )
    print(
        f"  inside_window  (encounter <{args.threshold} AU, INSIDE DR3 window, "
        f"anomaly): {counts['inside_window']} "
        f"({100*counts['inside_window']/n:.1f}%)"
    )
    print(
        f"  no_encounter_in_range (never <{args.threshold} AU in baseline): "
        f"{counts['no_encounter_in_range']} "
        f"({100*counts['no_encounter_in_range']/n:.1f}%)"
    )
    print("=" * 78)
    logger.info("Wrote outputs to %s", out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
