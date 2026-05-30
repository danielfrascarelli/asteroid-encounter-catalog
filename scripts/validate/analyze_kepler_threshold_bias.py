"""Characterise the Kepler-vs-N-body threshold crossings — Track B Stage 3.

``FROZEN_RUN.md`` reported that 25,283 pairs cross the 0.05 AU detection
threshold when refined from Kepler two-body to N-body — *all in the same
direction* (Kepler ``< 0.05``, N-body ``>= 0.05``), with zero crossings the
other way, and read this as "a small systematic over-detection bias near the
threshold, no false negatives".

This script measures that claim on ``encounters_catalog_hybrid_stageb.parquet``
and adds the crucial methodological caveat the freeze note glossed over: the
frozen catalogue **only contains pairs Kepler already placed below 0.05 AU**
(that is the detection threshold), and the N-body subset was selected by an
*orbital* criterion (``q_min < 1.8 AU OR e_max > 0.3``), not by distance. So
the "0 reverse crossings" is largely **censoring**, not a measurement: pairs
Kepler put *above* 0.05 AU were never written and therefore cannot be observed
crossing downward. The genuinely measurable quantities are:

1. Among refined pairs, how many cross upward and what is the Δdist
   distribution near the threshold (the real over-detection rate).
2. How the crossing rate depends on proximity to the threshold, relative
   velocity, and orbital band (q, e) — i.e. where Kepler is least reliable.
3. An honest statement of what *cannot* be concluded (the false-negative rate)
   and the experiment that would measure it.

Everything is computed with DuckDB so the ~1 GB / multi-million-row catalogue
never has to be materialised in memory.

Usage
-----
    docker compose run --rm pipeline python -m scripts.validate.analyze_kepler_threshold_bias \\
        --out-prefix data/output/kepler_bias/threshold
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import duckdb

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
logger = logging.getLogger(__name__)

_THRESHOLD_AU = 0.05


def _sql_path(path: Path) -> str:
    """Return a DuckDB-safe single-quoted path literal."""
    return "'" + str(path).replace("'", "''") + "'"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument(
        "--catalog",
        type=Path,
        default=Path("data/output/encounters_catalog_hybrid_stageb.parquet"),
    )
    parser.add_argument(
        "--elements",
        type=Path,
        default=Path("data/cache/nbody_validation/mpcorb_stageb_elements.parquet"),
        help="Cached MPCORB elements (number,a_au,e) for the orbital-band join.",
    )
    parser.add_argument(
        "--out-prefix", type=Path, default=Path("data/output/kepler_bias/threshold")
    )
    args = parser.parse_args()

    if not args.catalog.exists():
        logger.error("Catalog not found: %s", args.catalog)
        return 1
    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    cat = _sql_path(args.catalog)
    con = duckdb.connect()

    # Restrict to N-body-refined rows: only there do we have both distances.
    refined = f"""
        SELECT number_1, number_2, jd_tdb, rel_vel_au_day,
               dist_au_kepler, dist_au_nbody,
               (dist_au_nbody - dist_au_kepler) AS delta_dist_au
        FROM read_parquet({cat})
        WHERE refinement_method = 'nbody'
          AND dist_au_nbody IS NOT NULL
          AND dist_au_kepler IS NOT NULL
    """
    con.execute(f"CREATE TEMP VIEW refined AS {refined}")

    # ---- headline counts -----------------------------------------------------
    counts = con.execute(f"""
        SELECT
            count(*) AS n_refined,
            sum(CASE WHEN dist_au_kepler < {_THRESHOLD_AU} THEN 1 ELSE 0 END) AS n_kepler_below,
            sum(CASE WHEN dist_au_kepler >= {_THRESHOLD_AU} THEN 1 ELSE 0 END) AS n_kepler_above,
            sum(CASE WHEN dist_au_kepler < {_THRESHOLD_AU}
                      AND dist_au_nbody >= {_THRESHOLD_AU} THEN 1 ELSE 0 END) AS n_cross_up,
            sum(CASE WHEN dist_au_kepler >= {_THRESHOLD_AU}
                      AND dist_au_nbody < {_THRESHOLD_AU} THEN 1 ELSE 0 END) AS n_cross_down,
            avg(delta_dist_au) AS mean_delta,
            median(delta_dist_au) AS median_delta,
            stddev_samp(delta_dist_au) AS std_delta,
            avg(abs(delta_dist_au)) AS mean_abs_delta
        FROM refined
        """).fetchone()
    keys = [
        "n_refined",
        "n_kepler_below",
        "n_kepler_above",
        "n_cross_up",
        "n_cross_down",
        "mean_delta",
        "median_delta",
        "std_delta",
        "mean_abs_delta",
    ]
    head = {k: (float(v) if v is not None else None) for k, v in zip(keys, counts)}
    n_below = head["n_kepler_below"] or 0
    head["cross_up_rate_of_kepler_below"] = head["n_cross_up"] / n_below if n_below else None

    # ---- Δdist by proximity-to-threshold band (Kepler distance) --------------
    band_rows = con.execute(f"""
        SELECT
            floor(dist_au_kepler / 0.005) * 0.005 AS kepler_band_lo,
            count(*) AS n,
            avg(delta_dist_au) AS mean_delta,
            median(delta_dist_au) AS median_delta,
            stddev_samp(delta_dist_au) AS std_delta,
            sum(CASE WHEN dist_au_nbody >= {_THRESHOLD_AU} THEN 1 ELSE 0 END) AS n_cross_up,
            avg(rel_vel_au_day) AS mean_rel_vel
        FROM refined
        WHERE dist_au_kepler < {_THRESHOLD_AU}
        GROUP BY kepler_band_lo
        ORDER BY kepler_band_lo
        """).fetchall()
    bands = [
        {
            "kepler_band_lo_au": float(r[0]),
            "kepler_band_hi_au": float(r[0]) + 0.005,
            "n": int(r[1]),
            "mean_delta_au": float(r[2]),
            "median_delta_au": float(r[3]),
            "std_delta_au": float(r[4]) if r[4] is not None else None,
            "n_cross_up": int(r[5]),
            "cross_up_rate": int(r[5]) / int(r[1]) if r[1] else 0.0,
            "mean_rel_vel_au_day": float(r[6]),
        }
        for r in band_rows
    ]

    # ---- relative-velocity dependence of crossings (Kepler-below pairs) ------
    vel_rows = con.execute(f"""
        WITH b AS (
            SELECT delta_dist_au, rel_vel_au_day,
                   (dist_au_nbody >= {_THRESHOLD_AU}) AS crossed,
                   ntile(5) OVER (ORDER BY rel_vel_au_day) AS vq
            FROM refined
            WHERE dist_au_kepler < {_THRESHOLD_AU}
        )
        SELECT vq, count(*) AS n,
               min(rel_vel_au_day) AS vlo, max(rel_vel_au_day) AS vhi,
               avg(delta_dist_au) AS mean_delta,
               sum(CASE WHEN crossed THEN 1 ELSE 0 END) AS n_cross_up
        FROM b GROUP BY vq ORDER BY vq
        """).fetchall()
    vel_quintiles = [
        {
            "quintile": int(r[0]),
            "n": int(r[1]),
            "rel_vel_lo_au_day": float(r[2]),
            "rel_vel_hi_au_day": float(r[3]),
            "mean_delta_au": float(r[4]),
            "n_cross_up": int(r[5]),
            "cross_up_rate": int(r[5]) / int(r[1]) if r[1] else 0.0,
        }
        for r in vel_rows
    ]

    # ---- orbital-band dependence: join crossing pairs to MPCORB elements -----
    orbital: dict = {"available": False}
    if args.elements.exists():
        el = _sql_path(args.elements)
        orbital_rows = con.execute(f"""
            WITH crossings AS (
                SELECT number_1, number_2
                FROM refined
                WHERE dist_au_kepler < {_THRESHOLD_AU} AND dist_au_nbody >= {_THRESHOLD_AU}
            ),
            el AS (SELECT number, a_au, e FROM read_parquet({el})),
            joined AS (
                SELECT
                    least(e1.a_au * (1 - e1.e), e2.a_au * (1 - e2.e)) AS q_min,
                    greatest(e1.e, e2.e) AS e_max
                FROM crossings c
                JOIN el e1 ON c.number_1 = e1.number
                JOIN el e2 ON c.number_2 = e2.number
            )
            SELECT
                count(*) AS n_joined,
                sum(CASE WHEN q_min < 1.8 THEN 1 ELSE 0 END) AS n_q_min_lt_1p8,
                sum(CASE WHEN e_max > 0.3 THEN 1 ELSE 0 END) AS n_e_max_gt_0p3,
                avg(q_min) AS mean_q_min,
                median(q_min) AS median_q_min,
                avg(e_max) AS mean_e_max,
                median(e_max) AS median_e_max
            FROM joined
            """).fetchone()
        ok = [
            "n_joined",
            "n_q_min_lt_1p8",
            "n_e_max_gt_0p3",
            "mean_q_min",
            "median_q_min",
            "mean_e_max",
            "median_e_max",
        ]
        orbital = {"available": True}
        orbital.update({k: (float(v) if v is not None else None) for k, v in zip(ok, orbital_rows)})

    result = {
        "catalog": str(args.catalog),
        "threshold_au": _THRESHOLD_AU,
        "headline": head,
        "kepler_distance_bands": bands,
        "rel_vel_quintiles": vel_quintiles,
        "orbital_band_of_crossings": orbital,
        "caveat": (
            "The catalog contains only Kepler<0.05 AU pairs and the N-body subset "
            "was selected by q_min<1.8 OR e_max>0.3, not by distance. The zero "
            "downward crossings is therefore censoring, not a measured false-negative "
            "rate; measuring the latter requires N-body-refining a sample of pairs "
            "with Kepler distance in [0.05, ~0.06] AU, which this catalog excludes."
        ),
    }

    summary_json = args.out_prefix.with_name(args.out_prefix.name + "_summary.json")
    summary_json.write_text(json.dumps(result, indent=2))

    import csv as _csv

    bands_csv = args.out_prefix.with_name(args.out_prefix.name + "_bands.csv")
    if bands:
        with bands_csv.open("w", newline="") as fh:
            w = _csv.DictWriter(fh, fieldnames=list(bands[0].keys()))
            w.writeheader()
            w.writerows(bands)

    # ---- console report ------------------------------------------------------
    print("\n=== Kepler vs N-body threshold-crossing analysis ===")
    print(f"refined (N-body) pairs : {head['n_refined']:,.0f}")
    print(f"  Kepler < 0.05 AU     : {head['n_kepler_below']:,.0f}")
    print(f"  Kepler >= 0.05 AU    : {head['n_kepler_above']:,.0f}  (in-catalog by construction)")
    print(f"cross up   (K<0.05 -> N>=0.05): {head['n_cross_up']:,.0f}")
    print(f"cross down (K>=0.05 -> N<0.05): {head['n_cross_down']:,.0f}")
    if head["cross_up_rate_of_kepler_below"] is not None:
        print(f"over-detection rate (of Kepler<0.05): {head['cross_up_rate_of_kepler_below']:.4%}")
    print(
        f"Δdist (N-K): mean={head['mean_delta']:.2e}  median={head['median_delta']:.2e}  "
        f"std={head['std_delta']:.2e} AU"
    )
    print("\nKepler-distance band  |     n      | mean Δ (AU) | cross-up rate")
    print("-" * 64)
    for b in bands:
        print(
            f"  [{b['kepler_band_lo_au']:.3f},{b['kepler_band_hi_au']:.3f})  "
            f"{b['n']:>11,}  {b['mean_delta_au']:>11.2e}  {b['cross_up_rate']:>12.4%}"
        )
    if orbital.get("available"):
        print(
            f"\ncrossings orbital band: {orbital['n_q_min_lt_1p8']:,.0f}/{orbital['n_joined']:,.0f} "
            f"have q_min<1.8 AU; median q_min={orbital['median_q_min']:.2f} AU, "
            f"median e_max={orbital['median_e_max']:.3f}"
        )
    print(f"\nWrote:\n  {summary_json}\n  {bands_csv}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
