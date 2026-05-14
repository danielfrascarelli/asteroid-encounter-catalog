"""Full encounter-detection pipeline using Gaia DR3 orbital elements.

Loads gaiadr3.sso_orbits, applies subset filters from config, supplements with
JPL Horizons elements for major bodies absent from Gaia DR3 (too bright to be
observed by Gaia: Ceres, Vesta, Pallas, Hygiea), builds the temporal grid over
the Gaia DR3 observation window, runs the detection (prefilter → KD-tree scan →
refinement), and writes the encounter catalog.

For N > 5 000 asteroids the orbital pair prefilter is skipped automatically;
the cKDTree spatial query at the configured threshold_au provides equivalent
filtering without O(N²) memory cost.

Usage:
    docker compose run --rm pipeline python -m scripts.run_pipeline
    docker compose run --rm pipeline python -m scripts.run_pipeline --config config.yaml
"""

from __future__ import annotations

import argparse
import logging
import sys
import time
from pathlib import Path

import polars as pl
from astropy.time import Time

from src.detect.pipeline import detect_encounters
from src.ingest.gaia_orbits import load_gaia_orbits
from src.propagate.grid import make_time_grid
from src.utils.config import load_config

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)

# Bodies that MUST appear in the catalog (too bright for Gaia → not in sso_orbits).
# Queried from JPL Horizons at an epoch within the Gaia DR3 window.
_REQUIRED_BODIES = [1, 2, 4, 10]  # Ceres, Pallas, Vesta, Hygiea
_BODY_NAMES = {1: "Ceres", 2: "Pallas", 4: "Vesta", 10: "Hygiea"}
_SUPPLEMENT_EPOCH_JD = 2457200.5  # 2015-07-01 TDB — mid-Gaia window


def _fetch_horizons_elements(numbers: list[int], epoch_jd: float) -> pl.DataFrame | None:
    """Query JPL Horizons for heliocentric osculating elements at *epoch_jd*.

    Returns a DataFrame with the same schema as gaia_orbits, or None on failure.
    """
    try:
        from astroquery.jplhorizons import Horizons
    except ImportError:
        logger.warning("astroquery not available — skipping Horizons supplement")
        return None

    rows = []
    for num in numbers:
        try:
            name = _BODY_NAMES.get(num, str(num))
            h = Horizons(id=name, id_type="smallbody", epochs=epoch_jd)
            el = h.elements()
            rows.append(
                {
                    "number": int(num),
                    "designation": str(el["targetname"][0]).split("(")[0].strip(),
                    "a_au": float(el["a"][0]),
                    "e": float(el["e"][0]),
                    "i_deg": float(el["incl"][0]),
                    "Omega_deg": float(el["Omega"][0]),
                    "omega_deg": float(el["w"][0]),
                    "M_deg": float(el["M"][0]),
                    "epoch_jd": float(el["datetime_jd"][0]),
                }
            )
            logger.info(
                "Horizons supplement: (%d) %s  a=%.4f AU  epoch=%s",
                num,
                rows[-1]["designation"],
                rows[-1]["a_au"],
                Time(rows[-1]["epoch_jd"], format="jd", scale="tdb").utc.iso[:10],
            )
        except Exception as exc:
            logger.warning("Failed to fetch Horizons elements for (%d): %s", num, exc)

    if not rows:
        return None
    return pl.DataFrame(rows).cast(
        {
            "number": pl.Int32,
            "a_au": pl.Float64,
            "e": pl.Float64,
            "i_deg": pl.Float64,
            "Omega_deg": pl.Float64,
            "omega_deg": pl.Float64,
            "M_deg": pl.Float64,
            "epoch_jd": pl.Float64,
        }
    )


def _supplement_major_bodies(elements: pl.DataFrame) -> pl.DataFrame:
    """Add required bodies missing from *elements* via JPL Horizons."""
    present = set(elements["number"].to_list())
    missing = [n for n in _REQUIRED_BODIES if n not in present]
    if not missing:
        return elements

    logger.info(
        "Bodies absent from gaia_orbits (not observed by Gaia — too bright): %s. "
        "Fetching elements from JPL Horizons at epoch JD %.1f…",
        missing,
        _SUPPLEMENT_EPOCH_JD,
    )
    supplement = _fetch_horizons_elements(missing, _SUPPLEMENT_EPOCH_JD)
    if supplement is None or len(supplement) == 0:
        logger.warning("Could not supplement major bodies — they will be absent from the catalog.")
        return elements

    # Keep only columns present in elements (drop any extras from Horizons)
    cols = elements.columns
    supplement = supplement.select([c for c in cols if c in supplement.columns])
    return pl.concat([supplement, elements])


def _apply_subset(df: pl.DataFrame, cfg) -> pl.DataFrame:
    sub = cfg.subset
    df = df.filter(
        (pl.col("a_au") >= sub.semimajor_axis_au.min)
        & (pl.col("a_au") <= sub.semimajor_axis_au.max)
    )
    if sub.max_asteroids is not None:
        df = df.head(sub.max_asteroids)
    return df


def _verify_major_bodies(results: pl.DataFrame) -> None:
    """Log which required bodies appear in the catalog."""
    for n in _REQUIRED_BODIES:
        hits = results.filter(
            (pl.col("number_1") == n) | (pl.col("number_2") == n)
        )
        if len(hits) == 0:
            logger.warning("Gate check FAILED: (%d) has no encounters in catalog.", n)
        else:
            closest = hits["dist_au"].min()
            logger.info("Gate check OK: (%d) — %d encounters, closest %.6f AU", n, len(hits), closest)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run Gaia asteroid encounter detection")
    parser.add_argument("--config", default="config.yaml")
    args = parser.parse_args()

    cfg = load_config(args.config)

    # --- Load ---
    orbits_path = Path(cfg.paths.raw) / "gaia_orbits.parquet"
    logger.info("Loading orbital elements from %s", orbits_path)
    elements = load_gaia_orbits(orbits_path)
    logger.info("Loaded %d asteroids (Gaia DR3 sso_orbits)", len(elements))

    # --- Supplement major bodies absent from Gaia DR3 ---
    elements = _supplement_major_bodies(elements)

    # --- Subset ---
    elements = _apply_subset(elements, cfg)
    sub = cfg.subset
    logger.info(
        "After subset filter: %d asteroids  (a=[%.2f, %.2f] AU%s)",
        len(elements),
        sub.semimajor_axis_au.min,
        sub.semimajor_axis_au.max,
        f", max_asteroids={sub.max_asteroids}" if sub.max_asteroids else "",
    )

    if len(elements) < 2:
        logger.error("Fewer than 2 asteroids after filtering — nothing to detect.")
        return 1

    # --- Time grid ---
    tw = cfg.time_window
    t_start = Time(tw.start, scale=tw.scale).tdb.jd
    t_end = Time(tw.end, scale=tw.scale).tdb.jd
    grid = make_time_grid(t_start, t_end, step_hours=cfg.propagation.time_step_hours)
    logger.info(
        "Time grid: %s → %s  (%d steps, Δt=%.1fh)",
        tw.start,
        tw.end,
        len(grid),
        cfg.propagation.time_step_hours,
    )

    # --- Detection ---
    det = cfg.detection
    par = cfg.parallel
    n_workers = par.n_workers if par.enabled else 1

    logger.info(
        "Starting encounter detection  (workers=%s, chunk_size=%.0f days)…",
        n_workers,
        par.chunk_size_days,
    )
    t0 = time.monotonic()

    results = detect_encounters(
        elements,
        grid,
        threshold_au=det.threshold_au,
        semimajor_diff_max_au=det.prefilter.semimajor_diff_max_au,
        inclination_diff_max_deg=det.prefilter.inclination_diff_max_deg,
        leaf_size=det.kdtree.leaf_size,
        fine_step_seconds=det.refinement.fine_time_step_seconds,
        window_hours=det.refinement.window_hours,
        prefilter_enabled=det.prefilter.enabled,
        refinement_enabled=det.refinement.enabled,
        n_workers=n_workers,
        chunk_size_days=par.chunk_size_days,
    )

    elapsed = time.monotonic() - t0
    logger.info(
        "Detection complete in %.1fs — %d encounters ≤ %.4f AU",
        elapsed,
        len(results),
        det.threshold_au,
    )

    # --- Gate check: major bodies ---
    _verify_major_bodies(results)

    # --- Save ---
    out_dir = Path(cfg.paths.output)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{cfg.output.filename}.{cfg.output.format}"

    if cfg.output.format == "parquet":
        results.write_parquet(out_path, compression=cfg.output.compression)
    else:
        results.write_csv(out_path)

    logger.info("Catalog saved → %s  (%d rows)", out_path, len(results))

    # --- Top encounters ---
    if len(results) > 0:
        logger.info("Top 10 closest encounters:")
        for row in results.head(10).iter_rows(named=True):
            t = Time(row["jd_tdb"], format="jd", scale="tdb")
            logger.info(
                "  (%d) %-20s — (%d) %-20s  %.6f AU  %s",
                row["number_1"],
                row["designation_1"],
                row["number_2"],
                row["designation_2"],
                row["dist_au"],
                t.utc.iso[:10],
            )

    return 0


if __name__ == "__main__":
    sys.exit(main())
