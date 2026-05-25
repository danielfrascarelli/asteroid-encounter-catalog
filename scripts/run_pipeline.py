"""Full encounter-detection pipeline using MPCORB orbital elements.

Loads all asteroids from MPCORB.DAT (not limited to those observed by Gaia),
applies subset filters from config, builds the temporal grid over the Gaia DR3
observation window, runs the detection (prefilter → KD-tree scan → refinement),
and writes the encounter catalog.

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

from src.catalog.writer import write_detection_sidecar
from src.detect.pipeline import detect_encounters
from src.ingest.mpcorb import parse_mpcorb
from src.ingest.mpcorb_archive import discover_snapshots, select_for_window
from src.propagate.grid import make_time_grid, propagate_full_grid
from src.utils.config import load_config

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)

# Bodies whose presence in the output catalog is verified as a gate check.
_REQUIRED_BODIES = [1, 2, 4, 10]  # Ceres, Pallas, Vesta, Hygiea


def _apply_subset(df: pl.DataFrame, cfg) -> pl.DataFrame:
    """Truncate elements DataFrame to cfg.subset.max_asteroids if set."""
    if cfg.subset.max_asteroids is not None:
        df = df.head(cfg.subset.max_asteroids)
    return df


def _verify_major_bodies(results: pl.DataFrame) -> None:
    """Log which required bodies appear in the catalog."""
    for n in _REQUIRED_BODIES:
        hits = results.filter((pl.col("number_1") == n) | (pl.col("number_2") == n))
        if len(hits) == 0:
            logger.warning("Gate check FAILED: (%d) has no encounters in catalog.", n)
        else:
            closest = hits["dist_au"].min()
            logger.info(
                "Gate check OK: (%d) — %d encounters, closest %.6f AU", n, len(hits), closest
            )


def main() -> int:
    parser = argparse.ArgumentParser(description="Run Gaia asteroid encounter detection")
    parser.add_argument("--config", default="config.yaml")
    args = parser.parse_args()

    cfg = load_config(args.config)

    # --- Time window (needed to pick the right MPCORB snapshot) ---
    tw = cfg.time_window
    t_start = Time(tw.start, scale=tw.scale).tdb.jd
    t_end = Time(tw.end, scale=tw.scale).tdb.jd

    # --- Select MPCORB snapshot whose epoch is closest to the window centre ---
    snapshots = discover_snapshots(Path(cfg.paths.raw))
    if not snapshots:
        logger.error(
            "No MPCORB snapshots found under %s. Run scripts.download_mpcorb or "
            "scripts.download_mpcorb_historical first.",
            cfg.paths.raw,
        )
        return 1
    snap = select_for_window(snapshots, t_start, t_end)
    centre_jd = 0.5 * (t_start + t_end)
    centre_iso = Time(centre_jd, format="jd", scale="tdb").utc.iso[:10]
    snap_epoch_iso = Time(snap.epoch_jd, format="jd", scale="tdb").utc.iso[:10]
    offset_years = (snap.epoch_jd - centre_jd) / 365.25
    logger.info(
        "Window centre %s → selected snapshot %s (epoch %s, |Δt|=%.2f yr)",
        centre_iso,
        snap.path.name,
        snap_epoch_iso,
        abs(offset_years),
    )
    if len(snapshots) > 1:
        for s in snapshots:
            mark = " ← selected" if s.path == snap.path else ""
            logger.info(
                "  available: %-48s epoch=%s%s",
                s.path.name,
                Time(s.epoch_jd, format="jd", scale="tdb").utc.iso[:10],
                mark,
            )

    # --- Load ---
    sub = cfg.subset
    elements = parse_mpcorb(
        snap.path,
        only_numbered=sub.only_numbered,
        semimajor_min_au=sub.semimajor_axis_au.min,
        semimajor_max_au=sub.semimajor_axis_au.max,
    )
    logger.info("Loaded %d asteroids from MPCORB", len(elements))

    # --- Subset: max_asteroids cap ---
    elements = _apply_subset(elements, cfg)
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
    # Strategy A: when ``coarse_step_hours`` is set and exceeds the fine step,
    # the bulk N-body cache is built (and the KD-tree scan walks) at the coarse
    # step. Refinement still uses the fine ``time_step_hours`` indirectly via
    # the Kepler propagator. When ``coarse_step_hours`` is unset, the bulk grid
    # falls back to ``time_step_hours`` (original behaviour).
    fine_step_hours = cfg.propagation.time_step_hours
    coarse_step_hours = cfg.propagation.coarse_step_hours or fine_step_hours
    use_tiered = coarse_step_hours > fine_step_hours
    grid = make_time_grid(t_start, t_end, step_hours=coarse_step_hours)
    if use_tiered:
        logger.info(
            "Tiered cache: bulk grid Δt=%.1fh (%d steps); refinement falls back "
            "to Kepler on the fine grid Δt=%.1fh",
            coarse_step_hours,
            len(grid),
            fine_step_hours,
        )
    else:
        logger.info(
            "Time grid: %s → %s  (%d steps, Δt=%.1fh)",
            tw.start,
            tw.end,
            len(grid),
            coarse_step_hours,
        )

    par = cfg.parallel

    # --- Propagation (N-body branch precomputes the trajectory) ---
    positions = None
    if cfg.propagation.method.lower() == "rebound":
        logger.info(
            "Propagation method: rebound  (integrator=%s, planets=%s, major_asteroids=%s)",
            cfg.propagation.rebound.integrator,
            cfg.propagation.rebound.include_planets,
            cfg.propagation.rebound.include_major_asteroids,
        )
        rebound_kwargs = {
            "include_planets": cfg.propagation.rebound.include_planets,
            "include_major_asteroids": cfg.propagation.rebound.include_major_asteroids,
            "integrator": cfg.propagation.rebound.integrator,
            "n_workers": par.n_workers if par.enabled else 1,
            # The integrator step stays at the fine resolution to preserve
            # accuracy; only the snapshot cadence (the time_grid spacing) is
            # coarsened. This is the standard trick — full N-body precision
            # with 12× less disk and 12× less work in the KD-tree scan.
            "dt_days": fine_step_hours / 24.0,
        }
        cache_dir = cfg.paths.cache if cfg.propagation.cache_results else None
        cache_key = None
        if cache_dir is not None:
            from src.propagate.cache import build_cache_key

            cache_key = build_cache_key(
                snapshot_sha=snap.path,
                time_grid=grid,
                method="rebound",
                rebound_kwargs=rebound_kwargs,
                n_asteroids=len(elements),
            )
        t_prop = time.monotonic()
        positions = propagate_full_grid(
            elements,
            grid,
            method="rebound",
            rebound_kwargs=rebound_kwargs,
            cache_dir=cache_dir,
            cache_key=cache_key,
            cache_format=cfg.propagation.cache_format,
        )
        logger.info(
            "Propagation done in %.1fs — trajectory shape %s",
            time.monotonic() - t_prop,
            positions.shape if positions is not None else None,
        )

    # --- Detection ---
    det = cfg.detection
    n_workers = par.n_workers if par.enabled else 1

    # Strategy A: widen the KD-tree query radius to cover encounters whose
    # closest approach falls between two coarse samples. Half-step displacement
    # at v_max sets the bound; we use the full step for extra safety margin.
    query_radius_au = det.threshold_au
    if use_tiered:
        v_max_au_per_day = det.max_relative_velocity_km_s * 86_400.0 / 1.495_978_707e8
        widen_au = v_max_au_per_day * (coarse_step_hours / 24.0)
        query_radius_au = det.threshold_au + widen_au
        logger.info(
            "Coarse Δt=%.1fh × v_max=%.1f km/s → widening +%.5f AU " "(query radius %.5f AU)",
            coarse_step_hours,
            det.max_relative_velocity_km_s,
            widen_au,
            query_radius_au,
        )

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
        positions=positions,
        query_radius_au=query_radius_au,
        force_kepler_refine=use_tiered,
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
        results.write_parquet(out_path, compression=cfg.output.compression)  # type: ignore[arg-type]
    else:
        results.write_csv(out_path)

    logger.info("Catalog saved → %s  (%d rows)", out_path, len(results))

    # --- Provenance sidecar ---
    # Records scan/refine method, time-grid params, prefilter, MPCORB hash, and
    # git commit so the catalog can be audited without relying on the filename.
    # The filename may say "rebound" while refine actually ran Kepler when tiered.
    run_id = Time.now().utc.iso.replace(" ", "T").replace(":", "").replace("-", "")[:15] + "Z"
    write_detection_sidecar(
        out_path,
        run_id=run_id,
        n_encounters=len(results),
        cfg=cfg,
        mpcorb_path=Path(snap.path),
        coarse_step_hours=coarse_step_hours,
        fine_step_hours=fine_step_hours,
        use_tiered=use_tiered,
        force_kepler_refine=use_tiered,
    )

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
