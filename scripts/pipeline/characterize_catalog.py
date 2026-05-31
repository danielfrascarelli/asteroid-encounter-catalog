"""Enrich the encounter catalog with physical and observational properties.

Reads the detection output (encounters_catalog.parquet), characterizes each
encounter, and writes an enriched catalog (encounters_characterized.parquet)
along with a JSON metadata sidecar.

MPCORB snapshot selection: same logic as ``scripts.pipeline.run_pipeline`` —
the snapshot whose epoch is closest to the observation-window centre is used,
not the current ``data/raw/MPCORB.DAT``.  This keeps provenance consistent
between detection and characterisation.  ``--mpcorb`` overrides the
auto-selection if needed.

Usage:
    docker compose run --rm pipeline python -m scripts.pipeline.characterize_catalog
    docker compose run --rm pipeline python -m scripts.pipeline.characterize_catalog --config config.yaml
"""

from __future__ import annotations

import argparse
import dataclasses
import logging
import time
from datetime import UTC, datetime
from pathlib import Path

import polars as pl
import pyarrow.parquet as pq
from astropy.time import Time

from src.catalog.writer import write_catalog
from src.characterize.encounter import characterize_catalog, characterize_catalog_streaming
from src.ingest.gaia_orbits import load_gaia_orbits
from src.ingest.mpcorb import parse_mpcorb
from src.ingest.mpcorb_archive import discover_snapshots, select_for_window
from src.utils.config import load_config

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)

_REQUIRED_BODIES = [1, 2, 4, 10]
_BODY_NAMES = {1: "Ceres", 2: "Pallas", 4: "Vesta", 10: "Hygiea"}


def _supplement_elements(elements: pl.DataFrame, mpcorb: pl.DataFrame) -> pl.DataFrame:
    """Add major bodies missing from gaia_orbits using MPCORB elements."""
    present = set(elements["number"].to_list())
    missing = [n for n in _REQUIRED_BODIES if n not in present]
    if not missing:
        return elements
    logger.info("Supplementing elements from MPCORB for bodies: %s", missing)
    supplement = mpcorb.filter(pl.col("number").is_in(missing)).select(elements.columns)
    return pl.concat([supplement, elements])


def main() -> int:
    parser = argparse.ArgumentParser(description="Characterize encounter catalog")
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument(
        "--input",
        default=None,
        help="Override input catalog path (default: data/output/encounters_catalog.parquet)",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Override output catalog path (default: data/output/encounters_characterized.parquet)",
    )
    parser.add_argument(
        "--mpcorb",
        default=None,
        help=(
            "Explicit MPCORB.DAT path to use.  Default: auto-select the "
            "archived snapshot whose epoch is closest to the window centre, "
            "matching what scripts.pipeline.run_pipeline does."
        ),
    )
    parser.add_argument(
        "--streaming",
        choices=["auto", "on", "off"],
        default="auto",
        help=(
            "Chunked streaming characterisation to bound peak RAM. 'auto' "
            "(default) streams when the input exceeds --streaming-threshold rows; "
            "'on'/'off' force the mode. Streaming output is NOT globally sorted "
            "by dist_au (preserves input order)."
        ),
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=1_000_000,
        help="Rows per chunk in streaming mode (default: 1,000,000).",
    )
    parser.add_argument(
        "--streaming-threshold",
        type=int,
        default=2_000_000,
        help="Row count above which 'auto' uses streaming (default: 2,000,000).",
    )
    args = parser.parse_args()

    cfg = load_config(args.config)
    out_dir = Path(cfg.paths.output)

    in_path = Path(args.input) if args.input else out_dir / "encounters_catalog.parquet"

    # --- Decide streaming vs in-memory (cheap: read parquet row count metadata) ---
    n_rows = pq.ParquetFile(in_path).metadata.num_rows
    if args.streaming == "on":
        use_streaming = True
    elif args.streaming == "off":
        use_streaming = False
    else:
        use_streaming = n_rows > args.streaming_threshold
    # Streaming defaults to a distinct filename so it never clobbers the small
    # 158k in-memory run that downstream tooling/tests still point at.
    default_out = (
        "encounters_characterized_full.parquet"
        if use_streaming
        else "encounters_characterized.parquet"
    )
    out_path = Path(args.output) if args.output else out_dir / default_out
    logger.info(
        "Input %s has %d rows → %s characterisation (chunk_size=%d)",
        in_path,
        n_rows,
        "STREAMING" if use_streaming else "in-memory",
        args.chunk_size,
    )

    # --- Load orbital elements ---
    orbits_path = Path(cfg.paths.raw) / "gaia_orbits.parquet"
    elements = load_gaia_orbits(orbits_path)

    # --- Resolve MPCORB snapshot ---
    # If --mpcorb wasn't passed, mirror run_pipeline's logic: pick the archived
    # snapshot whose epoch is closest to the observation-window centre.  This
    # avoids the silent inconsistency where characterisation uses a newer
    # MPCORB.DAT than the one that produced the detection catalog.
    if args.mpcorb is not None:
        mpcorb_path = Path(args.mpcorb)
    else:
        tw = cfg.time_window
        t_start = Time(tw.start, scale=tw.scale).tdb.jd
        t_end = Time(tw.end, scale=tw.scale).tdb.jd
        snapshots = discover_snapshots(Path(cfg.paths.raw))
        if snapshots:
            snap = select_for_window(snapshots, t_start, t_end)
            mpcorb_path = Path(snap.path)
            logger.info(
                "Auto-selected MPCORB snapshot %s (epoch %s)",
                mpcorb_path.name,
                Time(snap.epoch_jd, format="jd", scale="tdb").utc.iso[:10],
            )
        else:
            mpcorb_path = Path(cfg.paths.raw) / "MPCORB.DAT"
            logger.warning(
                "No archived snapshots found under %s; falling back to MPCORB.DAT",
                cfg.paths.raw,
            )
    logger.info("Parsing MPCORB: %s", mpcorb_path)
    mpcorb = parse_mpcorb(str(mpcorb_path))

    # --- Supplement elements for major bodies not in gaia_orbits ---
    elements = _supplement_elements(elements, mpcorb)

    run_id = datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- Streaming path: bounded RAM, disk-to-disk, no global sort ---
    if use_streaming:
        t0 = time.monotonic()
        summary = characterize_catalog_streaming(
            str(in_path),
            elements,
            mpcorb,
            str(out_path),
            run_id,
            albedo=cfg.characterize.default_albedo,
            chunk_size=args.chunk_size,
            mpcorb_path=mpcorb_path,
            config_dict=dataclasses.asdict(cfg),
        )
        logger.info(
            "Streaming characterization complete in %.1fs — %d rows, %d Gaia-observable, %d chunks",
            time.monotonic() - t0,
            summary["n_encounters"],
            summary["n_gaia_observable"],
            summary["n_chunks"],
        )
        for n in _REQUIRED_BODIES:
            g = summary["gate"][str(n)]
            if g["present"]:
                logger.info(
                    "Gate OK: (%d) %s — %d encounters, closest %.6f AU",
                    n,
                    _BODY_NAMES[n],
                    g["n_encounters"],
                    g["closest_au"],
                )
            else:
                logger.warning("Gate: (%d) %s not in enriched catalog", n, _BODY_NAMES[n])
        return 0

    # --- In-memory path (small catalogs, e.g. the 158k detection run) ---
    logger.info("Loading detection catalog: %s", in_path)
    encounters = pl.read_parquet(in_path)
    logger.info("%d encounters loaded", len(encounters))

    t0 = time.monotonic()
    enriched = characterize_catalog(
        encounters,
        elements,
        mpcorb,
        albedo=cfg.characterize.default_albedo,
    )
    elapsed = time.monotonic() - t0
    logger.info("Characterization complete in %.1fs", elapsed)

    # --- Save catalog + metadata sidecar ---
    write_catalog(
        enriched,
        out_path,
        run_id=run_id,
        mpcorb_path=mpcorb_path,
        config_dict=dataclasses.asdict(cfg),
    )

    # --- Summary stats ---
    observable = enriched.filter(pl.col("gaia_observable"))
    logger.info("Gaia-observable encounters: %d / %d", len(observable), len(enriched))
    logger.info(
        "Velocity range (km/s): %.3f – %.3f",
        float(enriched["rel_vel_km_s"].min()),  # type: ignore[arg-type]
        float(enriched["rel_vel_km_s"].max()),  # type: ignore[arg-type]
    )
    logger.info(
        "Diameter range body 1 (km): %.1f – %.1f",
        float(enriched["diameter_1_km"].drop_nulls().min()),  # type: ignore[arg-type]
        float(enriched["diameter_1_km"].drop_nulls().max()),  # type: ignore[arg-type]
    )

    # --- Gate check: major bodies ---
    for n in _REQUIRED_BODIES:
        hits = enriched.filter((pl.col("number_1") == n) | (pl.col("number_2") == n))
        if len(hits) == 0:
            logger.warning("Gate: (%d) not in enriched catalog", n)
        else:
            row = hits.head(1).row(0, named=True)
            logger.info(
                "Gate OK: (%d) — closest %.6f AU  D=%.0f km  class=%s",
                n,
                float(hits["dist_au"].min()),  # type: ignore[arg-type]
                (
                    float(hits.filter(pl.col("number_1") == n).head(1)["diameter_1_km"][0])
                    if len(hits.filter(pl.col("number_1") == n)) > 0
                    else float("nan")
                ),
                row["class_1"] if row["number_1"] == n else row["class_2"],
            )

    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main())
