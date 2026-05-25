"""Catalog writer — enforces typed schema and writes reproducibility sidecar."""

from __future__ import annotations

import dataclasses
import hashlib
import importlib.metadata
import json
import logging
import subprocess
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import polars as pl

from src.catalog.schema import CATALOG_SCHEMA

logger = logging.getLogger(__name__)


def _hash_file(path: Path) -> str:
    """Return a 16-char SHA-256 prefix for a file."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:16]


def _dep_versions() -> dict[str, str]:
    """Collect installed versions of key pipeline dependencies."""
    deps = ["astropy", "scipy", "polars", "numpy"]
    out: dict[str, str] = {}
    for dep in deps:
        try:
            out[dep] = importlib.metadata.version(dep)
        except importlib.metadata.PackageNotFoundError:
            out[dep] = "unknown"
    return out


def _git_commit_info() -> dict[str, str]:
    """Best-effort git commit + dirty flag for provenance. Empty if unavailable."""
    out: dict[str, str] = {}
    try:
        sha = subprocess.run(
            ["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=True
        ).stdout.strip()
        out["commit"] = sha
        status = subprocess.run(
            ["git", "status", "--porcelain"], capture_output=True, text=True, check=True
        ).stdout.strip()
        out["dirty"] = "true" if status else "false"
    except (subprocess.CalledProcessError, FileNotFoundError):
        pass
    return out


def write_detection_sidecar(
    catalog_path: Path,
    run_id: str,
    n_encounters: int,
    cfg: Any,
    *,
    mpcorb_path: Path | None = None,
    coarse_step_hours: float | None = None,
    fine_step_hours: float | None = None,
    use_tiered: bool = False,
    force_kepler_refine: bool = False,
) -> Path:
    """Write a provenance sidecar JSON next to a *detection* catalog parquet.

    Unlike :func:`write_catalog` (which writes a sidecar for the *characterised*
    catalog), this captures the parameters that produced the geometric encounter
    list: scan method, refine method, time-grid resolution, prefilter, and the
    full propagation config.  This makes it auditable whether a file labelled
    ``..._rebound_...`` actually used N-body for the final reported distance
    (it usually does NOT when ``use_tiered`` is True: the scan uses the cached
    N-body trajectory, but the refinement falls back to a Kepler propagator).

    Parameters
    ----------
    catalog_path:
        Path of the detection parquet that was just written.
    run_id:
        Unique run identifier (e.g. ISO timestamp).
    n_encounters:
        Number of rows in the catalog.
    cfg:
        Pipeline config object (must support ``dataclasses.asdict``).
    mpcorb_path:
        MPCORB snapshot used; hashed for provenance.
    coarse_step_hours, fine_step_hours:
        Resolved time-grid parameters (after ``coarse_step_hours or fine`` fallback).
    use_tiered:
        Whether the bulk grid ran at ``coarse_step_hours`` while refinement used Kepler.
    force_kepler_refine:
        Whether refinement was forced to Kepler regardless of scan method.

    Returns
    -------
    Path
        Where the sidecar JSON was written (``<stem>_provenance.json`` next to
        the catalog).
    """
    scan_method = cfg.propagation.method.lower()
    refine_method = "kepler" if (force_kepler_refine or scan_method == "kepler") else scan_method

    meta: dict[str, Any] = {
        "run_id": run_id,
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "catalog_path": str(catalog_path),
        "n_encounters": n_encounters,
        "threshold_au": cfg.detection.threshold_au,
        "scan": {
            "method": scan_method,
            "grid_step_hours": coarse_step_hours,
            "rebound": (
                dataclasses.asdict(cfg.propagation.rebound) if scan_method == "rebound" else None
            ),
        },
        "refine": {
            "method": refine_method,
            "enabled": cfg.detection.refinement.enabled,
            "fine_step_seconds": cfg.detection.refinement.fine_time_step_seconds,
            "window_hours": cfg.detection.refinement.window_hours,
            "forced_by_tiered": force_kepler_refine and scan_method != "kepler",
        },
        "tiered_mode": use_tiered,
        "fine_step_hours": fine_step_hours,
        "prefilter": {
            "enabled": cfg.detection.prefilter.enabled,
            "semimajor_diff_max_au": cfg.detection.prefilter.semimajor_diff_max_au,
            "inclination_diff_max_deg": cfg.detection.prefilter.inclination_diff_max_deg,
        },
        "config": dataclasses.asdict(cfg),
        "dependencies": _dep_versions(),
        "git": _git_commit_info(),
    }
    if mpcorb_path and mpcorb_path.exists():
        meta["mpcorb"] = {
            "path": str(mpcorb_path),
            "sha256_prefix": _hash_file(mpcorb_path),
            "size_bytes": mpcorb_path.stat().st_size,
        }

    sidecar = catalog_path.parent / (catalog_path.stem + "_provenance.json")
    sidecar.write_text(json.dumps(meta, indent=2, default=str))
    logger.info("Detection provenance sidecar → %s", sidecar)
    return sidecar


def write_catalog(
    df: pl.DataFrame,
    out_path: Path,
    run_id: str,
    mpcorb_path: Path | None = None,
    config_dict: dict[str, Any] | None = None,
) -> None:
    """Write enriched catalog with typed schema + metadata sidecar.

    Parameters
    ----------
    df:
        Enriched DataFrame from ``characterize_catalog()``.  Must NOT already
        contain a ``run_id`` column; it is added here.
    out_path:
        Destination ``.parquet`` path.
    run_id:
        Unique run identifier written to every row and the sidecar.
    mpcorb_path:
        Path to ``MPCORB.DAT`` for SHA-256 provenance hashing.  Optional.
    config_dict:
        Serialisable pipeline config for the sidecar.  Optional.
    """
    enriched = df.with_columns(pl.lit(run_id).alias("run_id"))

    # Select and cast to schema column order (skip columns absent from df)
    present = [c for c in CATALOG_SCHEMA if c in enriched.columns]
    enriched = enriched.select(present).with_columns(
        [pl.col(c).cast(CATALOG_SCHEMA[c]) for c in present]
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    enriched.write_parquet(out_path, compression="zstd")
    logger.info("Catalog written → %s  (%d rows)", out_path, len(enriched))

    # Build sidecar metadata
    meta: dict[str, Any] = {
        "run_id": run_id,
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "n_encounters": len(enriched),
        "n_gaia_observable": int(enriched["gaia_observable"].sum()),
        "catalog_path": str(out_path),
        "schema_columns": present,
        "dependencies": _dep_versions(),
    }
    if mpcorb_path and mpcorb_path.exists():
        meta["mpcorb"] = {
            "path": str(mpcorb_path),
            "sha256_prefix": _hash_file(mpcorb_path),
            "size_bytes": mpcorb_path.stat().st_size,
        }
    if config_dict is not None:
        meta["config"] = config_dict

    sidecar = out_path.parent / (out_path.stem + "_metadata.json")
    sidecar.write_text(json.dumps(meta, indent=2, default=str))
    logger.info("Metadata sidecar → %s", sidecar)
