"""Catalog writer — enforces typed schema and writes reproducibility sidecar."""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import logging
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
