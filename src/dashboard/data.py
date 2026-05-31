"""Dashboard data layer — memory-safe access to the encounter catalog.

Separated from ``app.py`` so it is pure (no Streamlit) and unit-testable. The
characterised catalog can be either the small 158 k in-memory run or the full
72 M-row streaming run (``encounters_characterized_full.parquet``, ~5.8 GB).
Loading the latter eagerly would exhaust RAM, so:

- global headline stats (total, Gaia-observable, closest approach) come from a
  cheap lazy aggregation that never materialises the frame;
- interactive views load only the **closest N** encounters (a lazy top-k on
  ``dist_au``) — exploration is about close approaches, so this is the useful
  slice and its memory is bounded by ``cap``.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl

# Prefer the full streaming-characterised catalog; fall back to the 158 k run.
_FULL_CATALOG = Path("data/output/encounters_characterized_full.parquet")
_SMALL_CATALOG = Path("data/output/encounters_characterized.parquet")

# Above this row count, interactive views load only the closest-N subset.
DISPLAY_CAP = 300_000


def resolve_catalog_path() -> Path | None:
    """Return the best available characterised catalog, or None if absent."""
    if _FULL_CATALOG.exists():
        return _FULL_CATALOG
    if _SMALL_CATALOG.exists():
        return _SMALL_CATALOG
    return None


def catalog_stats(path: Path) -> dict:
    """Global headline stats via a streaming lazy aggregation (no full load).

    Returns ``{n_total, n_gaia_observable, dist_min_au}``.
    """
    agg = (
        pl.scan_parquet(path)
        .select(
            pl.len().alias("n_total"),
            pl.col("gaia_observable").sum().alias("n_gaia"),
            pl.col("dist_au").min().alias("d_min"),
        )
        .collect(engine="streaming")
    )
    row = agg.row(0, named=True)
    return {
        "n_total": int(row["n_total"]),
        "n_gaia_observable": int(row["n_gaia"] or 0),
        "dist_min_au": float(row["d_min"]),
    }


def load_catalog_display(path: Path, cap: int = DISPLAY_CAP) -> tuple[pl.DataFrame, bool]:
    """Load a bounded, display-ready slice of the catalog.

    If the catalog has more than ``cap`` rows, returns the ``cap`` closest
    encounters (lazy top-k on ``dist_au``) and ``capped=True``; otherwise the
    whole catalog and ``capped=False``.
    """
    n_total = pl.scan_parquet(path).select(pl.len()).collect().row(0)[0]
    if n_total > cap:
        df = pl.scan_parquet(path).sort("dist_au").head(cap).collect(engine="streaming")
        return df, True
    return pl.read_parquet(path), False
