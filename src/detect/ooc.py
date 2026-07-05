"""Out-of-core detection: scan → disk shards → DuckDB dedup → batched refine.

The in-memory detection path holds the full set of unique candidate pairs in a
Python dict and then the refined catalogue as a DataFrame. At the scale of the
numbered main-belt population (~449 k bodies → tens of millions of unique
candidate pairs) this dict alone is ~15--20 GB and, together with the trajectory
working set, overflows a 24 GB Docker Desktop VM (OOM-kill, ExitCode 137, observed
at ~40/104 scan chunks regardless of worker count).

This module keeps parent memory bounded:

1. **Scan → shards.** The parallel scan streams each chunk's candidate list to a
   Parquet shard on disk (via the ``on_chunk`` sink of :func:`scan_parallel`);
   nothing accumulates in RAM.
2. **Dedup out-of-core.** DuckDB reads all shards and keeps, per ``(idx_i,
   idx_j)`` pair, the row with the minimum ``(d_coarse_au, t_coarse_jd)`` — the
   same lexicographic tie-break as the in-memory ``_merge_into`` — writing a
   single deduped candidates Parquet. DuckDB spills to disk, so this never
   materialises the full set in RAM.
3. **Batched refine + streaming write.** The deduped candidates are read in
   batches; each batch is refined with the ordinary Kepler refiner and appended
   to the output Parquet via a streaming writer. Peak memory is one batch.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import polars as pl
import pyarrow.parquet as pq

from src.catalog.schema import CATALOG_SCHEMA  # noqa: F401  (kept for callers)
from src.detect.refine import refine_candidates

logger = logging.getLogger(__name__)

_SHARD_SCHEMA = {
    "idx_i": pl.Int64,
    "idx_j": pl.Int64,
    "t_coarse_jd": pl.Float64,
    "d_coarse_au": pl.Float64,
}


def make_shard_sink(shard_dir: Path):
    """Return an ``on_chunk`` callback that writes each scan chunk to a shard.

    The returned callable is passed as ``scan_parallel(on_chunk=...)``; it keeps a
    monotonically increasing shard index and returns nothing. Empty chunks are
    skipped (no zero-row shard files).
    """
    shard_dir = Path(shard_dir)
    shard_dir.mkdir(parents=True, exist_ok=True)
    state = {"i": 0}

    def sink(chunk_result: list[tuple[int, int, float, float]]) -> None:
        if not chunk_result:
            return
        idx_i, idx_j, t_c, d_c = zip(*chunk_result)
        pl.DataFrame(
            {
                "idx_i": pl.Series(idx_i, dtype=pl.Int64),
                "idx_j": pl.Series(idx_j, dtype=pl.Int64),
                "t_coarse_jd": pl.Series(t_c, dtype=pl.Float64),
                "d_coarse_au": pl.Series(d_c, dtype=pl.Float64),
            },
            schema=_SHARD_SCHEMA,
        ).write_parquet(shard_dir / f"chunk_{state['i']:05d}.parquet", compression="zstd")
        state["i"] += 1

    return sink


def dedup_shards(shard_dir: Path, out_parquet: Path, memory_limit: str = "4GB") -> int:
    """Dedup scan shards out-of-core with DuckDB → single candidates Parquet.

    Keeps, per ``(idx_i, idx_j)``, the row with minimum ``(d_coarse_au,
    t_coarse_jd)`` — identical to the in-memory ``_merge_into`` tie-break.

    Returns the number of unique candidate pairs written. Returns 0 (and writes an
    empty, correctly-typed Parquet) when there are no shards.
    """
    import duckdb

    shard_dir = Path(shard_dir)
    out_parquet = Path(out_parquet)
    out_parquet.parent.mkdir(parents=True, exist_ok=True)
    shards = sorted(shard_dir.glob("chunk_*.parquet"))
    if not shards:
        pl.DataFrame(schema=_SHARD_SCHEMA).write_parquet(out_parquet)
        return 0

    con = duckdb.connect()
    con.execute(f"PRAGMA memory_limit='{memory_limit}'")
    con.execute(f"PRAGMA temp_directory='{shard_dir / '_duckdb_tmp'}'")
    glob = str(shard_dir / "chunk_*.parquet")
    con.execute(
        """
        COPY (
            SELECT idx_i, idx_j, t_coarse_jd, d_coarse_au
            FROM read_parquet($glob)
            QUALIFY row_number() OVER (
                PARTITION BY idx_i, idx_j ORDER BY d_coarse_au, t_coarse_jd
            ) = 1
        ) TO $out (FORMAT PARQUET, COMPRESSION ZSTD)
        """,
        {"glob": glob, "out": str(out_parquet)},
    )
    row = con.execute(f"SELECT count(*) FROM read_parquet('{out_parquet}')").fetchone()
    con.close()
    return int(row[0]) if row else 0


def refine_streaming(
    candidates_parquet: Path,
    elements: pl.DataFrame,
    out_path: Path,
    *,
    threshold_au: float,
    fine_step_seconds: float,
    window_hours: float,
    n_workers: int,
    batch_size: int = 2_000_000,
) -> int:
    """Refine deduped candidates in batches, appending to *out_path* incrementally.

    Reads ``candidates_parquet`` in ``batch_size`` row batches (bounding memory),
    refines each with the Kepler refiner, and streams the result to a single
    output Parquet via a :class:`pyarrow.parquet.ParquetWriter`. The final catalog
    is NOT globally sorted by distance (unlike the in-memory path) — sort
    downstream if needed. Returns the total number of refined encounters written.
    """
    from src.detect.pipeline import _SCHEMA

    candidates_parquet = Path(candidates_parquet)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    pf = pq.ParquetFile(candidates_parquet)
    writer: pq.ParquetWriter | None = None
    n_total = 0
    empty = pl.DataFrame(schema=_SCHEMA)
    try:
        for batch in pf.iter_batches(batch_size=batch_size):
            tbl = pl.from_arrow(batch)
            if isinstance(tbl, pl.Series):  # single-column safety
                tbl = tbl.to_frame()
            if tbl.is_empty():
                continue
            cands = list(
                zip(
                    tbl["idx_i"].to_list(),
                    tbl["idx_j"].to_list(),
                    tbl["t_coarse_jd"].to_list(),
                    tbl["d_coarse_au"].to_list(),
                )
            )
            refined = refine_candidates(
                elements,
                cands,
                threshold_au,
                fine_step_seconds,
                window_hours,
                n_workers=n_workers,
            )
            if refined.is_empty():
                continue
            arrow = refined.to_arrow()
            if writer is None:
                writer = pq.ParquetWriter(out_path, arrow.schema, compression="zstd")
            writer.write_table(arrow)
            n_total += len(refined)
    finally:
        if writer is not None:
            writer.close()
    if writer is None:  # no encounter survived refinement → write empty typed file
        empty.write_parquet(out_path, compression="zstd")
    logger.info("Streaming refine complete: %d encounters → %s", n_total, out_path)
    return n_total


def detect_encounters_ooc(
    elements: pl.DataFrame,
    time_grid: np.ndarray,
    *,
    threshold_au: float,
    leaf_size: int,
    fine_step_seconds: float,
    window_hours: float,
    n_workers: int | str,
    chunk_size_days: float,
    out_path: Path,
    spill_dir: Path,
    positions: np.ndarray | None = None,
    query_radius_au: float | None = None,
    refine_batch_size: int = 2_000_000,
    keep_scan_checkpoint: bool = True,
) -> int:
    """Full out-of-core detection: scan→shards→dedup→batched refine→*out_path*.

    Writes the final catalogue to *out_path* and returns the row count. Parent
    memory stays bounded (one chunk during scan, one batch during refine); DuckDB
    handles the dedup out-of-core. The deduped candidates Parquet is retained at
    ``spill_dir/scan_candidates.parquet`` (resumable refine) unless
    ``keep_scan_checkpoint`` is False.

    Prefilter is intentionally not applied here — the out-of-core path targets the
    large-N production run where the pair-list prefilter is skipped anyway.
    """
    from src.detect.parallel import resolve_n_workers, scan_parallel

    spill_dir = Path(spill_dir)
    shard_dir = spill_dir / "scan_shards"
    candidates_parquet = spill_dir / "scan_candidates.parquet"
    nw = resolve_n_workers(n_workers)

    logger.info("OOC detection: scan → shards in %s", shard_dir)
    sink = make_shard_sink(shard_dir)
    scan_parallel(
        elements,
        time_grid,
        None,  # no prefilter pair-list at production scale
        threshold_au,
        leaf_size,
        n_workers,
        chunk_size_days,
        positions=positions,
        query_radius_au=query_radius_au,
        on_chunk=sink,
    )

    logger.info("OOC detection: DuckDB dedup of shards → %s", candidates_parquet)
    n_unique = dedup_shards(shard_dir, candidates_parquet)
    logger.info("OOC detection: %d unique candidate pairs", n_unique)

    n_out = refine_streaming(
        candidates_parquet,
        elements,
        out_path,
        threshold_au=threshold_au,
        fine_step_seconds=fine_step_seconds,
        window_hours=window_hours,
        n_workers=nw,
        batch_size=refine_batch_size,
    )

    # Shards are large and redundant once deduped; drop them, keep the deduped
    # candidates parquet as the resumable checkpoint.
    for s in shard_dir.glob("chunk_*.parquet"):
        s.unlink()
    tmp = shard_dir / "_duckdb_tmp"
    if tmp.exists():
        for f in tmp.glob("*"):
            f.unlink()
        tmp.rmdir()
    try:
        shard_dir.rmdir()
    except OSError:
        pass
    if not keep_scan_checkpoint and candidates_parquet.exists():
        candidates_parquet.unlink()
    return n_out
