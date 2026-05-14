"""Download and cache Gaia DR3 SSO observations via TAP.

The Gaia archive table ``gaiadr3.sso_observation`` contains one row per
asteroid transit. The ``epoch`` column is in TCB (Barycentric Coordinate Time)
as a Julian Date. All other processing converts to TDB via
``src.utils.time_utils.tcb_to_tdb`` before use.

Download strategy
-----------------
The ``number_mp`` range is split into fixed-size batches. Each batch is
fetched as a synchronous TAP query (/tap/sync) and written to its own Parquet
file in *cache_dir* immediately upon completion (atomic rename from ``.part``).

On rerun, completed chunks are skipped. If *batch_size* changed between runs,
existing chunks are reused when their union fully covers a new range — no
re-download needed. Failed batches are retried with exponential backoff
(30 s, 60 s, 120 s ± 10 % jitter) before being marked as permanently failed.

Usage (via download script):
    docker compose run --rm pipeline python -m scripts.download_gaia_sso
"""

from __future__ import annotations

import logging
import os
import random
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import polars as pl
from astroquery.utils.tap.core import TapPlus

logger = logging.getLogger(__name__)

_REQUIRED_COLUMNS = {
    "solution_id",
    "source_id",
    "denomination",
    "number_mp",
    "epoch",
}

# Gaia DR3: numbered SSO objects go up to ~158 000
_DEFAULT_MP_MAX = 160_000
_DEFAULT_BATCH_SIZE = 5_000
_RETRY_BASE_SECONDS = 30.0


def _chunk_path(cache_dir: Path, mp_start: int | None, mp_end: int | None) -> Path:
    if mp_start is None:
        return cache_dir / "unnumbered.parquet"
    return cache_dir / f"mp_{mp_start:07d}_{mp_end:07d}.parquet"


def _build_chunk_from_cache(cache_dir: Path, mp_start: int, mp_end: int, dest: Path) -> bool:
    """Try to build chunk [mp_start, mp_end] from existing cached files.

    Scans cache_dir for mp_*.parquet files that overlap [mp_start, mp_end],
    checks if their union fully covers the range (no gaps), and if so merges
    them into dest (filtered to the exact range). Returns True if built from
    cache, False if a TAP fetch is needed.
    """
    existing: list[tuple[int, int, Path]] = []
    for p in cache_dir.glob("mp_*.parquet"):
        parts = p.stem.split("_")
        if len(parts) != 3:
            continue
        try:
            cs, ce = int(parts[1]), int(parts[2])
        except ValueError:
            continue
        if cs <= mp_end and ce >= mp_start:
            existing.append((cs, ce, p))

    if not existing:
        return False

    existing.sort()
    current = mp_start - 1
    covering: list[Path] = []
    for cs, ce, p in existing:
        if cs > current + 1:
            return False  # gap in coverage
        covering.append(p)
        current = max(current, ce)
    if current < mp_end:
        return False  # range extends past last chunk

    frames = [pl.read_parquet(p) for p in covering]
    merged = pl.concat(frames).filter(
        (pl.col("number_mp") >= mp_start) & (pl.col("number_mp") <= mp_end)
    )
    part = dest.with_suffix(".part")
    merged.write_parquet(part, compression="zstd")
    part.rename(dest)
    return True


def _fetch_range(
    archive_url: str,
    col_list: str,
    mp_start: int | None,
    mp_end: int | None,
    cache_path: Path,
    max_retries: int = 3,
) -> tuple[pl.DataFrame, float]:
    """Fetch one number_mp sub-range (or NULL) via a synchronous TAP query.

    Uses /tap/sync — the result comes back in the HTTP response body, so there
    is no job polling and no "Cannot find result" race condition.

    Retries with exponential backoff (±10 % jitter) on failure. Writes result
    atomically to cache_path before returning.

    Creates its own TapPlus instance so it is safe to call from threads.
    """
    if mp_start is None:
        where = "number_mp IS NULL"
    else:
        where = f"number_mp BETWEEN {mp_start} AND {mp_end}"
    adql = f"SELECT {col_list} FROM gaiadr3.sso_observation WHERE {where}"
    label = "unnumbered" if mp_start is None else f"mp {mp_start}–{mp_end}"

    last_exc: Exception | None = None
    for attempt in range(max_retries + 1):
        if attempt > 0:
            wait = _RETRY_BASE_SECONDS * (2 ** (attempt - 1))
            wait *= 1 + random.uniform(-0.1, 0.1)
            logger.warning(
                "  ↻ retry %d/%d for %s in %.0fs: %s",
                attempt,
                max_retries,
                label,
                wait,
                last_exc,
            )
            time.sleep(wait)
        try:
            t0 = time.monotonic()
            tap = TapPlus(url=archive_url)
            job = tap.launch_job(adql)
            table = job.get_results()
            elapsed = time.monotonic() - t0
            df = pl.DataFrame() if len(table) == 0 else pl.from_pandas(table.to_pandas())
            df = df.rename({c: c.lower() for c in df.columns if c != c.lower()})
            part = cache_path.with_suffix(".part")
            df.write_parquet(part, compression="zstd")
            part.rename(cache_path)
            return df, elapsed
        except Exception as exc:
            last_exc = exc
            is_last = attempt == max_retries
            logger.warning(
                "  %s attempt %d/%d failed: %s",
                label,
                attempt + 1,
                max_retries + 1,
                exc,
                exc_info=is_last,
            )

    raise RuntimeError(f"All {max_retries + 1} attempts failed for {label}") from last_exc


def _query_mp_max(archive_url: str) -> int:
    """Return MAX(number_mp) from the table (fast single-row query)."""
    tap = TapPlus(url=archive_url)
    job = tap.launch_job("SELECT MAX(number_mp) AS mp_max FROM gaiadr3.sso_observation")
    table = job.get_results()
    val = table["mp_max"][0]
    return int(val) if val is not None else _DEFAULT_MP_MAX


def download_gaia_sso(
    archive_url: str,
    columns: list[str],
    dest: str | Path,
    *,
    n_workers: int | str = "auto",
    mp_max: int | None = None,
    batch_size: int = _DEFAULT_BATCH_SIZE,
    cache_dir: Path | None = None,
    max_retries: int = 3,
) -> pl.DataFrame:
    """Download ``gaiadr3.sso_observation`` via parallel TAP jobs.

    Each batch is saved to *cache_dir* as soon as it completes. On rerun,
    completed chunks are skipped. If *batch_size* differs from a previous run,
    existing chunks are reused when they fully cover a new range (no re-download).

    Parameters
    ----------
    archive_url:
        Base URL of the Gaia TAP server.
    columns:
        Column names to retrieve. Must include all of ``_REQUIRED_COLUMNS``.
    dest:
        Output path for the merged Parquet file.
    n_workers:
        Number of parallel TAP connections. ``"auto"`` uses
        ``min(os.cpu_count(), 8)``.
    mp_max:
        Upper bound of ``number_mp`` range. If ``None``, queried first.
    batch_size:
        Number of ``number_mp`` values per TAP job.
    cache_dir:
        Directory for per-chunk Parquet files. Defaults to
        ``<dest.parent.parent>/cache/gaia_sso_chunks``.
    max_retries:
        Number of retry attempts per failed TAP job (exponential backoff).

    Returns
    -------
    polars.DataFrame
        All downloaded observations, concatenated and sorted by ``source_id``.

    Raises
    ------
    ValueError
        If *columns* is missing any required column.
    """
    missing = _REQUIRED_COLUMNS - set(columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)

    if cache_dir is None:
        cache_dir = dest.parent.parent / "cache" / "gaia_sso_chunks"
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    if n_workers == "auto":
        n_workers = min(os.cpu_count() or 4, 8)

    col_list = ", ".join(columns)

    if mp_max is None:
        logger.info("Querying MAX(number_mp) from %s…", archive_url)
        mp_max = _query_mp_max(archive_url)

    step = max(1, batch_size)
    ranges: list[tuple[int | None, int | None]] = [
        (start, min(start + step - 1, mp_max)) for start in range(1, mp_max + 1, step)
    ]
    ranges.append((None, None))

    # Exact-name resume: chunks already present for this batch_size
    cached = [(s, e) for s, e in ranges if _chunk_path(cache_dir, s, e).exists()]
    pending = [(s, e) for s, e in ranges if not _chunk_path(cache_dir, s, e).exists()]

    # Cross-batch reuse: build pending numbered chunks from existing chunks of any batch_size
    rebuilt: list[tuple[int | None, int | None]] = []
    still_pending: list[tuple[int | None, int | None]] = []
    for s, e in pending:
        if (
            s is not None
            and e is not None
            and _build_chunk_from_cache(cache_dir, s, e, _chunk_path(cache_dir, s, e))
        ):
            rebuilt.append((s, e))
        else:
            still_pending.append((s, e))

    if rebuilt:
        logger.info("  Rebuilt %d chunk(s) from prior-run cache (no TAP needed)", len(rebuilt))
    cached = cached + rebuilt
    pending = still_pending

    logger.info(
        "gaiadr3.sso_observation — %d cached | %d pending | "
        "%d workers | batch_size %d | number_mp 1–%d + unnumbered",
        len(cached),
        len(pending),
        n_workers,
        step,
        mp_max,
    )

    failed: list[tuple[int | None, int | None]] = []
    completed = len(cached)

    if pending:
        with ThreadPoolExecutor(max_workers=int(n_workers)) as pool:
            futures = {
                pool.submit(
                    _fetch_range,
                    archive_url,
                    col_list,
                    s,
                    e,
                    _chunk_path(cache_dir, s, e),
                    max_retries,
                ): (s, e)
                for s, e in pending
            }
            for fut in as_completed(futures):
                s, e = futures[fut]
                completed += 1
                pct = 100 * completed / len(ranges)
                label = "unnumbered" if s is None else f"mp {s}–{e}"
                try:
                    _df, elapsed = fut.result()
                    logger.info(
                        "  [%d/%d | %5.1f%%] %-22s %.1fs",
                        completed,
                        len(ranges),
                        pct,
                        label,
                        elapsed,
                    )
                except Exception:
                    logger.exception(
                        "  [%d/%d | %5.1f%%] %-22s FAILED (all retries exhausted)",
                        completed,
                        len(ranges),
                        pct,
                        label,
                    )
                    failed.append((s, e))

    if failed:
        logger.warning(
            "%d ranges failed permanently: %s",
            len(failed),
            [f"mp {s}–{e}" if s else "unnumbered" for s, e in failed],
        )

    frames: list[pl.DataFrame] = []
    missing_chunks = []
    for s, e in ranges:
        cp = _chunk_path(cache_dir, s, e)
        if cp.exists():
            chunk = pl.read_parquet(cp)
            if len(chunk) > 0:
                frames.append(chunk)
        else:
            missing_chunks.append((s, e))

    if missing_chunks:
        logger.warning(
            "No chunk file for %d ranges — those observations are absent.", len(missing_chunks)
        )

    if not frames:
        logger.warning("No rows returned from Gaia TAP — writing empty file")
        result = pl.DataFrame()
    else:
        result = pl.concat(frames).sort("source_id")

    result.write_parquet(dest, compression="zstd")
    logger.info("Saved %d observations to %s", len(result), dest)
    return result


def load_gaia_sso(path: str | Path) -> pl.DataFrame:
    """Load a previously downloaded Gaia SSO Parquet file.

    Parameters
    ----------
    path:
        Path to the Parquet file written by :func:`download_gaia_sso`.

    Returns
    -------
    polars.DataFrame

    Notes
    -----
    The ``epoch`` column is in TCB. Use
    ``src.utils.time_utils.tcb_to_tdb`` to convert to TDB before
    propagation.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Gaia SSO file not found: {path}")
    return pl.read_parquet(path)
