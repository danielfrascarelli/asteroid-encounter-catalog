"""Download and cache Gaia DR3 SSO observations via TAP.

The Gaia archive table ``gaiadr3.sso_observation`` contains one row per
asteroid transit. The ``epoch`` column is in TCB (Barycentric Coordinate Time)
as a Julian Date. All other processing converts to TDB via
``src.utils.time_utils.tcb_to_tdb`` before use.

Download strategy
-----------------
Instead of sequential OFFSET pagination (O(n²) server-side scans), the
download splits the ``number_mp`` range into fixed-size batches and fetches
up to *n_workers* batches concurrently via ``ThreadPoolExecutor``.  Small
batches (default 5 000) keep each async TAP job under ~3 minutes, avoiding
the server-side connection resets that occur on long-running jobs.
Unnumbered objects (``number_mp IS NULL``) are fetched as a separate job.

Usage (via download script):
    docker compose run --rm pipeline python -m scripts.download_gaia_sso
"""

from __future__ import annotations

import logging
import os
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


def _fetch_range(
    archive_url: str,
    col_list: str,
    mp_start: int | None,
    mp_end: int | None,
) -> tuple[pl.DataFrame, float]:
    """Fetch one number_mp sub-range (or NULL) via a single async TAP job.

    Creates its own TapPlus instance so it is safe to call from threads.

    Returns
    -------
    tuple of (DataFrame, elapsed_seconds)
    """
    if mp_start is None:
        where = "number_mp IS NULL"
    else:
        where = f"number_mp BETWEEN {mp_start} AND {mp_end}"

    adql = f"SELECT {col_list} FROM gaiadr3.sso_observation WHERE {where}"
    logger.debug("TAP async job: %s", adql)

    t0 = time.monotonic()
    tap = TapPlus(url=archive_url)
    job = tap.launch_job_async(adql)
    table = job.get_results()
    elapsed = time.monotonic() - t0

    if len(table) == 0:
        return pl.DataFrame(), elapsed
    return pl.from_pandas(table.to_pandas()), elapsed


def _query_mp_max(archive_url: str) -> int:
    """Return MAX(number_mp) from the table (fast single-row query)."""
    tap = TapPlus(url=archive_url)
    job = tap.launch_job("SELECT MAX(number_mp) AS mp_max FROM gaiadr3.sso_observation")
    table = job.get_results()
    val = table["mp_max"][0]
    return int(val) if val is not None else _DEFAULT_MP_MAX


_DEFAULT_BATCH_SIZE = 5_000


def download_gaia_sso(
    archive_url: str,
    columns: list[str],
    dest: str | Path,
    *,
    n_workers: int | str = "auto",
    mp_max: int | None = None,
    batch_size: int = _DEFAULT_BATCH_SIZE,
) -> pl.DataFrame:
    """Download ``gaiadr3.sso_observation`` via parallel TAP jobs.

    Splits the ``number_mp`` range into batches of *batch_size* and fetches
    up to *n_workers* batches concurrently.  Using small batches keeps each
    TAP job short (< 5 min), avoiding server-side connection resets that occur
    on long-running async jobs.

    Parameters
    ----------
    archive_url:
        Base URL of the Gaia TAP server.
    columns:
        Column names to retrieve. Must include all of ``_REQUIRED_COLUMNS``.
    dest:
        Output path for the Parquet file (written once, at the end).
    n_workers:
        Number of parallel TAP connections. ``"auto"`` uses
        ``min(os.cpu_count(), 8)`` (Gaia TAP handles ~8 concurrent jobs
        without throttling).
    mp_max:
        Upper bound of ``number_mp`` range. If ``None``, queried first with
        ``MAX(number_mp)``.
    batch_size:
        Number of ``number_mp`` values per TAP job. Smaller values reduce
        per-job duration and the risk of server-side connection resets.
        Default: 5 000 (≈ 2–3 min per job on the Gaia archive).

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

    if n_workers == "auto":
        n_workers = min(os.cpu_count() or 4, 8)

    col_list = ", ".join(columns)

    if mp_max is None:
        logger.info("Querying MAX(number_mp) from %s…", archive_url)
        mp_max = _query_mp_max(archive_url)

    step = max(1, batch_size)
    ranges: list[tuple[int | None, int | None]] = [
        (start, min(start + step - 1, mp_max))
        for start in range(1, mp_max + 1, step)
    ]
    ranges.append((None, None))  # unnumbered objects

    logger.info(
        "gaiadr3.sso_observation — %d parallel workers | %d ranges | "
        "batch_size %d | number_mp 1–%d + unnumbered",
        n_workers, len(ranges), step, mp_max,
    )

    frames: list[pl.DataFrame] = []
    total_ranges = len(ranges)
    completed = 0

    with ThreadPoolExecutor(max_workers=n_workers) as pool:
        futures = {
            pool.submit(_fetch_range, archive_url, col_list, mp_start, mp_end): (mp_start, mp_end)
            for mp_start, mp_end in ranges
        }
        for fut in as_completed(futures):
            mp_start, mp_end = futures[fut]
            completed += 1
            pct = 100 * completed / total_ranges
            label = "unnumbered" if mp_start is None else f"mp {mp_start}–{mp_end}"
            try:
                df, elapsed = fut.result()
                if len(df) > 0:
                    frames.append(df)
                logger.info("  [%d/%d | %5.1f%%] %-22s %.1fs",
                            completed, total_ranges, pct, label, elapsed)
            except Exception:
                logger.exception("  [%d/%d | %5.1f%%] %-22s FAILED",
                                 completed, total_ranges, pct, label)

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
