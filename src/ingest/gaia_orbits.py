"""Download and cache Gaia DR3 SSO orbital elements via TAP.

``gaiadr3.sso_orbits`` has one row per numbered asteroid with osculating
orbital elements derived from Gaia observations.  The osculating epoch
(``osc_epoch``) falls within the Gaia DR3 observation window (2014–2017),
so propagation to that window requires at most a few months rather than the
~10-year backward extrapolation that MPCORB (epoch ~2025) requires.

Column mapping (Gaia → pipeline schema)
-----------------------------------------
``gaiadr3.sso_orbits``          pipeline column       unit
------------------------------  --------------------  ---------------
number_mp      (Int64)       →  number    (Int32)
denomination   (String)      →  designation (Utf8)
semi_major_axis (AU)         →  a_au       (Float64)
eccentricity                 →  e          (Float64)
inclination    (rad)         →  i_deg      (Float64)  × 180/π
arg_perihelion (rad)         →  omega_deg  (Float64)  × 180/π
long_asc_node  (rad)         →  Omega_deg  (Float64)  × 180/π
mean_anomaly   (rad)         →  M_deg      (Float64)  × 180/π
osc_epoch      (JD, TDB)     →  epoch_jd   (Float64)

Download strategy
-----------------
Identical to ``gaia_sso``: number_mp range split into fixed-size batches,
each fetched as a synchronous TAP query (/tap/sync) and cached atomically
as Parquet.  Completed chunks are skipped on rerun.  Failed batches are
retried with exponential backoff (30 s, 60 s, 120 s ± 10 % jitter).
"""

from __future__ import annotations

import logging
import math
import os
import random
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import polars as pl
from astroquery.utils.tap.core import TapPlus

logger = logging.getLogger(__name__)

_ORBITS_COLUMNS = [
    "number_mp",
    "denomination",
    "semi_major_axis",
    "eccentricity",
    "inclination",
    "arg_perihelion",
    "long_asc_node",
    "mean_anomaly",
    "osc_epoch",
]

_DEFAULT_MP_MAX = 160_000
_DEFAULT_BATCH_SIZE = 5_000
_RETRY_BASE_SECONDS = 30.0
_RAD_TO_DEG = 180.0 / math.pi


def _chunk_path(cache_dir: Path, mp_start: int, mp_end: int) -> Path:
    return cache_dir / f"orb_{mp_start:07d}_{mp_end:07d}.parquet"


def _fetch_chunk(
    archive_url: str,
    mp_start: int,
    mp_end: int,
    cache_path: Path,
    max_retries: int = 3,
) -> tuple[pl.DataFrame, float]:
    """Fetch orbital elements for number_mp in [mp_start, mp_end] via synchronous TAP.

    Retries with exponential backoff (±10 % jitter) on failure.  Writes
    result atomically to cache_path before returning.
    """
    col_list = ", ".join(_ORBITS_COLUMNS)
    adql = (
        f"SELECT {col_list} FROM gaiadr3.sso_orbits "
        f"WHERE number_mp BETWEEN {mp_start} AND {mp_end}"
    )
    label = f"orb {mp_start}–{mp_end}"

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


def _to_pipeline_schema(df: pl.DataFrame) -> pl.DataFrame:
    """Rename Gaia columns and convert angular elements from radians to degrees."""
    return (
        df.rename(
            {
                "number_mp": "number",
                "denomination": "designation",
                "semi_major_axis": "a_au",
                "eccentricity": "e",
                "osc_epoch": "epoch_jd",
            }
        )
        .with_columns(
            [
                (pl.col("inclination") * _RAD_TO_DEG).alias("i_deg"),
                (pl.col("arg_perihelion") * _RAD_TO_DEG).alias("omega_deg"),
                (pl.col("long_asc_node") * _RAD_TO_DEG).alias("Omega_deg"),
                (pl.col("mean_anomaly") * _RAD_TO_DEG).alias("M_deg"),
            ]
        )
        .drop(["inclination", "arg_perihelion", "long_asc_node", "mean_anomaly"])
        .with_columns([pl.col("number").cast(pl.Int32)])
    )


def _query_mp_max(archive_url: str) -> int:
    """Return MAX(number_mp) from the orbits table (fast single-row query)."""
    tap = TapPlus(url=archive_url)
    job = tap.launch_job("SELECT MAX(number_mp) AS mp_max FROM gaiadr3.sso_orbits")
    table = job.get_results()
    val = table["mp_max"][0]
    return int(val) if val is not None else _DEFAULT_MP_MAX


def download_gaia_orbits(
    archive_url: str,
    dest: str | Path,
    *,
    mp_max: int | None = None,
    batch_size: int = _DEFAULT_BATCH_SIZE,
    cache_dir: Path | None = None,
    n_workers: int | str = "auto",
    max_retries: int = 3,
) -> pl.DataFrame:
    """Download ``gaiadr3.sso_orbits`` via parallel TAP jobs.

    Parameters
    ----------
    archive_url:
        Base URL of the Gaia TAP server.
    dest:
        Output path for the merged Parquet file (pipeline schema).
    mp_max:
        Upper bound of ``number_mp`` range. If ``None``, queried first.
    batch_size:
        Number of ``number_mp`` values per TAP job.
    cache_dir:
        Directory for per-chunk Parquet files. Defaults to
        ``<dest.parent.parent>/cache/gaia_orbits_chunks``.
    n_workers:
        Number of parallel TAP connections. ``"auto"`` uses
        ``min(os.cpu_count(), 8)``.
    max_retries:
        Number of retry attempts per failed TAP job (exponential backoff).

    Returns
    -------
    polars.DataFrame
        Orbital elements in pipeline schema (number, designation, a_au, e,
        i_deg, Omega_deg, omega_deg, M_deg, epoch_jd), sorted by number.
    """
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)

    if cache_dir is None:
        cache_dir = dest.parent.parent / "cache" / "gaia_orbits_chunks"
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    if n_workers == "auto":
        n_workers = min(os.cpu_count() or 4, 8)

    if mp_max is None:
        logger.info("Querying MAX(number_mp) from %s…", archive_url)
        mp_max = _query_mp_max(archive_url)

    step = max(1, batch_size)
    ranges = [
        (start, min(start + step - 1, mp_max))
        for start in range(1, mp_max + 1, step)
    ]

    cached = [(s, e) for s, e in ranges if _chunk_path(cache_dir, s, e).exists()]
    pending = [(s, e) for s, e in ranges if not _chunk_path(cache_dir, s, e).exists()]

    logger.info(
        "gaiadr3.sso_orbits — %d cached | %d pending | %d workers | batch_size %d | mp 1–%d",
        len(cached),
        len(pending),
        n_workers,
        step,
        mp_max,
    )

    failed: list[tuple[int, int]] = []
    completed = len(cached)

    if pending:
        with ThreadPoolExecutor(max_workers=int(n_workers)) as pool:
            futures = {
                pool.submit(
                    _fetch_chunk,
                    archive_url,
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
                label = f"orb {s}–{e}"
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
            [f"orb {s}–{e}" for s, e in failed],
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
            "No chunk file for %d ranges — those orbits are absent.", len(missing_chunks)
        )

    if not frames:
        logger.warning("No rows returned from Gaia TAP — writing empty file")
        result = pl.DataFrame()
    else:
        result = _to_pipeline_schema(pl.concat(frames)).sort("number")

    result.write_parquet(dest, compression="zstd")
    logger.info("Saved %d orbital element sets to %s", len(result), dest)
    return result


def load_gaia_orbits(path: str | Path) -> pl.DataFrame:
    """Load a previously downloaded Gaia orbital elements Parquet file.

    Parameters
    ----------
    path:
        Path to the Parquet file written by :func:`download_gaia_orbits`.

    Returns
    -------
    polars.DataFrame
        Orbital elements in pipeline schema: ``number``, ``designation``,
        ``a_au``, ``e``, ``i_deg``, ``Omega_deg``, ``omega_deg``,
        ``M_deg``, ``epoch_jd`` — ready for :func:`detect_encounters`.

    Notes
    -----
    ``epoch_jd`` is in TDB (Barycentric Dynamical Time), consistent with
    the internal convention of the propagation and detection modules.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Gaia orbits file not found: {path}")
    return pl.read_parquet(path)
