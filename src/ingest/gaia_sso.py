"""Download and cache Gaia DR3 SSO observations via TAP.

The Gaia archive table ``gaiadr3.sso_observation`` contains one row per
asteroid transit. The ``epoch`` column is in TCB (Barycentric Coordinate Time)
as a Julian Date. All other processing converts to TDB via
``src.utils.time_utils.tcb_to_tdb`` before use.

Usage (via download script):
    docker compose run --rm pipeline python -m scripts.download_gaia_sso
"""

from __future__ import annotations

import logging
from pathlib import Path

import polars as pl
from astroquery.utils.tap.core import TapPlus

logger = logging.getLogger(__name__)

# Minimum columns required by downstream pipeline stages
_REQUIRED_COLUMNS = {
    "solution_id",
    "source_id",
    "denomination",
    "number_mp",
    "epoch",
}


def download_gaia_sso(
    archive_url: str,
    columns: list[str],
    dest: str | Path,
    *,
    chunk_size: int = 50_000,
) -> pl.DataFrame:
    """Download ``gaiadr3.sso_observation`` via TAP and save to Parquet.

    The download uses chunked ADQL queries (ORDER BY + OFFSET) to avoid
    TAP server row-count limits. The ``epoch`` column is preserved as-is
    in TCB; conversion to TDB happens at the propagation stage.

    Parameters
    ----------
    archive_url:
        Base URL of the Gaia TAP server.
    columns:
        List of column names to retrieve. Must include all of
        ``_REQUIRED_COLUMNS``.
    dest:
        Output path for the Parquet file.
    chunk_size:
        Rows per TAP query (default 50 000, well within the 2M row limit).

    Returns
    -------
    polars.DataFrame
        All downloaded observations.

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

    col_list = ", ".join(columns)
    tap = TapPlus(url=archive_url)

    frames: list[pl.DataFrame] = []
    offset = 0

    logger.info("Downloading gaiadr3.sso_observation from %s", archive_url)

    while True:
        adql = (
            f"SELECT TOP {chunk_size} {col_list} "
            f"FROM gaiadr3.sso_observation "
            f"ORDER BY source_id "
            f"OFFSET {offset}"
        )
        logger.debug("TAP query: %s", adql)
        job = tap.launch_job(adql)
        table = job.get_results()

        if len(table) == 0:
            break

        df = pl.from_pandas(table.to_pandas())
        frames.append(df)
        logger.info("  fetched %d rows (total so far: %d)", len(df),
                    sum(len(f) for f in frames))

        if len(table) < chunk_size:
            break

        offset += chunk_size

    if not frames:
        logger.warning("No rows returned from Gaia TAP — writing empty file")
        result = pl.DataFrame()
    else:
        result = pl.concat(frames)

    # Annotate epoch scale in metadata via schema comment (Parquet metadata)
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
