"""Query API for the close-encounter catalog."""

from __future__ import annotations

from pathlib import Path

import polars as pl


def load_catalog(path: str | Path) -> pl.DataFrame:
    """Load the encounter catalog from a parquet file."""
    return pl.read_parquet(path)


def filter_encounters(
    df: pl.DataFrame,
    min_dist_au: float | None = None,
    max_dist_au: float | None = None,
    date_start: str | None = None,
    date_end: str | None = None,
    body_ids: list[int] | None = None,
    gaia_observable_only: bool = False,
) -> pl.DataFrame:
    """Filter the encounter catalog by various criteria.

    Parameters
    ----------
    df:
        Catalog DataFrame (from :func:`load_catalog`).
    min_dist_au:
        Keep only encounters with ``dist_au >= min_dist_au``.
    max_dist_au:
        Keep only encounters with ``dist_au <= max_dist_au``.
    date_start:
        Keep only encounters on or after this ISO date string (``"yyyy-mm-dd"``).
    date_end:
        Keep only encounters on or before this ISO date string.
    body_ids:
        Keep only encounters involving at least one of these MPC numbers.
    gaia_observable_only:
        If True, keep only Gaia-observable encounters.

    Returns
    -------
    pl.DataFrame
        Filtered catalog.
    """
    mask = pl.lit(True)
    if min_dist_au is not None:
        mask = mask & (pl.col("dist_au") >= min_dist_au)
    if max_dist_au is not None:
        mask = mask & (pl.col("dist_au") <= max_dist_au)
    if date_start is not None:
        mask = mask & (pl.col("date_utc") >= date_start)
    if date_end is not None:
        mask = mask & (pl.col("date_utc") <= date_end)
    if body_ids is not None:
        mask = mask & (pl.col("number_1").is_in(body_ids) | pl.col("number_2").is_in(body_ids))
    if gaia_observable_only:
        mask = mask & pl.col("gaia_observable")
    return df.filter(mask)


def top_encounters(
    df: pl.DataFrame,
    n: int = 10,
    by: str = "dist_au",
    ascending: bool = True,
) -> pl.DataFrame:
    """Return the top-N encounters sorted by a given column.

    Parameters
    ----------
    df:
        Catalog DataFrame.
    n:
        Number of rows to return.
    by:
        Column to sort by.  Defaults to ``"dist_au"`` (closest first).
    ascending:
        Sort direction.  ``True`` = smallest value first.

    Returns
    -------
    pl.DataFrame
        Top-N rows sorted by ``by``.
    """
    return df.sort(by, descending=not ascending).head(n)
