"""Generate null perturbers for the specificity test (Stage 3 of the mass layer).

For a given (target, real_perturber, encounter_date) triple, a *null* is an
asteroid that
  1. lives in the same orbital band as ``real_perturber``
     (``|a_null - a_real| <= a_window``), so any astrometric effect on
     ``target`` is geometrically plausible;
  2. does **not** have a close encounter with ``target`` within the Gaia DR3
     window. We use the precomputed encounter catalog as the source of truth
     for "no close approach" -- a pair that does not appear with
     ``dist_au < min_separation_au`` is treated as having no encounter.

The motivation: if the mass fit on a null returns the same mass distribution
as on the real perturber, the pipeline is not detecting a perturbation; it is
absorbing systematic noise.
"""

from __future__ import annotations

import logging

import numpy as np
import polars as pl

logger = logging.getLogger(__name__)


def asteroids_close_to_target(
    catalog: pl.DataFrame,
    target_number: int,
    min_separation_au: float = 0.1,
) -> set[int]:
    """Return MPCORB numbers that come within ``min_separation_au`` of the target.

    Uses the precomputed encounter catalog (any row with the target on either
    side of the pair counts).
    """
    close = catalog.filter(
        ((pl.col("number_1") == target_number) | (pl.col("number_2") == target_number))
        & (pl.col("dist_au") <= min_separation_au)
    )
    others = set()
    for row in close.iter_rows(named=True):
        n1, n2 = int(row["number_1"]), int(row["number_2"])
        if n1 == target_number:
            others.add(n2)
        else:
            others.add(n1)
    return others


def sample_null_perturbers(
    target_number: int,
    real_perturber_number: int,
    real_perturber_a_au: float,
    mpcorb: pl.DataFrame,
    encounters: pl.DataFrame,
    n_nulls: int,
    a_window_au: float = 0.5,
    min_separation_au: float = 0.1,
    h_window_mag: float | None = 1.5,
    seed: int = 42,
) -> list[int]:
    """Pick ``n_nulls`` MPCORB numbers eligible as nulls for the given target.

    Eligibility rules
    -----------------
    * Numbered asteroid (already enforced by the MPCORB loader caller).
    * ``|a - real_perturber_a| <= a_window_au``.
    * Not the real perturber and not the target itself.
    * Does not appear within ``min_separation_au`` of ``target`` in the
      precomputed encounter catalog.
    * Optionally within ``h_window_mag`` of the real perturber's H magnitude
      (size-matched, so the prior on the fit -- which initialises from H --
      starts in a comparable place).

    The function picks deterministically via ``seed`` so a re-run produces the
    same nulls. If fewer than ``n_nulls`` candidates remain after filters, all
    of them are returned.
    """
    rng = np.random.default_rng(seed)
    df = mpcorb
    if "number" not in df.columns:
        raise ValueError("mpcorb DataFrame must contain a 'number' column")
    df = df.filter(
        pl.col("number").is_not_null()
        & pl.col("a_au").is_not_null()
        & (pl.col("a_au") >= max(0.0, real_perturber_a_au - a_window_au))
        & (pl.col("a_au") <= real_perturber_a_au + a_window_au)
        & (pl.col("number") != real_perturber_number)
        & (pl.col("number") != target_number)
    )
    if h_window_mag is not None and "H" in df.columns:
        real_row = mpcorb.filter(pl.col("number") == real_perturber_number)
        if real_row.height == 1 and real_row["H"][0] is not None:
            h_real = float(real_row["H"][0])
            df = df.filter(
                pl.col("H").is_not_null()
                & (pl.col("H") >= h_real - h_window_mag)
                & (pl.col("H") <= h_real + h_window_mag)
            )

    blocked = asteroids_close_to_target(
        encounters, target_number, min_separation_au=min_separation_au
    )
    if blocked:
        df = df.filter(~pl.col("number").is_in(list(blocked)))

    available = df["number"].to_numpy()
    n_available = len(available)
    if n_available == 0:
        logger.warning(
            "No null perturbers eligible for target %d (real perturber %d).",
            target_number,
            real_perturber_number,
        )
        return []
    if n_available <= n_nulls:
        logger.info(
            "Only %d nulls eligible for target %d; returning all.",
            n_available,
            target_number,
        )
        return [int(x) for x in available]
    chosen = rng.choice(available, size=n_nulls, replace=False)
    return sorted(int(x) for x in chosen)


__all__ = [
    "asteroids_close_to_target",
    "sample_null_perturbers",
]
