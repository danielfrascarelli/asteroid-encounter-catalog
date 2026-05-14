"""KD-tree scan — identifies candidate close-encounter epochs.

For each step in *time_grid*:

1. Propagate all asteroids to that time.
2. Build a ``cKDTree`` over the (N, 3) position array.
3. Query all pairs within *threshold_au*.
4. Intersect with the set of compatible pairs from the prefilter.
5. Track the best (minimum) distance seen per surviving pair.

Returns one ``(idx_i, idx_j, best_t_jd, best_dist_au)`` tuple per pair
that had at least one epoch within *threshold_au*.
"""

from __future__ import annotations

import logging

import numpy as np
import polars as pl
from scipy.spatial import cKDTree
from tqdm import tqdm

from src.propagate.grid import propagate_grid

logger = logging.getLogger(__name__)


def scan_time_grid(
    elements: pl.DataFrame,
    time_grid: np.ndarray,
    pairs: np.ndarray | None,
    threshold_au: float,
    leaf_size: int = 30,
) -> list[tuple[int, int, float, float]]:
    """Scan *time_grid* for compatible pairs closer than *threshold_au*.

    Parameters
    ----------
    elements:
        Orbital elements DataFrame (columns required by
        :func:`src.propagate.kepler.propagate_df`).
    time_grid:
        Array of JD TDB values to evaluate.
    pairs:
        ``(M, 2)`` int32 array of prefilter-compatible pair indices.
    threshold_au:
        Distance threshold in AU.
    leaf_size:
        ``cKDTree`` leaf_size parameter.

    Returns
    -------
    list of (idx_i, idx_j, best_t_jd, best_dist_au)
        One entry per pair where the minimum observed distance is below
        *threshold_au*.  Pairs never seen within threshold are omitted.
    """
    if pairs is not None and len(pairs) == 0:
        return []

    # Build a fast set for O(1) membership tests (query_pairs guarantees i<j).
    # When pairs is None the orbital prefilter was skipped; all spatially close
    # pairs returned by query_pairs are accepted.
    compatible_set: set[tuple[int, int]] | None = (
        None if pairs is None else {(int(p[0]), int(p[1])) for p in pairs}
    )

    best: dict[tuple[int, int], tuple[float, float]] = {}

    for step_idx, (t_jd, pos) in enumerate(
        tqdm(propagate_grid(elements, time_grid), total=len(time_grid), desc="KD-tree scan", unit="step", leave=False)
    ):
        tree = cKDTree(pos, leafsize=leaf_size)
        raw: np.ndarray = tree.query_pairs(threshold_au, output_type="ndarray")

        if len(raw) == 0:
            continue

        # Vectorised distance computation for all pairs within threshold
        diffs = pos[raw[:, 0]] - pos[raw[:, 1]]
        dists = np.linalg.norm(diffs, axis=1)

        for (a_idx, b_idx), d in zip(raw, dists):
            key = (int(a_idx), int(b_idx))  # query_pairs guarantees i < j
            if compatible_set is not None and key not in compatible_set:
                continue
            d_f = float(d)
            if key not in best or d_f < best[key][1]:
                best[key] = (t_jd, d_f)

        if (step_idx + 1) % 100 == 0:
            logger.debug(
                "Scanned %d/%d steps; %d candidate pairs so far",
                step_idx + 1,
                len(time_grid),
                len(best),
            )

    logger.info(
        "KD-tree scan complete: %d candidate pairs in %d steps",
        len(best),
        len(time_grid),
    )
    return [(k[0], k[1], v[0], v[1]) for k, v in best.items()]
