"""Orbital prefilter — heuristic pair reduction before the KD-tree scan.

Two asteroids are dropped from the candidate list if their semimajor axes
differ by more than *semimajor_diff_max_au* OR their inclinations differ by
more than *inclination_diff_max_deg*.  These are cheap heuristic criteria —
they are **not** a proof of geometric impossibility: high-eccentricity or
high-inclination crossing orbits can come within the encounter threshold
even when Δa or Δi is large, so a real encounter can be missed by this
filter.

Empirically the criteria eliminate the vast majority of N*(N-1)/2 pairs
without affecting recall on the MBA population that dominates the catalog,
but completeness on the high-e tail has not been quantified (audit blocker
#2).  For N > 5000 the caller in :mod:`src.detect.pipeline` skips the
prefilter entirely and relies on the cKDTree spatial query at the
configured threshold radius — which IS exact under the propagation model.
"""

from __future__ import annotations

import logging

import numpy as np
import polars as pl

logger = logging.getLogger(__name__)


def compatible_pairs(
    elements: pl.DataFrame,
    semimajor_diff_max_au: float = 0.5,
    inclination_diff_max_deg: float = 30.0,
) -> np.ndarray:
    """Return row-index pairs that pass the orbital prefilter.

    Parameters
    ----------
    elements:
        DataFrame with at least columns ``a_au`` (Float64) and ``i_deg``
        (Float64).  Row order determines the indices in the output.
    semimajor_diff_max_au:
        Maximum allowed |a₁ - a₂| in AU.
    inclination_diff_max_deg:
        Maximum allowed |i₁ - i₂| in degrees.

    Returns
    -------
    np.ndarray of shape (M, 2), dtype int32
        Each row is a pair ``(idx_i, idx_j)`` with ``idx_i < idx_j`` that
        passed both criteria.
    """
    a = elements["a_au"].to_numpy()
    i = elements["i_deg"].to_numpy()
    n = len(a)

    if n < 2:
        return np.empty((0, 2), dtype=np.int32)

    ii, jj = np.triu_indices(n, k=1)

    da = np.abs(a[ii] - a[jj])
    di = np.abs(i[ii] - i[jj])

    mask = (da <= semimajor_diff_max_au) & (di <= inclination_diff_max_deg)
    pairs = np.stack([ii[mask], jj[mask]], axis=1).astype(np.int32)

    total = n * (n - 1) // 2
    logger.info(
        "Prefilter: %d / %d pairs survive (%.1f%%)",
        len(pairs),
        total,
        100.0 * len(pairs) / total if total > 0 else 0.0,
    )
    return pairs
