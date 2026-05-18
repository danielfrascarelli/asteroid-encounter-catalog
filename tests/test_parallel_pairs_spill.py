"""Test that the parallel scan spills large ``pairs`` arrays to a tempfile.

When prefilter survives many pairs (e.g. 2000 asteroids → ~1.3 M pairs ≈ 20 MB)
the previous code path passed the array through ``Pool`` initargs, which
pickled it once per worker and could deadlock the forkserver. The fix in
``src/detect/parallel.py`` spills arrays >1 MB to a temp directory and lets
each worker mmap them.
"""

from __future__ import annotations

import numpy as np
import polars as pl

from src.detect.parallel import scan_parallel


def _make_elements(n: int) -> pl.DataFrame:
    rng = np.random.default_rng(seed=7)
    return pl.DataFrame(
        {
            "number": np.arange(1, n + 1, dtype=np.int32),
            "designation": [f"({i+1})" for i in range(n)],
            "a_au": 2.7 + rng.normal(0.0, 0.001, size=n),
            "e": np.zeros(n),
            "i_deg": np.zeros(n),
            "Omega_deg": np.zeros(n),
            "omega_deg": np.zeros(n),
            "M_deg": rng.uniform(0.0, 360.0, size=n),
            "epoch_jd": np.full(n, 2457200.5),
        },
        schema_overrides={"number": pl.Int32, "designation": pl.Utf8},
    )


def test_scan_parallel_with_spilled_pairs() -> None:
    """A pairs array >1 MB must trigger the tempfile spill path and complete."""
    n_ast = 50
    elements = _make_elements(n_ast)
    time_grid = np.linspace(2457200.5, 2457220.5, 21)  # 20 days, daily

    # Synthesise a pair list larger than the 1 MB spill threshold by repeating
    # the (0, 1) pair many times — the scan tolerates duplicate pairs.
    n_pairs = 200_000  # 200k pairs × 2 int64 = 3.2 MB → spill triggers
    pairs = np.zeros((n_pairs, 2), dtype=np.int64)
    pairs[:, 1] = 1
    assert pairs.nbytes > 1_000_000

    # Two workers exercise the worker-init path.
    result = scan_parallel(
        elements,
        time_grid,
        pairs=pairs,
        threshold_au=0.5,
        leaf_size=10,
        n_workers=2,
        chunk_size_days=5.0,
        positions=None,
    )
    # We don't assert specific contents — just that the call returned without
    # hanging on the pickle-via-initargs deadlock the spill path is meant to
    # avoid.
    assert isinstance(result, list)
