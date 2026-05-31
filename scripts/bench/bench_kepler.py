"""Benchmark the Kepler solver / propagator throughput (Track C1 evaluation).

Backs the numba decision in ``docs/perf_evaluation.md``: ``solve_kepler`` and
``kepler_to_cartesian`` are already numpy-vectorised (the only loop is the fixed
Newton-iteration count over the whole array), so they run at several M
elements/s single-threaded — and the detection/refinement pipeline already
parallelises across worker *processes*, multiplying that by ``n_workers``.

Usage:
    docker compose run --rm pipeline python -m scripts.bench.bench_kepler
"""

from __future__ import annotations

import os
import time

import numpy as np

from src.propagate.kepler import kepler_to_cartesian, solve_kepler


def main() -> int:
    rng = np.random.default_rng(42)
    print(f"cpu_count = {os.cpu_count()}")
    print("--- solve_kepler (Newton, vectorised) ---")
    for n in (100_000, 1_000_000, 5_000_000):
        m = rng.uniform(0, 2 * np.pi, n)
        e = rng.uniform(0, 0.4, n)
        solve_kepler(m[:1000], e[:1000])  # warm
        t = time.perf_counter()
        solve_kepler(m, e)
        dt = time.perf_counter() - t
        print(f"  N={n:>9,}: {dt * 1000:8.1f} ms  ->  {n / dt / 1e6:6.1f} M elem/s")

    print("--- kepler_to_cartesian (full propagation, the refinement hot path) ---")
    n = 1_000_000
    a = rng.uniform(1.5, 4, n)
    e = rng.uniform(0, 0.4, n)
    i = rng.uniform(0, 0.5, n)
    om_big = rng.uniform(0, 2 * np.pi, n)
    om_small = rng.uniform(0, 2 * np.pi, n)
    m0 = rng.uniform(0, 2 * np.pi, n)
    epoch = np.full(n, 2457400.5)
    t_jd = np.full(n, 2457000.5)
    t = time.perf_counter()
    kepler_to_cartesian(a, e, i, om_big, om_small, m0, epoch, t_jd)
    dt = time.perf_counter() - t
    print(f"  N={n:,}: {dt * 1000:.1f} ms  ->  {n / dt / 1e6:.1f} M elem/s single-thread")
    print(
        f"  × n_workers (process-parallel scan/refine) ≈ {n / dt / 1e6 * (os.cpu_count() or 1):.0f} M elem/s aggregate"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
