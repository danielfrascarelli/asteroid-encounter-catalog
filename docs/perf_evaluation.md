# Performance evaluation — numba & persistent fine-grid cache (Track C1)

> Evaluation of the two optional perf follow-ups in `current_working_plan.md`
> Track C Stage 1. **Verdict: neither is warranted now.** The Kepler hot path is
> already numpy-vectorised and the pipeline is already process-parallel; the
> expensive coarse trajectory is already cached. Documented so the decision is
> not re-litigated without a new driver.
>
> Reproduce: `docker compose run --rm pipeline python -m scripts.bench.bench_kepler`

## 1. numba on `solve_kepler` / the refinement hot path — NOT warranted

`src/propagate/kepler.py::solve_kepler` solves Kepler's equation by
Newton-Raphson with a **fixed iteration loop over the whole array** — every
operation (`np.sin`, `np.cos`, the update) is a vectorised numpy call over all
elements at once. There is **no per-element Python loop**, which is the only
thing numba accelerates. `kepler_to_cartesian` is the same: pure broadcast numpy.

Measured throughput (single thread, 28-core host):

| function | throughput |
|---|---|
| `solve_kepler` (Newton, vectorised) | **~7.3 M elem/s** (flat from 1e5 to 5e6) |
| `kepler_to_cartesian` (full propagation) | **~3.7 M elem/s** |

Why numba would not help — and could hurt:

- The array path is already near-optimal numpy; numba's win is over scalar Python
  loops, which this code does not have. A numba scalar-Newton kernel would have to
  *beat* vectorised numpy, not just match it.
- The detection scan and the refinement **already parallelise across worker
  processes** (`src/detect/parallel.py`, `refine_candidates(n_workers=…)`), so the
  effective throughput is ~7 M elem/s × `n_workers` (~28× here). A numba
  `parallel=True` kernel would oversubscribe against that process pool.
- numba is **not** in the stack (`ModuleNotFoundError`); adding it means a new
  dependency + image rebuild + JIT warmup — CLAUDE.md discourages new deps without
  justification, and there is no measured win to justify it.

**Recommendation:** do not add numba. Revisit only if a future run profiles the
Kepler path as a dominant, non-parallelisable bottleneck (not the case today —
refinement is at its ~4.5–10× plateau).

## 2. Persistent fine-grid cache between runs — NOT warranted

The refinement re-propagates Kepler on a fine sub-grid per candidate each run. A
persistent on-disk cache of those fine-grid positions would only pay off for
**identical re-runs** (same pairs, same window, same epochs). In practice runs
differ (snapshot, threshold, body subset), so cache hits would be rare, while the
cache itself adds keying/invalidation complexity and disk.

Moreover the **expensive** part — the coarse N-body bulk trajectory — is
*already* cached between runs (`src/propagate/cache.py`, zarr/memmap, cache hit
< 1 s). The fine-grid Kepler refinement is cheap (≈ 7 M elem/s) by comparison, so
caching it has a poor cost/benefit.

**Recommendation:** do not add a persistent fine-grid cache. The existing coarse
trajectory cache already covers the costly re-computation.

## Summary

Both Track C1 items are **evaluated and declined** with measurements. The
pipeline's performance posture (vectorised numpy + process-level parallelism +
cached coarse trajectory) is already appropriate for the catalog scale; neither
optional optimisation clears the cost/benefit bar without a new, profiled driver.
