# Kepler-refine vs N-body refinement — error characterisation

**Track 1 Stage A** of [DEEPWORK_PLAN.md](../DEEPWORK_PLAN.md). Quantifies the
error introduced by the Kepler-2-body sub-grid refinement that produces the
final `dist_au`, `t_min`, `rel_vel_au_day` columns of the frozen catalog
([FROZEN_RUN.md](../FROZEN_RUN.md)).

## TL;DR

Over a stratified sample of **796 pairs** drawn from the 72 M-row catalog,
each refined under full N-body (rebound: WHFast warmup + IAS15 across a ±12 h
window with 60 s sampling, Sun + Jupiter + Saturn + Ceres/Pallas/Vesta/Hygiea
as massive bodies):

| Metric | Median | p95 | p99 | Max |
|---|---|---|---|---|
| `|Δdist_au|` | **1.6e-5** | 6.5e-4 | 2.2e-3 | 5.6e-3 |
| `|Δt_min|` (h) | 1.6 | 8.0 | 12.0¹ | 12.0¹ |
| `|Δrel_vel|` (AU/d) | 1.5e-7 | 4.0e-6 | 1.1e-5 | 2.6e-5 |

¹ Clipped at the ±12 h integration window for 26 / 796 pairs (3.3 %). For
these pairs the true N-body minimum lies outside the window, so the reported
`dist_au_nbody` is the value at the window boundary and is an *upper bound*
on the true Kepler-vs-N-body disagreement.

**Bottom line**: for 99 % of the catalog the Kepler-refined `dist_au` agrees
with the full N-body minimum to within **~2 mAU**, well below the 0.05 AU
encounter threshold. The refinement does not flip detection status for any
pair in the sample. It does, however, shift the ranking of *which* pair is
closest within ~mAU bands, and is materially relevant for any downstream task
that consumes sub-mAU geometry (e.g. mass determinations).

## Setup

- **Sample**: `data/cache/nbody_validation/sample_1000.parquet` — 796 pairs
  stratified across 5 axes (`a_1`, `e_1`, `i_1`, `dist_au`, `|Δa|`),
  intentionally over-sampling the high-e / high-i / small-q tails where
  Kepler ↔ N-body disagreement is expected to be largest. Generator:
  [scripts/validate/sample_for_nbody_check.py](../scripts/validate/sample_for_nbody_check.py).
- **N-body refiner**: per-pair simulation in REBOUND.
  - Massive bodies: Sun + Jupiter + Saturn + Ceres + Pallas + Vesta + Hygiea.
  - WHFast warmup at dt = 600 s from MPCORB epoch (2016-02-17) to the window
    start.
  - IAS15 across ±12 h around the Kepler-refined `t_min`, sampling
    `dist(t)` every 60 s; parabolic fit at the discrete minimum to recover
    sub-step precision.
  - Heliocentric ecliptic J2000 frame, consistent with the rest of the
    pipeline.
  - Code: [scripts/validate/refine_pair_nbody.py](../scripts/validate/refine_pair_nbody.py),
    tests in [tests/test_refine_pair_nbody.py](../tests/test_refine_pair_nbody.py).
- **Comparator**: parallel pool over the 796 pairs.
  - 24 workers × spawn context (REBOUND-safe), chunksize 4.
  - Wall-clock: **~9 s** for the full sample on a 28-core machine.
  - Energy drift: max **3.4e-14** — symplectic warmup + IAS15 well within
    spec.
  - Output: [`data/output/kepler_vs_nbody_comparison.parquet`](../data/output/kepler_vs_nbody_comparison.parquet).
  - Code: [scripts/validate/compare_kepler_vs_nbody.py](../scripts/validate/compare_kepler_vs_nbody.py).

## Stratified results

### Error grows with eccentricity

| `e_max` bin | n | med `|Δdist|` | p95 | max |
|---|---|---|---|---|
| < 0.10 | 78 | 6e-6 | 2.9e-4 | 9.0e-4 |
| 0.10 – 0.20 | 266 | 1.3e-5 | 4.6e-4 | 2.3e-3 |
| 0.20 – 0.30 | 265 | 1.6e-5 | 5.8e-4 | 3.6e-3 |
| 0.30 – 0.45 | 117 | 2.8e-5 | 1.4e-3 | 3.1e-3 |
| 0.45 – 0.70 | 70 | **7.0e-5** | **1.3e-3** | **5.6e-3** |

A factor ~12× rise in the median as `e_max` goes from quiescent to
high-eccentricity. Expected: higher e ⇒ stronger non-Keplerian acceleration
near perihelion ⇒ Kepler 2-body diverges from N-body faster.

### Error grows as perihelion drops

| `q_min` bin (AU) | n | med `|Δdist|` | p95 | max |
|---|---|---|---|---|
| < 1.3 (NEA-like) | 72 | **6.5e-5** | 1.1e-3 | **5.6e-3** |
| 1.3 – 1.8 | 134 | 2.8e-5 | 7.6e-4 | 2.6e-3 |
| 1.8 – 2.2 | 333 | 1.8e-5 | 5.8e-4 | 3.6e-3 |
| 2.2 – 2.6 | 178 | 1.0e-5 | 6.5e-4 | 2.3e-3 |
| 2.6 – 3.0 | 75 | 7e-6 | 2.9e-4 | 9.0e-4 |
| > 3.0 | 4 | 2e-6 | 3.1e-5 | 3.1e-5 |

Same story from the other side: pairs with at least one Earth-crosser-like
member have the largest Kepler/N-body disagreements.

### Inclination has no clean trend

Median `|Δdist|` is essentially flat in `i_max` across 0 – 25° (1e-5 to 2e-5).
The lone bin at i > 25° has only 9 samples, so the apparent rise there is
noise. Take-away: the dynamical sensitivity here is e and q, not i.

## Subset analysis

| Subset | n | Frac | med `|Δdist|` | p95 | p99 | max |
|---|---|---|---|---|---|---|
| **Full sample** | 796 | 100 % | 1.6e-5 | 6.5e-4 | 2.2e-3 | 5.6e-3 |
| **Main-belt** (q∈[1.8,3.0] AU, e<0.3) | 526 | 66 % | 1.3e-5 | 5.7e-4 | 1.4e-3 | 3.6e-3 |
| **NEA-like** (q<1.3 AU) | 72 | 9 % | 6.5e-5 | 1.1e-3 | — | 5.6e-3 |

The catalog is overwhelmingly main-belt (≈66 % of the *stratified* sample;
the unstratified catalog skews even more heavily main-belt because the
high-e and low-q corners we deliberately over-sampled are relatively rare).
For the typical main-belt pair, Kepler-refine agrees with N-body to better
than 1.4 mAU at p99.

## How many pairs exceed each threshold

| Threshold | Pairs over | % |
|---|---|---|
| `|Δdist|` > 1×10⁻⁴ AU | 190 / 796 | 23.9 % |
| `|Δdist|` > 5×10⁻⁴ AU | 54 / 796 | 6.8 % |
| `|Δdist|` > 1×10⁻³ AU | 25 / 796 | 3.1 % |
| `|Δdist|` > 2×10⁻³ AU | 9 / 796 | 1.1 % |
| `|Δdist|` > 5×10⁻³ AU | 1 / 796 | 0.1 % |

For comparison, the catalog encounter threshold is 0.05 AU = 5×10⁻² AU.
None of the 796 pairs in the sample flips status (close vs not-close) when
re-refined under N-body. But the **ranking** of the closest 1 % of pairs in
the catalog could shift on the order of mAU when re-refined.

## Recommendation for Stage B

The plan ([DEEPWORK_PLAN.md](../DEEPWORK_PLAN.md) § Stage B) defines a
selective N-body refinement on a subset where Kepler is "not defensible".
Based on the measured error distribution, **Stage B is warranted but
narrowly scoped**:

- **Defensible Kepler-refine** (no Stage B needed): `q_min ≥ 1.8 AU` and
  `e_max ≤ 0.3`. 66 % of the catalog. p99 error 1.4 mAU — well below any
  scientifically interesting threshold for *detection*.
- **Refine under N-body** (Stage B target): `q_min < 1.8 AU` **or**
  `e_max > 0.3`. ~34 % of pairs in the stratified sample, but likely
  smaller in the full catalog because we over-sampled the tail. A
  back-of-the-envelope using the unstratified fractions (e>0.3: ~10 %,
  q<1.8: ~15 % with overlap) suggests **~20 % of catalog pairs**
  (~14 M out of 72 M). At ~9 s / 796 pairs measured here ⇒ extrapolated
  cost ~44 h CPU on 24 cores. Feasible.

The harder claim (Stage B should chase **mAU precision** for mass-determination
candidates) is essentially independent of the headline characterisation: any
mass-fitting subset can be re-refined per-pair under N-body at negligible
marginal cost (~10 ms / pair).

## What this report does NOT say

1. **It does not say the catalog is wrong.** Detection (is dist < 0.05 AU?)
   is unaffected; only sub-mAU geometry is.
2. **It does not validate against truth (JPL Horizons).** The N-body
   refiner is checked against Horizons in one canonical pair
   (test_refine_pair_matches_jpl_horizons, currently `@pytest.mark.skip` for
   network reasons) — a broader cross-check should be done before publishing.
3. **It does not address the mass-layer problems** (Track 2 of the deepwork
   plan). Those are orthogonal: even with Stage B applied, the χ²_red ≈ 425
   issue and the 0/41 specificity result remain open.

## Artifacts

- Sample: `data/cache/nbody_validation/sample_1000.parquet` (796 rows)
- Per-pair comparison: `data/output/kepler_vs_nbody_comparison.parquet`
  (796 rows × 23 cols)
- Notebook: [notebooks/nbody_error_characterization.ipynb](../notebooks/nbody_error_characterization.ipynb)

## Reproducibility

```bash
# Re-sample (deterministic with --seed 42)
docker compose run --rm pipeline python -m scripts.validate.sample_for_nbody_check

# Re-run the comparator (deterministic given the same MPCORB snapshot)
docker compose run --rm pipeline python -m scripts.validate.compare_kepler_vs_nbody \\
    --workers 24 --window-hours 12.0

# Run tests
docker compose run --rm test pytest tests/test_refine_pair_nbody.py -v
```
