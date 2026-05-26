# Kepler-refine vs N-body refinement — error characterisation

**Track 1 Stage A** of [DEEPWORK_PLAN.md](../DEEPWORK_PLAN.md). Quantifies the
error introduced by the Kepler-2-body sub-grid refinement that produces the
final `dist_au`, `t_min`, `rel_vel_au_day` columns of the frozen catalog
([FROZEN_RUN.md](../FROZEN_RUN.md)).

## TL;DR

Over a stratified sample of **964 pairs** drawn from the 72 M-row catalog,
each refined under full N-body (rebound: WHFast warmup + IAS15 across a ±12 h
window with 60 s sampling, Sun + Jupiter + Saturn + Ceres/Pallas/Vesta/Hygiea
as massive bodies):

| Metric | Median | p95 | p99 | Max |
|---|---|---|---|---|
| `|Δdist_au|` | **1.2e-5** | 6.8e-4 | 2.5e-3 | 1.1e-2 |
| `|Δt_min|` (h) | 1.7 | 9.4 | 12.0¹ | 12.0¹ |
| `|Δrel_vel|` (AU/d) | 1.3e-7 | 5.1e-6 | 1.2e-5 | 3.8e-5 |

¹ Clipped at the ±12 h integration window for 33 / 964 pairs (3.4 %). For
these pairs the true N-body minimum lies outside the window, so the reported
`dist_au_nbody` is the value at the window boundary and is an *upper bound*
on the true Kepler-vs-N-body disagreement.

**Bottom line**: for 99 % of the catalog the Kepler-refined `dist_au` agrees
with the full N-body minimum to within **~2.5 mAU**, well below the 0.05 AU
encounter threshold. The refinement does not flip detection status for any
pair in the sample. It does, however, shift the ranking of *which* pair is
closest within ~mAU bands, and is materially relevant for any downstream task
that consumes sub-mAU geometry (e.g. mass determinations).

The largest single disagreement we observed was **11.3 mAU** for a pair with
`e_max ≈ 0.49, q_min ≈ 1.36, i_max ≈ 9°` — a high-eccentricity Mars-grazing
pair whose Kepler trajectory clearly diverges from N-body across the ±12 h
encounter window.

## Setup

- **Sample**: `data/cache/nbody_validation/sample_1000.parquet` — 964 pairs
  stratified across 5 axes (`a_mid`, `e_max`, `i_max`, `q_min`, `dist_au`).
  Stratification variables are **symmetric** in body 1 / body 2, matching
  the way the analysis bins the error: previous versions of this sampler
  used `e_1`/`i_1`/`q_1`, which under-sampled pairs where body 2 was the
  high-e member. Generator:
  [scripts/validate/sample_for_nbody_check.py](../scripts/validate/sample_for_nbody_check.py).
  All 200 occupied bins satisfy the plan's `≥3 candidates per bin`
  acceptance criterion (enforced by the `--min-per-bin` flag).
- **N-body refiner**: per-pair simulation in REBOUND.
  - Massive bodies: Sun + Jupiter + Saturn + Ceres + Pallas + Vesta + Hygiea.
    If either target asteroid is itself one of Ceres/Pallas/Vesta/Hygiea, the
    refiner excludes that body from the perturber list to avoid the
    double-counting bug fixed in
    `tests/test_refine_pair_nbody.py::test_refine_pair_target_is_major_asteroid_no_duplication`.
  - WHFast warmup at dt = 600 s from MPCORB epoch (2016-02-17) to the window
    start.
  - IAS15 across ±12 h around the Kepler-refined `t_min`, sampling
    `dist(t)` every 60 s; parabolic fit at the discrete minimum to recover
    sub-step precision.
  - Heliocentric ecliptic J2000 frame, consistent with the rest of the
    pipeline.
  - Code: [scripts/validate/refine_pair_nbody.py](../scripts/validate/refine_pair_nbody.py),
    tests in [tests/test_refine_pair_nbody.py](../tests/test_refine_pair_nbody.py).
- **Comparator**: parallel pool over the 964 pairs.
  - 24 workers × spawn context (REBOUND-safe), chunksize 4.
  - Wall-clock: **~7 s** for the full sample on a 28-core machine.
  - Energy drift: max **3.5e-14** — symplectic warmup + IAS15 well within
    spec.
  - Output: [`data/output/kepler_vs_nbody_comparison.parquet`](../data/output/kepler_vs_nbody_comparison.parquet).
  - Code: [scripts/validate/compare_kepler_vs_nbody.py](../scripts/validate/compare_kepler_vs_nbody.py).

## Stratified results

### Error grows with eccentricity

| `e_max` bin | n | med `|Δdist|` | p95 | max |
|---|---|---|---|---|
| < 0.10 | 208 | 1.0e-5 | 2.9e-4 | 7.8e-4 |
| 0.10 – 0.20 | 257 | 0.7e-5 | 3.2e-4 | 1.3e-3 |
| 0.20 – 0.30 | 229 | 1.4e-5 | 8.0e-4 | 3.6e-3 |
| 0.30 – 0.45 | 186 | 3.0e-5 | 1.0e-3 | 5.6e-3 |
| 0.45 – 0.70 | 84 | **5.5e-5** | **2.8e-3** | **1.1e-2** |

A factor ~6× rise in the median as `e_max` goes from quiescent to
high-eccentricity, and a factor ~10× rise in p95. Expected: higher e ⇒
stronger non-Keplerian acceleration near perihelion ⇒ Kepler 2-body
diverges from N-body faster.

### Error grows as perihelion drops

| `q_min` bin (AU) | n | med `|Δdist|` | p95 | max |
|---|---|---|---|---|
| < 1.3 (NEA-like) | 106 | **5.5e-5** | 7.5e-4 | 2.9e-3 |
| 1.3 – 1.8 | 183 | 1.9e-5 | 8.9e-4 | **1.1e-2** |
| 1.8 – 2.2 | 265 | 1.4e-5 | 7.8e-4 | 5.6e-3 |
| 2.2 – 2.6 | 276 | 0.9e-5 | 3.6e-4 | 2.5e-3 |
| 2.6 – 3.0 | 134 | 0.6e-5 | 5.1e-4 | 2.4e-3 |

Same story from the other side: pairs with at least one Earth/Mars-crosser-
like member have the largest Kepler/N-body disagreements. The single worst
pair lives in the `1.3-1.8 AU` bin — high-eccentricity Mars-grazing.

### Inclination has no clean trend in the median

Median `|Δdist|` is flat in `i_max` across 0 – 15° (≈1e-5) and rises only
mildly above 15° (to ~2e-5). The p95 wanders without monotone trend. The
dynamical sensitivity here is e and q, not i.

## Subset analysis

| Subset | n | Frac | med `|Δdist|` | p95 | p99 | max |
|---|---|---|---|---|---|---|
| **Full sample** | 964 | 100 % | 1.2e-5 | 6.8e-4 | 2.5e-3 | 1.1e-2 |
| **Main-belt** (q∈[1.8,3.0] AU, e<0.3) | 576 | 60 % | 0.9e-5 | 4.2e-4 | 1.3e-3 | 3.6e-3 |
| **NEA-like** (q<1.3 AU) | 106 | 11 % | 5.5e-5 | 7.2e-4 | — | 2.9e-3 |

The catalog is overwhelmingly main-belt (60 % of the *stratified* sample;
the unstratified catalog skews even more heavily main-belt because the
high-e and low-q corners we deliberately over-sampled are relatively rare).
For the typical main-belt pair, Kepler-refine agrees with N-body to better
than 1.3 mAU at p99.

## How many pairs exceed each threshold

| Threshold | Pairs over | % |
|---|---|---|
| `|Δdist|` > 1×10⁻⁴ AU | 213 / 964 | 22.1 % |
| `|Δdist|` > 5×10⁻⁴ AU | 69 / 964 | 7.2 % |
| `|Δdist|` > 1×10⁻³ AU | 27 / 964 | 2.8 % |
| `|Δdist|` > 2×10⁻³ AU | 11 / 964 | 1.1 % |
| `|Δdist|` > 5×10⁻³ AU | 2 / 964 | 0.2 % |
| `|Δdist|` > 1×10⁻² AU | 1 / 964 | 0.1 % |

For comparison, the catalog encounter threshold is 0.05 AU = 5×10⁻² AU.
None of the 964 pairs in the sample flips status (close vs not-close) when
re-refined under N-body. But the **ranking** of the closest 1 % of pairs in
the catalog could shift on the order of mAU when re-refined.

## Recommendation for Stage B

The plan ([DEEPWORK_PLAN.md](../DEEPWORK_PLAN.md) § Stage B) defines a
selective N-body refinement on a subset where Kepler is "not defensible".
Based on the measured error distribution, **Stage B is warranted but
narrowly scoped**:

- **Defensible Kepler-refine** (no Stage B needed): `q_min ≥ 1.8 AU` and
  `e_max ≤ 0.3`. 60 % of the catalog. p99 error 1.3 mAU — well below any
  scientifically interesting threshold for *detection*.
- **Refine under N-body** (Stage B target): `q_min < 1.8 AU` **or**
  `e_max > 0.3`. ~40 % of pairs in the stratified sample, but probably
  smaller in the full catalog because we over-sampled the tail. A
  back-of-the-envelope using the unstratified fractions (e>0.3: ~10 %,
  q<1.8: ~15 % with overlap) suggests **~20 % of catalog pairs**
  (~14 M out of 72 M). At ~7 ms / pair measured here ⇒ extrapolated
  cost ~27 h CPU on 24 cores. Feasible.

The harder claim (Stage B should chase **mAU precision** for mass-determination
candidates) is essentially independent of the headline characterisation: any
mass-fitting subset can be re-refined per-pair under N-body at negligible
marginal cost (~10 ms / pair).

## What this report does NOT say

1. **It does not say the catalog is wrong.** Detection (is dist < 0.05 AU?)
   is unaffected; only sub-mAU geometry is.
2. **It does not validate against truth (JPL Horizons) at scale.** The
   N-body refiner is cross-checked against Horizons on a small fixture of
   reference pairs (`tests/fixtures/jpl_horizons_pairs.json`,
   tests/test_refine_pair_nbody.py marker `horizons`). A broader 1000-pair
   cross-check is future work — current evidence is consistent with the
   refiner being accurate but the breadth of validation is limited.
3. **It does not address the mass-layer problems** (Track 2 of the deepwork
   plan). Those are orthogonal: even with Stage B applied, the χ²_red ≈ 425
   issue and the 0/41 specificity result remain open.

## Artifacts

- Sample: `data/cache/nbody_validation/sample_1000.parquet` (964 rows)
- Per-pair comparison: `data/output/kepler_vs_nbody_comparison.parquet`
  (964 rows × 24 cols)
- Notebook: [notebooks/nbody_error_characterization.ipynb](../notebooks/nbody_error_characterization.ipynb)

## Reproducibility

```bash
# Re-sample (deterministic with --seed 42)
docker compose run --rm pipeline python -m scripts.validate.sample_for_nbody_check

# Re-run the comparator (deterministic given the same MPCORB snapshot)
docker compose run --rm pipeline python -m scripts.validate.compare_kepler_vs_nbody \\
    --workers 24 --window-hours 12.0

# Run tests (offline; horizons marker deselected by default)
docker compose run --rm test pytest tests/test_refine_pair_nbody.py -v

# Run horizons cross-check (requires network)
docker compose run --rm test pytest tests/test_refine_pair_nbody.py -v -m horizons
```
