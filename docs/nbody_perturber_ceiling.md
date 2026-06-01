# N-body perturber set — systematic-error ceiling

> FOLLOWUP_PLAN item 2. Quantifies the precision ceiling imposed by the frozen
> run's truncated N-body perturber set.

## The problem

The frozen catalog's coarse scan used a REBOUND integration with massive bodies
= **Sun + Jupiter + Saturn only** (`include_planets: [sun, jupiter, saturn]`,
`include_major_asteroids: false` — see the provenance sidecar and
[FROZEN_RUN.md](../FROZEN_RUN.md)). Planetary states come from astropy's
**builtin** JPL ephemeris (`solar_system_ephemeris.set('builtin')`), a low-
precision offline bundle, **not** DE440/SPICE
([src/propagate/nbody.py:12](../src/propagate/nbody.py#L12)).

Two consequences for how the error budget must be read:

1. **The Stage A/B error budget is *internal*.** The Kepler-refine validation in
   FROZEN_RUN.md compares the reported Kepler-2-body distances against this same
   3-body N-body model. It measures the Kepler-vs-N-body gap, **not** how far the
   truncated perturber set is from the true dynamics.
2. **No DE440-grade external truth.** The only external check is the JPL Horizons
   cross-validation (which uses a DE44x ephemeris), and that is sampling-cadence-
   limited over a handful of literature pairs ([VALIDATION_SUMMARY.md](../VALIDATION_SUMMARY.md)).

So the missing perturbers — Uranus, Neptune, the terrestrials — are an
**unquantified systematic** until measured directly. This document measures it.

## Method

[scripts/validate/measure_nbody_perturber_ceiling.py](../scripts/validate/measure_nbody_perturber_ceiling.py)
draws a stratified sample from the frozen catalog and re-refines each pair under
N-body **twice**, changing only the perturber set:

* **baseline** — `(Sun, Jupiter, Saturn)`, the frozen configuration;
* **full** — all eight planets `(Sun … Neptune)`.

Major asteroids are excluded from both (matching the frozen scan), so the
comparison isolates the *planetary* perturber set. The statistic of interest is
`Δdist = dist_full − dist_baseline` at closest approach. Strata mirror the
FROZEN_RUN high-error cut: `cold` (`q_min ≥ 1.8 AU ∧ e_max ≤ 0.3`), `high_e`
(`e_max > 0.3`), `low_q` (`q_min < 1.8 AU`).

Run (seed 42, 30 pairs/stratum, ±6 h IAS15 window, 60 s sampling):

```
docker compose run --rm pipeline python -m scripts.validate.measure_nbody_perturber_ceiling \
    --sample-per-stratum 30 --subsample-every 1500
```

## Result

90 pairs, `|Δdist|` (distance shift from adding Uranus/Neptune/terrestrials):

| stratum | n | median | p95 | p99 | max |
|---------|---|--------|-----|-----|-----|
| **all** | 90 | **1.3 μAU** | **36 μAU** | **67 μAU** | **80 μAU** |
| cold    | 30 | 1.0 μAU | 14 μAU | 32 μAU | 40 μAU |
| low_q   | 30 | 1.5 μAU | 31 μAU | 41 μAU | 45 μAU |
| high_e  | 30 | 2.5 μAU | 57 μAU | 76 μAU | 80 μAU |

Signed median ≈ 0 (−5 × 10⁻⁹ AU): the effect is scatter, not a one-directional
bias. The shift grows on the high-eccentricity / low-perihelion tail, as
expected (those orbits sample the inner system where the terrestrials matter).

## Interpretation

The perturber-set ceiling is **~1 μAU median, ≤ 80 μAU max** on this stratified
sample. For context:

- The detection threshold is 0.05 AU = **50,000 μAU**. The worst-case perturber
  shift (80 μAU) is **0.16 %** of the threshold.
- It is **two to three orders of magnitude smaller** than the Kepler-2-body
  refinement error already documented in FROZEN_RUN.md (Stage A p99 = 2.5 mAU,
  Stage B max = 15.2 mAU = 15,200 μAU).

**Conclusion.** The truncated `Sun + Jupiter + Saturn` perturber set is **not**
the dominant error term in the frozen catalog — the Kepler-2-body refinement is,
by ~100×. Adding the remaining planets would change closest-approach distances
by ~1 μAU typically and ≤ 80 μAU in the worst (high-e) case, well within the
existing distance error bars. The builtin-ephemeris / 3-body model is adequate
for a *candidate* catalog at the 0.05 AU threshold; it would need revisiting only
for a precision (sub-mAU, DE440-referenced) catalog, which this freeze explicitly
does not claim to be.

A full DE440-referenced bound (vs. *true* dynamics rather than the 8-planet
model) still requires a Horizons/SPICE comparison at finer cadence — out of scope
here, but the 8-planet result caps the *planetary truncation* contribution.
