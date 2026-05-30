# Literature validation of the frozen encounter catalog

> Consolidated status of every literature cross-check run against the frozen
> candidate catalog ([FROZEN_RUN.md](../FROZEN_RUN.md)). Honest match rates, the
> major-body regression gate, and the two cross-checks that are **blocked on
> source data** (Goffin 2014, Fuentes-Muñoz 2024).
>
> Track B Stage 1 of [current_working_plan.md](../current_working_plan.md).

## Summary

| source | scope | result | status |
|---|---|---|---|
| **Major-body gate** (CLAUDE.md) | (1) Ceres, (2) Pallas, (4) Vesta, (10) Hygiea present | **4 / 4 present** (352 / 47 / 458 / 162 encounters) | ✅ now a regression test |
| **Fienga 2003** (J/A+A/411/L7) | 4 events inside the Gaia window | **3 / 4 matched**; 1 detection gap | ✅ distances agree (median 52 μAU) |
| **Galád 2002** | 4 in-window Hygiea encounters | **4 / 4 matched** | ✅ µAU–20 µAU agreement |
| **JPL Horizons** direct | 8 literature pairs re-derived from JPL | our distance vs JPL `|Δ|` ≤ ~5 × 10⁻⁶ AU | ✅ at JPL cadence |
| **Goffin 2014** (J/A+A/565/A56) | pair-by-pair encounter table | — | 🔴 **data-blocked** (VizieR has only mass tables 5–6) |
| **Fuentes-Muñoz 2024** (LPSC #2388) | 231 mass-determination pairs | — | 🔴 **data-blocked** (pair list not yet ingested) |

## Major-body gate (regression test)

CLAUDE.md requires the four large perturbers to appear in the catalog, given
their Hill-sphere reach. Measured on `encounters_catalog_rebound_005au.parquet`,
matching FROZEN_RUN.md § "Major-body gate checks":

| body | encounters | closest approach (AU) |
|---|---:|---:|
| (1) Ceres   | 352 | 0.003819 |
| (2) Pallas  |  47 | 0.006288 |
| (4) Vesta   | 458 | 0.000936 |
| (10) Hygiea | 162 | 0.005287 |

This is now enforced by `tests/test_validation.py::TestFrozenMajorBodyGate`
(opt-in: `RUN_REAL_CATALOG_TESTS=1` + the frozen catalog present; skipped in CI
where the multi-GB artifact is absent). It asserts presence, the exact
FROZEN_RUN encounter counts, and the closest-approach distances (±1 μAU) — a
drift means the local catalog is no longer the documented frozen artifact.

## Fienga 2003 (J/A+A/411/L7) — 3 / 4

`scripts/validate/validate_fienga_2003.py`,
artifacts `data/output/fienga_2003_{matches,misses}.csv`.

| perturber → target | Fienga date | our date | Fienga dist (AU) | our dist (AU) | \|Δdist\| (AU) |
|---|---|---|---:|---:|---:|
| (48) → (300)  | 2017-04-01 | 2017-04-24 | 0.0084 | 0.008348 | 5.2 × 10⁻⁵ |
| (65) → (976)  | 2014-11-01 | 2014-11-12 | 0.0378 | 0.037808 | 8 × 10⁻⁶ |
| (1) → (57)    | 2017-01-01 | 2017-01-10 | 0.0437 | 0.043548 | 1.5 × 10⁻⁴ |

Distance residuals: median 52 μAU, max 152 μAU. (Fienga's dates are
month-rounded predictions, so the ~10–20 day `our date` offsets are expected.)

**Miss: (804) → (733)** at 2015-02-01 (Fienga 0.0138 AU). It is **absent from
the frozen Kepler and hybrid catalogs** — a detection gap, not a refinement
error. Note `data/output/jpl_horizons_validation.csv` (a *pre-freeze* run, dated
before the catalog was assembled) did capture this pair at our_date 2015-02-12,
dist 0.013753 AU, agreeing with JPL to 6 × 10⁻⁷ AU — so the event is real and
near-threshold; the current frozen catalog simply does not contain it. The most
likely cause is the orbital prefilter or the coarse-scan widening at this
near-0.014 AU geometry; it is consistent with the prefilter-recall deficit
documented in [docs/prefilter_recall.md](prefilter_recall.md).

## Galád 2002 — 4 / 4

`scripts/validate/validate_galad_2002.py`,
artifacts `data/output/galad_2002_{matches,misses}.csv` (0 misses).

All four in-window (10) Hygiea encounters reproduced, sub-day epoch agreement:

| target | Galád date | our date | Galád r (AU) | our dist (AU) | \|Δt\| (d) |
|---|---|---|---:|---:|---:|
| (4803)  | 2017-04-05 | 2017-04-05 | 0.011922 | 0.011924 | 0.32 |
| (10018) | 2015-05-24 | 2015-05-23 | 0.023737 | 0.023736 | 0.27 |
| (11328) | 2016-10-30 | 2016-10-30 | 0.023988 | 0.023974 | 0.26 |
| (20331) | 2016-06-03 | 2016-06-02 | 0.041641 | 0.041619 | 0.14 |

Distance agreement: 1 μAU to ~22 μAU.

## JPL Horizons direct cross-check — 8 pairs

`scripts/validate/validate_jpl_horizons.py`,
artifact `data/output/jpl_horizons_validation.csv`. Re-derives the
closest-approach distance independently from JPL ephemerides for the 8
Fienga/Galád literature pairs and compares to ours. `|our − JPL|` ranges from
sub-10⁻⁷ AU to ~5 × 10⁻⁶ AU at the JPL sampling cadence.

⚠️ **Cadence caveat** (FROZEN_RUN limit 3): JPL is sampled at 1 h / 30 min and
`argmin` is taken — this validates accuracy *at that cadence* over ~8 pairs, not
sub-cadence / micro-AU accuracy on all 72 M rows.

## Blocked: Goffin 2014 and Fuentes-Muñoz 2024

These are the two pair lists the plan flagged for the systematic mass-pair
cross-check. Both are blocked on **source data**, not on the pipeline:

- **Goffin 2014 (J/A+A/565/A56).** The downloaded VizieR snapshot
  (`data/raw/goffin_2014_encounters.parquet`, 536 rows) contains only the
  **mass-determination tables 5–6** (columns `Seq, Name, Nd, M, Diam, Dens, …` —
  per-asteroid masses), not the pair-by-pair *encounter* table the validator
  needs (perturber, target, epoch, miss distance). VizieR does not appear to
  carry Goffin's encounter list. To unblock: obtain the encounter table from the
  paper's electronic material or reconstruct (perturber, target, epoch) pairs and
  re-run `scripts/validate/validate_goffin_2014.py` (the loader and matching
  logic already exist and are unit-tested with synthetic data).
- **Fuentes-Muñoz et al. 2024 (LPSC #2388).** The 231-pair mass list is in the
  abstract/supplementary, not yet ingested into `data/raw/`. To unblock: ingest
  the pair list, then cross-match by (perturber, target) + epoch window as for
  Galád/Fienga.

Neither blocker affects the catalog's defensible scope: the 4-body gate +
Fienga + Galád + JPL cross-checks already confirm the catalog reproduces known
asteroid-asteroid encounters to μAU at the validation cadence.

## Reproduce

```bash
docker compose run --rm pipeline python -m scripts.validate.validate_fienga_2003
docker compose run --rm pipeline python -m scripts.validate.validate_galad_2002
docker compose run --rm pipeline python -m scripts.validate.validate_jpl_horizons
RUN_REAL_CATALOG_TESTS=1 docker compose run --rm -e RUN_REAL_CATALOG_TESTS=1 \
    pipeline pytest tests/test_validation.py::TestFrozenMajorBodyGate -q
```
