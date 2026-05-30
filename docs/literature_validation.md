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
| **Fuentes-Muñoz 2025** (AJ 170, 353) | 40,004 perturber→target pairs (Gaia FPR) | **11,804 (29.5 %) present** in the DR3 catalog | ✅ overlap is a lower bound (see below) |
| **Goffin 2014** (J/A+A/565/A56) | pair-by-pair encounter table | — | 🔴 **no machine-readable encounter list exists** (VizieR confirmed: masses only) |

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

## Fuentes-Muñoz et al. 2025 (AJ 170, 353) — 11,804 confirmations

`scripts/ingest/download_fuentes_munoz.py` →
`scripts/validate/validate_fuentes_munoz_2025.py`,
artifacts `data/output/literature_validation/fuentes_munoz_2025_{matches.parquet,summary.json}`.

The published version of the Gaia-FPR mass study (first shown as LPSC 2024 #2388)
ships a machine-readable **Table 5**: each perturber asteroid with a
pipe-delimited list of the (≤100 highest-signal) test asteroids whose astrometry
showed a measurable mass signal — i.e. a dynamically significant close encounter.
Parsing it (numbered perturbers only; provisional designations such as
`2013 KY18` and provisional targets dropped) yields **40,004 unique numbered
(perturber, target) pairs** from 1,645 numbered perturbers.

Cross-matching by unordered (number, number) against the frozen catalog:

- **11,804 / 40,004 (29.5 %) of the Fuentes-Muñoz pairs are present** in our
  DR3-window catalog as < 0.05 AU encounters. Matched separations span
  0.0002–0.050 AU (median 0.020 AU).

⚠️ **This overlap is a LOWER BOUND, not a recall.** Fuentes-Muñoz fit orbits over
the full Gaia FPR + archival astrometry baseline, so a pair's signal-producing
encounter can fall at *any* epoch — usually **outside** our DR3 window
(2014-07-25 → 2017-05-28), or beyond 0.05 AU, or outside our `a∈[1.5,4.0]`
numbered-MBA scope. A Fuentes-Muñoz pair absent from our catalog is therefore
*expected*, not a miss. The defensible statement is the positive one: **11,804
pairs flagged by an independent Gaia-FPR mass study have their close approach
reproduced inside our window**, a large-scale cross-confirmation. The parser is
unit-tested (`tests/test_validation.py::TestFuentesMunozParse`).

## Blocked: Goffin 2014 — no machine-readable encounter list

`scripts/validate/validate_goffin_2014.py` (loader + matcher exist, unit-tested
with synthetic data) cannot run because the source data does not exist in
machine-readable form. The VizieR catalog `J/A+A/565/A56` was inspected directly
(ReadMe + the downloaded `data/raw/goffin_2014_encounters.parquet`, 536 rows): it
contains **only** table5 ("Asteroid masses obtained"), table6 ("overview of mass
determinations" — a per-perturber literature mass comparison: `Seq, Name, M, e_M,
Type, Ref`), and `refs`. **There is no pair-by-pair encounter table** (perturber,
target, epoch, miss distance) — Goffin's underlying encounter list was never
published as a VizieR table. To unblock: extract the encounter list from the
paper's body / electronic material, or treat Fuentes-Muñoz 2025 (above) as the
machine-readable successor for the mass-pair cross-check. This is a source-data
limitation, not a pipeline one.

None of this affects the catalog's defensible scope: the 4-body gate + Fienga +
Galád + JPL + 11,804 Fuentes-Muñoz confirmations already show the catalog
reproduces known asteroid-asteroid encounters to μAU at the validation cadence.

## Reproduce

```bash
docker compose run --rm pipeline python -m scripts.validate.validate_fienga_2003
docker compose run --rm pipeline python -m scripts.validate.validate_galad_2002
docker compose run --rm pipeline python -m scripts.validate.validate_jpl_horizons
docker compose run --rm pipeline python -m scripts.ingest.download_fuentes_munoz
docker compose run --rm pipeline python -m scripts.validate.validate_fuentes_munoz_2025
RUN_REAL_CATALOG_TESTS=1 docker compose run --rm -e RUN_REAL_CATALOG_TESTS=1 \
    pipeline pytest tests/test_validation.py::TestFrozenMajorBodyGate -q
```
