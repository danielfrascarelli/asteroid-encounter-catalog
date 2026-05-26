# Frozen run — Kepler-refined geometric candidate catalog

> Canonical reference for the geometric close-encounter **candidate** catalog
> generated under the Gaia DR3 observation window (2014-07-25 → 2017-05-28).
> Anything claimed in follow-up analyses (figures, tables, papers, dashboards)
> must refer to the values *and the scope* in this document.

## Scope and limits (read this first)

This freeze documents a **candidate catalog under frozen assumptions**, not a
complete or fully N-body-refined catalog and not a mass catalog. Three hard
limits constrain what can be defensibly claimed from it:

1. **Final distances are Kepler-2-body, not N-body.** When tiered mode is on
   (it is — `forced_by_tiered = true` in the sidecar), the rebound trajectory
   is used only for the coarse KD-tree scan. The sub-grid refinement that
   produces the reported minimum distance, encounter epoch, and relative
   velocity runs through `kepler_to_cartesian`. So the closest-approach
   numbers are geometric values under a two-body model, not a full
   gravitational solution. **Measured error budget**
   ([docs/kepler_refine_error_report.md](docs/kepler_refine_error_report.md),
   964 stratified pairs re-refined under full N-body):
   median `|Δdist|` = **12 μAU**, p95 = **678 μAU**, p99 = **2.5 mAU**,
   max = **11.3 mAU**. Error correlates with e_max and inversely with
   q_min — the largest disagreements (>1 mAU) are concentrated in
   high-eccentricity pairs (e_max > 0.45) and NEA-like pairs
   (q_min < 1.3 AU). None of the 964 sample pairs changes detection status
   when re-refined, but ranking *within* the closest mAU band is not
   stable.
2. **Not complete.** The orbital prefilter (\|Δa\| ≤ 0.5 AU, \|Δi\| ≤ 30°)
   is a heuristic that can drop real high-eccentricity / high-inclination
   crossing orbits. Recall on the high-e/i tail has not been quantified.
   Audit blocker #2 is still open. The word "complete" must not be used.
3. **Validation precision is sampling-cadence-limited.** The Horizons cross-
   checks (`scripts/validate/validate_jpl_horizons.py`,
   `scripts/validate/validate_novel_a.py`) sample JPL at 1 h or 30 min and
   take `argmin` — sub-cadence precision is *not* validated. The "0 μAU MAE"
   headline in `VALIDATION_SUMMARY.md` is at that cadence over ~8 literature
   pairs, not a global proof of micro-AU accuracy on 72 M rows.

What this freeze **does** support: claims about the candidate list under the
exact configuration recorded in the provenance sidecar — "pairs whose
Kepler-refined minimum distance was ≤ 0.05 AU under the prefilter that was
applied". Anything stronger (completeness, sub-km accuracy, mass detection)
needs separate validation work that is **not** in this freeze.

## TL;DR

| field | value |
|---|---|
| catalog file | `data/output/encounters_catalog_rebound_005au.parquet` |
| catalog SHA-256 | `b0272be7aab649b4f01d85f79011c72074f3c01b7695613ce033cd16bb0fb5e6` |
| size on disk | 2,806,825,159 B (2.6 GiB) |
| rows (Kepler-refined candidates) | **72,236,904** |
| threshold | 0.05 AU |
| provenance sidecar | `data/output/encounters_catalog_rebound_005au_provenance.json` |
| MPCORB snapshot used | `MPCORB_20160217.DAT` (SHA-256 prefix `3e44e7d36b59a7ff`) |
| scan method | rebound (whfast, Sun + Jupiter + Saturn, dt = 1 h) |
| coarse grid | Δt = 12 h |
| refine method | **Kepler 2-body** (forced by tiered mode) on Δt = 120 s window of ±2 h |
| prefilter | enabled — \|Δa\| ≤ 0.5 AU, \|Δi\| ≤ 30° (heuristic; recall not quantified) |
| pipeline code | `main` at commit `b1c4d9a` (audit rounds 1+2 merged) plus the `fine_time_step_seconds=120` setting backported to `config.yaml` so this catalog is reproducible from current main with the same config |

## What this run **is**

A list of asteroid pairs whose Kepler-refined heliocentric ecliptic 3-D
separation came within 0.05 AU during the Gaia DR3 observation window,
under the prefilter and propagation configuration recorded above.

## What this run is **not**

- **Not a mass catalog.** 41 pairs in
  `data/output/mass_followup_candidates.csv` are *follow-up targets* for
  potential mass-fitting work, not mass measurements. The mass-fitting
  layer (`scripts/mass/fit_mass_gaia_loo.py` and friends) is exploratory
  with ~arcsecond systematic residuals. Diagnostic batch fits in
  `data/output/loo_batch_results.csv` show median `chi²_red_window ≈ 425`
  and a maximum of ~7.2 × 10⁵ — clear evidence that the forward model is
  mis-specified and is absorbing orbital drift, not isolating gravitational
  perturbation. The specificity test (`data/output/specificity_ranking.csv`)
  returns **0 / 41** encounter-specific detections. The mass layer cannot
  be cited as a result.
- **Not a fully N-body-refined catalog.** See limit 1 above. The filename
  contains `rebound` because the coarse scan used N-body; the *final*
  reported distance for each row came out of the Kepler refiner.
- **Not strictly complete.** See limit 2 above. Audit blocker #2.
- **Not validated at sub-cadence precision.** See limit 3 above. The
  validation MAE only constrains accuracy at the JPL sampling cadence
  (1 h / 30 min) over a small literature sample, not at micro-AU.

## Inputs

| input | path | SHA-256 prefix | size |
|---|---|---|---|
| MPCORB snapshot | `data/raw/mpcorb_archive/MPCORB_20160217.DAT` | `3e44e7d36b59a7ff` | 305 MiB |
| Gaia DR3 SSO observations | `data/raw/gaia_sso.parquet` | (16,000 rows used for validation only) | 905 KiB |
| Gaia DR3 SSO orbits | `data/raw/gaia_orbits.parquet` | (used for downstream characterization) | 4.9 MiB |
| MPCORB metadata sidecar | `data/raw/MPCORB.json` | `47d55f34c28b448f` | — |

The MPCORB snapshot from 2016-02-17 was selected because its osculating
epoch sits at the centre of the Gaia DR3 observation window, minimising
the propagation distance for Kepler refinement.

## Configuration (frozen)

Full machine-readable copy lives in
`data/output/encounters_catalog_rebound_005au_provenance.json` under the
`config` key.  Key parameters:

```yaml
observation_window:    2014-07-25 → 2017-05-28  (Gaia DR3 SSO)
threshold_au:          0.05
time_step_hours:       1.0          # fine grid for refinement
coarse_step_hours:     12.0         # bulk N-body cache
fine_step_seconds:     60.0         # Kepler refinement sub-grid
window_hours:          2.0          # ±2 h around each scan minimum
propagation.method:    rebound      # whfast, Sun+Jupiter+Saturn
prefilter:             enabled, Δa≤0.5 AU, Δi≤30°
subset:                only_numbered = true
```

## Claims (numbers from this run)

### Encounter counts by separation

| distance bound | encounters | fraction |
|---|---:|---:|
| d < 0.001 AU |     26,038 | 0.04 % |
| d < 0.005 AU |    704,413 | 0.98 % |
| d < 0.010 AU |  2,833,425 | 3.92 % |
| d < 0.020 AU | 11,403,496 | 15.79 % |
| d < 0.030 AU | 25,779,325 | 35.69 % |
| d < 0.040 AU | 46,035,452 | 63.73 % |
| d < 0.050 AU | 72,236,904 | 100.00 % |

### Major-body gate checks (must be present)

| asteroid | encounters | closest approach (AU) |
|---|---:|---:|
| (1) Ceres   | 352 | 0.003819 |
| (2) Pallas  |  47 | 0.006288 |
| (4) Vesta   | 458 | 0.000936 |
| (10) Hygiea | 162 | 0.005287 |

All four required major bodies are present, satisfying the regression gate.

### Top-10 closest encounters

| body 1 | body 2 | JD (TDB) | d (AU) |
|---|---|---:|---:|
| (153222) 2000 YD43 | (238587) 2004 YX3   | 2457500 | 6.6 × 10⁻⁶ |
| (15072) Landolt    | (387599) 2001 XF180 | 2457000 | 1.2 × 10⁻⁵ |
| (270730) 2002 QE130| (366918) 2005 UC211 | 2457500 | 1.5 × 10⁻⁵ |
| (161150) 2002 SL25 | (412792) 2014 PU21  | 2457600 | 1.6 × 10⁻⁵ |
| (117065) 2004 KD9  | (439086) 2011 QP5   | 2457200 | 1.7 × 10⁻⁵ |
| (52249) 1981 EK21  | (408138) 2013 CL75  | 2457700 | 1.7 × 10⁻⁵ |
| (17067) 1999 GF19  | (236737) 2007 JC18  | 2457100 | 2.0 × 10⁻⁵ |
| (435807) 2008 VV60 | (436353) 2010 JC112 | 2457400 | 2.4 × 10⁻⁵ |
| (110273) 2001 SX251| (221390) 2005 YS34  | 2457100 | 2.5 × 10⁻⁵ |
| (209619) 2005 AT19 | (304025) 2006 DR59  | 2457000 | 2.6 × 10⁻⁵ |

These figures come from the Kepler refiner; positions from the rebound scan
agreed to within ~10⁻⁷ AU (≈ 12 km) on the regression benchmark — see
`monitoring/2026-05-23_pipeline_run.md`.

## Caveats and known limitations (freeze-aware)

- **Observability columns are computed elsewhere.**  `solar_elongation_deg`,
  `gaia_observable`, and apparent magnitude are *not* on this parquet —
  they live on `data/output/encounters_characterized.parquet`, which is a
  characterisation of the smaller 158 k-row detection run, not this 72 M
  catalog.  Re-characterising 72 M rows requires a streaming refactor that
  is **not** part of this freeze.
- **Frame fix landed after the catalog was written.**  PR #21 (audit
  blocker #1) corrected the Earth-position frame in
  `src.characterize.observability`.  The data in this parquet is unaffected
  because the frame bug was downstream of detection; any future
  characterisation of this catalog will use the corrected frame.
- **Prefilter recall is unverified.**  Audit blocker #2 — high-eccentricity
  pairs with \|Δa\| > 0.5 AU could be missing.  Do not claim "complete".
- **MPCORB.DAT in `data/raw/` is the *current* download**, not the
  snapshot used for this run.  Use `data/raw/mpcorb_archive/MPCORB_20160217.*`
  when reproducing.

## Reproducing this run

```bash
# 1. Restore the exact MPCORB snapshot
cp data/raw/mpcorb_archive/MPCORB_20160217.DAT data/raw/MPCORB.DAT

# 2. Use the local config that selects rebound mode + 0.05 AU threshold
docker compose run --rm pipeline python -m scripts.pipeline.run_pipeline \
    --config config.local.yaml

# 3. Bit-identical replay requires the same code commit and dependency
#    versions listed in the provenance sidecar.
```

## Code commit

This freeze documents output produced under the `perf/refine-kepler-cache`
branch (HEAD at `06de6d0` when this catalog was written).  Subsequent
audit fixes merged to `main` do **not** change the data in this parquet:

| PR | What changed | Effect on this catalog |
|---|---|---|
| #21 | Earth frame in `src/characterize/observability.py` (ICRS → helio ecliptic) | None — bug was downstream of detection |
| #22 | Gaia SSO epoch convention in docstrings only | None — code already correct |
| #23 | `write_detection_sidecar` and backfill script | None — adds the sidecar that documents this catalog |
| #24 | This document | — |
| #25 | Scripts moved into `ingest/`, `pipeline/`, etc. | None — same code, different paths |
| #26 | Audit round-2 cleanups (docstrings, MPCORB-snapshot selection in characterize) | None — affects future runs |

To reproduce the catalog from current `main`:

1. Restore the MPCORB snapshot:
   `cp data/raw/mpcorb_archive/MPCORB_20160217.DAT data/raw/MPCORB.DAT`
2. Use the default `config.yaml` (which now has
   `fine_time_step_seconds: 120`, matching what produced this catalog).
3. Bit-identical replay still requires the same dependency versions
   listed in the provenance sidecar.

Anything generated from this catalog (figures, derived tables) must cite
both this freeze (catalog SHA `b0272be7…`) and the code commit at *its own*
generation time.
