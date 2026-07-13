# Orbital prefilter recall on the high-e / high-i tail (audit blocker #2)

> Quantifies what fraction of real close encounters the orbital prefilter
> drops, in the adverse regime where it is most likely to fail. Closes the
> *measurement* side of audit blocker #2 / FROZEN_RUN.md caveat #2: the recall
> is now a number, not an unknown.
>
> **Reproduce**: `docker compose run --rm pipeline python -m
> scripts.validate.benchmark_prefilter_recall`
> (artifacts under `data/output/prefilter_recall/`, gitignored).

## TL;DR

| quantity | value |
|---|---|
| adverse subset | numbered, `a∈[1.5,4.0] AU`, `e>0.3 ∨ i>15°` — **52,411 bodies** |
| MPCORB snapshot | `MPCORB_20160217.DAT` (the frozen-run snapshot) |
| real encounters `<0.05 AU` (no prefilter, Kepler) | **606,393** |
| kept by the current prefilter | 463,164 |
| **missed by the current prefilter** | **143,229** |
| **recall (current prefilter)** | **76.38 %** [95 % Wilson CI 76.27–76.49 %] |
| dominant cause | the `|Δa| ≤ 0.5 AU` cut (142,645 of 143,229 misses; 99.6 %) |
| recommended fix | replace `|Δa|≤0.5` with a **threshold-padded radial-overlap test** → **100 % recall**, same cost |

**Verdict: the current prefilter is NOT safe for a completeness claim on the
high-e/high-i tail.** Recall there is ~76 %, far below the 99 % bar. Caveat #2
stays open as a *quantified* caveat for the frozen catalog, with a concrete,
zero-cost fix recommended for any future run.

## What the prefilter does, and why it can fail

[src/detect/prefilter.py](../src/detect/prefilter.py) keeps a pair only if

```
|a₁ − a₂| ≤ 0.5 AU   AND   |i₁ − i₂| ≤ 30°
```

before the KD-tree spatial scan, to avoid the O(N²) pair explosion. It is a
heuristic on osculating elements, **not** a proof of geometric impossibility: a
high-eccentricity orbit reaches well inside the belt at perihelion and well
outside at aphelion, so it can physically cross — and come within 0.05 AU of —
an orbit whose semimajor axis differs by far more than 0.5 AU. The `|Δa|≤0.5`
cut is blind to eccentricity, so it discards exactly those crossings.

## Method

The prefilter is a **pure deterministic mask on orbital elements**. Running
detection *with* the prefilter is therefore mathematically identical to running
it *without* and intersecting the result with the mask — the KD-tree scan and
the Kepler refinement are byte-for-byte the same on the surviving pairs. We
exploit this:

1. **Adverse subset** — numbered MBAs in the frozen-catalog scope
   (`a∈[1.5,4.0]`, `only_numbered`) with `e>0.3 ∨ i>15°`: 52,411 bodies, the
   regime where the Δa/Δi heuristic is most likely to fail.
2. **Ground truth `F`** — detection **without** the prefilter over the full
   Gaia DR3 window (2014-07-25 → 2017-05-28). With `pairs=None` the scan is
   KD-tree-spatial-only, O(N log N) per step, so the full 52 k-body population
   is tractable without materialising O(N²) pairs. `F` = 606,393 encounters
   with Kepler-refined minimum distance ≤ 0.05 AU.
3. **Apply the mask** to `F` → `P` (what the prefilter keeps).
   `recall = |P| / |F|`.
4. **Cross-check** — on a 2,000-body sub-sample we additionally ran the *real*
   pipeline with the prefilter enabled and confirmed it reproduces the analytic
   mask exactly (`pipeline_equals_analytic_mask = true`, 0 pairs in either
   symmetric difference, `with_prefilter ⊆ no_prefilter`). The analytic shortcut
   is therefore exact, not an approximation.

**Propagation is Kepler 2-body throughout** (`method=kepler`). This matches the
model that produces the *final reported distances* in the frozen catalog (there
the rebound trajectory drives only the coarse scan; every reported minimum
distance comes out of the Kepler refiner). The recall *ratio* is robust to the
propagation model regardless, since both the with- and without-prefilter sets
use the same model.

Scope note: this measures recall for **adverse–adverse** encounters (both
bodies in the 52 k subset). It does not cover adverse–normal pairs, so the
143,229 missed encounters are a **lower bound** on what the prefilter removes
from the high-e/i tail of the full catalog.

## Results

### Headline

- **Recall = 76.38 %** [95 % CI 76.27–76.49 %]. 143,229 of 606,393 real
  adverse-tail encounters are absent from a prefiltered catalog.

### Recall by band

Recall vs `|Δa|` is a cliff at the 0.5 AU cut — ~99.9 % below it, 0 % above
(by construction, every `|Δa|>0.5` pair is discarded):

| `|Δa|` (AU) | n encounters | recall | missed |
|---|---:|---:|---:|
| [0, 0.1)   | 220,968 | 99.93 % | 164 |
| [0.1, 0.25)|  87,767 | 99.77 % | 198 |
| [0.25, 0.5)| 155,013 | 99.86 % | 222 |
| [0.5, 0.75)| 101,692 |  0.00 % | 101,692 |
| [0.75, 1)  |  33,283 |  0.00 % | 33,283 |
| [1, 2)     |   7,660 |  0.00 % | 7,660 |
| [2, 10)    |      10 |  0.00 % | 10 |

So **142,645 of 143,229 misses (99.6 %) are the `|Δa|>0.5` cut**; only 584 are
the `|Δi|>30°` cut. Max `|Δa|` among real encounters: **2.33 AU**.

Recall degrades smoothly with eccentricity — the heuristic is worst exactly
where it matters:

| `e_max` of the pair | n | recall |
|---|---:|---:|
| [0, 0.1)   |  76,887 | 99.86 % |
| [0.1, 0.2) | 173,098 | 87.45 % |
| [0.2, 0.3) | 172,248 | 73.47 % |
| [0.3, 0.5) | 168,517 | 60.16 % |
| [0.5, 0.7) |  12,056 | 47.07 % |
| [0.7, 1)   |   3,587 | 39.17 % |

(Full per-band tables incl. `Δi` and `i_max` in
`data/output/prefilter_recall/summary.json`.)

## Recommended fix: eccentricity-aware radial-overlap prefilter

Two bodies within `D` AU in 3-D necessarily have heliocentric distances within
`D`. So their radial ranges `[q,Q] = [a(1−e), a(1+e)]`, **padded by the
threshold**, must overlap. This is a *provable necessary condition* for a `<D`
encounter — and just as cheap as the current cut (vectorised min/max
comparisons, O(N²) but no propagation):

```
keep iff  max(q₁, q₂) − D ≤ min(Q₁, Q₂)
```

Measured recall on the adverse subset (no Δi cut needed):

| prefilter criterion | recall | missed |
|---|---:|---:|
| current: `|Δa|≤0.5 ∧ |Δi|≤30°` | **76.38 %** | 143,229 |
| radial-overlap, no pad | 99.62 % | 2,332 |
| **radial-overlap, pad = 0.05 AU (= threshold)** | **100.0000 %** | **0** |
| widen `|Δa|≤1.5 ∧ |Δi|≤30°` | 99.79 % | 1,290 |
| widen `|Δa|≤2.0 ∧ |Δi|≤30°` | 99.83 % | 1,054 |

The threshold-padded radial-overlap test recovers **every** missed encounter
(0 misses, provably complete under the propagation model) at the same cost. The
residual 2,332 misses of the un-padded version are pairs whose ranges sit
within 0.05 AU of overlapping (one near aphelion, the other near perihelion) —
exactly what the pad accounts for.

If a minimal change is preferred over swapping the criterion, **widening `|Δa|`
to ≥1.5 AU** brings recall to 99.8 %; the inclination cut should also be relaxed
(to ≥45°, or dropped — it costs ~1,000 high-Δi crossings) since radial overlap
already encodes the physics.

## Impact on the frozen catalog and on caveat #2

> **Corrección (2026-05-31, reafirmada 2026-07-04):** el pair-list prefilter
> solo se construye para N ≤ 5.000; el freeze corrió a escala de cinturón
> (449.454 cuerpos), así que `compatible_pairs` fue **saltado**
> (`prefilter.effective = "skipped_large_n"` en el sidecar) y solo se aplicó la
> query espacial del KD-tree. Las cifras de recall de este documento son
> **contrafactuales para el freeze**: miden el daño que el prefiltro causaría
> *si se aplicara a N chico*, no un déficit sufrido por el catálogo congelado.

- El freeze **no** perdió los 143.229 encuentros adversos por el prefiltro (el
  corte `|Δa|≤0.5` no le dio forma a ese catálogo). Su incompletitud real viene
  de otras fuentes: cadencia gruesa de 12 h, censura del umbral Kepler 0.05 AU
  (0.70 % medido) y — hasta la regeneración post-B1 — el recorte de ventana de
  refinamiento (`docs/tribunal_cientifico_2026-07-04.md`, B1).
- Las cifras de esta medición aplican a cualquier corrida futura con
  N ≤ 5.000 y prefiltro activo (p. ej. subsets de validación): ahí el
  radial-overlap con pad es el criterio recomendado (recall 100 % medido).
- **Caveat #2 sigue abierto pero por las causas de arriba**, no por el
  prefiltro. Cerrarlo = presupuesto de completitud de la Tarea 5 del plan de
  remediación (`planning/TRIBUNAL_REMEDIATION_PLAN.md`).

## Artifacts

- [scripts/validate/benchmark_prefilter_recall.py](../scripts/validate/benchmark_prefilter_recall.py)
- `data/output/prefilter_recall/summary.json` — recall + all per-band tables + cross-check
- `data/output/prefilter_recall/adverse_no_prefilter_encounters.parquet` — ground-truth `F` with per-pair (Δa, Δi, e_max, i_max)
- `data/output/prefilter_recall/missed_pairs.parquet` — the 143,229 dropped encounters
