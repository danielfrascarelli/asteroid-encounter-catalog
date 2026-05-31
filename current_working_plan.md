# Current working plan — trabajo pendiente post-cierre de la capa de masas

> Plan vivo del trabajo pendiente tras disolver `FOLLOWUP_PLAN.md` (Track A capa
> de masas cerrada, Track B B1/B2/B3 completos — todo consolidado en
> [ROADMAP.md](ROADMAP.md) § "Estado actual" y los `docs/mass_layer_*.md` /
> [docs/kepler_threshold_bias_paper.md](docs/kepler_threshold_bias_paper.md)).
>
> Cubre tres bloques, en orden de valor:
>
> - **Track A — Completitud y rigor del catálogo**: lo que separa al catálogo
>   (el producto **publicable**, Track 1) de poder reclamar completitud y llevar
>   columnas de observabilidad. Camino crítico.
> - **Track B — Cierre del proyecto (Fase 7)**: validación contra literatura,
>   dashboard y README final / reproducibilidad.
> - **Track C — Opcionales / deuda técnica**: no bloqueantes (perf, experimento
>   de falsos negativos del threshold, bookkeeping).
>
> Este archivo es el **único registro vivo** del avance pendiente.

---

## Cómo usar este documento

**Para arrancar / retomar cold**:

1. Leer "Estado global" — qué etapa está activa.
2. Saltar a la etapa `🟡 IN PROGRESS` o la siguiente `⚪ PENDING`.
3. Leer "Cómo retomar" de esa etapa.
4. Confirmar contra el filesystem los entregables esperados.
5. Avanzar.

**Para pausar**:

1. Commitear el WIP en su branch (`wip(track-X-stageY): <qué quedó>`).
2. Actualizar "Progreso" de la etapa: fecha, qué quedó, qué bloquea, próximo paso.
3. Mover estado a `🟠 PAUSED` si se deja por >1 semana, sino `🟡 IN PROGRESS`.
4. Push de la branch (no merge a `main` hasta etapa `🟢 DONE`).

**Convenciones de estado**: ⚪ `PENDING` · 🟡 `IN PROGRESS` · 🟠 `PAUSED` ·
🔴 `BLOCKED` · 🟢 `DONE` · ⚫ `ABANDONED`.

**Branches**: una por etapa, `track{A|B|C}/stage{N}-{slug}`. **Nunca trabajar en `main`**.

---

## Estado global

**Fecha de creación**: 2026-05-30
**Última actualización**: 2026-05-31
**Estado**: 🟢 **PLAN COMPLETO**. Las 7 etapas cerradas: **A1 (PR #53)**, **B1
(#54)**, **A2 (#55)**, **B2 (#56)**, **C3 (#57)**, **C2 (#58)** y **C1**
(perf, evaluado y descartado). Único residual no-autónomo: **QA visual humana del
dashboard** (B2). La capa de masas sigue cerrada (no determinable en DR3; sólo
reabre con DR4/FPR).

| Track | Etapa | Estado | Branch | PR | Notas |
|-------|-------|--------|--------|----|-------|
| A | 1: recall del prefiltro (audit #2) | 🟢 DONE | `trackA/stage1-prefilter-recall` | (pendiente) | Recall=76.4% en cola adversa; fix radial-overlap → 100%. Cuantificado, no cerrado |
| A | 2: caracterización catálogo 72M (streaming) | 🟢 DONE | `trackA/stage2-characterize-bigcatalog` | (pendiente) | 72.2M caracterizadas (18.9% observables), 0 OOM; streaming chunked + tests de paridad |
| B | 1: validación literatura completa (Fase 7) | 🟢 DONE | `trackB/stage1-literature-validation` | (pendiente) | Gate 4 cuerpos + Fuentes-Muñoz 2025 (11.8k confirmaciones) + consolidación; Goffin sin datos en VizieR |
| B | 2: dashboard + README final + reproducibilidad | 🟡 IN PROGRESS | `trackB/stage2-dashboard-readme` | (pendiente) | Data layer memory-safe + README final + reproducibilidad DONE; falta QA visual humana del dashboard |
| C | 1: perf followups (numba / cache persistente) | 🟢 DONE | `trackC/stage1-perf` | (pendiente) | Evaluado → numba y cache persistente descartados con mediciones (docs/perf_evaluation.md) |
| C | 2: experimento falsos negativos threshold Kepler | 🟢 DONE | `trackC/stage2-kepler-false-negatives` | (pendiente) | Tasa medida 0.70 % en [0.05,0.06); matriz de confusión cerrada; RNAAS standalone viable (opcional) |
| C | 3: bookkeeping (audit #6, refs stale) | 🟢 DONE | `trackC/stage3-bookkeeping` | (pendiente) | #6 resuelto por cierre; ref stale podado (perf/refine-kepler-cache conservado); ROADMAP actualizado |

**Recomendación de orden**:

1. **A1 primero** (prefilter recall). Es el único caveat científico que hoy
   impide reclamar completitud del catálogo, y es acotado (~2-3 h).
2. **B1 y B2** (Fase 7) llevan el proyecto a "cerrado/publicable".
3. **A2** (caracterización 72M) si se necesita observabilidad sobre el catálogo
   completo; es un refactor mayor, evaluable después de A1.
4. **Track C** oportunista.

---

## Track A — Completitud y rigor del catálogo

### Stage 1 — Cuantificar el recall del prefiltro orbital (audit blocker #2)

**Estado**: 🟢 DONE
**Estimación**: ~2-3 h (mayormente cómputo)
**Branch**: `trackA/stage1-prefilter-recall`
**Depende de**: nada.

#### Objetivo

El prefiltro orbital ([src/detect/prefilter.py](src/detect/prefilter.py)) descarta
pares con `|Δa| > 0.5 AU` o `|Δi| > 30°` antes del KD-tree, por costo. Es un
heurístico cuyo **recall en la cola de alta e / alta i no está medido** — hoy el
catálogo **no puede reclamar completitud** ([FROZEN_RUN.md](FROZEN_RUN.md) caveat
#2). Medir qué fracción de encuentros reales (<0.05 AU) pierde el prefiltro.

#### Plan técnico

1. **Seleccionar un subset adverso** (~½ h): asteroides de alta e (>0.3) y/o
   alta i (>15°) donde el prefiltro Δa/Δi es más probable que falle. ~2000-5000
   cuerpos del MPCORB snapshot del run congelado.
2. **Correr detección SIN prefiltro** sobre ese subset (cómputo N², acotado por
   el tamaño del subset) y CON prefiltro, mismos parámetros que el run congelado.
3. **Comparar conjuntos de encuentros** (<0.05 AU): recall = |prefiltered ∩
   no-prefilter| / |no-prefilter|. Caracterizar los **perdidos** por (Δa, Δi, e, i).
4. **Documento** [docs/prefilter_recall.md]: recall global + por banda, y una
   recomendación (ensanchar Δa/Δi, o cuantificar el caveat para el paper).

#### Entregables

- [x] `scripts/validate/benchmark_prefilter_recall.py`
- [x] `data/output/prefilter_recall/` (subset, encuentros sin prefiltro, missed, summary.json)
- [x] `docs/prefilter_recall.md`
- [x] Actualizar el caveat #2 en [FROZEN_RUN.md](FROZEN_RUN.md) con el número medido

#### Criterios de aceptación

- [x] Recall medido sobre el subset adverso con barra de incertidumbre.
- [x] Veredicto emitido: recall = **76.38 %** < 99 % → el prefiltro **no es seguro**
  en la cola; sesgo cuantificado (99.6 % de las pérdidas son el corte `|Δa|≤0.5`)
  y fix recomendado (prefiltro de **solapamiento radial** con pad = threshold →
  **100 % recall**, mismo costo). Caveat #2 queda **cuantificado, no cerrado**
  (cerrarlo exige re-correr el catálogo con el nuevo prefiltro).

#### Cómo retomar

```bash
git checkout trackA/stage1-prefilter-recall
cat docs/prefilter_recall.md          # qué se concluyó
ls data/output/prefilter_recall/      # qué corrió
```

#### Progreso

**2026-05-30 — DONE.** Medido sobre el subset adverso completo (52,411 cuerpos
numerados, `a∈[1.5,4.0]`, `e>0.3 ∨ i>15°`, snapshot `MPCORB_20160217.DAT`):
606,393 encuentros reales `<0.05 AU` sin prefiltro (Kepler); **recall = 76.38 %**
[IC95 76.27–76.49 %], 143,229 perdidos. El cross-check (2,000 cuerpos) prueba que
la máscara analítica == pipeline-con-prefiltro byte-a-byte (0 diferencias), así
que el atajo analítico es exacto. 99.6 % de las pérdidas son el corte `|Δa|≤0.5`
(ciego a la excentricidad). Fix recomendado y verificado: prefiltro de
solapamiento radial `max(q₁,q₂) − D ≤ min(Q₁,Q₂)` → **100.0000 % recall, 0
perdidos**, mismo O(N²). Detalle en [docs/prefilter_recall.md](docs/prefilter_recall.md).
**Próximo**: el fix NO se aplicó al catálogo congelado (requiere re-corrida full);
queda como recomendación para una futura corrida o para DR4/FPR.

---

### Stage 2 — Caracterización del catálogo completo (72M filas, streaming)

**Estado**: 🟢 DONE (refactor + tests + corrida 72M completa, 0 OOM)
**Estimación**: ~2-4 días (refactor + corrida)
**Branch**: `trackA/stage2-characterize-bigcatalog`
**Depende de**: nada (pero correr después de A1).

#### Objetivo

`characterize_catalog` ([src/characterize/](src/characterize/)) materializa el
catálogo entero en memoria → **OOMea a 31 GB** sobre los 72M. El catálogo
caracterizado actual (`encounters_characterized.parquet`) es solo el run de
**158k** filas, no los 72M. Refactorizar a streaming para que el catálogo
completo lleve observabilidad Gaia + magnitudes/diámetros estimados.

#### Plan técnico (preliminar — completar al arrancar)

1. Refactor de `characterize_catalog` a `pl.LazyFrame` + `sink_parquet` por chunks
   (no materializar; respetar la corrección de frame de PR #21 y el reordenamiento
   `_1`=cuerpo mayor del audit round 5).
2. Verificar paridad con el run de 158k (mismas columnas, mismos valores sobre el
   solapamiento).
3. Correr sobre `encounters_catalog_hybrid_stageb.parquet` (o el 72M Kepler).
4. Tests de no-regresión + sidecar de provenance.

#### Entregables

- [x] Refactor streaming en `src/characterize/` (`characterize_catalog_streaming`, chunked, RAM acotada por chunk)
- [x] `data/output/encounters_characterized_full.parquet` (72,236,904 filas, 5.77 GB) + sidecar `_metadata.json`
- [x] Tests de paridad (`tests/test_characterize.py::TestStreamingParity`: chunked==in-memory, archivo streamed==in-memory, sidecar)
- [x] Nota en [FROZEN_RUN.md](FROZEN_RUN.md) (observabilidad ahora cubre el catálogo completo)

#### Criterios de aceptación

- [x] Corre sobre 72M sin OOM: 73 chunks de 1M en 8635 s (~2.4 h), exit 0, pico acotado por chunk (~<2 GB)
  en una máquina que antes OOMeaba a 31 GB. 13,640,870 (18.9 %) Gaia-observables; gate 4/4 idéntico (352/47/458/162).
- [x] Paridad verificada sobre input idéntico: in-memory `characterize_catalog(sort=False)` == primeras 50k filas
  del streamed 72M (exacto), más los unit tests. La "paridad vs 158k" del plan original no aplica: ese run es una
  **detección distinta**, no un subset del 72M; el test correcto es paridad sobre el mismo input.

#### Cómo retomar

```bash
git checkout trackA/stage2-characterize-bigcatalog
cat logs/characterize_full.log            # progreso de la corrida 72M
ls -la data/output/encounters_characterized_full.parquet  # sidecar _metadata.json al lado
# Re-correr: docker compose run --rm pipeline python -m scripts.pipeline.characterize_catalog \
#   --input data/output/encounters_catalog_hybrid_stageb.parquet --streaming on
```

#### Progreso

**2026-05-30 — refactor DONE, corrida en curso.** `characterize_catalog` ahora acepta
`sort=False`; nueva `characterize_catalog_streaming(input, …)` lee el parquet por
chunks (`pyarrow.iter_batches`, sólo las 7 columnas de detección), caracteriza cada
chunk y lo escribe incrementalmente con `pq.ParquetWriter` (RAM acotada por chunk,
no por catálogo). El script `characterize_catalog.py` decide streaming por conteo de
filas (`--streaming auto|on|off`, `--chunk-size`), y el output de streaming va a
`encounters_characterized_full.parquet` (no pisa el run de 158k). Output **no**
ordenado globalmente por dist (preserva orden de input; caracterización es
row-independiente). Tests de paridad verdes. **Corrida 72M completa** (2026-05-31 01:50, exit 0):
72,236,904 filas en 73 chunks, 8635 s; 13,640,870 (18.9 %) Gaia-observables; gate
4/4 idéntico al FROZEN_RUN (352/47/458/162). Output 5.77 GB +
`encounters_characterized_full_metadata.json`. Paridad sobre input idéntico
verificada (primeras 50k == in-memory). Sin OOM (máquina que antes OOMeaba a 31 GB).
FROZEN_RUN actualizado. **Etapa cerrada.**

---

## Track B — Cierre del proyecto (Fase 7 del roadmap)

### Stage 1 — Validación contra literatura completa

**Estado**: 🟢 DONE (gate + Fuentes-Muñoz 2025 + consolidación; Goffin sin datos en VizieR)
**Estimación**: ~2-3 días
**Branch**: `trackB/stage1-literature-validation`
**Depende de**: nada.

#### Objetivo

Cerrar el cruce con literatura sobre el catálogo congelado. Ya existen
`validate_fienga_2003.py` (4/4), `validate_galad_2002.py` (4/4),
`validate_jpl_horizons.py`. Falta el cruce sistemático con **Goffin (2014)** y
**Fuentes-Muñoz et al. (2024)** (los pares usados para determinación de masas) y
confirmar que los encuentros con (1) Ceres, (4) Vesta, (2) Pallas, (10) Hygiea
aparecen (gate de regresión del CLAUDE.md).

#### Plan técnico

1. Parsear/obtener las listas de pares de Goffin 2014 y Fuentes-Muñoz 2024.
2. Cruzar contra `encounters_catalog_hybrid_stageb.parquet` (matching por par +
   ventana temporal); reportar match rate y discrepancias.
3. Verificar el gate de los 4 cuerpos grandes (ya OK en FROZEN_RUN, formalizar
   en test de regresión si no lo está).
4. Documento [docs/literature_validation.md] consolidado.

#### Entregables

- [x] `scripts/ingest/download_fuentes_munoz.py` + `scripts/validate/validate_fuentes_munoz_2025.py` (parser MRT + cross-match)
- [x] `docs/literature_validation.md` — consolidado (gate + Fienga 3/4 + Galád 4/4 + JPL + Fuentes-Muñoz + Goffin)
- [x] Test de regresión `tests/test_validation.py` (gate de 4 cuerpos `TestFrozenMajorBodyGate` + parser `TestFuentesMunozParse`)
- [~] `validate_goffin_2014.py`: imposible — VizieR `J/A+A/565/A56` no tiene tabla de encuentros (solo masas). No es blocker del pipeline.

#### Criterios de aceptación

- [x] Los 4 cuerpos grandes presentes (gate del CLAUDE.md) → **4/4**, ahora test de regresión.
- [x] Match rate reportado honestamente: Fienga **3/4**, Galád **4/4**, JPL 8 pares ≤~5e-6 AU,
  **Fuentes-Muñoz 2025: 11,804/40,004 (29.5 %) presentes** (lower bound — FPR cubre baseline
  fuera de la ventana DR3). Goffin: sin lista de encuentros publicada en VizieR (documentado).

#### Cómo retomar

```bash
git checkout trackB/stage1-literature-validation
cat docs/literature_validation.md     # estado consolidado
# Goffin sólo se puede cerrar si aparece la tabla de ENCUENTROS (no las de masas
# 5-6) en el material electrónico del paper; Fuentes-Muñoz 2025 ya cubre el cruce
# sistemático de pares de masas como sucesor machine-readable.
```

#### Progreso

**2026-05-30 — DONE.** (1) **Gate de los 4 cuerpos** como test de regresión
`TestFrozenMajorBodyGate` contra el catálogo congelado (352/47/458/162, closest a
±1 µAU; opt-in `RUN_REAL_CATALOG_TESTS=1`, skip en CI). (2) **Fuentes-Muñoz et al.
2025 (AJ 170, 353)** — bajé la Tabla 5 machine-readable de IOP (fuente oficial),
parseé 40,004 pares numerados (perturber→target) de 1,645 perturbadores numerados
(provisionales tipo `2013 KY18` y targets provisionales descartados; parser
unit-tested `TestFuentesMunozParse`), y crucé contra el catálogo: **11,804/40,004
(29.5 %) presentes** como encuentros <0.05 AU en la ventana DR3 (lower bound, no
recall: FPR ajusta sobre baseline completo → la mayoría de los encuentros caen
fuera de la ventana/scope; los presentes son confirmaciones independientes).
(3) **docs/literature_validation.md** consolida todo: Fienga 3/4, Galád 4/4, JPL 8
pares, Fuentes-Muñoz 11,804 confirmaciones, y Goffin documentado como
**imposible** (VizieR sólo trae masas, confirmado vía ReadMe). Suite: 30 passed /
26 skipped por defecto; lint+black+mypy verdes.

---

### Stage 2 — Dashboard, README final y reproducibilidad

**Estado**: 🟡 IN PROGRESS (data layer + README + reproducibilidad DONE; falta QA visual humana del dashboard)
**Estimación**: ~2-3 días
**Branch**: `trackB/stage2-dashboard-readme`
**Depende de**: B1 (idealmente, para mostrar validación en el dashboard).

#### Objetivo

`src/dashboard/app.py` ya existe (Streamlit). Pulirlo contra el catálogo
congelado/híbrido, escribir el **README final** para humanos, y verificar la
**reproducibilidad** end-to-end desde cero (todos los comandos del README
funcionan en Docker limpio).

#### Plan técnico

1. Revisar/pulir el dashboard: que levante contra `encounters_catalog_hybrid_stageb.parquet`,
   filtros por distancia/fecha/cuerpo, vistas de los encuentros notables.
2. README final: descripción, quickstart Docker, estructura, alcance defensible
   (citar FROZEN_RUN), limitaciones (capa de masas cerrada, completitud del
   prefiltro — resultado de A1).
3. Smoke test de reproducibilidad: `docker compose build` + comandos clave desde cero.

#### Entregables

- [x] Data layer memory-safe (`src/dashboard/data.py`): stats globales por agregación lazy + display capado a los N más cercanos; el dashboard ahora usa `encounters_characterized_full.parquet` (72M) si existe, con fallback al de 158k. Tests `tests/test_dashboard_data.py`.
- [x] `README.md` final: validación corregida (Fienga 3/4, Galád 4/4, Fuentes-Muñoz 11.8k, gate 4 cuerpos), recall del prefiltro, catálogo caracterizado full, paths de scripts corregidos.
- [x] Checklist de reproducibilidad (sección "🔁 Reproducibilidad" en README, end-to-end desde `docker compose build`).
- [~] QA visual del dashboard (boot headless OK, health 200; falta revisión humana en navegador — Streamlit no es verificable sin cliente).

#### Criterios de aceptación

- [x] `docker compose up dashboard` levanta (boot headless verificado: health HTTP 200, sin traceback); carga memory-safe verificada contra el 72M real (stats 0.7s, closest-300k 4.5s).
- [x] README permite reproducir desde cero (checklist con todos los comandos en Docker).

#### Cómo retomar

```bash
git checkout trackB/stage2-dashboard-readme
docker compose up dashboard   # QA visual en http://localhost:8501 (revisión humana)
```

#### Progreso

**2026-05-31 — data layer + README + reproducibilidad DONE.** (1) Dashboard:
extraje el acceso a datos a `src/dashboard/data.py` (puro, testeable): stats
globales por agregación lazy (no materializa el frame) y display capado a los
300k encuentros más cercanos (top-k lazy) → maneja el 72M sin OOM; prefiere
`encounters_characterized_full.parquet` con fallback al de 158k. `app.py`
consume eso; corregí el header/docstring que sobre-afirmaban (la capa de masas
está cerrada). Boot headless verificado (health 200), `tests/test_dashboard_data.py`
verde, y verificado contra el 72M real. Monté `./src` en el servicio dashboard
del compose. (2) README final: corregí la validación (era "4/4 Fienga / 4M
encuentros" → Fienga 3/4, Galád 4/4, **Fuentes-Muñoz 11.8k**, gate 4 cuerpos),
añadí recall del prefiltro (76.4 %), el catálogo caracterizado full, paths de
scripts (`scripts.ingest.*`/`scripts.validate.*`), y una sección de
**Reproducibilidad** end-to-end. **Falta**: QA visual humana del dashboard en el
navegador (no verificable de forma autónoma).

---

## Track C — Opcionales / deuda técnica

### Stage 1 — Perf followups (no bloqueante)

**Estado**: 🟢 DONE (evaluado → ambas descartadas con mediciones) · **Branch**: `trackC/stage1-perf`

- [x] **numba en `solve_kepler` / hot-path Kepler**: **NO se justifica.** El solver
      ya es numpy-vectorizado (el único loop es la iteración Newton sobre todo el
      array) → ~7.3 M elem/s `solve_kepler`, ~3.7 M elem/s `kepler_to_cartesian`
      single-thread, y el pipeline ya paraleliza por procesos (×n_workers ≈ 100 M
      elem/s agregado). numba acelera loops escalares (no hay), añade dependencia
      (CLAUDE.md la desaconseja) y oversubscribiría contra el Pool de procesos.
- [x] **Cache persistente de fine-grid**: **NO se justifica.** Sólo ayudaría en
      re-corridas idénticas (raras); la parte cara (trayectoria N-body coarse) ya
      está cacheada (`src/propagate/cache.py`, hit <1 s). El refinamiento Kepler
      fino es barato.
- Bench reproducible: `scripts/bench/bench_kepler.py`. Detalle: [docs/perf_evaluation.md](docs/perf_evaluation.md).

---

### Stage 2 — Experimento de falsos negativos del threshold Kepler

**Estado**: 🟢 DONE · **Branch**: `trackC/stage2-kepler-false-negatives`
**Depende de**: nada.

**2026-05-31 — DONE.** `scripts/validate/measure_threshold_false_negatives.py`:
detección Kepler a 0.06 AU sobre 10.000 cuerpos → **17.469 pares en [0.05,0.06)**
re-refinados con N-body (±12 h, IAS15; `spawn` para evitar el deadlock BLAS+fork).
**122 cruzan hacia abajo → tasa de falsos negativos = 0.70 %** [IC95 0.59–0.83 %]
(~0.42 % excluyendo near_boundary). Δdist mediana −3×10⁻⁷, σ=3.7×10⁻⁴ AU:
scatter simétrico, NO sesgo. Matriz de confusión cerca de 0.05 AU cerrada
(~1.5 % arriba / ~0.4–0.7 % abajo). Extrapolación: ~10⁵ encuentros <0.05 censurados
catalog-wide. Actualicé [docs/kepler_threshold_bias_paper.md](docs/kepler_threshold_bias_paper.md)
(sección medida + decisión: **standalone RNAAS ahora es viable, opcional — llamada del autor**)
y FROZEN_RUN límite 1. Artefactos en `data/output/kepler_false_negatives/`.

La nota del sesgo Kepler ([docs/kepler_threshold_bias_paper.md](docs/kepler_threshold_bias_paper.md))
mide la sobre-detección pero **no** la tasa de falsos negativos (el catálogo no
tiene pares con Kepler ≥ 0.05 AU). Para medirla: re-refinar con N-body una muestra
de pares con `d_Kepler ∈ [0.05, 0.06] AU` y contar cuántos bajan de 0.05.
**Gatillo**: si se completa, el conjunto (matriz de confusión del prefiltro
cerca del threshold) justificaría una **nota técnica standalone** (p. ej. RNAAS).

- [ ] Re-detección con threshold 0.06 sobre muestra representativa (solo prefiltro Kepler)
- [ ] N-body refine de los pares en [0.05, 0.06] AU
- [ ] Matriz de confusión completa → actualizar la nota; decidir standalone

---

### Stage 3 — Bookkeeping

**Estado**: 🟢 DONE · **Branch**: `trackC/stage3-bookkeeping`

- [x] **Audit blocker #6 (mass fitting)** marcado **resuelto por cierre** en ROADMAP
      (la capa de masas no es determinable en DR3). Los otros dos abiertos también
      quedan cerrados: **#2** (recall del prefiltro) cuantificado por A1, y
      **big-catalog characterize** hecho por A2.
- [x] Podado el ref stale `origin/docs/followup-pause-a2.6` (borrado del remoto).
      **Conservado** `perf/refine-kepler-cache` (FROZEN_RUN lo cita como provenance
      `06de6d0`, verificado presente).
- [x] `SCIENTIFIC_AUDIT.md` ya no existe en el repo (fue removido; sólo aparece en
      mensajes de commit). El tracking vivo de bloqueantes está en ROADMAP.md /
      FROZEN_RUN.md, que quedaron actualizados (bug #2 → cuantificado; nota de cierre
      del current_working_plan con PRs #53–#56).

---

## Apéndices

### A. Dependencias entre etapas

```
Track A (camino crítico del catálogo publicable):
  A1 prefilter recall ──► cierra el caveat de completitud
  A2 characterize 72M ──► observabilidad sobre el catálogo completo (refactor)

Track B (cierre del proyecto, Fase 7):
  B1 validación literatura ──► B2 dashboard + README final

Track C (opcionales, independientes):
  C1 perf · C2 falsos negativos Kepler · C3 bookkeeping
```

### B. Qué NO está acá (cerrado / fuera de alcance)

- **Capa de masas (Track A del follow-up anterior)**: CERRADA. No determinable en
  DR3. Sólo reabre con DR4/FPR. Ver [docs/mass_layer_track_a_closure.md](docs/mass_layer_track_a_closure.md).
- **Refinamiento N-body universal del catálogo (Stage C)**: descartado por costo
  (~400 días-CPU). El híbrido sobre el 12 % crítico es la solución adoptada.

### C. Referencias para arranque

- [ROADMAP.md](ROADMAP.md) § "Estado actual" y "Bugs / limitaciones" — verdad viva del proyecto.
- [FROZEN_RUN.md](FROZEN_RUN.md) — alcance defensible del catálogo congelado.
- [src/detect/prefilter.py](src/detect/prefilter.py) — Track A1.
- [src/characterize/](src/characterize/) — Track A2.
- [docs/kepler_threshold_bias_paper.md](docs/kepler_threshold_bias_paper.md) — Track C2.

### D. Bitácora de cambios al plan

| Fecha | Cambio | Autor |
|-------|--------|-------|
| 2026-05-30 | Plan creado al disolver `FOLLOWUP_PLAN.md` (capa de masas cerrada; info consolidada en ROADMAP + docs). | DF |
| 2026-05-30 | **Track A Stage 1 DONE**: recall del prefiltro medido (76.4 % en cola adversa, 143k encuentros perdidos por `|Δa|≤0.5`); fix radial-overlap → 100 %. Caveat #2 cuantificado en FROZEN_RUN. | DF |
| 2026-05-30 | **Track B Stage 1 DONE**: gate 4 cuerpos como test de regresión + Fuentes-Muñoz 2025 (AJ 170,353) Tabla 5 de fuente oficial → 11,804/40,004 pares (29.5 %) confirmados en el catálogo DR3. Goffin documentado como imposible (VizieR sin tabla de encuentros). | DF |
| 2026-05-31 | **Track A Stage 2 DONE**: refactor streaming de caracterización (`characterize_catalog_streaming`); corrida sobre el híbrido 72.2M sin OOM (73 chunks, 18.9 % observables, gate 4/4); tests de paridad. | DF |
| 2026-05-31 | **Track B Stage 2 (parcial)**: data layer memory-safe del dashboard (`src/dashboard/data.py`, usa el 72M), README final (validación corregida + recall prefiltro + reproducibilidad). Falta QA visual humana del dashboard. | DF |
| 2026-05-31 | **Track C Stage 3 DONE**: bookkeeping — audit #6 resuelto por cierre en ROADMAP, bug #2 actualizado a "cuantificado", nota de cierre PRs #53–#56; podado `origin/docs/followup-pause-a2.6` (conservado `perf/refine-kepler-cache`). | DF |
| 2026-05-31 | **Track C Stage 2 DONE**: medida la tasa de falsos negativos del threshold Kepler (0.70 % en [0.05,0.06) sobre 17.469 pares N-body); matriz de confusión cerrada; nota actualizada (RNAAS standalone viable, opcional). | DF |
| 2026-05-31 | **Track C Stage 1 DONE**: perf followups evaluados y descartados con mediciones (numba innecesario sobre solver vectorizado + paralelo por procesos; cache fine-grid de bajo valor). `docs/perf_evaluation.md` + `scripts/bench/bench_kepler.py`. **Plan completo.** | DF |
