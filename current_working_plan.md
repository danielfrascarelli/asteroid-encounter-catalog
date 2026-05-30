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
**Última actualización**: 2026-05-30
**Estado**: 🟡 **EN PROGRESO**. **A1 (recall del prefiltro) cerrada** en
`trackA/stage1-prefilter-recall` (recall=76.4 % en la cola adversa; fix
radial-overlap → 100 %; caveat #2 cuantificado). Siguiente recomendado: B1/B2
(Fase 7) o A2 (caracterización 72M). Toda la capa de masas sigue cerrada (no
determinable en DR3); este plan es sobre el **catálogo** y el cierre del proyecto.

| Track | Etapa | Estado | Branch | PR | Notas |
|-------|-------|--------|--------|----|-------|
| A | 1: recall del prefiltro (audit #2) | 🟢 DONE | `trackA/stage1-prefilter-recall` | (pendiente) | Recall=76.4% en cola adversa; fix radial-overlap → 100%. Cuantificado, no cerrado |
| A | 2: caracterización catálogo 72M (streaming) | ⚪ PENDING | — | — | Refactor; el catálogo caracterizado actual es solo 158k filas |
| B | 1: validación literatura completa (Fase 7) | ⚪ PENDING | — | — | Goffin 2014, Fuentes-Muñoz 2024 sobre el catálogo congelado |
| B | 2: dashboard + README final + reproducibilidad | ⚪ PENDING | — | — | `src/dashboard/app.py` ya existe; falta pulido + README |
| C | 1: perf followups (numba / cache persistente) | ⚪ PENDING | — | — | No bloqueante; refinement ya en meseta ~4.5-10× |
| C | 2: experimento falsos negativos threshold Kepler | ⚪ PENDING | — | — | Gatillo para una nota standalone (decisión B3) |
| C | 3: bookkeeping (audit #6, refs stale) | ⚪ PENDING | — | — | Cerrar #6 (respondido por cierre Track A); podar refs |

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

**Estado**: ⚪ PENDING
**Estimación**: ~2-4 días (refactor + corrida)
**Branch propuesta**: `trackA/stage2-characterize-bigcatalog`
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

- [ ] Refactor streaming en `src/characterize/`
- [ ] `data/output/encounters_characterized_full.parquet` (72M)
- [ ] Tests de paridad
- [ ] Nota en [FROZEN_RUN.md](FROZEN_RUN.md) (las columnas de observabilidad ahora cubren el catálogo completo)

#### Criterios de aceptación

- Corre sobre 72M sin OOM (pico de RAM acotado y documentado).
- Paridad con el run de 158k sobre el solapamiento.

#### Cómo retomar

— No arrancada —

#### Progreso

— No arrancada —

---

## Track B — Cierre del proyecto (Fase 7 del roadmap)

### Stage 1 — Validación contra literatura completa

**Estado**: ⚪ PENDING
**Estimación**: ~2-3 días
**Branch propuesta**: `trackB/stage1-literature-validation`
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

- [ ] `scripts/validate/validate_goffin_2014.py` extendido / `validate_fuentes_munoz_2024.py`
- [ ] `data/output/literature_validation/` (matches, misses)
- [ ] `docs/literature_validation.md`
- [ ] Test de regresión `tests/test_validation.py` (gate de 4 cuerpos)

#### Criterios de aceptación

- Match rate reportado honestamente para Goffin/Fuentes-Muñoz.
- Los 4 cuerpos grandes presentes (gate del CLAUDE.md).

#### Cómo retomar

— No arrancada —

#### Progreso

— No arrancada —

---

### Stage 2 — Dashboard, README final y reproducibilidad

**Estado**: ⚪ PENDING
**Estimación**: ~2-3 días
**Branch propuesta**: `trackB/stage2-dashboard-readme`
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

- [ ] Dashboard pulido (`src/dashboard/app.py`) levantando en `http://localhost:8501`
- [ ] `README.md` final
- [ ] Checklist de reproducibilidad verificado

#### Criterios de aceptación

- `docker compose up dashboard` levanta y muestra el catálogo.
- README permite reproducir desde cero.

#### Cómo retomar

— No arrancada —

#### Progreso

— No arrancada —

---

## Track C — Opcionales / deuda técnica

### Stage 1 — Perf followups (no bloqueante)

**Estado**: ⚪ PENDING · **Branch**: `trackC/stage1-perf`

El refinement ya está en su meseta (~4.5-10× sobre baseline; ver memoria
`project_perf_followups`). Quedan optimizaciones opcionales: **numba** en el
hot-path de Kepler (`src/propagate/kepler.py`), **cache persistente** de
propagaciones entre runs. Sólo si una corrida futura grande lo justifica.

- [ ] Evaluar numba en `solve_kepler` / `_refine_chunk_arr`
- [ ] Cache persistente de fine-grid entre runs

---

### Stage 2 — Experimento de falsos negativos del threshold Kepler

**Estado**: ⚪ PENDING · **Branch**: `trackC/stage2-kepler-false-negatives`
**Depende de**: nada.

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

**Estado**: ⚪ PENDING · **Branch**: `trackC/stage3-bookkeeping`

- [ ] **Cerrar audit blocker #6 (mass fitting)**: ya quedó **respondido** por el
      cierre de Track A (capa no viable en DR3). Marcarlo resuelto en el tracking,
      no como abierto. (Memoria `audit_followups`.)
- [ ] Podar refs/branches stale: `perf/refine-kepler-cache` (su contenido ya está
      en `main`; **conservar el ref** porque FROZEN_RUN lo cita como provenance
      `06de6d0`), `origin/docs/followup-pause-a2.6`.
- [ ] Revisar `SCIENTIFIC_AUDIT.md` y cerrar formalmente los bloqueantes resueltos.

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
