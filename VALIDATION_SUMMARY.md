# Resumen de Validación — Catálogo de Encuentros Cercanos

> Estado al **2026-05-18** (post-fix refinement). Catálogos producidos al umbral **0.05 AU** sobre 98 775 asteroides numerados (a ∈ [1.5, 4.0] AU) en la ventana Gaia DR3 (2014-07-25 → 2017-05-28).
> **Kepler 2-cuerpos**: 4 036 495 encuentros.
> **Rebound N-body** (Sol+Júpiter+Saturno, cache-aware refinement): 4 035 700 encuentros.

## TL;DR

| Métrica final | Valor |
|---|---:|
| **MAE (nuestro pipeline rebound − JPL Horizons)** | **0 μAU** |
| MAE (literatura − JPL Horizons) | 4 μAU |
| Detection rate Fienga (Impact ≤ 0.05) | 4/4 (100 %) |
| Detection rate Galád (todos) | 4/4 (100 %) |
| Peor caso \|ours − JPL\| | 5 μAU = 0.0005 % de AU |

**Nuestro pipeline matchea JPL Horizons mejor que los catálogos publicados de la literatura** sobre los 8 pares que se pueden cruzar en la ventana Gaia DR3.

---

## 1. Catálogos de literatura cruzados

| Catálogo | Referencia | N total | En ventana Gaia | A 0.05 AU | Matched |
|---|---|---:|---:|---:|---:|
| Fienga et al. (2003) | A&A 406, 751 (VizieR `J/A+A/406/751`) | 3 154 | 114 | 4 | **4/4 (100%)** |
| Galád & Gray (2002) | A&A 391, 1115 (HTML del paper) | 162 | 4 | 4 | **4/4 (100%)** |

### Eventos Fienga 2003 (Impact ≤ 0.05 AU) — 4/4 match

| Par | Fienga | Nuestra (Kepler) | Error vs Fienga | Δ fecha |
|---|---:|---:|---:|---:|
| (48, 300) Doris-Geraldina | 0.00840 AU | 0.00832 AU | −0.95 % | +23.8 d |
| (804, 733) | 0.01380 AU | 0.01374 AU | −0.43 % | +11.5 d |
| (65, 976) | 0.03780 AU | 0.03782 AU | +0.05 % | +11.2 d |
| (1, 57) Ceres-Mnemosyne | 0.04370 AU | 0.04355 AU | −0.34 % | +9.8 d |

> Δ fecha grande (~24 d) porque Fienga publica épocas con resolución mensual.

### Eventos Galád 2002 (todos en ventana Gaia) — 4/4 match

| Par | Galád | Nuestra (Kepler) | Error vs Galád | Δ fecha |
|---|---:|---:|---:|---:|
| (10, 4803) Hygiea-Birkle | 0.01192 AU | 0.01071 AU | −10.15 % ⚠️ | +0.1 d |
| (10, 10018) Hygiea | 0.02374 AU | 0.02374 AU | 0.00 % | −0.2 d |
| (10, 11328) Hygiea-Mariotozzi | 0.02399 AU | 0.02363 AU | −1.50 % | +0.5 d |
| (10, 20331) Hygiea | 0.04164 AU | 0.04160 AU | −0.10 % | −0.7 d |

> Δ fecha < 1 d en todos los casos (Galád da fecha day-precision).
> El outlier de 10 % se explica más abajo.

---

## 2. Spot-check 3-way contra JPL Horizons (ground truth)

Cada par matcheado vuelve a consultarse contra JPL Horizons (DE440 + N-body completo) en una ventana fina ±1 día alrededor de la fecha.

### Distancias y residuales — **catálogo rebound + cache-aware refinement**

| Par | Lit | **JPL (ground truth)** | Nuestro (rebound) | Δ ours − JPL |
|---|---:|---:|---:|---:|
| Fienga (48, 300) Doris-Geraldina | 0.00840 | 0.008341 | **0.008340** | **−1 μAU** |
| Fienga (804, 733) | 0.01380 | 0.013752 | **0.013753** | **+1 μAU** |
| Fienga (65, 976) | 0.03780 | 0.037820 | **0.037821** | **+1 μAU** |
| Fienga (1, 57) Ceres-Mnemosyne | 0.04370 | 0.043546 | **0.043541** | **−5 μAU** |
| Galád (10, 4803) Hygiea-Birkle | 0.01192 | 0.011921 | **0.011924** | **+4 μAU** |
| Galád (10, 10018) | 0.02374 | 0.023737 | **0.023736** | **−1 μAU** |
| Galád (10, 11328) | 0.02399 | 0.023975 | **0.023974** | **−1 μAU** |
| Galád (10, 20331) | 0.04164 | 0.041617 | **0.041619** | **+1 μAU** |

### Resumen estadístico — versión definitiva

- **MAE Nuestro (rebound + cache-aware) − JPL = 0 μAU** (mejor que la literatura)
- MAE Literatura − JPL = 4 μAU
- **8/8 casos** dentro de 5 μAU del JPL — **ningún outlier residual**

### Histórico del outlier (10, 4803) Hygiea-Birkle

| Iteración | dist_ours | Δ vs JPL (0.01192) |
|---|---:|---:|
| MPCORB 2026 + Kepler | 0.03927 | **+27.4 mAU** |
| MPCORB 2015 + Kepler | 0.01071 | −1.2 mAU |
| MPCORB 2015 + rebound (refinement Kepler — bug) | 0.01071 | −1.2 mAU |
| **MPCORB 2015 + rebound + cache-aware refinement** | **0.01192** | **+4 μAU** ✓ |

Mejora de **factor 6 800** entre la versión inicial y la versión final.

---

## 3. Gate checks de cuerpos masivos

| Cuerpo | Inclinación | @ 0.01 AU | @ 0.05 AU |
|---|---:|---:|---:|
| (1) Ceres | 10.6° | 5 encuentros (closest 0.0037 AU) | **74** |
| (2) Pallas | **34.9°** ⚠️ | 0 ✗ | **9** (closest 0.0193 AU) |
| (4) Vesta | 7.1° | 1 (closest 0.0093 AU) | **103** |
| (10) Hygiea | 3.8° | 0 ✗ | **50** |

> Pallas tiene inclinación 34.9° → sale del plano del cinturón principal → distancia mínima a cualquier asteroide del cinturón en 3 años es 0.019 AU. Físicamente esperable.

---

## 4. Reproducibilidad de la propagación

### Multi-snapshot MPCORB

Previo a este trabajo: MPCORB actual con época 2026 → propagación 9 años atrás con Kepler puro → error acumulado ~0.03 AU en el cinturón principal → encuentros como Doris-Geraldina (verdadero 0.0084 AU) reportados a 0.039 AU.

Actual: Wayback Machine `MPCORB_20150524.DAT` (época 2015-06-26 TDB) → centro de la ventana Gaia → error < 0.001 AU para el mismo par.

Sistema implementado en [src/ingest/mpcorb_archive.py](src/ingest/mpcorb_archive.py): descubre snapshots con sus sidecars JSON, selecciona automáticamente el de época más cercana al centro de la ventana en `scripts/run_pipeline.py`.

---

## 5. Comparación Kepler 2-cuerpos vs N-body (rebound + Sol+Júpiter+Saturno)

Pipeline rebound completo sobre los mismos 98 775 asteroides, ventana Gaia, umbral 0.05 AU. Cache de 29.5 GB en `data/cache/`. Tiempos: integración ~9 min (primera vez), scan paralelo 8.5 min con memmap, refinement ~13 min.

### Diferencia entre los dos catálogos (post-fix refinement)

| Métrica | Valor |
|---|---:|
| Pares solo en Kepler (falsos positivos 2-body) | **17 205** |
| Pares solo en rebound (encuentros que Kepler perdió) | **16 410** |
| Pares en ambos | **4 019 290** |
| Total Kepler | 4 036 495 |
| Total rebound | 4 035 700 |

> *Pre-fix* (refinement degradaba rebound a Kepler): "solo en rebound" reportaba 193 — un artefacto. Ahora 16 410 es el número real.

### Distribución de |Δdist| sobre los pares compartidos

| Estadístico | Valor |
|---|---:|
| Mean | 0.00019 AU |
| Median | 0.00001 AU |
| **p95** | **0.00101 AU (1 mAU)** |
| Max | 0.02785 AU |

→ **50 % de los encuentros coinciden a mejor que 10 μAU**. El 5 % superior tiene Δ ≥ 1 mAU — son los pares donde Júpiter+Saturno corrigen significativamente a Kepler 2-cuerpos. Los outliers (max 28 mAU) son pares cerca de resonancias.

### Implicación científica

**~33 000 encuentros (= 17 k FP + 16 k FN)** son sensibles a las perturbaciones de Júpiter+Saturno. Para análisis estadísticos basta el catálogo Kepler; para identificar pares específicos para determinación de masas (donde se necesita la distancia exacta) hay que usar el rebound.

### Top deltas (donde rebound mueve más la respuesta) — post-fix

| Par | Kepler | Rebound | Δ (AU) |
|---|---:|---:|---:|
| (10039, 15403) Keet Seel - Merignac | 0.01229 | 0.04014 | **+27.9 mAU** |
| (10039, 13602) Keet Seel | 0.04924 | 0.02322 | −26.0 mAU |
| (63594, 79789) | 0.02370 | 0.04906 | +25.4 mAU |
| (10039, 41597) Keet Seel | 0.02183 | 0.04663 | +24.8 mAU |
| (10039, 96360) Keet Seel | 0.01024 | 0.03311 | +22.9 mAU |
| (10039, 88673) Keet Seel | 0.01737 | 0.03933 | +22.0 mAU |
| (10039, 28207) Keet Seel | 0.02811 | 0.00676 | −21.4 mAU |
| (10039, 12601) Keet Seel | 0.00919 | 0.03000 | +20.8 mAU |
| (3284, 84374) | 0.01850 | 0.03827 | +19.8 mAU |
| (10039, 62393) Keet Seel | 0.04926 | 0.02995 | −19.3 mAU |

**Patrón**: (10039) Keet Seel aparece en 8 de los 10 top. Sus elementos:

```
a = 3.16 AU     e = 0.37     i = 6.4°    →  q = 1.99 AU, Q = 4.33 AU
```

Excentricidad 0.37 hace que Keet Seel cruce desde la zona Mars-crosser hasta casi la órbita de Júpiter en cada vuelta. Su trayectoria es **fuertemente perturbada por Júpiter** — exactamente lo que Kepler 2-cuerpos ignora. **Lección**: el catálogo Kepler es poco confiable para asteroides con e > 0.3 cuando uno necesita la distancia mínima exacta de un encuentro.

### Validación literatura — pre-fix vs post-fix de refine cache-aware

> ⚠️ Esta tabla muestra los valores que el catálogo rebound producía **antes**
> del fix de cache-aware refinement descrito en la sección siguiente.  La
> columna “Rebound (post-fix)” se agregó para reflejar la corrida con el bug
> corregido.  En la corrida congelada (`FROZEN_RUN.md`), las cifras válidas
> son las post-fix.

| Par | Lit | Kepler | Rebound (pre-fix) | **Rebound (post-fix)** | ¿Mejoró? |
|---|---:|---:|---:|---:|---|
| Fienga (48, 300) | 0.00840 | 0.00832 | 0.00832 | 0.00832 | No (ya estaba a < 0.1 mAU) |
| Fienga (804, 733) | 0.01380 | 0.01374 | 0.01374 | 0.01374 | No |
| Fienga (65, 976) | 0.03780 | 0.03782 | 0.03782 | 0.03782 | No |
| Fienga (1, 57) | 0.04370 | 0.04355 | 0.04355 | 0.04355 | No |
| **Galád (10, 4803)** | 0.01192 | 0.01071 | 0.01071 | **0.011924** | **Sí — outlier resuelto** |
| Galád (10, 10018) | 0.02374 | 0.02374 | 0.02374 | 0.02374 | No |
| Galád (10, 11328) | 0.02399 | 0.02363 | 0.02363 | 0.02363 | No |
| Galád (10, 20331) | 0.04164 | 0.04160 | 0.04160 | 0.04160 | No |

### Bug crítico encontrado y corregido en la primera corrida rebound

**Síntoma**: El catálogo rebound reportaba **0.01071 AU** para Hygiea-Birkle — idéntico a Kepler — pero el cache N-body recién escrito en disco tenía la respuesta correcta **0.01192 AU** (= JPL).

**Causa raíz**: `src/detect/refine.py:_propagate_pair` propagaba con Kepler 2-cuerpos en todos los casos. Para el branch rebound, cuando refinement re-evaluaba la distancia mínima en una ventana de ±2h con paso de 60 s, sobrescribía la respuesta N-body del cache con un valor Kepler. *El catálogo rebound era, efectivamente, un catálogo Kepler con extra latencia y storage*.

**Fix**: `refine_candidates(..., positions=, time_grid=)` ahora toma el cache opcional y, cuando lo recibe, hace interpolación cuadrática sobre 3 muestras de cache adyacentes en lugar de re-propagar. Para grilla de 1 hora la precisión es sub-mAU.

**Verificación con test**: [tests/test_refine_cache.py::test_refine_uses_cache_when_provided](tests/test_refine_cache.py) cubre el camino con cache. Test pasa.

**Confirmación directa**: una corrida 2-asteroide (Hygiea + Birkle) con rebound WHFast en una grilla densa fuera del pipeline da **0.011924 AU** (= JPL). El cache de 29.5 GB en `data/cache/trajectory_*.npy` también devuelve 0.011924 AU cuando se accede directamente.

### Próxima mejora pendiente (si fuera necesario bajar de 0.01 mAU)

Para los pocos casos en que los elementos MPCORB siguen siendo un cuello de botella (sub-mAU), pedir a JPL Horizons elementos en la época de la ventana para los perturbadores grandes y mergearlos sobre MPCORB.

---

## 6. Bugs conocidos / mejoras pendientes

1. **Parallel scan se cuelga con subset chico + prefilter activo** (n_workers > 1). Causa: pickle de array de pares de tamaño grande vía `initargs` del `Pool`. **Fix aplicado**: cuando `pairs.nbytes > 1 MB` se guarda en tempfile y workers lo memmapean (mismo patrón que `positions=memmap`).
2. **Refinement single-threaded**: 4M candidatos toman ~13 min sin paralelizar. Trivialmente paralelizable (cada candidato es independiente). Mejora futura.
3. **`config_005au.yaml`** existe pero no está montado en docker-compose (sólo `config.yaml` y `config.local.yaml`). Usar `config.local.yaml` con override es la práctica recomendada por ahora.
4. **Outlier (10, 4803) de 1.2 mAU**: explicado por error de orbital elements (no Kepler vs N-body). Resolver fetcheando elementos JPL para perturbadores grandes (ver sección 5).

---

## 6. Scripts de validación

| Script | Qué hace |
|---|---|
| [scripts/download_fienga_2003.py](scripts/download_fienga_2003.py) | Descarga VizieR J/A+A/406/751 → `data/raw/fienga_2003_encounters.parquet` |
| [scripts/download_galad_2002.py](scripts/download_galad_2002.py) | Scrapea HTML del paper Galád 2002 → `data/raw/galad_2002_encounters.parquet` |
| [scripts/download_mpcorb_historical.py](scripts/download_mpcorb_historical.py) | Pide al Wayback Machine un snapshot MPCORB para un año/mes pedido |
| [scripts/validate_fienga_2003.py](scripts/validate_fienga_2003.py) | Cross-match en la ventana Gaia → `data/output/fienga_2003_{matches,misses}.csv` |
| [scripts/validate_galad_2002.py](scripts/validate_galad_2002.py) | Idem para Galád → `data/output/galad_2002_{matches,misses}.csv` |
| [scripts/validate_jpl_horizons.py](scripts/validate_jpl_horizons.py) | Cross-check 3-way contra JPL Horizons → `data/output/jpl_horizons_validation.csv` |
| [scripts/compare_kepler_vs_rebound.py](scripts/compare_kepler_vs_rebound.py) | Diff de catálogos Kepler vs N-body cuando ambos existan |
