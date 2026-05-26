# Deepwork plan — N-body refinement + mass layer redesign

> Plan de trabajo profundo para resolver las dos limitaciones científicas
> identificadas en el audit round 5 que NO se pueden cerrar con un commit:
>
> - **Track 1**: el catálogo final es Kepler 2-cuerpos, no N-body.
> - **Track 2**: la capa de masas no es publicable (χ²_red ≈ 425, specificity 0/41).
>
> Este documento es el **único registro vivo** del avance. Cualquier
> retoma de trabajo debe empezar leyendo este archivo, no la conversación
> de Claude.

---

## Cómo usar este documento

**Para arrancar / retomar trabajo cold**:

1. Leer la sección "Estado global" (abajo) → te dice qué etapa está activa.
2. Saltar a la etapa marcada `🟡 IN PROGRESS` o la siguiente `⚪ PENDING`.
3. Leer la sección "Cómo retomar" de esa etapa específica.
4. Confirmar el estado contra el filesystem (los entregables esperados).
5. Avanzar.

**Para pausar**:

1. Commitear el WIP en su branch (`wip(track-N-stageX): <qué quedó hecho>`).
2. Actualizar el "Progreso" de la etapa activa: fecha, qué quedó hecho,
   qué bloquea (si bloquea), próximo paso concreto.
3. Mover estado a `🟠 PAUSED` si se va a dejar por >1 semana, sino `🟡 IN PROGRESS`.
4. Push de la branch (no merge a main hasta que la etapa esté `🟢 DONE`).

**Convenciones de estado**:

- ⚪ `PENDING` — no arrancada
- 🟡 `IN PROGRESS` — trabajo activo, branch abierta
- 🟠 `PAUSED` — branch abierta, sin actividad reciente, hay un blocker o cambio de prioridad
- 🔴 `BLOCKED` — no se puede avanzar sin input externo (datos, decisión)
- 🟢 `DONE` — PR mergeado a main, criterios de aceptación verificados

**Branches**:

Cada etapa va en su propia branch con el formato `track{N}/stage{X}-{slug}`.
Ejemplo: `track1/stageA-kepler-vs-nbody-error`. **Nunca trabajar en `main`**.

---

## Estado global

**Fecha de creación**: 2026-05-26
**Última actualización**: 2026-05-26
**Etapa activa**: Track 1 / Stage A (🟡 IN PROGRESS — A.1–A.4 listos, falta commit + PR)
**Próxima etapa a arrancar**: revisar PR de Stage A → arrancar Stage B (refinamiento selectivo) con el criterio `q_min < 1.8 ∨ e_max > 0.3`

| Track | Etapa | Estado | Branch | PR | Inicio | Fin | Notas |
|-------|-------|--------|--------|----|--------|-----|-------|
| 1     | A: caracterizar error Kepler vs N-body | 🟡 IN PROGRESS | `track1/stageA-kepler-vs-nbody-error` | #31 | 2026-05-26 | — | A.1–A.4 completos. p99 \|Δdist\| = 2.5 mAU sobre 964 pares; el error escala con e y 1/q. Falta merge. |
| 1     | B: refinamiento N-body selectivo | ⚪ PENDING | — | — | — | — | Depende de A |
| 1     | C: refinamiento N-body universal | ⚪ PENDING | — | — | — | — | Probablemente no se hace; depende de B |
| 2     | 1: joint fit órbita + masa | ⚪ PENDING | — | — | — | — | — |
| 2     | 2: covarianza Gaia AL | ⚪ PENDING | — | — | — | — | Puede hacerse en paralelo con 2.1 |
| 2     | 3: specificity test riguroso | ⚪ PENDING | — | — | — | — | Depende de 2.1 |
| 2     | 4: validación contra masas conocidas | ⚪ PENDING | — | — | — | — | Gate antes de publicar |

**Recomendación de orden**: Track 1 Stage A primero (es la más barata y decide
si vale la pena seguir con Track 1). Después Track 2 en paralelo si A muestra
que Kepler-refine es defendible.

---

## Track 1 — Refinamiento N-body

### Contexto del problema

El catálogo congelado (`encounters_catalog_rebound_005au.parquet`, 72 M rows,
SHA-256 `b0272be7…`) usa rebound (whfast, Sol + Júpiter + Saturno, dt = 1 h)
para el **coarse scan** del KD-tree, pero el **refinamiento sub-grid** que
produce la distancia mínima reportada corre `kepler_to_cartesian` 2-cuerpos
en una ventana ±2 h con dt = 120 s. Ver
[src/detect/pipeline.py:166-183](src/detect/pipeline.py#L166-L183) y
[FROZEN_RUN.md:14-22](FROZEN_RUN.md#L14-L22).

Consecuencia honesta: las distancias mínimas, épocas y velocidades relativas
del catálogo son **valores geométricos bajo Kepler 2-cuerpos**, no soluciones
gravitacionales completas. Pares en cruzadores de órbitas cerca de resonancias
o con perihelios bajos pueden discrepar de N-body en mAU.

### Stage A — Caracterizar el error Kepler vs N-body

**Estado**: 🟡 IN PROGRESS (A.1–A.4 done, falta merge)
**Branch**: `track1/stageA-kepler-vs-nbody-error`
**PR**: pendiente
**Estimación original**: 1 semana — **real**: ~1 día (sample muy paralelizable)
**Bloquea a**: Stage B, decisión de seguir con Track 1

#### Objetivo

Cuantificar, sobre una muestra representativa, cuánto se mueven
`dist_min`, `t_min` y `rel_vel` al pasar de refinamiento Kepler a
refinamiento N-body. Esto convierte la afirmación cualitativa
"Kepler puede discrepar de N-body en mAU" en una **distribución de error
medida** que se puede citar.

#### Por qué arrancar por acá

Antes de gastar semanas en Stage B (refinamiento N-body selectivo) o
descartar Stage C (universal), necesitamos saber:

1. ¿Cuál es el error mediano? ¿Es <0.0005 AU? Si sí, Kepler-refine es
   defendible para la mayoría del catálogo.
2. ¿Cuál es la cola peor? ¿Hay un subset claramente identificable
   (e.g., e > 0.3, |Δi| > 15°, perihelio < 1.3 AU) donde Kepler es
   inaceptable?
3. Si la cola peor es chica (<1% del catálogo), Stage B es suficiente.
4. Si la cola peor es grande, hay que considerar Stage C o descartar
   la pretensión de precisión N-body global.

#### Plan técnico

**A.1 — Muestrear ~1000 candidatos del catálogo congelado** (½ día)

- Script nuevo: `scripts/validate/sample_for_nbody_check.py`.
- Lee `data/output/encounters_catalog_rebound_005au.parquet`.
- Sampleo estratificado por (a_1, e_1, i_1, dist_au, |Δa|):
  - 200 bins, ~5 candidatos por bin.
  - Buscar cobertura uniforme del espacio orbital, no del espacio
    "más probable" — queremos detectar la cola peor.
- Output: `data/cache/nbody_validation/sample_1000.parquet` con columnas
  del catálogo original más `bin_id`.
- **Criterio**: cada bin orbital con ≥3 candidatos; rango de e cubierto
  hasta 0.7; rango de i hasta 25°.

**A.2 — Refinador N-body por par** (1.5 días)

- Script nuevo: `scripts/validate/refine_pair_nbody.py`.
- Para un par dado `(number_1, number_2, jd_tdb_kepler, dist_au_kepler)`:
  1. Cargar elementos MPCORB de ambos cuerpos en su época.
  2. Crear simulación rebound: Sol + Júpiter + Saturno como cuerpos
     masivos, asteroide_1 y asteroide_2 como test particles.
     Considerar incluir Ceres/Vesta/Pallas/Hygiea como cuerpos masivos
     si el par está en el cinturón principal (controla con flag).
  3. Integrar desde la época MPCORB hasta `t_min_kepler − 6 h` con whfast,
     dt = 600 s (warmup).
  4. Cambiar a `IAS15` (adaptativo, alto orden) y muestrear posiciones
     cada 60 s desde `t_min_kepler − 6 h` hasta `t_min_kepler + 6 h`.
  5. Calcular `dist(t)` entre los dos asteroides, ajustar parábola en el
     mínimo para sub-paso, devolver `dist_min_nbody`, `t_min_nbody`,
     `rel_vel_nbody`.
- **Frame**: heliocéntrico eclíptico J2000, igual que el resto del pipeline.
- **Tests** en `tests/test_refine_pair_nbody.py`:
  - Caso conocido: (1) Ceres vs un par cercano — `dist_min_nbody` dentro
    de 0.0001 AU de JPL Horizons.
  - Conservación de energía relativa <1e-9 después de la integración.

**A.3 — Correr el comparador sobre los 1000 candidatos** (1 día)

- Script: `scripts/validate/compare_kepler_vs_nbody.py`.
- Lee `sample_1000.parquet`.
- Por cada par: llama al refinador N-body (A.2), captura
  `(dist_min_nbody, t_min_nbody, rel_vel_nbody)`.
- Paralelizar con `multiprocessing.Pool` (CPU-bound, rebound es C bajo
  el capó; 28 workers debería bajar 1000 pares a <30 min).
- Output: `data/output/kepler_vs_nbody_comparison.parquet` con
  columnas:
  - `number_1`, `number_2`, `bin_id`
  - `dist_au_kepler`, `dist_au_nbody`, `delta_dist_au`
  - `t_min_kepler_jd`, `t_min_nbody_jd`, `delta_t_min_hours`
  - `rel_vel_kepler`, `rel_vel_nbody`, `delta_rel_vel_au_day`
  - `a_1`, `e_1`, `i_1`, `q_1` (perihelio), idem `_2`
  - `nbody_converged` (bool), `nbody_energy_drift`

**A.4 — Análisis y reporte** (1 día)

- Notebook nuevo: `notebooks/nbody_error_characterization.ipynb`.
- Histogramas:
  - `delta_dist_au` global + por bin de e, i, |Δa|, q_min.
  - `delta_t_min_hours` global + por bin.
- Mapa 2D `(e_max, i_max)` → `delta_dist_au_p95` para identificar la
  cola peor.
- Tabla de candidatos con `|delta_dist_au| > 0.001` listados con sus
  elementos orbitales — buscar patrón.
- **Documento**: `docs/kepler_refine_error_report.md` con:
  - Mediana, p95, p99 del error en distancia.
  - Mediana, p95, p99 del error en tiempo.
  - Identificación del subset peor (criterio orbital explícito).
  - Recomendación: ¿hace falta Stage B? ¿Sobre qué subset?

#### Entregables

- [ ] `scripts/validate/sample_for_nbody_check.py`
- [ ] `scripts/validate/refine_pair_nbody.py`
- [ ] `scripts/validate/compare_kepler_vs_nbody.py`
- [ ] `tests/test_refine_pair_nbody.py` (≥3 tests, incl. cross-check JPL Horizons)
- [ ] `data/cache/nbody_validation/sample_1000.parquet`
- [ ] `data/output/kepler_vs_nbody_comparison.parquet`
- [ ] `notebooks/nbody_error_characterization.ipynb`
- [ ] `docs/kepler_refine_error_report.md`
- [ ] Update a `FROZEN_RUN.md` § "Scope and limits" con el error medido (no más "puede discrepar en mAU" cualitativo).

#### Criterios de aceptación

- Comparación corre en <1 h sobre 1000 pares en una máquina con 28 cores.
- Cross-check con JPL Horizons para 5 pares al azar: `|dist_min_nbody − dist_horizons| < 0.0001 AU`.
- El reporte responde explícitamente: ¿Stage B es necesario? ¿Sobre qué subset?
- PR mergeada a main.

#### Cómo retomar Stage A

```bash
git checkout track1/stageA-kepler-vs-nbody-error
ls scripts/validate/                     # qué scripts ya existen
ls data/cache/nbody_validation/          # qué sample ya está
ls data/output/kepler_vs_nbody*.parquet  # qué corridas ya están
cat docs/kepler_refine_error_report.md   # qué se concluyó
```

Luego seguir el sub-paso A.X que esté incompleto.

#### Progreso

| Sub-paso | Estado | Fecha | Comentario |
|---|---|---|---|
| A.1 sample | 🟢 done | 2026-05-26 | `scripts/validate/sample_for_nbody_check.py`; output `data/cache/nbody_validation/sample_1000.parquet` con 964 pares estratificados sobre `(a_mid, e_max, i_max, q_min, dist_au)` (binning simétrico; 200 bins, mínimo 3/bin enforced). Pool intermedio de 1.5M (--pool-size) para que las colas simétricas queden pobladas. |
| A.2 refiner | 🟢 done | 2026-05-26 | `scripts/validate/refine_pair_nbody.py` + `tests/test_refine_pair_nbody.py`. WHFast warmup desde MPCORB epoch → IAS15 en ±12h con muestreo de 60s, ajuste parabólico del mínimo. Bug crítico arreglado: el target ya no se duplica como perturber si coincide con Ceres/Pallas/Vesta/Hygiea. Pytest 5/5 passed (+ fixture Horizons offline). |
| A.3 comparison run | 🟢 done | 2026-05-26 | `scripts/validate/compare_kepler_vs_nbody.py`. 964/964 ok en ~7 s con 24 workers; ventana ±12h (la inicial ±6h truncaba la cola); max energy_drift 3.5e-14. Output: `data/output/kepler_vs_nbody_comparison.parquet`. |
| A.4 report | 🟢 done | 2026-05-26 | `docs/kepler_refine_error_report.md` + `notebooks/nbody_error_characterization.ipynb`. p99 \|Δdist\| = 2.5 mAU; recomendación Stage B = subset `q_min < 1.8 ∨ e_max > 0.3` (~20% del catálogo). FROZEN_RUN.md actualizado. |

**Resultados clave** (para informar Stage B):
- Mediana `|Δdist|` = 12 μAU, p95 = 678 μAU, p99 = 2.5 mAU, max = 11.3 mAU.
- El error escala con `e_max` (factor 12× entre e<0.10 y e>0.45) e inversamente con `q_min` (factor 10× entre q>2.6 y q<1.3).
- 3.4% near-boundary: el verdadero mínimo N-body podría estar fuera de ±12h; subestimación del error en esos casos.
- Ninguno de los 964 pares cambia status de detección al re-refinarse.

**Cómo retomar**:
1. Verde merge del PR de Stage A.
2. Arrancar Stage B usando el criterio `q_min < 1.8 ∨ e_max > 0.3` derivado en A.4.

---

### Stage B — Refinamiento N-body selectivo

**Estado**: ⚪ PENDING (depende de Stage A)
**Branch**: (no creada)
**Estimación**: 2 semanas
**Bloquea a**: nada estrictamente; mejora del catálogo

#### Objetivo

Refinar con N-body **solo** el subset identificado en Stage A como
"Kepler-refine no defendible". El resto del catálogo conserva Kepler.
Output: catálogo híbrido con columna `refinement_method ∈ {kepler, nbody}`
y, donde corresponda, `dist_au_nbody`, `t_min_nbody_jd`, `rel_vel_nbody`.

#### Plan técnico (preliminar — refinar al iniciar)

- Reusar `refine_pair_nbody.py` de Stage A, optimizado para batch.
- Pipeline: filtrar el catálogo por el criterio orbital de Stage A,
  correr N-body en paralelo, joinear de vuelta al parquet principal.
- Estimar volumen: si la cola peor es 1% del catálogo, son 720k pares.
  A 0.5 s por par × 28 workers ≈ 3.6 h. Si es 10%, son 36 h.
- Actualizar schema con las columnas nuevas, marcadas como opcionales.
- Update a `FROZEN_RUN.md`: nueva tabla TL;DR con el catálogo híbrido.

#### Criterios de aceptación (preliminar)

- Cada candidato del subset crítico tiene `dist_au_nbody` y
  `refinement_method = "nbody"`.
- Los demás conservan `dist_au` Kepler y `refinement_method = "kepler"`.
- Conservación de energía relativa <1e-9 en cada integración.
- Validación contra Goffin (2014) y Fienga (2003): los pares de literatura
  caen dentro de 1×10⁻⁴ AU del valor publicado.

#### Cómo retomar

(Se completa cuando Stage A termine y el plan se concrete con los
parámetros reales del subset crítico.)

#### Progreso

— No arrancada —

---

### Stage C — Refinamiento N-body universal

**Estado**: ⚪ PENDING (probablemente **no se hace**)
**Estimación**: 2-3 meses
**Decisión**: posponer salvo que Stage A/B muestren que es necesario para una claim que vale el costo.

#### Por qué probablemente no

Refinar N-body los 72 M pares del catálogo requiere o bien una integración
rebound global con dt fino y dump de posiciones (≈80 TB inviable), o bien
72 M integraciones individuales (≈400 días-CPU). Ninguna es razonable
salvo que el target sea un paper sobre **el pipeline mismo**.

Más realista: si Stage B no alcanza, considerar adoptar un dataset externo
(JPL SBDB close-approaches API, NEODyS-2) en lugar de regenerar.

#### Cómo retomar

Solo si Stage B se completa y el resultado no satisface — entonces
re-discutir scope con stakeholders.

---

## Track 2 — Capa de masas

### Contexto del problema

`scripts/mass/fit_mass_gaia_loo.py` y compañía:
- χ²_red mediano ≈ 425 (un fit correcto debería estar ≈ 1).
- Specificity test: 0/41 detecciones específicas vs nulls.
- Masas estimadas 100–10⁴× sobreestimadas vs literatura.

Diagnóstico (ver [encounter_analysis/DETECTIONS.md](encounter_analysis/DETECTIONS.md)):
el forward model actual ajusta **una sola masa M_perturber** a partir
de residuos astrométricos del target. Pero los residuos están dominados
por **drift orbital** del target (errores acumulados de `a`, `e`, `i`
desde MPCORB) y **sistemática observacional** Gaia (along-scan vs
across-scan), no por la deflección gravitacional del encuentro.

El modelo absorbe drift y sistemática como "masa" → señal espuria,
specificity nula.

### Stage 1 — Joint fit órbita + masa

**Estado**: ⚪ PENDING
**Branch**: (no creada)
**Estimación**: 3-4 semanas
**Depende de**: ninguna (puede empezar en paralelo con Track 1 Stage A)

#### Objetivo

Cambiar el forward model de **1 parámetro (M)** a **7 parámetros
(M, Δa, Δe, Δi, ΔΩ, Δω, ΔM₀)** del target, donde los seis Δ-elementos
absorben el drift orbital y solo lo correlacionado con la geometría del
encuentro va a M.

#### Plan técnico

**1.1 — Inventario del código actual** (½ día)

- Leer `scripts/mass/fit_mass_gaia_loo.py` y dependencias en `src/mass/`
  (si existe) o donde estén las funciones del forward model.
- Documentar en `docs/mass_layer_audit.md`:
  - Qué se está fittando (parámetros, datos, likelihood).
  - Qué supuestos hace (errores diagonales, sistema de coordenadas, etc.).
  - Punto exacto de inserción de los 6 parámetros extra.

**1.2 — Diseño del forward model extendido** (2 días)

- Definir cómo aplicar los Δ-elementos: perturbar los elementos MPCORB
  del target en su época, propagar, calcular residuos vs observaciones
  Gaia DR3.
- Decidir parametrización: ¿Δ absolutos o relativos (Δa/a)?
  Recomendación: relativos para mejor condicionamiento numérico.
- Priors: gaussianos centrados en 0 con σ = incertidumbre típica de
  MPCORB. Esto los hace una regularización Bayesiana, no parámetros
  libres puros.
- Documentar en `docs/mass_layer_design.md`.

**1.3 — Implementación** (1.5 semanas)

- Nuevo módulo: `src/mass/forward_model_joint.py`.
- Función `predict_residuals(params, obs, perturber_elements) -> ndarray`
  donde `params = (M, da, de, di, dOmega, domega, dM0)`.
- Optimizer: `scipy.optimize.least_squares` con `method='trf'` y
  `jac='2-point'` para Jacobiano numérico (analítico es overkill por ahora).
- Bounds: M ∈ [0, 10× masa Hill plausible]; Δ-elementos dentro de
  3σ del prior.
- Output: `(M_fit, M_err, da_fit, ..., chi2_red, n_obs)`.

**1.4 — Wrapper LOO actualizado** (3 días)

- Refactor de `scripts/mass/fit_mass_gaia_loo.py` para usar el forward
  model joint.
- Output schema extendido: además de `M_fit`, agregar `da_fit`, `de_fit`,
  etc. (los seis Δ).
- Re-correr sobre los 41 candidatos.
- Output: `data/output/loo_batch_results_joint.csv`.

**1.5 — Diagnóstico** (2 días)

- Comparar χ²_red joint vs simple sobre los 41 candidatos.
- Si χ²_red joint ≈ 1 (o al menos < 10), la hipótesis "el drift se
  comía la señal" se confirma.
- Si sigue siendo >100, hay otro problema (covarianza, sistemática
  no modelada) — bloquear Stage 1 hasta entender.

#### Entregables

- [ ] `docs/mass_layer_audit.md`
- [ ] `docs/mass_layer_design.md`
- [ ] `src/mass/forward_model_joint.py` (con tests unitarios)
- [ ] `tests/test_forward_model_joint.py`
- [ ] `scripts/mass/fit_mass_gaia_loo_joint.py` (o flag `--joint` en el script existente)
- [ ] `data/output/loo_batch_results_joint.csv`
- [ ] Sección de comparación en `encounter_analysis/DETECTIONS.md`.

#### Criterios de aceptación

- χ²_red mediano del fit joint < 10 (idealmente ~1) sobre los 41 candidatos.
- Tests unitarios pasan: el forward model joint, dado parámetros
  conocidos, reproduce residuos sintéticos dentro de error numérico.
- PR mergeada con review del diagnóstico.

#### Cómo retomar

```bash
git checkout track2/stage1-joint-fit
cat docs/mass_layer_audit.md           # qué entendiste del modelo viejo
cat docs/mass_layer_design.md          # qué diseñaste
ls src/mass/forward_model_joint.py     # qué está implementado
docker compose run --rm test pytest tests/test_forward_model_joint.py -v
```

#### Progreso

| Sub-paso | Estado | Fecha | Comentario |
|---|---|---|---|
| 1.1 audit | ⚪ | — | — |
| 1.2 design | ⚪ | — | — |
| 1.3 implementation | ⚪ | — | — |
| 1.4 LOO wrapper | ⚪ | — | — |
| 1.5 diagnostic | ⚪ | — | — |

---

### Stage 2 — Covarianza Gaia AL

**Estado**: ⚪ PENDING
**Estimación**: 1 semana
**Depende de**: Stage 1 (puede arrancar en paralelo desde la sub-tarea de inventario)

#### Objetivo

Gaia tiene incertidumbres muy distintas según la dirección de escaneo:
along-scan (AL) ≪ across-scan (AC), tipicamente 10× más preciso.
El likelihood actual probablemente usa σ uniforme. Hay que extraer la
covarianza por observación y usar χ² Mahalanobis correcto.

#### Plan técnico (preliminar)

- Inventariar qué columnas de Gaia DR3 (o FPR) tienen errores AL/AC.
- Construir matriz de covarianza 2×2 por observación (rotada según
  el ángulo de escaneo).
- Reescribir el likelihood: `(r - r_pred)ᵀ Σ⁻¹ (r - r_pred)` por
  observación, sumar.
- Verificar que sobre datos sintéticos sin masa, el χ²_red baja a ~1
  cuando la covarianza es correcta.

#### Criterios de aceptación

- Likelihood Mahalanobis implementado en `src/mass/likelihood_al.py`.
- Test sintético: dado un dataset con errores AL/AC conocidos y sin
  perturbación, fit devuelve M ≈ 0 con χ²_red ≈ 1.
- Re-corrida de los 41 candidatos con joint fit + Mahalanobis: χ²_red
  no debería empeorar.

#### Cómo retomar / Progreso

— No arrancada —

---

### Stage 3 — Specificity test riguroso

**Estado**: ⚪ PENDING
**Estimación**: 1 semana
**Depende de**: Stage 1 (y preferiblemente Stage 2)

#### Objetivo

Para cada candidato real, generar N=100 "null encounters" geométricamente
compatibles (mismo target, mismo período Gaia, perturber elegido random
de la población compatible pero **sin** pasaje cercano real). Correr el
mismo fit y mostrar que la distribución de M_fit / χ²_real está
significativamente sesgada vs la distribución de nulls.

#### Plan técnico (preliminar)

- Función `generate_null_perturbers(target, n=100) -> list[Asteroid]`:
  - Mismo `a` ± 0.5 AU, misma clase orbital, sin acercamiento <0.1 AU al target en la ventana Gaia.
- Correr `fit_mass_gaia_loo_joint` sobre cada null.
- Métrica: p-value de M_fit_real vs distribución de M_fit_null.
- Output: `data/output/specificity_test_v2.csv` con columnas
  `target, perturber_real, p_value_M, p_value_chi2, n_nulls`.

#### Criterios de aceptación

- Para ≥3 candidatos de mayor score, p_value < 0.05.
- Si no, conclusión honesta: no hay señal de masa en el dataset actual.
  Documentar y decidir si seguir.

#### Cómo retomar / Progreso

— No arrancada —

---

### Stage 4 — Validación contra masas conocidas

**Estado**: ⚪ PENDING
**Estimación**: 1 semana
**Depende de**: Stage 1, 2, 3
**Gate**: antes de publicar cualquier masa nueva.

#### Objetivo

Re-fittear (1) Ceres, (4) Vesta, (2) Pallas, (10) Hygiea — masas con
1% de error en literatura — con el pipeline completo (joint fit +
Mahalanobis). Si el pipeline las reproduce dentro de 3σ, podemos
citar masas nuevas.

#### Plan técnico (preliminar)

- Identificar los target asteroides perturbados por cada uno en la
  ventana Gaia DR3 (ya hay candidatos en `mass_followup_candidates.csv`).
- Fittear y comparar contra:
  - Ceres: M = (4.71 ± 0.04)×10²⁰ kg (DAWN, 2015)
  - Vesta: M = (2.59 ± 0.01)×10²⁰ kg (DAWN, 2011)
  - Pallas: M = (2.05 ± 0.05)×10²⁰ kg (Goffin 2014)
  - Hygiea: M = (8.3 ± 0.4)×10¹⁹ kg (Vernazza+ 2020)
- Métrica: `(M_fit - M_lit) / sqrt(σ_fit² + σ_lit²)` < 3 → ✓.

#### Criterios de aceptación

- Al menos 3 de los 4 calibradores reproducidos dentro de 3σ.
- Reporte en `docs/mass_layer_validation.md`.

#### Cómo retomar / Progreso

— No arrancada —

---

## Apéndices

### A. Dependencias entre etapas

```
Track 1:
  Stage A ──► (decisión: vale Stage B?) ──► Stage B ──► (rara vez) Stage C

Track 2:
  Stage 1 ──┬──► Stage 3 ──┐
            │              ├──► Stage 4 (gate)
  Stage 2 ──┴──────────────┘

Cross-track:
  Stage A es independiente de todo Track 2 — pueden correr en paralelo.
  Si Stage B se hace, podría reemplazar los candidatos de Track 2 con
    distancias más precisas, lo cual revaluaría los inputs de Stage 1.4.
```

### B. Costo computacional estimado

| Etapa | Compute | Wall-clock (1 máquina 28 cores) |
|-------|---------|----------------------------------|
| T1 A  | bajo    | 1 día de trabajo + 30 min cómputo |
| T1 B  | medio   | 2 semanas trabajo + 4-36 h cómputo |
| T1 C  | enorme  | meses; descartado |
| T2 1  | medio   | 3-4 semanas trabajo + horas cómputo |
| T2 2  | bajo    | 1 semana trabajo |
| T2 3  | medio   | 1 semana trabajo + 1-2 días cómputo (41×100 fits) |
| T2 4  | bajo    | 1 semana trabajo |

### C. Decisiones aún no tomadas

- **Incluir Ceres/Vesta/Pallas/Hygiea como cuerpos masivos en el N-body
  de Stage A.2?** Default actual del pipeline: no. Para validación
  contra Goffin, sí debería incluirlos (Goffin los usa). Decidir al
  arrancar A.2.
- **Joint fit: priors gaussianos vs uniformes?** Default propuesto:
  gaussianos con σ del catálogo MPCORB. Decidir en Stage 1.2.
- **Mahalanobis: usar covarianza por observación o por sub-grupo de
  observaciones?** Decidir en Stage 2.

### D. Glosario rápido (para retomar cold)

- **Refinamiento Kepler/N-body**: paso final del detector que computa la
  distancia mínima exacta entre dos asteroides en una ventana ±2 h
  alrededor del candidato del KD-tree.
- **LOO (Leave-One-Out)**: técnica de fit donde se ajusta el modelo
  excluyendo una observación a la vez para estimar la robustez.
- **AL/AC (along-scan / across-scan)**: ejes de la geometría de escaneo
  de Gaia; AL es ~10× más preciso que AC.
- **Specificity**: fracción de candidatos donde el fit detecta la masa
  con significancia estadística vs la distribución de nulls.
- **χ²_red**: chi² reducido = χ² / (n_obs − n_params). Modelo bien
  especificado: ≈ 1. >> 1: modelo mal especificado o errores subestimados.

### E. Referencias para arranque

- [FROZEN_RUN.md](FROZEN_RUN.md) — alcance defensible del catálogo congelado.
- [ROADMAP.md](ROADMAP.md) § "Estado actual" — limitaciones listadas.
- [encounter_analysis/DETECTIONS.md](encounter_analysis/DETECTIONS.md) — diagnóstico actual de la capa de masas.
- [src/detect/pipeline.py](src/detect/pipeline.py) — donde Kepler-refine vive.
- [src/propagate/nbody.py](src/propagate/nbody.py) — propagador rebound ya implementado.
- Goffin (2014), A&A 565, A56 — masas de literatura para validación.
- DAWN papers — Ceres, Vesta ground truth.

---

## Bitácora de cambios al plan

| Fecha | Cambio | Autor |
|-------|--------|-------|
| 2026-05-26 | Plan creado tras audit round 5. | DF |
| 2026-05-26 | Stage A completa (A.1–A.4); PR #31. p99 \|Δdist\| = 2.5 mAU sobre 964 pares (estratificación simétrica); recomendación Stage B = subset (e_max>0.3 ∨ q_min<1.8). | DF |
