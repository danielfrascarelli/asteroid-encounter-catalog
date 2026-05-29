# Follow-up plan — Post-deepwork follow-on work

> Plan vivo para trabajo posterior al cierre del deepwork (PRs #31-#38).
> Cubre dos bloques:
>
> - **Track A**: tres direcciones para destrabar la capa de masas tras el
>   gate FAIL de Stage 4. Secuenciales y condicionales — el resultado de
>   cada una decide si arrancar la siguiente.
> - **Track B**: tres trabajos ortogonales pendientes (outliers Stage 2,
>   specificity completo, side-paper sobre el sesgo Kepler). Independientes
>   entre sí y de Track A.
>
> Este archivo es el **único registro vivo** del avance.
>
> Diagnósticos congelados que motivan este plan:
> [docs/mass_layer_validation.md](docs/mass_layer_validation.md),
> [docs/mass_layer_stage3_diagnostic.md](docs/mass_layer_stage3_diagnostic.md),
> [docs/mass_layer_stage2_diagnostic.md](docs/mass_layer_stage2_diagnostic.md),
> [FROZEN_RUN.md](FROZEN_RUN.md).

---

## Cómo usar este documento

**Para arrancar / retomar cold**:

1. Leer la sección "Estado global" — qué etapa está activa.
2. Saltar a la etapa marcada `🟡 IN PROGRESS` o la siguiente `⚪ PENDING`.
3. Leer la sección "Cómo retomar" de esa etapa específica.
4. Confirmar contra el filesystem (los entregables esperados).
5. Avanzar.

**Para pausar**:

1. Commitear el WIP en su branch (`wip(track-X-stageY): <qué quedó hecho>`).
2. Actualizar "Progreso" de la etapa activa: fecha, qué quedó, qué bloquea,
   próximo paso concreto.
3. Mover estado a `🟠 PAUSED` si se va a dejar por >1 semana, sino
   `🟡 IN PROGRESS`.
4. Push de la branch (no merge a `main` hasta etapa `🟢 DONE`).

**Convenciones de estado**:

- ⚪ `PENDING` — no arrancada
- 🟡 `IN PROGRESS` — trabajo activo, branch abierta
- 🟠 `PAUSED` — branch abierta, sin actividad reciente, blocker o cambio
  de prioridad
- 🔴 `BLOCKED` — no se puede avanzar sin input externo (datos, decisión)
- 🟢 `DONE` — PR mergeada a main, criterios verificados
- ⚫ `ABANDONED` — etapa descartada (decisión científica), con justificación

**Branches**: una por etapa, formato `track{A|B}/stage{N}-{slug}`.
Ejemplo: `trackA/stage1-tighten-priors`. **Nunca trabajar en `main`**.

---

## Estado global

**Fecha de creación**: 2026-05-29
**Última actualización**: 2026-05-29
**Etapa activa**: ninguna — plan recién creado.
**Próxima etapa a arrancar**: Track A Stage 1 (tighten priors) por ser el
experimento más barato y el que decide si Track 2 tiene futuro estructural.

| Track | Etapa | Estado | Branch | PR | Inicio | Fin | Notas |
|-------|-------|--------|--------|----|--------|-----|-------|
| A | 1: tighten priors | ⚪ PENDING | — | — | — | — | Gate decisión Track A |
| A | 2: multi-target joint fit | ⚪ PENDING | — | — | — | — | Solo si A1 no resuelve |
| A | 3: OU forward model para drift | ⚪ PENDING | — | — | — | — | Solo si A2 todavía deja bias |
| B | 1: investigar outliers Stage 2 (Alkeste/57942, Eros/176865) | ⚪ PENDING | — | — | — | — | Independiente; científico |
| B | 2: specificity sobre 22/27 restantes | ⚪ PENDING | — | — | — | — | Completar Stage 3 |
| B | 3: side-paper 25,283 false-positives Kepler | ⚪ PENDING | — | — | — | — | Posible mini-paper |

**Recomendación de orden**:

1. **Track A Stage 1 primero**. Es 1 día de trabajo y decide si el bias
   de Stage 4 es estructural (entonces Track A1 falla y hay que ir a A2/A3)
   o sólo es overfitting (entonces A1 lo arregla y Stage 4 pasa).
2. **Track B en paralelo si hay capacidad**. B1 y B2 son scope acotado y
   no dependen de Track A. B3 es opcional (side-paper, no urgente).

---

## Track A — Destrabar la capa de masas

### Contexto

El gate Stage 4 del deepwork ([docs/mass_layer_validation.md](docs/mass_layer_validation.md))
mostró que el pipeline joint+Mahalanobis subestima sistemáticamente las
masas literatura: ratios fit/lit = 0.77 (Ceres) / 0.57 (Pallas) / 0.24
(Hygiea), todos con |z| > 6. El bias escala con `1/M_real`. Mecanismo
hipotético: los 6 deltas orbitales del forward model joint absorben
parcialmente la deflección gravitatoria del encuentro, dejando una
masa fitted sesgada hacia abajo.

Las tres etapas de este track atacan el problema en niveles crecientes
de complejidad.

### Stage 1 — Tighten priors orbitales

**Estado**: ⚪ PENDING
**Estimación**: ~1 día
**Branch propuesta**: `trackA/stage1-tighten-priors`
**Depende de**: nada — arrancar cuando se reactive el track.

#### Objetivo

Apretar los priors gaussianos sobre los 6 Δ-elementos del joint fit por
un factor 5–10× y re-correr Stage 4. Si las masas calibradoras convergen
hacia el valor literatura (|z| < 3 sobre al menos 3 de 4), el bias es de
overfitting y Track A está resuelto. Si no, el bias es estructural y hay
que ir a Stage 2.

#### Por qué arrancar acá

Los priors actuales (en [src/mass/forward_model_joint.py](src/mass/forward_model_joint.py)::`JointFitPriors`):

```
sigma_da_rel = 2e-4
sigma_de     = 5e-4
sigma_di_deg = 0.05
sigma_dOmega_deg = 0.2
sigma_domega_deg = 0.2
sigma_dM_deg     = 0.5
```

fueron elegidos anchos para permitir absorber drift orbital MPCORB durante
Stage 1. Pero son ~10× más anchos que la incertidumbre típica reportada
por MPCORB para asteroides numerados bien observados. Apretarlos al
nivel real de MPCORB reduce el espacio donde los deltas pueden absorber
señal de masa.

#### Plan técnico

1. **Medir σ reales de MPCORB** (½ día). Para los targets de los 4
   calibradores (`stage4_validation_summary.csv`), extraer las
   incertidumbres formales por elemento del catálogo MPCORB y/o JPL SBDB.
   Construir distribución σ_a/a, σ_e, σ_i, σ_Ω, σ_ω, σ_M. Tomar mediana
   por elemento como nuevo prior.
2. **Variante `JointFitPriors`** (½ día). Agregar `TightPriors`
   instanciable que herede de `JointFitPriors` con los nuevos σ. Exponer
   flag `--priors {default,tight}` en `fit_mass_gaia_joint.py`,
   `run_joint_batch.py`, `run_stage4_validation.py`.
3. **Re-correr Stage 4 calibradores con priors tight**. Output:
   `data/output/stage4_validation_tight_summary.csv`.
4. **Re-correr batch 27 con priors tight** (verificar que χ²_red no se
   degrada). Output: `data/output/loo_batch_results_joint_mahal_tight.csv`.
5. **Diagnóstico** [docs/mass_layer_stage_a1_tight_priors.md]: tabla
   comparativa default vs tight para los 11 calibradores + los 27 fits.
   Veredicto: ¿|z| < 3 para ≥ 3 calibradores?

#### Entregables

- [ ] `data/output/mpcorb_uncertainties_per_element.csv` (medición)
- [ ] Nueva `TightPriors` en [src/mass/forward_model_joint.py](src/mass/forward_model_joint.py)
- [ ] Flag `--priors` en los 3 scripts
- [ ] `data/output/stage4_validation_tight_summary.csv`
- [ ] `docs/mass_layer_stage_a1_tight_priors.md`

#### Criterios de aceptación

- **Veredicto positivo**: ≥ 3/4 calibradores dentro de |z| < 3, χ²_red
  mediano del batch 27 no empeora más de 2×. Conclusión: bias era
  overfitting. Stage A2/A3 no necesarias.
- **Veredicto negativo**: |z| > 3 en ≥ 2 calibradores. Conclusión: bias
  estructural. Pasar a Stage A2.
- **Veredicto intermedio**: |z| ~ 3-5, ratios mejoran pero no pasan
  estricto. Decisión humana: ¿reportar como "bias parcialmente
  corregido" y abandonar, o seguir con A2?

#### Cómo retomar

```bash
git checkout trackA/stage1-tighten-priors
cat docs/mass_layer_stage_a1_tight_priors.md      # qué se concluyó
ls data/output/stage4_validation_tight_*.csv      # qué corrió ya
```

Si no hay archivos: arrancar por el sub-paso 1 (medir σ MPCORB).

#### Progreso

— No arrancada —

---

### Stage 2 — Multi-target joint fit

**Estado**: ⚪ PENDING (sólo si Stage 1 falla el gate)
**Estimación**: ~1 semana
**Branch propuesta**: `trackA/stage2-multitarget-joint`
**Depende de**: Stage 1 con veredicto negativo o intermedio.

#### Objetivo

Romper la degeneración M ↔ deltas haciendo un fit conjunto que comparte
`M_perturber` entre todos los targets de un mismo cuerpo grande, con
deltas por-target libres. La consistencia inter-target restringe la masa
mucho más que cualquier target individual.

#### Plan técnico

1. **Forward model multi-target** (3 días). Nueva función
   `residuals_joint_multitarget(params, target_bundles, perturber_elements)`
   donde `params = (log10_M, da_1, de_1, ..., dM_1, da_2, de_2, ...)`
   con N_targets×6 deltas + 1 masa compartida. Vectoriza el cálculo
   de residuos sobre todos los targets.
2. **Optimizer y bounds** (1 día). `scipy.optimize.least_squares`,
   bounds por bloque (mismo prior por target). Jacobian numérico
   (analítico es overkill).
3. **Tests sintéticos** (1 día). Sintetizar 3 targets con masa real
   conocida + ruido AL/AC realista; fit conjunto debe recuperar la
   masa dentro de 1% sobre toda la suite de prueba.
4. **Aplicar a calibradores** (1 día). Para Pallas (5 targets de Stage 4):
   un solo fit conjunto. Resultado esperado: σ_M se reduce por sqrt(5)
   vs single-target si los deltas son ortogonales a la masa; si la
   masa fit converge a 2.05×10²⁰ kg dentro de 3σ, problema resuelto
   para Pallas. Repetir para Hygiea (5 targets), Ceres (1 target — no
   aplica). Vesta sigue afuera (no hay fits exitosos).
5. **Diagnóstico** (1 día). Tabla comparativa single-target vs
   multi-target sobre los calibradores.

#### Entregables

- [ ] `src/mass/forward_model_joint_multitarget.py`
- [ ] `tests/test_forward_model_joint_multitarget.py`
- [ ] `scripts/mass/fit_mass_gaia_multitarget.py`
- [ ] `data/output/stage_a2_multitarget_validation.csv`
- [ ] `docs/mass_layer_stage_a2_multitarget.md`

#### Criterios de aceptación

- Tests sintéticos pasan dentro de 1% en M_fit.
- ≥ 2 de 3 calibradores (Pallas, Hygiea, Ceres en la medida que aplique)
  dentro de |z| < 3.
- Si pasa: Track A resuelto a este nivel. Si no: pasar a Stage 3.

#### Cómo retomar

```bash
git checkout trackA/stage2-multitarget-joint
ls src/mass/forward_model_joint_multitarget.py
docker compose run --rm test pytest tests/test_forward_model_joint_multitarget.py -v
```

#### Progreso

— No arrancada —

---

### Stage 3 — Forward model físico para drift (OU)

**Estado**: ⚪ PENDING (sólo si Stage 2 todavía deja bias)
**Estimación**: 2-3 semanas
**Branch propuesta**: `trackA/stage3-ou-drift`
**Depende de**: Stage 2 con veredicto negativo.

#### Objetivo

Reemplazar los 6 Δ-elementos absolutos del joint fit por un proceso
estocástico de Ornstein-Uhlenbeck (OU) con varianza calibrada a partir
del catálogo MPCORB. Reduce los grados de libertad efectivos del modelo
de 6 (por target) a 2 (escala temporal τ + amplitud σ), forzando que
el drift se comporte de manera coherente con la dinámica conocida y
no absorba señal de masa.

#### Plan técnico (preliminar — completar al arrancar)

1. Caracterizar σ y τ del drift residual MPCORB-vs-Gaia sobre una
   muestra de asteroides numerados sin perturbers grandes cerca.
2. Implementar evaluación de la log-likelihood OU sobre la trayectoria
   del target durante la ventana de fit.
3. Reemplazar el bloque de 6 deltas en `forward_model_joint.py` por
   parámetros OU `(τ, σ)` con priors empíricos del paso 1.
4. Tests sintéticos + re-corrida calibradores.
5. Si pasa: Track A resuelto. Si no: la limitación es del dataset DR3
   (no hay leverage suficiente) y se cierra Track A esperando DR4.

#### Entregables

- [ ] `src/mass/drift_ou.py`
- [ ] `tests/test_drift_ou.py`
- [ ] `docs/mass_layer_stage_a3_ou_drift.md`
- [ ] `data/output/stage_a3_ou_validation.csv`

#### Criterios de aceptación

- ≥ 2/3 calibradores dentro de |z| < 3.
- O conclusión honesta: el dataset DR3 no soporta más; esperar DR4.

#### Cómo retomar

— No arrancada —

#### Progreso

— No arrancada —

---

## Track B — Trabajos ortogonales pendientes

### Stage 1 — Investigar outliers Stage 2 (Alkeste / 57942 y Eros / 176865)

**Estado**: ⚪ PENDING
**Estimación**: 3-5 días
**Branch propuesta**: `trackB/stage1-stage2-outliers`
**Depende de**: nada.

#### Objetivo

Diagnosticar por qué (124) Alkeste → 57942 quedó en χ²_red = 84.45 y
(389) Industria → 176865 saltó de 4.86 (AL) a 31.87 (Mahalanobis 2D)
tras Stage 2. Hipótesis a testear:

- **Perturber secundario no modelado** durante la ventana de fit.
- **Sistemática AC del catálogo Gaia** que el likelihood Mahalanobis
  expone pero que el AL-only ocultaba.
- **Outlier en algún transit individual** (cosmic ray, edge of FOV) que
  sesga el fit.

#### Plan técnico

1. **Per-observation residuals** (1 día). Modificar
   `fit_mass_gaia_joint.py --debug-residuals` para escribir un CSV con
   `(jd_tdb, r_AL, r_AC, σ_AL, σ_AC, chi2_per_obs)` por observación.
   Identificar transits con `chi2_per_obs > 10` y verificar.
2. **Búsqueda de perturber secundario** (1 día). Para cada outlier,
   consultar el catálogo híbrido por encuentros del target con
   *otros* asteroides en la ventana ±90 días, dist < 0.3 AU. Si hay
   un perturber secundario, incluirlo como cuerpo masivo en el fit
   y ver si χ²_red baja.
3. **Diagnóstico AC vs AL por observación** (1 día). Visualizar r_AL
   vs r_AC para ver si la sistemática es 1D (perturber missing) o
   2D (Gaia scan systematic). Plot en notebook.
4. **Documento** (1 día): `docs/mass_layer_stage2_outliers.md` con
   conclusiones, fixes (si los hay), y nota en el diagnóstico
   original Stage 2.

#### Entregables

- [ ] `data/output/stage2_outliers/<perturber>_<target>_per_obs.csv` × 2
- [ ] `notebooks/stage2_outliers_diagnosis.ipynb`
- [ ] `docs/mass_layer_stage2_outliers.md`
- [ ] Posible fix en `forward_model_joint.py` si se identifica algo
  arreglable (e.g. perturber secundario configurable)

#### Cómo retomar

```bash
git checkout trackB/stage1-stage2-outliers
ls data/output/stage2_outliers/
cat docs/mass_layer_stage2_outliers.md
```

#### Progreso

— No arrancada —

---

### Stage 2 — Completar specificity sobre los 22/27 restantes

**Estado**: ⚪ PENDING
**Estimación**: 1-2 días (mostly cómputo)
**Branch propuesta**: `trackB/stage2-specificity-full`
**Depende de**: nada.

#### Objetivo

Stage 3 del deepwork sólo corrió specificity sobre 5 candidatos
seleccionados por χ²_red bajo. Completar los 22 restantes (todos los
que tuvieron fit exitoso en Stage 2: 27 - 5 = 22) para tener la
distribución completa de p-values y decidir cuáles candidatos son
estadísticamente reales.

#### Plan técnico

1. **Configurar batch sobre los 22** (½ día). Modificar (o reusar)
   `scripts/mass/run_specificity_test.py` con la lista completa de
   targets de `loo_batch_results_joint_mahal.csv`.
2. **Correr** (½ día cómputo, ~30 min con 24 workers para 22 × 50 nulls).
3. **Análisis y reporte** (½ día):
   - Histograma de p_χ² y p_mass sobre los 22 candidatos.
   - Tabla: cuántos pasan p_χ² ≤ 0.05; cuántos pasan p_mass ≤ 0.05.
   - Actualizar `docs/mass_layer_stage3_diagnostic.md` con apéndice
     "Specificity sobre los 22 restantes".

#### Entregables

- [ ] `data/output/specificity_v2/` extendido con los 22 candidatos
- [ ] `data/output/specificity_test_v2_full.csv`
- [ ] Apéndice en `docs/mass_layer_stage3_diagnostic.md`

#### Criterios de aceptación

- 22 fits corridos sin error.
- Reporte honesto del número de "detecciones" específicas (esperable: 0-3).

#### Cómo retomar

```bash
git checkout trackB/stage2-specificity-full
ls data/output/specificity_v2/ | wc -l       # cuántos fits hay
cat data/output/specificity_test_v2_full.csv
```

#### Progreso

— No arrancada —

---

### Stage 3 — Side-paper: 25,283 false positives Kepler (sesgo cerca del threshold)

**Estado**: ⚪ PENDING (opcional, side-project)
**Estimación**: 2-3 semanas (es un mini-paper)
**Branch propuesta**: `trackB/stage3-kepler-bias-paper`
**Depende de**: nada.

#### Objetivo

El Stage B del deepwork reveló que **25,283 pares cruzan el threshold
0.05 AU al refinar con N-body**, todos en la misma dirección (Kepler
reporta `<0.05`, N-body recomputa `≥0.05`); 0 cruzan en la dirección
contraria. Esto es un sesgo sistemático del prefiltro Kepler cerca del
threshold. Caracterizarlo cuantitativamente y publicarlo como side-paper
o nota técnica.

#### Plan técnico (preliminar)

1. **Magnitud y distribución del sesgo**:
   - Histograma de Δdist (Kepler - N-body) para los pares cerca del
     threshold (0.04 ≤ dist_kepler ≤ 0.05).
   - Tasa de cruce vs banda orbital, e, i, q.
2. **Mecanismo**:
   - ¿Es un efecto del prefiltro? ¿Del refinamiento Kepler sub-grid?
   - ¿Está correlacionado con dist_au, e, q, dist a planetas grandes?
3. **Comparación con literatura**: ¿se ha reportado este sesgo? Buscar
   en papers de close-approach catalogs (Fienga, JPL CNEOS).
4. **Implicación práctica para usuarios del catálogo**: corregir el
   threshold para reducir false-positive rate o usar el catálogo
   híbrido siempre.
5. **Borrador**: ~10-15 pp, figuras, tablas, código reproducible.

#### Entregables

- [ ] `notebooks/kepler_threshold_bias_analysis.ipynb`
- [ ] `docs/kepler_threshold_bias_paper.md` (borrador)
- [ ] Decisión: ¿submitirlo a Icarus/A&A o dejar como nota técnica?

#### Cómo retomar

— No arrancada —

#### Progreso

— No arrancada —

---

## Apéndices

### A. Dependencias entre etapas

```
Track A (secuencial, condicional):
  A1 tighten priors ──► (gate)
    ✓ pass → Track A resuelto
    ✗ fail → A2 multi-target ──► (gate)
                ✓ pass → Track A resuelto
                ✗ fail → A3 OU drift ──► (gate)
                            ✓ pass → Track A resuelto
                            ✗ fail → cerrar Track A (esperar DR4)

Track B (independientes entre sí, independientes de Track A):
  B1 outliers Stage 2 — pueden hacerse en cualquier momento
  B2 specificity completo — pueden hacerse en cualquier momento
  B3 side-paper Kepler — opcional, no urgente
```

### B. Costo computacional estimado

| Etapa | Trabajo | Wall-clock cómputo (1 máquina 24 cores) |
|-------|---------|------------------------------------------|
| A1    | 1 día   | ~10 min (re-correr Stage 4) |
| A2    | 1 semana| ~30 min (multi-target sobre calibradores) |
| A3    | 2-3 sem.| horas (modelo OU + re-corridas) |
| B1    | 3-5 días| ~5 min (debug fits) |
| B2    | 1-2 días| ~30 min (22 × 50 nulls) |
| B3    | 2-3 sem.| variable (análisis + drafting) |

### C. Referencias para arranque

- [docs/mass_layer_validation.md](docs/mass_layer_validation.md) — el gate FAIL de Stage 4 + sección "¿Se puede arreglar?" que originó este plan.
- [docs/mass_layer_stage3_diagnostic.md](docs/mass_layer_stage3_diagnostic.md) — specificity test, base para Track B Stage 2.
- [docs/mass_layer_stage2_diagnostic.md](docs/mass_layer_stage2_diagnostic.md) — outliers Alkeste y Eros, motivación de Track B Stage 1.
- [FROZEN_RUN.md](FROZEN_RUN.md) — el "25,283 false-positives" que motiva Track B Stage 3.
- [src/mass/forward_model_joint.py](src/mass/forward_model_joint.py) — `JointFitPriors`, lugar para Stage A1.
- [src/mass/likelihood_al.py](src/mass/likelihood_al.py) — Mahalanobis 2D base.

### D. Bitácora de cambios al plan

| Fecha | Cambio | Autor |
|-------|--------|-------|
| 2026-05-29 | Plan creado tras cierre del deepwork (PR #39). | DF |
