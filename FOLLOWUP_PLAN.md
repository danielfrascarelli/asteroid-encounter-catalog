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
**Última actualización**: 2026-05-30 (retomado)
**Estado**: ✅ **TRACK A CERRADO + TRACK B COMPLETO (B1/B2/B3)**. Tras retomar el
plan se ejecutó todo el trabajo pendiente ejecutable en paralelo. PRs #46–#49
abiertos a `main`:

- **#46 — Track A gate-check (cierre)**: el gate de A3 **falla**. Scans de χ²(masa)
  con ventana acotada ±60 d sobre los calibradores muestran que **ningún cuerpo
  ≠ Ceres tiene leverage ≳3σ** (Pallas σ≤1.2; Hygiea ≤2.1 de A2.6; Vesta sin
  targets). Un **experimento de sensibilidad a la ventana sobre Ceres** (PR #51)
  cierra el caso: aunque Ceres tiene el leverage más fuerte (σ~20), su mínimo de
  masa **vaga ~100× (0.04–4.37× lit) según la ventana, sin plateau, y la χ²(masa)
  es multimodal siempre** → **masa no identificable, incluido Ceres**. **A3 NO se
  arranca; Track A se cierra** con un negativo universal y limpio. Doc:
  [docs/mass_layer_track_a_closure.md](docs/mass_layer_track_a_closure.md).
- **#47 — Track B1 (outliers Stage 2)**: los dos χ² altos (57942=84.5, 176865=31.9)
  son **misfit along-scan pervasivo**, no sistemática AC (σ_AC≈612mas≫σ_AL≈3mas →
  χ²_red_AC≈0.1), no transit aislado, no perturbador secundario. Doc:
  [docs/mass_layer_stage2_outliers.md](docs/mass_layer_stage2_outliers.md).
- **#48 — Track B2 (specificity completo)**: **0/27 detecciones** específicas;
  `p_mass` nunca <0.08 (la masa real nunca sobresale del nulo). Apéndice en
  [docs/mass_layer_stage3_diagnostic.md](docs/mass_layer_stage3_diagnostic.md).
- **#49 — Track B3 (sesgo threshold Kepler)**: caracteriza los 25,283 cruces
  (0.29% de los pares refinados) y **corrige el framing de FROZEN_RUN**: el "0
  cruces hacia abajo" es **censura** (el catálogo no tiene pares Kepler≥0.05), y
  Δdist está dominado por scatter (mediana −1e-6 AU), no por un sesgo
  unidireccional. Nota técnica + notebook:
  [docs/kepler_threshold_bias_paper.md](docs/kepler_threshold_bias_paper.md).

**Conclusión global de la capa de masas**: convergente y honesta — la
determinación de masas **no es viable con Gaia DR3** en este pipeline (Track A
cerrado, specificity 0/27). El cuello de botella es el **leverage intrínseco**
del dataset, no el optimizador ni la parametrización del drift. El tooling
(joint multi-target, profiled, ventana acotada) queda listo para DR4/FPR.
**Historia de Track A** (resumen; detalle en los docs por etapa): A1 (tighten
priors) FAIL → A2 (multi-target) FAIL → closing-loop (el ratio era un bug del
optimizador, no físico) → A2.5 (profiled, bug arreglado; aflora sesgo real-data)
→ A2.6 ("9×" = basin espurio + drift; M↔deltas activa en datos reales;
no-identificable) → **gate-check A3 FAIL (#46) → Track A cerrado**. El leverage
intrínseco de DR3 es el límite; no hay nada que A3 (OU) pueda destrabar donde no
hay señal.

| Track | Etapa | Estado | Branch | PR | Inicio | Fin | Notas |
|-------|-------|--------|--------|----|--------|-----|-------|
| A | 1: tighten priors | 🟢 DONE | trackA/stage1-tighten-priors | #41 | 2026-05-29 | 2026-05-29 | Veredicto: bias estructural; pasar a A2 |
| A | 2: multi-target joint fit | 🟢 DONE | trackA/stage2-multitarget-joint | #42 | 2026-05-29 | 2026-05-29 | Gate FAIL 0/2; masa multi == single; refuta degeneración M↔deltas |
| A | 2.5: optimizador perfilado | 🟡 PARCIAL (mergeado) | trackA/stage2-multitarget-joint | #42 ✅ | 2026-05-29 | 2026-05-30 | Bug optimizador ARREGLADO (closing-loop ratio 1.000); pero Hygiea real → 9× lit. Queda sesgo real-data |
| A | 2.6: investigar sesgo real-data Hygiea | 🟢 DONE | trackA/stage2.6-realdata-bias | #44 | 2026-05-30 | 2026-05-30 | 9× = basin espurio + drift, NO sesgo físico. Masa no-identificable en DR3 (deltas absorben sobre datos reales). Solo Ceres defendible. Fix `--joint-window-days` incluido |
| A | gate-check A3 (ventana acotada) | 🟢 DONE | trackA/gate-a3-close | #46 | 2026-05-30 | 2026-05-30 | Gate FAIL: ningún calibrador ≠Ceres ≳3σ, Pallas σ≤1.2. A3 no se arranca |
| A | Ceres window-sensitivity | 🟢 DONE | trackA/ceres-window-sensitivity | #51 | 2026-05-30 | 2026-05-30 | Ceres no-identificable: mínimo barre 0.04–4.37× según ventana, multimodal. Cierre universal |
| A | 3: OU forward model para drift | ⚫ ABANDONED | — | — | — | — | No justificada: gate FAIL. El límite es leverage DR3, no la parametrización del drift |
| B | 1: outliers Stage 2 (Alkeste/57942, Industria/176865) | 🟢 DONE | trackB/stage1-stage2-outliers | #47 | 2026-05-30 | 2026-05-30 | Misfit along-scan pervasivo; refuta AC-systematic/transit-aislado/perturbador-secundario |
| B | 2: specificity sobre 22/27 restantes | 🟢 DONE | trackB/stage2-specificity-full | #48 | 2026-05-30 | 2026-05-30 | 0/27 detecciones; p_mass nunca <0.08 |
| B | 3: side-paper 25,283 false-positives Kepler | 🟢 DONE (draft) | trackB/stage3-kepler-bias-paper | #49 | 2026-05-30 | 2026-05-30 | "0 cruces abajo" = censura; Δdist scatter-dominado. Decisión publicar pendiente |

**Qué queda** (todo lo ejecutable ya se cerró 2026-05-30):

1. **Decisión humana sobre B3**: ¿publicar la nota del sesgo Kepler como apéndice
   del catálogo o como nota técnica standalone? (recomendación: apéndice).
2. **DR4 / FPR**: reabrir Track A (A3/multi-target, tooling listo) cuando arcos
   más largos den leverage a >1 calibrador. Medir falsos negativos del prefiltro
   (experimento en la nota B3).
3. **Resto del roadmap principal** (dashboard, caracterización big-catalog), fuera
   del scope de este plan de follow-up.

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

**Estado**: 🟡 IN PROGRESS
**Estimación**: ~1 día
**Branch**: `trackA/stage1-tighten-priors` (creada 2026-05-29)
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

- **2026-05-29**: Stage 1 completo. Veredicto **FAIL del gate**:
  - σ MPCORB medidos vía JPL SBDB (11 targets + 4 perturbers numbered):
    los priors `default` están 10^5–10^6× más sueltos que la
    incertidumbre formal — ratificó que apretar tiene sentido.
  - `TIGHT_PRIORS` definidos a ~10× σ p90 (2000–20000× más estrechos
    que default). Flag `--priors` agregado a los 4 scripts del pipeline.
  - Stage 4 re-corrido con tight priors (20 candidatos, 12 OK):
    todos los fits OK tienen |z| > 15. 0/4 calibradores pasan |z|<3.
  - En los 4 pares (perturber, target) matched contra el default
    (Pallas/28036/47563/73243, Hygiea/16772) las masas fit cambian
    < 0.05 %. Bias es estructural, no por overfitting de deltas.
  - Batch 27 (LOO Mahalanobis) re-corrido con tight: median
    mass_tight/mass_default = 0.93, pero 20/27 con |cambio|>5 %.
    Lectura: los Δ-elementos default absorbían ruido en candidatos
    de señal débil — coherente con specificity Stage 3.
  - Diagnóstico: [docs/mass_layer_stage_a1_tight_priors.md](docs/mass_layer_stage_a1_tight_priors.md).
  - Pendiente: PR a `main` con la branch `trackA/stage1-tighten-priors`.

---

### Stage 2 — Multi-target joint fit

**Estado**: 🟢 DONE — gate FAIL (2026-05-29)
**Estimación**: ~1 semana
**Branch**: `trackA/stage2-multitarget-joint` (creada 2026-05-29)
**Depende de**: Stage 1 (FAIL del gate, ver más arriba).
**Diagnóstico**: [docs/mass_layer_stage_a2_multitarget.md](docs/mass_layer_stage_a2_multitarget.md)

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

- [x] `src/mass/forward_model_joint_multitarget.py`
- [x] `tests/test_forward_model_joint_multitarget.py`
- [x] `scripts/mass/fit_mass_gaia_multitarget.py`
- [x] `data/output/stage_a2_multitarget_validation.csv`
- [x] `docs/mass_layer_stage_a2_multitarget.md`

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

Próximo sub-paso al retomar: **lanzar el fit calibrador**, p. ej. Pallas
(7 targets, perturber=2 en `data/output/stage4_validation_summary.csv`):

```bash
docker compose run --rm pipeline python -m scripts.mass.fit_mass_gaia_multitarget \
  --perturber 2 \
  --targets-csv data/output/stage4_validation_summary.csv \
  --likelihood mahalanobis2d --priors default \
  --output data/output/multitarget/fit_000002_pallas.json
```

#### Progreso

- **2026-05-29**: sub-pasos 1-3 listos. Pausa con WIP en la branch.
  - **Forward model**: `src/mass/forward_model_joint_multitarget.py` con
    `TargetBundle`, `residuals_joint_multitarget`, `fit_joint_multitarget`,
    `make_bounds`, `prior_residuals_multitarget`. Vector de parámetros
    `(1 + 6N)` (M compartida + 6 deltas por target). Likelihoods AL y
    Mahalanobis 2D heredados del joint single-target.
  - **Tests sintéticos**: `tests/test_forward_model_joint_multitarget.py`,
    6/6 pass. El test clave (`test_fit_joint_multitarget_recovers_mass_with_per_target_da`)
    verifica recuperación de M compartida + da_rel per-target usando
    una señal de masa time-varying (Heaviside en el encuentro) y un
    offset orbital constante por-target.
  - **CLI**: `scripts/mass/fit_mass_gaia_multitarget.py` lee
    `--targets-csv` (con columnas `perturber,target,encounter_date`)
    o `--targets-json`, hace LOO orbit fit por target, arma bundles
    y corre el joint multitarget. Output JSON con
    `mass_kg / mass_sigma_kg / per_target_deltas / chi2_red_joint`.
  - **Lint**: ruff limpio sobre los 3 archivos nuevos.
  - **Pendiente al retomar**: lanzar el fit sobre Pallas (7 targets) +
    Hygiea (8 targets), comparar contra literatura (|z| < 3 sobre ≥ 2/3)
    y, si pasa, escribir `docs/mass_layer_stage_a2_multitarget.md`.

- **2026-05-29 (retomado)**: fits calibradores corridos. **Gate FAIL 0/2**.
  - **Pallas** (perturber 2, 5 targets ok): M_fit = 1.171×10²⁰ kg,
    ratio = 0.5711, z = −17.0, χ²_red = 0.57.
  - **Hygiea** (perturber 10, 5 targets ok): M_fit = 1.919×10¹⁹ kg,
    ratio = 0.2313, z = −16.0, χ²_red = 74.0 (dominado por outlier 45989).
  - **Ceres**: solo 1 target ok → multi-target no aplica. **Vesta**: 0 fits.
  - **Hallazgo decisivo**: la masa multi-target es **idéntica** a la
    single-target (Pallas 0.5711=0.5711; Hygiea 0.2313 vs 0.234). Compartir
    M entre 5 targets no movió la masa → **refuta la degeneración M↔deltas**
    como mecanismo del bias. Los deltas ajustados son ~10⁻⁷ (no absorben
    señal). El bias es coherente por-target y escala con 1/M_real → apunta
    al modelo de deflección, no a la parametrización orbital.
  - Entregables escritos: CSV de validación + diagnóstico
    [docs/mass_layer_stage_a2_multitarget.md](docs/mass_layer_stage_a2_multitarget.md).
  - **Recomendación**: NO arrancar A3 a ciegas (su premisa quedó debilitada);
    primero un **closing-the-loop test** del forward model de deflección
    (inyectar deflección N-body con masa conocida sobre los mismos transits
    y verificar recuperación). ~1 día; discrimina bug vs límite del dataset.
  - Pendiente: PR a `main` con la branch `trackA/stage2-multitarget-joint`.

---

### Stage 2.5 — Optimizador perfilado (profiled likelihood)

**Estado**: 🟡 PARCIAL — bug del optimizador arreglado; sesgo real-data persiste
**Estimación**: ~3-4 días (implementación + validación calibradores)
**Branch**: `trackA/stage2-multitarget-joint` (PR #42)
**Depende de**: closing-the-loop test (hecho). Diagnósticos:
[docs/mass_layer_closing_loop_leverage.md](docs/mass_layer_closing_loop_leverage.md),
[docs/mass_layer_stage_a2_5_profiled.md](docs/mass_layer_stage_a2_5_profiled.md)

#### Motivación

El closing-the-loop test (inyectar masa conocida sobre transits reales sin
ruido) reveló que **el fit no mide masa — devuelve el prior fotométrico M_H**.
No es un bias físico ni la degeneración M↔deltas:

- El escaneo χ²(log10_M) con deltas congelados tiene un **mínimo profundo en la
  masa real**: Hygiea Δχ²~29 (~5σ), Ceres ~25σ. Solo Pallas es plano (~0.09σ).
- El gradiente de masa es real (dχ²/dlog10M≈−43) pero `least_squares` termina
  por `xtol` sin moverse de x0: trust-region mal condicionado (cond~10¹³, masa
  ~20 vs deltas ~10⁻⁴). Ni `diff_step` ni `x_scale` lo destraban.

#### Plan técnico

1. `fit_joint_multitarget_profiled`: outer 1-D (Brent acotado) sobre `log10_M`;
   inner least_squares solo sobre los 6N deltas (reparametrizados `u=δ/σ` para
   O(1), bien condicionado). σ_M de la curvatura del χ² perfilado en el mínimo.
2. Validar con el closing-loop (debe recuperar Hygiea 8.3e19 y Pallas
   leverage-limited con error enorme).
3. Re-correr calibradores reales (Ceres/Hygiea/Pallas) y batch 27.
4. Replicar el fix en el fit single-target `forward_model_joint.py`.

#### Criterios de aceptación

- Closing-loop recupera la masa inyectada para Hygiea/Ceres dentro de ~5%.
- Calibradores reales: Ceres y/o Hygiea dentro de |z|<3 contra literatura.
- Si pasa: **revierte el gate FAIL del deepwork**; la capa de masas es viable
  en DR3 para perturbadores con buena geometría. A3 (OU) deja de ser necesaria.

#### Entregables

- [x] `fit_joint_multitarget_profiled` en `src/mass/forward_model_joint_multitarget.py`
- [x] Flag `--optimizer {joint,profiled}` en `fit_mass_gaia_multitarget.py` + closing-loop
- [x] `data/output/stage_a2_5_profiled_validation.csv`
- [x] `docs/mass_layer_stage_a2_5_profiled.md`

#### Progreso

- **2026-05-29**: etapa creada a partir del closing-loop test. Tooling de
  diagnóstico ya commiteado (`closing_loop_test.py`, `probe_mass_sensitivity.py`).
- **2026-05-29 (implementado)**: `fit_joint_multitarget_profiled` (outer Brent
  1-D en log10_M, inner solo-deltas con `diff_step=0.1σ`). Resultados:
  - **Closing-loop sin ruido, profiled**: Hygiea recupera la masa inyectada
    **exacta** (ratio 1.000) vs joint viejo 0.234 → **bug del optimizador
    confirmado y arreglado**. Pallas no-identificable (χ²≈0 a masa 10¹⁵).
  - **Reales vs literatura**: Hygiea → 7.5×10²⁰ (**ratio 9.0**, z=167,
    χ²_red 1.19 buen ajuste); Pallas → 2.7×10¹⁶ (no-identificable, σ engañosa).
  - **Lectura**: el "ratio=M_H/M_lit" era el optimizador (no físico). Con el
    optimizador arreglado aflora un **sesgo real-data** (Hygiea 9× alto):
    no es leverage (~5σ) ni el optimizador → deriva orbital / sistemática Gaia.
  - **Limitación**: `--noise realistic` del closing-loop no es válido
    cuantitativamente (inyecta ruido AC ~433 mas inconsistente con whitening AL;
    χ²_red~4000). σ_M por curvatura no robusta en régimen plano.
  - Próximo: investigar el sesgo real-data de Hygiea (solapa Track B1).

---

### Stage 2.6 — Diagnóstico del sesgo real-data de Hygiea

**Estado**: 🟢 DONE — PR #44 mergeado a main (`93dc5b4`)
**Branch**: `trackA/stage2.6-realdata-bias` (eliminada tras merge)
**Depende de**: A2.5 (bug del optimizador arreglado).
**Diagnóstico**: [docs/mass_layer_stage_a2_6_realdata_bias.md](docs/mass_layer_stage_a2_6_realdata_bias.md)

**PR**: #44 (mergeado a main `93dc5b4`).

#### Objetivo

A2.5 dejó a Hygiea real en 9× literatura con χ²_red≈1.2 y lo llamó "sesgo
real-data". Stage 2.6 lo prueba escaneando la verosimilitud perfilada χ²(masa)
sobre datos REALES (por-target y joint), con per-observation residuals,
estabilidad target-a-target, y búsqueda de perturbador secundario.

#### Entregables

- [x] `scripts/mass/realdata_mass_scan.py` (scan profiled χ²(masa) real, per-target + joint, per-obs, ventana joint configurable)
- [x] `data/output/stage2_6/hygiea_*` (scan/perobs/summary; gitignored)
- [x] `docs/mass_layer_stage_a2_6_realdata_bias.md`

#### Hallazgos

1. **La χ²(masa) real es dentada y multimodal** (Δχ² salta 10⁴–10⁶ entre puntos
   de grilla adyacentes). El "9×" fue un **basin espurio** del Brent acotado, no
   un mínimo: re-escanear da el "mínimo" en 0.66× — resultado **no reproducible**.
   Causa: la **degeneración M↔deltas está activa sobre datos reales** (A2 la
   creyó descartada, pero eso valía solo para sintético sin ruido).
2. **Deriva orbital del arco joint unilateral**: la máscara joint de
   `_build_bundle` no tiene cota superior; para encuentros tempranos el arco
   post-encuentro corre ~2 años. 45989 (χ²_red@lit 1129; rampa 17→789 sobre
   ~814 d) lo evidencia.
3. **Ventana simétrica ±60 d cura 45989** (χ²_red@lit 1129→0.49) **pero no
   identifica la masa**: mínimos por-target 0.19–3.4× lit, curvas aún dentadas,
   sin leverage limpio ≳3σ.
4. **Perturbador secundario descartado**: ninguno de los 320 encuentros
   secundarios (<0.3 AU, ±90 d) involucra un cuerpo masivo.

#### Veredicto

El "9×" no era un sesgo físico coherente. **La masa de Hygiea no es determinable
con DR3 en este pipeline**; solo Ceres (25σ) tiene leverage defendible. El
cuello de botella es el **leverage**, no el optimizador.

#### Próximo paso

Decisión humana: (i) **cerrar Track A** con el veredicto "no determinable en DR3
salvo Ceres", o (ii) **A3 (OU)** — solo si antes, con ventana acotada, un
calibrador ≠ Ceres muestra ≳3σ. Fix barato recomendado: ventana joint
simétrica/acotada en `_build_bundle` (limpia el drift de 45989).

---

### Stage 3 — Forward model físico para drift (OU)

**Estado**: ⚫ ABANDONED — gate-check FAIL (2026-05-30, PR #46)
**Estimación**: 2-3 semanas (no realizada)
**Branch propuesta**: `trackA/stage3-ou-drift`
**Depende de**: Stage 2 con veredicto negativo (✓ cumplido, gate FAIL).

> **Cierre (2026-05-30).** El gate-check (escaneo χ²(masa) con ventana acotada
> ±60 d sobre los calibradores, [docs/mass_layer_track_a_closure.md](docs/mass_layer_track_a_closure.md))
> falla: ningún cuerpo ≠ Ceres alcanza ≳3σ de leverage (Pallas σ≤1.2, Hygiea
> ≤2.1, Vesta sin targets). A3 no se arranca: el cuello de botella es el leverage
> intrínseco de DR3, no la parametrización del drift. Track A cerrado.

> ⚠️ **Pre-requisito agregado tras A2 (2026-05-29)**: A2 mostró que los
> deltas orbitales son ~10⁻⁷ y NO absorben señal de masa — la premisa de A3
> (reducir dof de deltas evita la absorción) quedó debilitada. Antes de
> invertir 2-3 semanas en el modelo OU, correr un **closing-the-loop test**:
> inyectar una deflección N-body de masa conocida sobre los transits reales
> de un calibrador y verificar si el forward model joint la recupera. Si NO
> la recupera, el bias es un bug del modelo de deflección (ni A3 lo arregla);
> si SÍ la recupera, el límite es el dataset DR3 y A3/DR4 es el camino.
> Estimado: ~1 día. Es el primer sub-paso real de esta etapa.

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

### Stage 1 — Investigar outliers Stage 2 (Alkeste / 57942 e Industria / 176865)

**Estado**: 🟢 DONE — PR #47 (2026-05-30)
**Estimación**: 3-5 días (realizado en 1 sesión)
**Branch**: `trackB/stage1-stage2-outliers`
**Depende de**: nada.
**Doc**: [docs/mass_layer_stage2_outliers.md](docs/mass_layer_stage2_outliers.md)

> **Veredicto**: misfit **along-scan pervasivo** (σ_AC≈612mas≫σ_AL≈2-4mas →
> χ²_red_AC≈0.1, χ²_red_AL 66-349; |AL pull| 3-5σ en ~20-50% de los transits).
> **Refuta** las tres hipótesis de Stage 2: no es sistemática AC, no es un
> transit aislado (top-1 ≤8% del χ²), no hay perturbador secundario (0 masivos
> en 123/173 encuentros). No es señal de masa. Recomienda corte de calidad por
> residuo along-scan. Tool: `scripts/mass/diagnose_stage2_outliers.py`.

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

**Estado**: 🟢 DONE — PR #48 (2026-05-30)
**Estimación**: 1-2 días (realizado en 1 sesión)
**Branch**: `trackB/stage2-specificity-full`
**Depende de**: nada.
**Doc**: apéndice en [docs/mass_layer_stage3_diagnostic.md](docs/mass_layer_stage3_diagnostic.md)

> **Veredicto**: **0/27 detecciones** específicas (ninguna pasa p_chi2≤0.05 Y
> p_mass≤0.05). `p_mass` nunca <0.08 (mediana 0.26): la masa real nunca sobresale
> del nulo orbital/magnitud-matched. Confirma el cierre de Track A. Salida:
> `data/output/specificity_test_v2_full27.csv`.

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

**Estado**: 🟢 DONE (draft) — PR #49 (2026-05-30). Decisión de publicar pendiente.
**Estimación**: 2-3 semanas (análisis + draft realizados en 1 sesión)
**Branch**: `trackB/stage3-kepler-bias-paper`
**Depende de**: nada.
**Doc**: [docs/kepler_threshold_bias_paper.md](docs/kepler_threshold_bias_paper.md) + notebook

> **Veredicto**: 25,283/8.73M = 0.29% de cruces ascendentes, concentrados en
> [0.045,0.050) AU (1.51%) y en órbitas q_min<1.8 AU (90.9%) / encuentros
> rápidos. **Corrige FROZEN_RUN**: el "0 cruces hacia abajo" es **censura** (el
> catálogo no tiene pares Kepler≥0.05), no una tasa de falsos negativos medida;
> Δdist está dominado por scatter (mediana −1e-6 AU), no por un sesgo
> unidireccional. Propone el experimento (re-refinar [0.05,0.06] AU) para medir
> los falsos negativos. **Pendiente humano**: ¿apéndice del catálogo vs nota
> standalone? Tool: `scripts/validate/analyze_kepler_threshold_bias.py`.

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
Track A (RESUELTO — cerrado 2026-05-30):
  A1 tighten priors  ──► FAIL (bias estructural)
  A2 multi-target    ──► FAIL (masa multi==single)
  closing-loop       ──► bug del optimizador (no físico)
  A2.5 profiled      ──► bug arreglado; aflora sesgo real-data
  A2.6 real-data     ──► no-identificable (M↔deltas activa)
  A3 gate-check      ──► FAIL (ningún calibrador ≠Ceres ≳3σ)
  Ceres window-scan  ──► Ceres no-identificable (mínimo 0.04–4.37×, multimodal) ⇒ A3 ABANDONED, Track A CERRADO (negativo universal)

Track B (COMPLETO 2026-05-30):
  B1 outliers Stage 2     ──► DONE: misfit along-scan (#47)
  B2 specificity completo ──► DONE: 0/27 detecciones (#48)
  B3 side-paper Kepler    ──► DONE (draft): censura, no sesgo unidireccional (#49)
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
| 2026-05-29 | A1 gate FAIL (bias estructural). | DF |
| 2026-05-29 | A2 gate FAIL (masa multi==single; refuta degeneración M↔deltas). Próximo: closing-the-loop test antes de A3. | DF |
| 2026-05-29 | Closing-loop test: el ratio = M_H/M_lit es un **bug del optimizador** (no bias físico). χ²(masa) tiene mínimo real para Hygiea/Ceres; solo Pallas sin leverage. Nueva etapa A2.5 (optimizador perfilado). Revierte el veredicto "cerrar Track A". | DF |
| 2026-05-29 | A2.5 implementado (optimizador perfilado). Bug optimizador ARREGLADO (closing-loop ratio 1.000). Pero Hygiea real → 9× lit (χ²_red 1.19): aflora sesgo real-data, ahora aislado del optimizador y del leverage. Capa aún no validada. | DF |
| 2026-05-30 | PR #42 mergeado a main (`02b6946`); branch eliminada. Trabajo **PAUSADO**. Próxima al retomar: Stage 2.6 (investigar sesgo real-data Hygiea). | DF |
| 2026-05-30 | A2.6 cerrado. El "9×" de A2.5 = basin espurio sobre χ²(masa) dentada/multimodal (degeneración M↔deltas activa en datos reales) + deriva orbital del arco joint unilateral (45989). Ventana ±60 d cura 45989 pero no identifica la masa. Perturbador secundario descartado. Veredicto: masa de Hygiea no determinable en DR3; solo Ceres defendible. | DF |
| 2026-05-30 | PR #44 mergeado a main (`93dc5b4`); rama eliminada + podadas ramas stale. Fix `--joint-window-days` disponible en producción. Trabajo **PAUSADO**. Al retomar: decisión cerrar Track A vs A3, o pivotar a Track B. | DF |
| 2026-05-30 | **Retomado.** Decisión: ejecutar todo lo pendiente en paralelo. Gate-check A3 corrido (Pallas/Ceres ±60 d) → **FAIL**: solo Ceres con leverage, y su estimador es window-dependiente → **Track A cerrado**, A3 abandonada (PR #46). | DF |
| 2026-05-30 | Track B1 (PR #47): outliers Stage 2 = misfit along-scan pervasivo (σ_AC≫σ_AL), no AC-systematic/transit/perturbador-secundario. Track B2 (PR #48): specificity 0/27 detecciones. Track B3 (PR #49): los 25,283 cruces = censura + scatter, no sesgo unidireccional; nota técnica draft. | DF |
| 2026-05-30 | **Rigor pass (PR #51).** Experimento de sensibilidad a la ventana sobre Ceres: el mínimo de masa barre 0.04–4.37× lit sin plateau, multimodal en toda ventana → **Ceres tampoco es identificable**. Cierre de Track A endurecido a **negativo universal** (ninguna masa DR3, incluido Ceres). **Decisión B3**: nota va como **apéndice de métodos del catálogo**, no standalone (efecto de censura/selección estándar, FN sin medir); corregida la afirmación "no false negatives" en FROZEN_RUN.md. | DF |
