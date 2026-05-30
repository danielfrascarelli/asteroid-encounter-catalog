# Track A Stage 1 — Tighten priors orbitales

**Branch**: `trackA/stage1-tighten-priors`
**Fecha**: 2026-05-29
**Plan de origen**: follow-up post-deepwork (disuelto; ver [ROADMAP.md](../ROADMAP.md) § "Estado actual"), Track A Stage 1.
**Diagnóstico previo**: [docs/mass_layer_validation.md](mass_layer_validation.md) — gate FAIL del Stage 4 del deepwork (ratios fit/lit = 0.77 / 0.57 / 0.24 para Ceres / Pallas / Hygiea, |z| > 6 sistemático).

---

## TL;DR — Veredicto

**FAIL del gate.** Apretar los priors de los 6 Δ-elementos en el joint fit por
factores 2000–20 000× **no corrige el sesgo de masa**: en los pares
matched de Stage 4 (Pallas/28036, Pallas/47563, Pallas/73243, Hygiea/16772)
la masa fit con priors `tight` cambia <0.05 % vs `default`. Los 12 fits OK
de Stage 4 con priors tight tienen todos |z| > 15.

**Conclusión**: el bias **no** proviene de overfitting de los Δ-elementos.
Los deltas convergen a valores pequeños incluso con priors anchos — no
están absorbiendo señal de masa. El bias es **estructural** y el plan
debe pasar a **Track A Stage 2** (multi-target joint fit).

Hallazgo secundario relevante: en los 27 candidatos LOO del batch
(Mahalanobis joint fit), tight priors **sí** cambian las masas
significativamente (mediana mass_tight/mass_default = 0.93; 20 de 27 con
|cambio| > 5 %, hasta -97 % en un caso) y degradan χ²_red en 5 fits.
Esto sugiere que para los candidatos con señal débil los deltas estaban
absorbiendo **ruido**, no señal de masa — consistente con la conclusión
del specificity test (Stage 3 del deepwork) de que la mayoría no son
detecciones reales.

---

## Lo que hicimos

### 1. Medición empírica de σ orbitales

Script: [scripts/mass/measure_mpcorb_uncertainties.py](../scripts/mass/measure_mpcorb_uncertainties.py).
Salida: [data/output/mpcorb_uncertainties_per_element.csv](../data/output/mpcorb_uncertainties_per_element.csv).

Consultamos JPL SBDB (`https://ssd-api.jpl.nasa.gov/sbdb.api?full-prec=1`)
para los 4 perturbers (1, 2, 4, 10) y los 11 targets con fit `ok` en
[stage4_validation_summary.csv](../data/output/stage4_validation_summary.csv).
Cada respuesta incluye `orbit.elements[].sigma`: el σ formal del fit
heliocéntrico contra todas las observaciones MPC.

**Targets (11 numbered MBA, n=11):**

| elemento | mediana σ | p90 σ | máximo σ |
|---|---|---|---|
| a (AU)         | 1.6e-9 | 4.0e-9 | 7.0e-9 |
| e              | 8.8e-10| 2.0e-9 | 2.8e-9 |
| i (deg)        | 6.0e-8 | 8.2e-8 | 8.2e-8 |
| Ω (deg)        | 3.6e-7 | 7.3e-7 | 7.7e-7 |
| ω (deg)        | 1.1e-6 | 3.0e-6 | 3.3e-6 |
| M (deg)        | 1.2e-6 | 1.9e-6 | 4.0e-6 |

**Perturbers (1/2/4/10):** todos los σ son del mismo orden o más estrechos.

### 2. Definición de `TightPriors`

El plan literal pidió "5–10× más estrechos". Pero la medición mostró
que los priors `default` están 10^5–10^6× sueltos respecto al σ formal SBDB.
Para que el experimento sea informativo definimos `TightPriors` a un
nivel **principiado** (~10× el σ p90 de la población de targets), lo que
sigue siendo 2000–20 000× más estrecho que el default. Los bounds
acompañan a ~100× σ_tight, suficientes para que el optimizer pueda
moverse sin chocar contra hard cuts.

| parámetro | default σ | tight σ | ratio | tight bounds |
|---|---|---|---|---|
| da/a        | 2e-4 | 1e-7 | **2000×** | ±1e-5 |
| de          | 5e-4 | 1e-7 | **5000×** | ±1e-5 |
| di (deg)    | 5e-2 | 1e-5 | **5000×** | ±1e-3 |
| dΩ (deg)    | 2e-1 | 1e-5 | **20000×**| ±1e-3 |
| dω (deg)    | 2e-1 | 3e-5 | **6667×** | ±3e-3 |
| dM (deg)    | 5e-1 | 3e-5 | **16667×**| ±3e-3 |

Implementación: [src/mass/forward_model_joint.py](../src/mass/forward_model_joint.py)
agregó `TIGHT_PRIORS`, registro `PRIOR_PRESETS` y helper `resolve_priors()`.
Los 4 scripts del pipeline (`fit_mass_gaia_joint.py`, `run_joint_batch.py`,
`run_stage4_validation.py`, `summarize_joint_fits.py`) reciben ahora un
flag `--priors {default,tight}` y los outputs llevan sufijo `_tight` para
no colisionar con los outputs default.

### 3. Re-corrida Stage 4 con tight priors

```
docker compose run --rm pipeline python -m scripts.mass.run_stage4_validation \
  --priors tight --workers 12
```

Salida: [data/output/stage4_validation_summary_tight.csv](../data/output/stage4_validation_summary_tight.csv)
(20 candidatos, 12 con fit `ok`).

### 4. Re-corrida batch 27 con tight priors

Construimos `data/output/batch27_candidates.csv` desde los 27 (perturber,
target, date) ya en `loo_batch_results_joint_mahal.csv` y disparamos los
fits Mahalanobis-tight con `ProcessPoolExecutor(max_workers=12)`.

Salida: [data/output/loo_batch_results_joint_mahal_tight.csv](../data/output/loo_batch_results_joint_mahal_tight.csv)
(27/27 fits OK).

---

## Resultados Stage 4 — matched pairs

Las muestras de Stage 4 default y tight no son idénticas porque el
catálogo híbrido fue actualizado entre ambas corridas. Comparamos los
4 pares (perturber, target) que aparecen en ambas:

| perturber | target | χ²_red def | χ²_red tight | mass_fit def [kg]   | mass_fit tight [kg] | ratio def | ratio tight | Δratio    |
|-----------|--------|------------|--------------|---------------------|---------------------|-----------|-------------|-----------|
| Pallas    | 28036  | 0.360      | 0.366        | 1.171e+20           | 1.171e+20           | 0.5711    | 0.5710      | −7.2e-5   |
| Pallas    | 47563  | 0.698      | 0.698        | 1.168e+20           | 1.169e+20           | 0.5695    | 0.5700      | +4.8e-4   |
| Pallas    | 73243  | 0.750      | 0.751        | 1.171e+20           | 1.171e+20           | 0.5711    | 0.5710      | −8.1e-5   |
| Hygiea    | 16772  | 0.947      | 0.873        | 1.943e+19           | 1.943e+19           | 0.2341    | 0.2341      | +1.5e-5   |

**Cambio relativo en masa**: |Δratio| < 5e-4 en todos los casos. Las masas
de Pallas y Hygiea convergen al mismo valor (sesgado) sin importar el
prior. El bias **no** viene de los deltas absorbiendo señal.

Y los 12 fits OK del run tight (no sólo los 4 matched):

| perturber | n_ok | rango ratio_fit/lit | rango |z|       |
|-----------|------|---------------------|------------------|
| Ceres     | 3    | 0.72 – 0.74         | 17 – 33          |
| Pallas    | 3    | 0.570 – 0.571       | 17.5 – 17.6      |
| Vesta     | 1    | 1.90                | 62               |
| Hygiea    | 5    | 0.13 – 0.24         | 16 – 18          |

Ningún calibrador pasa |z| < 3. **0 de 4 perturbers pasa el gate.**

---

## Resultados batch 27 — diagnóstico colateral

Sobre las 27 (perturber, target) del batch LOO Mahalanobis matcheamos
default vs tight pair a pair:

- **mass_tight / mass_default**: mediana 0.93, p10 0.43, p90 1.11,
  min 0.025 (312/78961: 1.6e17 → 4.1e15), max 1.36 (235/118886).
- **|cambio| > 5 %**: 20 de 27 pares.
- **χ²_red_tight / χ²_red_default**: mediana 1.04, p90 3.53, max 45 (19/53467: 0.41 → 18.6).
- **χ²_red sube > 2×**: 5 pares (49/94474, 46/5998, 113/57806, 618/108638, 19/53467).

Lectura: para fits con señal débil (la mayoría de los 27) los Δ-elementos
del default fit absorben ruido astrométrico y producen masas espurias.
Los priors tight les quitan ese freedom, χ² empeora porque el modelo ya
no fita el ruido, y la masa migra (típicamente hacia abajo, lo que es
consistente con haber estado "infladas" por noise-absorption).

Este resultado refuerza la conclusión del Stage 3 del deepwork
([docs/mass_layer_stage3_diagnostic.md](mass_layer_stage3_diagnostic.md))
de que la mayoría de los 27 candidatos son nulls. **Pero no resuelve el
bias de los calibradores** — ese es estructural, no es un overfitting
de deltas.

---

## Implicaciones para el plan

- **Stage A1**: FAIL del gate (criterio: ≥3/4 calibradores con |z| < 3).
  Resultado: 0 de 4.
- **Stage A2** (multi-target joint fit) — recomendado. La hipótesis
  diagnóstica de A1 (deltas absorben señal de masa) queda descartada; el
  bias debe venir de otra parte del modelo. El multi-target fit comparte
  `M_perturber` entre N targets, lo que rompe la degeneración entre la
  masa y cualquier otro parámetro que NO se comparte entre targets
  (eg. errores sistemáticos de propagación del target, sesgos AC del
  catálogo Gaia por scan systematic, etc.).
- **Stage A3** (OU drift) — sigue siendo una segunda línea si A2 también
  deja bias residual.

### Sobre el batch 27 (side-benefit de A1)

Si bien tight priors NO ayudan para los calibradores, sí parecen ser un
filtro útil para descartar candidatos puramente ruidosos: un fit que
cambia su masa por factor 2× o su χ²_red por 10× cuando se aprietan
los priors es muy probablemente noise-driven. Vale la pena considerar
correr el specificity test (Track B Stage 2) **sobre los fits tight**
en vez de los default, o como cross-check, para ver si los fits que ya
parecían robustos en default siguen siendo consistentes en tight.

---

## Referencias

- [src/mass/forward_model_joint.py](../src/mass/forward_model_joint.py) — `JointFitPriors`, `TIGHT_PRIORS`, `PRIOR_PRESETS`, `resolve_priors`.
- [scripts/mass/measure_mpcorb_uncertainties.py](../scripts/mass/measure_mpcorb_uncertainties.py) — query SBDB.
- [data/output/mpcorb_uncertainties_per_element.csv](../data/output/mpcorb_uncertainties_per_element.csv) — input empírico para TightPriors.
- [data/output/stage4_validation_summary_tight.csv](../data/output/stage4_validation_summary_tight.csv) — re-corrida Stage 4 con tight.
- [data/output/loo_batch_results_joint_mahal_tight.csv](../data/output/loo_batch_results_joint_mahal_tight.csv) — re-corrida batch 27 con tight.
- [docs/mass_layer_validation.md](mass_layer_validation.md) — gate FAIL Stage 4 deepwork.
- [docs/mass_layer_stage3_diagnostic.md](mass_layer_stage3_diagnostic.md) — specificity test deepwork.
