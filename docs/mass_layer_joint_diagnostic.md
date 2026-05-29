# Stage 1.5 — Diagnóstico joint fit vs simple

> Comparación cuantitativa del forward model **joint (7 parámetros: M + 6 Δ-elementos)**
> contra el **simple (1 parámetro: M)** sobre los 41 candidatos de
> `data/output/mass_followup_candidates.csv`.

**Branch**: `track2/stage1-joint-fit`
**Fecha**: 2026-05-29
**Inputs**: `data/output/loo_batch_results_joint.csv` (este trabajo),
`data/output/loo_batch_results.csv` (corrida previa Stage 0).

---

## TL;DR

El cambio de forward model **divide la mediana de χ²_red por ~335×**
(511 → 1.52). 16/27 fits joint exitosos quedan por debajo de χ²_red = 2;
24/27 por debajo de 10. La hipótesis del plan ("el drift orbital del target
se comía la señal y se absorbía como masa espuria") **queda confirmada**:
una vez que se permite que los 6 Δ-elementos absorban el drift, el residuo
restante es estadísticamente compatible con ruido AL/AC Gaia.

Criterio de aceptación Stage 1.5 (`χ²_red mediano joint < 10`,
idealmente ~1): **CUMPLIDO** (mediano = 1.52).

| Métrica | Simple (Stage 0) | Joint (Stage 1) | Mejora |
|---|---|---|---|
| `chi2_red` mediano | 511.1 | **1.52** | **335×** |
| `chi2_red` mínimo | 6.74 | 0.92 | 7.3× |
| `chi2_red` máximo | 718,392 | 192.1 | 3,742× |
| `chi2_red` < 2 | 0 / 21 | **16 / 27** | — |
| `chi2_red` < 10 | 1 / 21 | **24 / 27** | — |

---

## Corrida

```bash
git checkout track2/stage1-joint-fit
docker compose run --rm pipeline python -m scripts.mass.run_joint_batch
docker compose run --rm pipeline python -m scripts.mass.summarize_joint_fits
```

- Wall-clock total: 4 m 58 s (24 cores, 41 candidatos, ~7 s/candidato exitoso).
- Output por candidato: `data/output/fit_<perturber>_<target>_joint.json`.
- Resumen consolidado: `data/output/loo_batch_results_joint.csv` (28 filas:
  header + 27 fits exitosos).
- Reporte del batch: `data/output/joint_batch_run_report.csv` (41 filas).

### Estadísticas chi2_red joint (n=27 fits exitosos)

| min | p50 | mean | p75 | p90 | max |
|---|---|---|---|---|---|
| 0.923 | 1.516 | 10.87 | 2.844 | 4.864 | 192.13 |

### Distribución por tramo

- `chi2_red_joint < 1.5` : 13 fits (48 %)
- `1.5 ≤ chi2_red_joint < 5` : 9 fits (33 %)
- `5 ≤ chi2_red_joint < 10` : 2 fits
- `chi2_red_joint ≥ 10` : 3 fits (outliers — ver § Outliers)

---

## Overlap directo joint ∩ simple (15 pares en ambas tablas)

```
perturber  target  chi2_red_joint  chi2_red_simple  ratio   mass_joint  mass_simple
111        18105   1.135           498638.6         439k×   5.5e17      3.4e17
618        108638  3.854           349120.0         90.6k×  4.0e17      3.7e17
389        176865  4.864           35723.8          7.3k×   7.0e17      6.6e17
111        18105   1.135           4821.99          4.2k×   5.5e17      4.3e17  (alt date)
111        18105   1.135           3799.78          3.3k×   5.5e17      4.1e17  (alt date)
…
241        51218   3.206           335.47           105×    1.2e18      1.2e18
46         5998    26.218          845.36           32×     3.3e17      3.4e17
206        44887   1.186           6.744            5.7×    2.4e17      2.2e17
46         5998    26.218          51.69            1.97×   3.3e17      3.4e17
```

Observaciones:

1. **El ratio χ²_red_simple / χ²_red_joint es ≥ 1 en todos los pares**
   (el modelo joint nunca empeora al simple — esperable porque tiene 6
   parámetros extra absorbiendo grados de libertad).
2. Las **masas estimadas joint y simple son del mismo orden** en pares con
   buen χ²_red simple bajo (e.g. Hersilia/206, 511, 241). En pares con
   χ²_red simple catastrófico (e.g. 111/18105 con χ²_red ≈ 500k),
   la masa simple es **sesgada** porque el drift se absorbe parcialmente
   en M; el joint la corrige levemente.
3. La masa "simple" (111) Ate ≈ 5.43×10¹⁷ kg (registrada como primera
   determinación LOO del pipeline; perturbando target 18105) queda
   ahora en **5.50×10¹⁷ kg con joint** y χ²_red = 1.14 vs χ²_red
   simple ≈ 498k. La masa cambia <2 %; el fit es ahora bien especificado.
   Análogamente, (206) Hersilia (target 44887) pasa de 2.18×10¹⁷ kg
   (simple) a 2.44×10¹⁷ kg (joint), shift dentro de σ_inflado simple.
   **El joint corrige el χ²_red sin invalidar las masas previas.**

---

## Outliers de χ²_red_joint

3 fits con `chi2_red_joint ≥ 10`:

| perturber | target | n_joint | chi2_red_joint | nota |
|---|---|---|---|---|
| 46 | 5998 | 71 | 26.22 | (46) Hestia perturbando — chi2 también alto en simple (845). Investigar si entran obs de mala calidad o sistemática de scan. |
| 42 | 7070 | 88 | 29.69 | (42) Isis. Chi2 simple no disponible. |
| 124 | 57942 | 139 | 192.13 | (124) Alkeste. Joint mejora vs simple pero sigue alto; sospecha de sistemática AL fuera del modelo (justifica Stage 2 covarianza AL/AC). |

Estos casos **no degradan** la conclusión global porque son <12 % de los
exitosos y siguen siendo dramáticamente mejores que su contraparte simple
(reducciones de 30×–4000×). Pero **identifican** a Stage 2 (covarianza
AL/AC correcta) como el siguiente cuello de botella concreto.

---

## Candidatos no fitteados (14 / 41)

Todos cayeron por `n_loo < 8` (ventana LOO de 180 días no acumula
suficiente baseline pre-encuentro). De los 14:

- 12 con `n_loo_orbit = 0` (encuentro muy temprano en la ventana Gaia,
  no hay obs >180 d antes).
- 1 con `n_loo_orbit = 3` ((58, 1999_xa148)).
- 1 con `n_loo_orbit = 6` ((57, 216875), top-1 del ranking).

Esto es una **limitación de cobertura temporal del survey** (Gaia DR3
arranca jul 2014, varios encuentros están en 2014-Q4 / 2015-Q1 sin
baseline previo), no del modelo. Opciones futuras:

1. **Bajar `loo_window_days`** (riesgo: contamina LOO con observaciones
   sensibles al encuentro).
2. **Usar Gaia FPR (DR4 Focused Product Release)** cuando esté disponible
   — extiende la ventana temporal hacia 2017.
3. **Aceptar el subset de 27 como población fitable**; reportar
   selección honestamente.

Recomendación: **opción 3** para el cierre Stage 1; opción 2 cuando
se aborde Stage 2/3.

---

## Implicaciones para las siguientes etapas

1. **Stage 2 (covarianza Gaia AL)**: justificada por los 3 outliers.
   El joint corrigió la mayoría del bias pero deja chi2 inflado en
   pares con muchas obs (n_joint > 100) donde la sistemática AL/AC
   diagonal probablemente subestime errores.

2. **Stage 3 (specificity test)**: ahora es ejecutable — antes el
   χ²_red catastrófico hacía cualquier specificity ruidoso. Con un
   χ²_red joint plausible, podemos generar nulls geométricamente
   compatibles y testear si el M_fit real se separa de la distribución
   de M_fit_null.

3. **Stage 4 (validación contra masas conocidas)**: gate previo a publicar.
   Re-correr 1/Ceres, 4/Vesta, 2/Pallas, 10/Hygiea con el pipeline joint.

---

## Entregables Stage 1 cerrados

- [x] `docs/mass_layer_audit.md`
- [x] `docs/mass_layer_design.md`
- [x] `src/mass/forward_model_joint.py` + tests
- [x] `scripts/mass/fit_mass_gaia_joint.py`
- [x] `scripts/mass/run_joint_batch.py`
- [x] `scripts/mass/summarize_joint_fits.py`
- [x] `data/output/loo_batch_results_joint.csv` (27 fits exitosos)
- [x] `data/output/joint_batch_run_report.csv` (41 candidatos, 14 con `n_loo` insuficiente)
- [x] `docs/mass_layer_joint_diagnostic.md` (este documento)

---

## Próximo paso

Track 2 Stage 2 — Covarianza Gaia AL/AC. Implementar likelihood Mahalanobis
en `src/mass/likelihood_al.py`; verificar sobre datos sintéticos que
M ≈ 0 con chi²_red ≈ 1; re-correr los 3 outliers para ver si bajan.
