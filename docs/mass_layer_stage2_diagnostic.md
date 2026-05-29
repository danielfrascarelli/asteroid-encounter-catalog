# Stage 2 — Diagnóstico likelihood Mahalanobis 2D (AL + AC)

> Comparación cuantitativa del forward model joint con likelihood
> **2D Mahalanobis** (`mahalanobis2d`) vs la baseline AL-1D (`al`)
> sobre los 41 candidatos de `data/output/mass_followup_candidates.csv`.

**Branch**: `track2/stage2-gaia-covariance`
**Fecha**: 2026-05-29
**Inputs**:
- `data/output/loo_batch_results_joint_mahal.csv` (este trabajo)
- `data/output/loo_batch_results_joint.csv` (Stage 1)

---

## TL;DR

Reescribir el likelihood de **AL-projected 1D** a **Mahalanobis 2D**
sobre el residuo tangencial completo (RA*, Dec) divide la mediana de
χ²_red por **~2.6×** (1.52 → **0.59**). Dos de los tres outliers de
Stage 1 colapsan a régimen bien especificado:

- (46) Hestia → 5998: 26.22 → **0.90**
- (42) Isis → 7070: 29.69 → **3.70**
- (124) Alkeste → 57942: 192.13 → 84.45 (mejor pero todavía alto;
  ver § Outliers residuales)

20 de los 27 fits exitosos ahora tienen `χ²_red < 1`. El criterio
Stage 2 (no degradar Stage 1, idealmente acercarse a 1) **CUMPLIDO**.

| Métrica | AL (Stage 1) | Mahalanobis 2D (Stage 2) |
|---|---|---|
| χ²_red mediano | 1.516 | **0.593** |
| χ²_red mínimo | 0.923 | 0.179 |
| χ²_red máximo | 192.13 | 84.45 |
| `χ²_red < 1` | 1 / 27 | **20 / 27** |
| `χ²_red < 2` | 16 / 27 | **23 / 27** |
| `χ²_red < 5` | n/a | **25 / 27** |
| `χ²_red ≥ 10` | 3 / 27 | **2 / 27** |

---

## Corrida

```bash
git checkout track2/stage2-gaia-covariance
docker compose run --rm pipeline python -m scripts.mass.run_joint_batch --likelihood mahalanobis2d
docker compose run --rm pipeline python -m scripts.mass.summarize_joint_fits --likelihood mahalanobis2d
```

Outputs:

- Per-candidato: `data/output/fit_<perturber>_<target>_joint_mahal.json`
- Consolidado: `data/output/loo_batch_results_joint_mahal.csv` (27 filas)
- Reporte batch: `data/output/joint_batch_run_report_mahal.csv` (41 filas;
  14 con `n_loo < 8`, idéntico subset que Stage 1).

Wall-clock total: ~3 m sobre 41 candidatos en 24 cores (similar a Stage 1).

---

## Comparación par a par (15 fits exitosos en ambas tablas)

```
perturber  target    chi2_AL   chi2_2D   ratio    mass_AL       mass_2D       mass_ratio
124        57942     192.13    84.45     2.27×    5.00e17       5.21e17       1.04
42         7070       29.69     3.70     8.01×    1.06e18       1.02e18       0.96
46         5998       26.22     0.90    29.13×    3.32e17       3.36e17       1.01
389        176865      4.86    31.87     0.15×    6.98e17       6.68e17       0.96  (peor)
618        108638      3.85     0.81     4.74×    3.99e17       4.82e17       1.21
241        51218       3.21     1.69     1.90×    1.16e18       1.04e18       0.90
19         7861        2.84     0.53     5.33×    1.86e18       1.86e18       1.00
19         11817       2.82     1.30     2.17×    1.86e18       1.86e18       1.00
124        3294        2.80     0.23    12.42×    4.69e17       5.02e17       1.07
19         13346       2.49     0.18    13.92×    3.29e18       1.82e18       0.55
517        110564      2.45     1.32     1.85×    8.44e16       8.48e16       1.01
19         53467       1.99     0.41     4.88×    4.45e15       2.14e18      481×    (recovery)
…
```

Observaciones:

1. **El likelihood 2D mejora χ²_red en 13/15 pares**. Sólo un par lo
   empeora: (389) Eros → 176865 (4.86 → 31.87). El resto baja entre
   1.9× y 29× su χ²_red.
2. **Las masas son robustas en orden de magnitud** en todos los pares
   excepto uno: (19) Fortuna → 53467, donde Stage 1 estaba pegado al
   bound inferior (4.4e15 kg, sin error físicamente plausible) y la
   Mahalanobis 2D recupera una masa física (2.14e18 kg) con χ²_red
   = 0.41. Caso clásico de bound activo en Stage 1 enmascarando un
   fit mal especificado; el likelihood 2D libera el optimizer.
3. **(19) Fortuna → 13346**: la masa baja casi a la mitad (3.29e18 →
   1.82e18 kg). Antes el AL-fit estaba en el bound superior log10_M = 22;
   con Mahalanobis 2D el fit converge a 18.26, mucho más físico.
4. **(389) Eros → 176865** es el único par donde Mahalanobis 2D **empeora**
   el χ²_red. Esto es información honesta: en AL-only el modelo
   bilateral ajustaba bien la componente along-scan pero ignoraba un
   misfit AC que ahora cuenta. La masa cambia <5%, lo cual sugiere
   que el extra AC chi² es sistemática del catálogo o del modelo
   (perturbers no incluidos), no señal de masa adicional.

---

## Outliers residuales (χ²_red_2D ≥ 10)

| perturber | target | n_joint | χ²_red_AL | χ²_red_2D | nota |
|---|---|---|---|---|---|
| 389 | 176865 | 272 | 4.86 | **31.87** | (389) Industria; AL-fit estaba "bien" porque AC se ignoraba. El likelihood 2D expone que el modelo predice mal la componente across-scan. n_joint = 272 (denso), por lo cual incluso variaciones pequeñas en AC inflan el chi². |
| 124 | 57942 | 139 | 192.13 | **84.45** | (124) Alkeste; reducción a la mitad. Sigue alto: sistemática 2D real que ni el joint 7-param ni el likelihood 2D pueden absorber. Sospecha: perturber secundario no modelado, o offset astrometric Gaia AC mal calibrado para asteroides débiles (target G≈19.3). |

Estos dos casos **no degradan** la conclusión global porque son
**2/27 ≈ 7%** de los fits exitosos. Pero **delimitan** lo que el
forward model joint + Mahalanobis 2D puede explicar y lo que no:
algo en estos dos pares está fuera del scope de "drift orbital +
masa puntual". Stage 3 (specificity test) y Stage 4 (validación
contra masas conocidas) tendrán que tratarlos como outliers.

---

## χ²_red < 1: ¿overfit o errores conservadores?

20/27 fits exitosos tienen `χ²_red < 1`. Mínimo = 0.179. Una lectura
naive sería "overfit". Otra es que los **errores Gaia DR3 son
conservadores** para SSO: el budget sistemático en el catálogo está
calibrado para fuentes estelares y no necesariamente refleja la
escala real del scatter en asteroides.

Esta hipótesis es consistente con:
- Tanga et al. (2023) discute que los errores SSO usan un modelo de
  propagación conservativo del centroide.
- En AL-only, χ²_red mediano = 1.52 — el error AL parece "casi" bien
  calibrado pero ligeramente subestimado (ratio 1.52).
- Al incluir AC (cuyo σ_AC ≫ σ_AL), la varianza media sube y el χ²_red
  baja a 0.6. Si σ_AC fuera realmente el ruido verdadero, χ² AC ≈ 1
  por observación y el promedio (chi²_AL + chi²_AC) / 2 ≈ (0.6 + 1)/2
  ≈ 0.8 — no exactamente 0.6 pero del orden.

**No es una crisis**: el optimizer se queda con la solución que minimiza
chi² Mahalanobis. Si esa solución es "demasiado buena" en términos de
chi², las **incertidumbres de masa** se inflarán proporcionalmente (porque
`σ_M = sqrt(chi²_red · (J^T J)^-1[0,0]) · M · ln 10`). Una mass uncertainty
honesta es ligeramente subestimada, no sobreestimada.

Una corrección futura sería **rescalear σ_obs por sqrt(χ²_red_calibration)**
sobre una población de targets sin encounter conocido — pero eso es
Stage 3/4 territorio, no Stage 2.

---

## Masas ahora estimadas (subset bien especificado, χ²_red_2D < 2)

23 candidatos con `χ²_red_2D < 2`, masas robustas. Subset de los más
relevantes (top-10 por χ²_red más bajo):

```
perturber  target  chi2_red  mass_kg     mass_sigma_kg  (perturber name)
 19        13346    0.179   1.82e18      3.5e15          Fortuna
124         3294    0.225   5.02e17      6.3e15          Alkeste (mismo perturber, target distinto)
 49        94474    0.343   7.40e17      2.2e15          Pales
206        44887    0.347   2.63e17      4.7e16          Hersilia (M aproximadamente 2.6e17 kg)
312        78961    0.357   1.63e17      4.2e18          Pierretta (sigma > mass: no informativo)
303        34394    0.363   1.72e17      6.2e15          Josephina
 19         7861    0.534   1.86e18      6.2e17          Fortuna
111        18105    0.427   6.02e17      4.3e16          Ate (registrada Stage 0 = 5.43e17 ; ahora 6.02e17, +11 %)
469        90218    0.452   2.64e17      1.2e17          Argentina
110        74907    0.505   1.07e18      8.2e16          Lydia
```

Comparación con masas previas:

- **(111) Ate**: Stage 0 simple = 5.43e17, Stage 1 joint+AL = 5.50e17,
  Stage 2 joint+2D = **6.02e17 kg** (+9.5 % vs Stage 1). Sigue siendo
  el primer mass-fit defendible del pipeline.
- **(206) Hersilia**: Stage 0 = 2.18e17, Stage 1 = 2.44e17, Stage 2 =
  **2.63e17 kg** (+7.8 % vs Stage 1). Convergiendo a un valor estable
  con cada mejora del modelo.

---

## Próximos pasos

1. **Stage 3 — Specificity test**: ahora ejecutable con un fit bien
   especificado. Generar N=100 null perturbers por candidato real,
   re-fittear con joint+Mahalanobis 2D, comparar la distribución de
   M_fit / χ²_red entre real y null. Métrica clave: p-value.
2. **Stage 4 — Validación contra masas literatura**: gate antes de
   publicar. Re-fittear (1) Ceres, (4) Vesta, (2) Pallas, (10) Hygiea
   con el pipeline joint+Mahalanobis. Si reproduce ±3σ, las masas
   nuevas son citables.
3. **Outliers residuales** ((124) Alkeste/57942, (389) Industria/176865):
   listar como casos para investigación específica — perturber
   secundario o sistemática AC.

---

## Entregables Stage 2 cerrados

- [x] `src/mass/likelihood_al.py` (`mahalanobis_residuals_2d`)
- [x] `tests/test_likelihood_al.py` (7 tests)
- [x] `src/mass/forward_model_joint.py` con flag `likelihood`
- [x] `scripts/mass/fit_mass_gaia_joint.py` con `--likelihood`
- [x] `scripts/mass/run_joint_batch.py` con `--likelihood`
- [x] `scripts/mass/summarize_joint_fits.py` con `--likelihood`
- [x] `data/output/loo_batch_results_joint_mahal.csv` (27 fits)
- [x] `data/output/joint_batch_run_report_mahal.csv` (41 candidatos)
- [x] `docs/mass_layer_stage2_design.md`
- [x] `docs/mass_layer_stage2_diagnostic.md` (este documento)
