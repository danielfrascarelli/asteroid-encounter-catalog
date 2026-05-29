# Stage 3 — Specificity test diagnóstico

> Test riguroso: ¿el fit joint+Mahalanobis es estadísticamente sensible
> al perturber real, o produce la misma "masa" para cualquier asteroide
> de tamaño similar en la misma banda orbital?

**Branch**: `track2/stage3-specificity`
**Fecha**: 2026-05-29
**Inputs**:
- 5 candidatos top del Stage 2 (`χ²_red_joint_mahal < 0.5`).
- 50 null perturbers por candidato, sampleados de MPCORB
  (`MPCORB_20160217.DAT`) con `|a - a_real| ≤ 0.5 AU`,
  `|H - H_real| ≤ 1.5 mag`, **sin** encuentro <0.1 AU contra el target
  en el catálogo `encounters_catalog_hybrid_stageb.parquet`.
**Output**: `data/output/specificity_test_v2.csv`,
`data/output/specificity_v2/specificity_<perturber>_<target>.json`.

---

## TL;DR

**2 / 5 candidatos pasan specificity en χ²_red** (p ≤ 0.05).
**0 / 5 candidatos pasan specificity en masa** (p > 0.15 para todos).

Para las masas que el pipeline reporta para (111) Ate, (206) Hersilia,
(124) Alkeste → 3294 — perturbers de las primeras tres determinaciones
del pipeline — el χ² del fit real es estadísticamente indistinguible
del χ² promedio que se obtiene fiteando un asteroide aleatorio del
mismo tamaño y banda orbital. La señal de masa **no es específica**.

Donde sí hay specificity en χ²: (19) Fortuna → 13346 (p_χ² = 0.04) y
(49) Pales → 94474 (p_χ² = 0.02). Ambos siguen mostrando masas
consistentes con los nulls (p_mass = 0.18 y 0.26) — el fit detecta
"algún tipo de señal" en la ventana pero no extrae una masa
distintiva.

Conclusión científica honesta: con la dataset actual y el modelo joint
7-param + Mahalanobis 2D, **no podemos atribuir una masa cuantitativa
a un perturber específico salvo, marginalmente, a Fortuna y Pales —
y aún así sólo en términos de χ², no de magnitud de masa**.

---

## Tabla principal (N=50 nulls/candidato)

```
target  real  chi2_real  chi2_null_p10  chi2_null_med  chi2_null_p90  mass_real     mass_null_med   p_chi2  p_mass
18105   111   0.427      0.419          0.421          0.439          6.02e17       2.52e17         0.80    0.26
44887   206   0.347      0.341          0.348          0.356          2.63e17       1.35e17         0.40    0.32
94474    49   0.343      0.344          0.363          0.460          7.40e17       3.46e17        *0.02*   0.18
13346    19   0.179      0.411          1.200          3.282          1.82e18       8.44e17        *0.04*   0.26
 3294   124   0.225      0.143          0.203          0.263          5.02e17       2.35e17         0.64    0.26
```

`p_chi2` = fracción de nulls con `χ²_red < χ²_real` (más bajo = mejor
specificity).
`p_mass` = fracción de nulls con `M_fit > M_real` (más bajo = real
está en el extremo alto de la distribución null).

---

## Lectura caso por caso

### (111) Ate → (18105): NO específico

`χ²_real = 0.427`, `null_med = 0.421` — virtualmente idéntico.
La distribución null está **muy concentrada** (p10 = 0.419, p90 = 0.439).
Esto significa que **cualquier asteroide en la banda orbital de Ate y
de tamaño similar produce un fit con χ²_red ≈ 0.42** sobre los obs
del target 18105. La masa real `6.02e17 kg` está apenas en el percentil
74 de las nulls (p_mass = 0.26) — masa real superior a la mediana pero
con 26% de nulls reportando una masa **mayor**.

**Implicación**: la determinación de masa de Ate **no se sostiene**
contra specificity. Hay que retractar la afirmación cuantitativa.

### (206) Hersilia → (44887): NO específico

`χ²_real = 0.347`, `null_med = 0.348` — peor que la mediana del null.
La masa real `2.63e17 kg` está en el percentil 68. Argumentos análogos
a Ate: el modelo absorbe astrometría comparable con cualquier perturber.

### (124) Alkeste → (3294): NO específico

`χ²_real = 0.225`, `null_med = 0.203` — el real es **peor** que la
mediana null. El mass_real `5.02e17 kg` está en el percentil 74. No hay
caso para una masa específica.

### (49) Pales → (94474): específico en χ²

`χ²_real = 0.343`, `null_med = 0.363`. El real cae en el **2 %** más
bajo de los nulls (p = 0.02 con N=50, resolución 0.02). Aunque la
diferencia con el null típico es pequeña en magnitud (~5 %), las
distribuciones null tienen una **cola larga** (p90 = 0.460) — el
modelo lucha con muchos nulls pero ajusta limpio con Pales. Esa
discriminación es señal real de geometría.

Sin embargo, `mass_real = 7.40e17 kg` está en el percentil 82 (p_mass
= 0.18) — la magnitud específica de la masa **todavía no se distingue**
del cluster null.

**Lectura**: Pales **sí** está deflectando al target (señal genuina en
χ²) pero la magnitud de masa absoluta no es extraíble con confianza.

### (19) Fortuna → (13346): específico en χ²

`χ²_real = 0.179`, `null_med = 1.200`. Diferencia dramática: la
distribución null es ancha (`p10 = 0.41, p90 = 3.28`) — la mayoría de
los nulls **no pueden ajustar** los obs del target 13346. Pero Fortuna
sí. p_χ² = 0.04 (2/50 nulls hicieron mejor).

`mass_real = 1.82e18 kg` está en el percentil 74 (p_mass = 0.26) — más
alto que la mediana null (`8.44e17 kg`) pero todavía 13/50 nulls
reportan más masa.

**Lectura**: Fortuna es el caso más fuerte. Hay deflección detectable
en los obs y el ajuste es claramente mejor con el perturber correcto.
La masa estimada `1.82e18 kg` no es estadísticamente única pero sí
está en el lado correcto de la distribución.

---

## Por qué los nulls también producen "masa"

El joint fit tiene **7 parámetros** (1 masa + 6 deltas orbitales).
Con 50–270 observaciones Gaia por candidato, el modelo tiene amplio
margen para encontrar una solución que minimice χ² **deformando la
órbita del target** y atribuyendo parte del residuo a una "masa"
del perturber, incluso si el perturber nunca estuvo cerca.

Los 6 deltas absorben drift orbital (la motivación original de Stage
1) y, en el proceso, también pueden absorber pequeños sesgos que
correlacionan con la geometría del null perturber. Si el null tiene
una órbita "razonablemente parecida" a la del real, el optimizer
encuentra una combinación (mass, deltas) que parece coherente.

Esto **no** es un bug del optimizer; es el comportamiento normal de
un modelo con 7 parámetros sobre datos ruidosos. Pero exige un test
de specificity como éste para no confundir ajuste con detección.

---

## Lo que el specificity test no captura

1. **Geometría 3D real**: dos asteroides pueden estar en la misma
   banda orbital pero a 1 AU de distancia 3D durante toda la ventana
   Gaia. El criterio actual (catálogo de encuentros, `dist < 0.1 AU`
   = "real encounter") excluye solamente lo que el catálogo ya
   identifica como cercano. Asteroides con distancia mínima 0.1–0.5
   AU son nulls "limpios" — y son los que producen `mass_null`
   intermedio.
2. **Forma de la cola de masa**: el test reporta percentiles, no la
   forma de la cola. Una distribución con masa típica `~2e17 kg`
   con cola pesada hasta `~1e18` puede tener percentil 26 idéntico
   al de Ate aunque tenga otra densidad cualitativa.
3. **Correlación de fits**: 50 fits independientes no son
   bootstrap — cada null es un asteroide real con su propia
   geometría. La distribución null es **una estimación empírica del
   nivel de absorción del modelo**, no un nulo gaussiano.

---

## Decisiones que esto fuerza

1. **No publicar** masas para (111) Ate, (206) Hersilia, (124)
   Alkeste/3294 sin specificity adicional. La memoria
   [project_ate_mass_result.md] (5.43e17 kg) y los registros previos
   son útiles internamente pero **no defendibles como detección**.
2. **Considerar publicables** los **χ²-significant**: (19) Fortuna
   y (49) Pales, con caveats:
   - La magnitud absoluta de masa requiere validación cross-track
     (Stage 4 vs literatura de DAWN / Goffin / Fienga).
   - Reportar `M_fit ± σ_fit` junto con `p_specificity_χ²` y la
     mediana null.
3. **Stage 4 (validación contra literatura)** es gate **obligatorio**
   antes de publicar cualquier masa. Si el pipeline ni siquiera
   reproduce Ceres/Vesta/Pallas/Hygiea dentro de 3σ, el problema
   metodológico es más profundo y posiblemente el catálogo Gaia DR3
   sin DR4 no es suficiente para mass-fits de calidad publicable.
4. **Investigar específicamente Fortuna y Pales**: ambos tienen
   `χ²_real < χ²_null_p10`, lo que significa que el optimizer encontró
   un mínimo que **no logra reproducir con perturbers comparables**.
   Eso es señal real — pero requiere refinamiento del null distribution
   (N=200+, banda orbital más estrecha) para reportar p-values
   convergidos.

---

## Próximos pasos

- Stage 4: re-fittear (1) Ceres, (4) Vesta, (2) Pallas, (10) Hygiea
  con el pipeline joint+Mahalanobis 2D y comparar contra masas DAWN /
  Goffin / Vernazza dentro de 3σ. Gate antes de cualquier
  publicación.
- Extender specificity al batch completo de 27 fits (no sólo los 5
  con mejor χ²) — calcular p-values para todos y reportar
  selectividad estadística honesta.
- Si Stage 4 pasa: re-correr specificity con `N=200` nulls sobre
  Fortuna y Pales para nivel de confianza publicable.
- Si Stage 4 falla: documentar como limitación fundamental del
  dataset DR3 y esperar DR4 / FPR.

---

## Entregables Stage 3

- [x] `src/mass/null_perturbers.py` + `tests/test_null_perturbers.py`
  (7 tests)
- [x] `scripts/mass/run_specificity_test.py`
- [x] `data/output/specificity_test_v2.csv`
- [x] `data/output/specificity_v2/` (per-candidate JSON + per-fit)
- [x] `docs/mass_layer_stage3_diagnostic.md` (este documento)
