# Track A Stage 2.6 — The Hygiea "9× real-data bias": diagnosis

> Estado: **diagnóstico cerrado**. Fecha: 2026-05-30.
> Herramienta: [scripts/mass/realdata_mass_scan.py](../scripts/mass/realdata_mass_scan.py).
> Datos: `data/output/stage2_6/hygiea_*` (scan, per-obs, summary).
> Antecedentes: [mass_layer_stage_a2_5_profiled.md](mass_layer_stage_a2_5_profiled.md),
> [mass_layer_closing_loop_leverage.md](mass_layer_closing_loop_leverage.md).

## Punto de partida

Stage 2.5 arregló el bug del optimizador (closing-loop sintético recupera la
masa inyectada con ratio 1.000) pero el fit **real** de Hygiea aterrizó en
**7.5×10²⁰ kg ≈ 9× la literatura** con un χ²_red≈1.19 (buen ajuste). El veredicto
de A2.5 fue: con el optimizador y la degeneración M↔deltas descartados, queda un
**"sesgo real-data"** coherente — sospechándose deriva orbital, sistemática
astrométrica Gaia, o un perturbador secundario.

Stage 2.6 prueba esas hipótesis directamente, escaneando la **verosimilitud
perfilada χ²(masa) sobre los datos reales** — por-target individual y joint —
con los 6 deltas re-optimizados en cada masa de prueba.

## Qué se midió

`realdata_mass_scan.py` construye los bundles reales de Hygiea (5 targets con
fit OK en Stage 4: 4803, 16772, 45989, 47605, 58775) y para cada grupo recorre
una grilla de `log10_M` de ±1.6 dex alrededor de la literatura (8.3×10¹⁹ kg),
profilando los deltas en cada punto. Reporta el mínimo, el Δχ²(lit−min), la
unimodalidad, y los residuos por-observación en el óptimo.

## Resultado 1 — la χ²(masa) real es dentada y multimodal (no hay "el" mínimo en 9×)

| grupo | M_min/lit | χ²r@min | χ²r@lit | Δχ²(lit−min) | σ | unimodal |
|-------|-----------|---------|---------|--------------|---|----------|
| target_4803  | 0.03 | 0.68 | 12.32 | 1035.8 | 32.2 | **no** |
| target_16772 | 0.19 | 1.12 | 1.27 | 18.9 | 4.3 | **no** |
| target_45989 | 0.29 | 1.51 | 1129.60 | 346322 | 588.5 | **no** |
| target_47605 | 1.51 | 1.28 | 1.84 | 121.0 | 11.0 | **no** |
| target_58775 | 0.13 | 0.52 | 0.55 | 1.8 | 1.3 | **no** |
| **joint** | **0.66** | 3.77 | 433.55 | 348555 | 590.4 | **no** |

Dos hechos demuelen la interpretación de A2.5:

1. **Ningún grupo es unimodal.** El Δχ² salta entre 10⁴ y 10⁶ entre puntos de
   grilla **adyacentes** (0.18 dex). Ejemplo del joint (ratio : Δχ² sobre el
   mínimo): `0.03:92926  0.04:82323  0.06:221922  0.09:137857  0.13:724355
   0.19:14038  0.29:5135  0.44:47147  0.66:0  1.00:348555  …`. Una verosimilitud
   física —aun multimodal— no oscila así a 0.18 dex de espaciado. Es el **fit
   interno de deltas cayendo en distintos basins** según la masa: la
   degeneración M↔deltas, *activa sobre datos reales*.

2. **El "9× lit" no es un mínimo.** En el joint, la región de 9× (ratio 7.74)
   tiene Δχ²≈514320 sobre el mínimo de grilla (0.66×). El 7.5×10²⁰ de A2.5 fue
   el **Brent acotado asentándose en un basin de ruido estrecho**, no el óptimo
   global. Re-escanear da un "mínimo" en 0.66× — es decir, el resultado **depende
   del camino del optimizador y no es reproducible**. No existe un "sesgo
   real-data de 9×"; existe **no-identificabilidad**.

## Resultado 2 — un driver del χ² salvaje es deriva orbital sobre la ventana joint unilateral

Los residuos por-observación a masa literatura exponen el target 45989 como
patológico y temporalmente estructurado:

| target | n_obs | mediana χ²/obs @lit | máx χ²/obs @lit | frac χ²>9 @lit |
|--------|-------|---------------------|------------------|----------------|
| 4803  | 48 | 6.21 | 93.7 | 0.33 |
| 16772 | 69 | 1.72 | 13.3 | 0.03 |
| **45989** | 157 | **639.0** | **38550** | **0.81** |
| 47605 | 111 | 2.06 | 33.0 | 0.07 |
| 58775 | 36 | 0.79 | 4.3 | 0.00 |

45989 muestra una **rampa temporal**: mediana χ²/obs de 17 en la primera mitad
de su arco a 789 en la segunda, sobre ~814 días. La causa es estructural: la
máscara joint de `_build_bundle` es **unilateral**
(`days_from_enc > -loo_window_days`, sin cota superior), así que para encuentros
tempranos (45989, 2015-05-22) el arco post-encuentro corre hasta el final de DR3
(~2 años), acumulando deriva orbital que la órbita LOO (ajustada >180 d antes del
encuentro) no captura. A masa literatura esa deriva no se puede absorber y los
residuos explotan; el fit entonces empuja la masa a un valor donde los 6 deltas
re-absorben la deriva.

## Resultado 3 — apretar la ventana arregla la patología pero NO identifica la masa

Re-corrido con una ventana **simétrica ±60 d** alrededor del encuentro:

| grupo | n_obs (±60d) | M_min/lit | χ²r@min | χ²r@lit | σ |
|-------|--------------|-----------|---------|---------|---|
| target_4803  | 30 | 0.19 | 0.29 | 5.97 | 17.4 |
| target_16772 | 25 | 0.19 | 0.90 | 1.00 | 2.1 |
| target_45989 | 13 | 0.19 | 0.41 | **0.49** | 1.3 |
| target_47605 | 43 | 3.41 | 0.38 | 0.42 | 1.8 |
| target_58775 | 18 | 2.27 | 0.47 | 0.48 | 0.5 |
| joint | 129 | 0.44 | 0.62 | 48.71 | 104.5 |

- **45989 se cura**: χ²_red@lit cae de **1129 → 0.49**. Confirma que su blow-up
  era deriva orbital del arco ancho, no masa.
- **Pero la masa sigue sin identificarse**: los mínimos por-target siguen
  **dispersos** (0.19, 0.19, 0.19, 3.41, 2.27 × lit), las curvas siguen
  **dentadas** (uni=no), y ningún encuentro tiene leverage limpio ≳3σ apuntando
  a literatura. χ²r@min < 1 a masas tan distintas como 0.19× y 3.4× ⇒ los deltas
  absorben la diferencia: **degeneración M↔deltas activa**.

## Resultado 4 — perturbador secundario descartado

Consulta sobre `encounters_catalog_hybrid_stageb.parquet`: para los 5 targets,
encuentros con **otros** cuerpos a <0.3 AU dentro de ±90 d del encuentro con
Hygiea. 320 candidatos, **ninguno** en el registro de asteroides masivos
(`_MAJOR_ASTEROIDS`). Los partners más cercanos (0.003–0.008 AU) son cuerpos
chicos (números >100k, masas ~10¹⁴–10¹⁵ kg, ≥4 órdenes bajo Hygiea): deflexión
despreciable. No hay perturbador secundario que explique el sesgo.

## Reconciliación con A2 y el closing-loop

A2 concluyó que "los deltas son ~10⁻⁷ y no absorben señal" — pero eso fue sobre
**datos sintéticos sin ruido** (deltas verdaderos = 0). Sobre **datos reales**,
con deriva orbital y sistemáticas astrométricas genuinas, los 6 deltas por target
**sí se mueven y absorben** estructura, y lo hacen de forma masa-dependiente. Por
eso la χ²(masa) real es dentada y la masa no se identifica, aun cuando el
closing-loop sintético la recuperaba exacta. El leverage intrínseco medido por
`probe_mass_sensitivity` (Ceres 25σ, Hygiea ~2σ, Pallas ~0.09σ) ya anticipaba
esto: Hygiea está en el régimen marginal donde el ruido/deriva real domina sobre
la señal de deflexión.

## Veredicto

El "sesgo real-data de 9×" de A2.5 **no era un sesgo físico coherente**. Era la
superposición de (a) un **basin espurio** del optimizador sobre una χ²(masa)
**dentada y multimodal** (no-identificabilidad: los deltas absorben masa sobre
datos reales), y (b) **contaminación por deriva orbital** del arco joint
unilateral, dominada por 45989.

Conclusión honesta: **la masa de Hygiea no es determinable con DR3 en este
pipeline.** Converge con el veredicto del deepwork ("capa no publicable") pero
con el mecanismo preciso. El único caso con leverage fuerte sigue siendo **Ceres**
(25σ, encuentro 18937); el joint multi-target **agrega encuentros débiles que
inyectan ruido, no restricción**.

## Implicaciones accionables

1. **La ventana joint debe ser simétrica y acotada** alrededor del encuentro
   (p. ej. ±60–90 d). La actual es unilateral y deja entrar deriva orbital. Es un
   fix barato que limpia 45989; aplicar en `_build_bundle` si se sigue usando el
   joint.
2. **No reportar masas multi-target de Hygiea/Pallas.** La σ por curvatura es
   engañosa sobre una superficie dentada. Solo perturbadores con leverage ≳10σ
   medido por `probe_mass_sensitivity` (hoy solo Ceres) son defendibles.
3. **El leverage, no el optimizador, es ahora el cuello de botella.** A3 (modelo
   OU para el drift) reduciría los dof de deltas y podría destrabar la
   degeneración M↔deltas sobre datos reales — pero solo vale la pena si antes se
   confirma, con la ventana acotada, que algún calibrador además de Ceres tiene
   señal ≳3σ. Hoy la evidencia dice que no en DR3.
