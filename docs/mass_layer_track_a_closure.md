# Track A — cierre: ninguna masa de perturbador es determinable en DR3 (incluido Ceres)

> ⚠️ **Veredicto acotado al MÉTODO de este doc (LOO secuencial), no al problema.**
> El motor `src/orbdet/` (joint orbit+mass sobre el arco completo, junio 2026)
> **sí determina masas** sobre Gaia FPR: 4/4 calibradores |z|<3 + (16) Psyche
> 2.43×10¹⁹ kg ±3.3%. Lo que fallaba era el ajuste secuencial órbita→masa, no el
> leverage de Gaia. Resultado actual:
> [`docs/mass_determination_results.md`](mass_determination_results.md). El análisis
> de abajo sigue siendo válido como diagnóstico de por qué el LOO no es identificable.

> Documento de cierre de **Track A** del follow-up post-deepwork (disuelto; ver [ROADMAP.md](../ROADMAP.md) § "Estado actual").
> Consolida el veredicto tras A1 → A2 → A2.5 → A2.6 y el **gate-check de A3**
> (este documento). Estado resultante de Track A: **⚫ CERRADO** — esperar DR4/FPR.

## TL;DR

El gate-check de A3 (escaneo de la χ²(masa) perfilada con **ventana joint
acotada ±60 d** sobre los calibradores ≠ Hygiea) **falla**: ningún calibrador
distinto de Ceres muestra leverage ≳3σ. Por lo tanto **no se arranca A3** (el
forward model OU del drift). El cuello de botella es el **leverage intrínseco
del dataset DR3**, no la parametrización del drift — y A3 no agrega leverage,
sólo restringe los grados de libertad del drift. La conclusión honesta:

> **Con DR3 y este pipeline, NINGUNA masa de perturbador es determinable —
> incluido Ceres.** Ceres tiene el leverage χ² más fuerte de la muestra (σ~20),
> pero un experimento de sensibilidad a la ventana (sección "Ceres tampoco se
> salva") muestra que su mínimo de masa perfilado **vaga ~100× (0.04–4.37×
> literatura) según la ventana de ajuste, sin plateau, y la χ²(masa) es
> multimodal en todas las ventanas**. El leverage es real pero no se traduce en
> una masa identificable. El negativo es universal y limpio.

## Por qué este gate-check decide A3

A3 propone reemplazar los 6 Δ-elementos absolutos por un proceso
Ornstein-Uhlenbeck `(τ, σ)` para reducir los dof del drift y evitar que absorba
señal de masa. Pero A2.6 ya mostró que el problema real es que **la χ²(masa)
sobre datos reales es dentada/multimodal** (degeneración M↔deltas activa) y que
**Hygiea no tiene leverage** (~2σ). A3 sólo vale la pena si **algún otro
calibrador**, con la ventana ya acotada (`--joint-window-days`, entregable de
A2.6), muestra una señal de masa limpia ≳3σ que A3 pudiera "limpiar". El
gate-check mide exactamente eso.

## Método

`scripts/mass/realdata_mass_scan.py` con `--joint-window-days 60` sobre los
calibradores con targets disponibles en `data/output/stage4_validation_summary.csv`:

```bash
# Pallas (el gate real: 5 targets en Stage 4)
docker compose run --rm pipeline python -m scripts.mass.realdata_mass_scan \
  --perturber 2 --targets-csv data/output/stage4_validation_summary.csv \
  --lit-mass-kg 2.05e20 --joint-window-days 60 \
  --out-prefix data/output/stage2_6/pallas_w60

# Ceres (referencia de leverage conocido)
docker compose run --rm pipeline python -m scripts.mass.realdata_mass_scan \
  --perturber 1 --targets-csv data/output/stage4_validation_summary.csv \
  --lit-mass-kg 4.71e20 --joint-window-days 60 \
  --out-prefix data/output/stage2_6/ceres_w60
```

- **Vesta** (perturber 4): 0 targets con fit exitoso en Stage 4 → no aplica.
- **Hygiea** (perturber 10): ya corrido en A2.6
  ([docs/mass_layer_stage_a2_6_realdata_bias.md](mass_layer_stage_a2_6_realdata_bias.md)).

`σ` = `sqrt(Δχ²(lit−min))` es la preferencia (en sigmas) del mínimo de grilla
frente a la masa de literatura; `local_min` cuenta mínimos locales interiores en
la grilla (un proxy de cuán dentada/multimodal es la curva).

## Resultados

Con la ventana ±60 d, muchos targets pierden observaciones y caen por debajo del
piso de 8 obs por lado; sobreviven los que tienen arco ddenso cerca del encuentro.

| Perturber | Grupo | n_obs | M_min/lit | χ²ᵣ@min | χ²ᵣ@lit | Δχ²(lit−min) | σ | min. locales |
|-----------|-------|------:|----------:|--------:|--------:|-------------:|----:|----:|
| **Ceres** (1) | target_18937 | 14 | **4.37** | 0.80 | 20.43 | 412.2 | **20.3** | 7 |
| **Pallas** (2) | target_28036 | 18 | 0.33 | 0.63 | 0.68 | 1.5 | 1.2 | 8 |
| Pallas (2) | target_59882 | 8 | 13.18 | 0.56 | 0.61 | 0.5 | 0.7 | 4 |
| Pallas (2) | joint | 26 | 0.08 | 0.60 | 0.62 | 1.0 | 1.0 | 6 |
| Hygiea (10) | (A2.6, ±60 d) | 13–43 | 0.19–3.41 | — | — | — | ≤2.1 | — |

Lectura:

1. **Pallas no tiene leverage** con ventana acotada: σ ≤ 1.2 en todos los grupos,
   Δχ² ≤ 1.5, y los mínimos por-target están **dispersos** (0.33×, 13.18×) sin
   coherencia. El joint cae en 0.08× lit con σ=1.0. Es estadísticamente plano.
2. **Hygiea tampoco** (A2.6): per-target σ ≤ 2.1, mínimos 0.19–3.4× lit, curvas
   dentadas.
3. **Ceres tiene leverage fuerte** (σ=20.3 ≫ 3σ) pero su mínimo cae en **4.37×
   lit** con ±60 d — ver el experimento dedicado abajo, que muestra que ese
   leverage **no** produce una masa identificable.
4. **Todas las curvas son multimodales** (7/8/4/6 mínimos locales interiores),
   incluida la de Ceres. Confirma el mecanismo de A2.6: sobre datos reales la
   degeneración M↔deltas vuelve la χ²(masa) dentada. Lo que distingue a Ceres no
   es una curva limpia sino un **basin profundo pero dentado** que hace al σ
   engañoso (ver abajo).

## Ceres tampoco se salva — sensibilidad a la ventana

El primer escaneo dejó la duda: Ceres tiene σ=20 pero su mínimo (4.37×) no
coincide con literatura. ¿Es un sesgo fijo corregible, o el estimador es
inestable? Para decidirlo se escaneó la χ²(masa) perfilada del único target de
Ceres (18937) variando la **media-ventana** de ajuste:

| ventana | n_obs | **M_min/lit** | χ²ᵣ@min | χ²ᵣ@lit | σ | mín. locales |
|---------|------:|--------------:|--------:|--------:|----:|----:|
| ±30 d | <8 | — (sin bundle) | — | — | — | — |
| ±45 d | 14 | **4.37** | 0.80 | 20.43 | 20.3 | 7 |
| ±60 d | 14 | **4.37** | 0.80 | 20.43 | 20.3 | 7 |
| ±90 d | 19 | **0.08** | 0.60 | 3.52 | 9.5 | 7 |
| ±120 d | 58 | **0.69** | 0.46 | 0.57 | 3.4 | 6 |
| ±180 d | 71 | **4.37** | 0.46 | 0.63 | 4.8 | 5 |
| one-sided | 234 | **0.04** | 1.35 | 2.35 | 21.5 | 6 |

El mínimo de masa de Ceres **barre ~100× (0.04–4.37× literatura) sin converger a
1.0** al crecer el arco, salta de forma **no monótona** (4.37→0.08→0.69→4.37→0.04)
y la curva tiene **5–7 mínimos locales en toda ventana**. El σ alto (3–21) es
**leverage ilusorio**: el basin es profundo pero dentado, así que el mínimo global
salta entre basins al cambiar la ventana. Mismo mecanismo que Hygiea (degeneración
M↔deltas activa sobre datos reales), solo que con basins más profundos.

**Conclusión del experimento**: el leverage de Ceres NO se traduce en una masa
identificable. No hay ventana con un plateau estable cerca de literatura. La masa
de Ceres **tampoco es determinable** con DR3 en este pipeline. (Datos:
`data/output/stage2_6/ceres_windows/`.)

## Veredicto del gate y decisión

**Gate FAIL**: ningún calibrador ≠ Ceres alcanza ≳3σ con la ventana acotada
(Pallas ≤1.2σ, Hygiea ≤2.1σ, Vesta sin datos), **y Ceres —el único con leverage
fuerte— da una masa no identificable** (window-inestable + multimodal). Por lo
tanto:

- **No se arranca A3.** Su premisa —que reducir los dof del drift destrabaría la
  masa— no se sostiene: donde no hay leverage (Pallas, Hygiea) no hay señal que
  A3 pueda limpiar; y donde sí lo hay (Ceres) el problema es la **no-
  identificabilidad** (χ²(masa) multimodal por la degeneración M↔deltas), que A3
  —que sigue dejando libres parámetros de drift por target— no resuelve.
- **Track A se cierra** con un veredicto **universal y limpio**: *ninguna masa de
  perturbador es determinable con Gaia DR3 en este pipeline, incluido Ceres*.
  Converge con el deepwork ("no publicable") con el mecanismo completamente
  caracterizado: **leverage insuficiente + no-identificabilidad M↔deltas sobre
  datos reales**. Esto es en sí un resultado metodológico citable (qué NO puede
  hacer la astrometría SSO de DR3 por deflexión), consistente con que la
  literatura de masas Gaia (p. ej. Fuentes-Muñoz 2024) use el FPR, no DR3.

## Qué reabriría Track A

- **Gaia DR4 / FPR**: arcos astrométricos más largos y mejor calibración AC
  elevan el leverage; los calibradores hoy planos (Pallas, Hygiea) podrían
  cruzar el umbral.
- **Si DR4 da leverage** para >1 calibrador, A3 (OU) y el fit multi-target joint
  (ya implementados, A2/A2.5) son el camino directo; el tooling queda listo.

## Linaje de Track A

| Etapa | Veredicto | Doc |
|-------|-----------|-----|
| A1 tighten priors | FAIL — bias estructural, no overfitting | [mass_layer_stage_a1_tight_priors.md](mass_layer_stage_a1_tight_priors.md) |
| A2 multi-target joint | FAIL — masa multi==single; refuta degeneración (sobre sintético) | [mass_layer_stage_a2_multitarget.md](mass_layer_stage_a2_multitarget.md) |
| closing-loop | bug del optimizador (no bias físico); mínimo real para Ceres/Hygiea | [mass_layer_closing_loop_leverage.md](mass_layer_closing_loop_leverage.md) |
| A2.5 profiled | bug ARREGLADO (ratio 1.000 sintético); aflora sesgo real-data | [mass_layer_stage_a2_5_profiled.md](mass_layer_stage_a2_5_profiled.md) |
| A2.6 real-data | "9×" = basin espurio + drift; M↔deltas activa en real; no-identificable | [mass_layer_stage_a2_6_realdata_bias.md](mass_layer_stage_a2_6_realdata_bias.md) |
| **A3 gate-check** | **FAIL — ningún calibrador ≠Ceres con ≳3σ** | este documento |
| **Ceres window-scan** | **Ceres no-identificable: mínimo barre 0.04–4.37× según ventana, multimodal. Cierre universal** | este documento |
