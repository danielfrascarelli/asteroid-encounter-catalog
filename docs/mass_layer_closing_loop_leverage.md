# Closing-the-loop test — why the mass fits collapse to the H-prior

> Estado: **hallazgo mayor**. Fecha: 2026-05-29.
> Revierte parcialmente la conclusión del deepwork ("capa de masas no publicable")
> y del veredicto preliminar de A2 ("cerrar Track A, esperar DR4").
> Herramientas: [scripts/mass/closing_loop_test.py](../scripts/mass/closing_loop_test.py),
> [scripts/dev/probe_mass_sensitivity.py](../scripts/dev/probe_mass_sensitivity.py).

## Punto de partida

A1 (tight priors) y A2 (multi-target) fallaron el gate, y A2 mostró que los
Δ-elementos ajustados son ~10⁻⁷ — no absorben señal. El ratio fit/lit de cada
calibrador resultó ser **exactamente** M_H/M_lit (la masa fotométrica derivada
de la magnitud absoluta H del perturbador), sin importar priors ni estimador.
Eso disparó la pregunta: ¿el forward model no tiene sensibilidad a la masa, o
hay un bug que impide al fit moverse del prior?

## Test 1 — closing-the-loop (inyección de masa conocida)

`closing_loop_test.py` reconstruye bundles reales (geometría Gaia + elementos
LOO), regenera RA/Dec sintéticos con `forward_model` a la **masa literatura**
(deltas = 0), y re-ajusta. Sin ruido, el fit debería recuperar la masa exacta.

| Perturber | inyectado | fit | ratio fit/iny | log10_fit vs H-init |
|-----------|-----------|-----|---------------|---------------------|
| Pallas (5 tgt) | 2.05×10²⁰ | 1.17×10²⁰ | 0.571 | ≡ H-init (Δ<10⁻¹³) |
| Hygiea (5 tgt) | 8.30×10¹⁹ | 1.94×10¹⁹ | 0.234 | ≡ H-init (Δ<10⁻¹³) |

El fit devuelve la masa inicial basada en H **exacta a 13 decimales**, con
nfev≈12–15. No mide masa: **echa el prior fotométrico**.

## Test 2 — sensibilidad astrométrica a la masa (`probe_mass_sensitivity.py`)

¿Es un bug (la masa no propaga) o falta de señal? El probe evalúa `forward_model`
sobre una grilla de masas y mide el desplazamiento along-scan (AL) inducido,
whitened con la covarianza real, **con los deltas congelados** (señal intrínseca,
sin degeneración). `chi_AL` = √Δχ² vs masa cero; una detección 3σ necesita ≳3.

| Perturber→target | N_obs | AL σ (mas) | chi_AL @ M_lit | régimen |
|------------------|-------|-----------|----------------|---------|
| Ceres→18937  | 234 | 3.8 | **25.2** | señal fuerte |
| Hygiea→58775 | 36  | 3.9 | **2.1**  | marginal |
| Pallas→28036 | 62  | 2.0 | **0.09** | sin señal |

La masa **sí** propaga (la deflexión escala ~linealmente con M). Pero el leverage
es radicalmente dependiente del encuentro: Ceres tiene 25σ de señal intrínseca,
Pallas prácticamente nada.

## Test 3 — escaneo de χ²(log10_M) con deltas congelados

`closing_loop_test.py --scan-mass` evalúa χ²(masa, δ=0) sobre los datos
sintéticos inyectados a la masa real:

- **Hygiea**: mínimo nítido y profundo **en la masa real** (χ²=0); en el prior-H
  χ²≈29. **Δχ² ≈ 29 (~5σ) a favor de la masa real.**
- **Pallas**: mínimo también en la masa real, pero la curva va de 0.036 a 0.000;
  **Δχ²(H→verdad) ≈ 0.008 (~0.09σ)** — plano.

## Diagnóstico — dos modos de falla distintos

| | Señal disponible | Causa del ratio = M_H/M_lit |
|--|------------------|------------------------------|
| **Pallas** | Δχ² ~ 0.008 (0.09σ) | **leverage genuinamente ausente** en DR3 (geometría débil / arco corto). Irrecuperable sin DR4. |
| **Hygiea, Ceres** | Δχ² ~ 29 / ~600 | **el optimizador no desciende**: la señal es real y profunda, pero el fit no se mueve del prior-H. |

Para Hygiea/Ceres **no es falta de leverage ni la degeneración M↔deltas** — es
una **patología numérica de `least_squares`**.

## Causa raíz del bug del optimizador

`fit_joint_multitarget` llama `least_squares(method="trf", ...)` sin `x_scale`
ni `diff_step`. Sobre los datos sintéticos de Hygiea:

- El gradiente de masa **existe y es grande**: con deltas congelados,
  dχ²/dlog10M ≈ **−42.8** al paso de 4×10⁻³ dex (sonda de 2 puntos). No es un
  Jacobiano cero.
- Aun así el fit termina por **`xtol`** con first-order optimality ~10¹¹, dejando
  la masa **exactamente** en x0.
- `jtj_condition` ~10¹³: log10_M (~20) y los deltas (~10⁻⁴) difieren 5 órdenes
  de magnitud. El trust-region, en coordenadas crudas, está dominado por la
  masa; un radio capaz de mover la masa lanzaría los deltas fuera de rango
  físico, así que el solver lo achica hasta que el paso de masa es ~10⁻⁴ dex →
  la masa se congela en el prior.

Se intentaron, sin éxito, las combinaciones de `diff_step` (2×10⁻⁴ para la masa;
0.1σ por delta) y `x_scale` (`'jac'`; array fijo masa=1/delta=σ): cada variante
congela una combinación distinta de parámetros, pero **ninguna mueve la masa**.
El problema es estructural del solver acoplado masa+deltas mal condicionado.

## Recomendación — optimizador perfilado (profiled likelihood)

El escaneo χ²(masa) es, de hecho, la verosimilitud perfilada con deltas=0, y su
mínimo cae en la masa real (en datos sin ruido). La solución robusta:

1. **Outer**: minimización 1-D de `log10_M` (Brent acotado) sobre el χ² perfilado.
2. **Inner**: para cada masa fija, ajustar solo los 6N deltas (problema homogéneo,
   bien condicionado; reparametrizar `u = δ/σ` para O(1)).
3. La curvatura del χ²(masa) en el mínimo da σ_M directamente.

Esto evita por completo el acoplamiento mal condicionado masa↔deltas. Sobre los
calibradores se espera: **Hygiea y Ceres recuperan la masa literatura** (señal
real ~5σ/25σ), **Pallas queda con error enorme** (sin leverage). Si se confirma,
**revierte la conclusión del deepwork** para perturbadores con buena geometría:
la capa de masas es viable en DR3 para Ceres/Hygiea, no para Pallas.

### Implicación

El "bias estructural que escala con 1/M_real" de Stage 4 no era físico: era el
fit devolviendo M_H (que casualmente da ratios 0.77/0.57/0.23 para
Ceres/Pallas/Hygiea porque M_H/M_lit escala así). El gate FAIL de A1/A2 y el
veredicto del deepwork se explican por este bug + la debilidad real de Pallas.
