# Track A Stage 2.5 — Optimizador perfilado: resultado

> Estado: **parcial — bug del optimizador arreglado, sesgo real-data persiste**.
> Fecha: 2026-05-29.
> Implementación: `fit_joint_multitarget_profiled` en
> [src/mass/forward_model_joint_multitarget.py](../src/mass/forward_model_joint_multitarget.py);
> flag `--optimizer {joint,profiled}` en
> [scripts/mass/fit_mass_gaia_multitarget.py](../scripts/mass/fit_mass_gaia_multitarget.py)
> y [scripts/mass/closing_loop_test.py](../scripts/mass/closing_loop_test.py).
> Datos: [data/output/stage_a2_5_profiled_validation.csv](../data/output/stage_a2_5_profiled_validation.csv).
> Antecedente: [mass_layer_closing_loop_leverage.md](mass_layer_closing_loop_leverage.md).

## Qué se hizo

El closing-loop test mostró que el fit conjunto `least_squares` no descendía en
la dirección de masa (trust-region mal condicionado, cond~10¹³) y devolvía el
prior fotométrico M_H. Stage 2.5 reemplaza ese fit por un **optimizador
perfilado**:

- **Outer**: minimización 1-D acotada (Brent) de la χ² perfilada sobre `log10_M`.
- **Inner**: para cada masa fija, `least_squares` sobre solo los 6N deltas
  (problema homogéneo, bien condicionado; `diff_step=0.1σ` por delta para no
  caer bajo el piso numérico del integrador).
- σ_M de la curvatura de χ²(log10_M) en el mínimo.

## Resultado

| Caso | Datos | optimizer | ref (kg) | fit (kg) | ratio | χ²_red |
|------|-------|-----------|----------|----------|-------|--------|
| Hygiea closing-loop | sintético sin ruido | **profiled** | 8.3×10¹⁹ (iny) | 8.30×10¹⁹ | **1.000** | 0.0 |
| Hygiea closing-loop | sintético sin ruido | joint (viejo) | 8.3×10¹⁹ (iny) | 1.94×10¹⁹ | 0.234 | 0.04 |
| Pallas closing-loop | sintético sin ruido | profiled | 2.05×10²⁰ (iny) | 2.1×10¹⁵ | ~0 | 0.0 |
| **Hygiea real** | Gaia DR3 | profiled | 8.3×10¹⁹ (lit) | **7.51×10²⁰** | **9.04** | 1.19 |
| **Pallas real** | Gaia DR3 | profiled | 2.05×10²⁰ (lit) | 2.66×10¹⁶ | ~0 | 0.60 |

## Lectura — tres problemas que estaban apilados

1. **Bug del optimizador (ARREGLADO).** Sobre datos sintéticos sin ruido
   inyectados a la masa real, el perfilado recupera la masa **exacta**
   (Hygiea ratio = 1.000), mientras el fit conjunto viejo se quedaba en
   M_H (ratio 0.234). Esto **confirma** que el "ratio = M_H/M_lit" de
   Stage 4 / A1 / A2 era el optimizador sin moverse, **no un bias físico**.

2. **No-identificabilidad de Pallas (REAL, no arreglable en DR3).** Con 0.09σ
   de señal intrínseca, los deltas absorben cualquier masa: el perfilado
   halla χ²≈0 en una masa absurda (10¹⁵). La σ por curvatura es engañosamente
   chica en este régimen plano — **no usar como detección**.

3. **Sesgo real-data (NUEVO, ahora expuesto).** Con el optimizador arreglado,
   Hygiea real ya no se queda en M_H: se va a **9× la literatura** con
   χ²_red=1.19 (buen ajuste). El sesgo cambió de dirección (0.23 → 9), lo que
   prueba que ahora la masa **sí** se ajusta a los datos — pero los datos
   reales contienen una sistemática (deriva orbital no capturada por los
   elementos LOO, sistemática astrométrica Gaia, o perturbador secundario)
   que la masa absorbe incorrectamente. No es leverage (Hygiea tiene ~5σ) ni
   el optimizador (arreglado): es modelado físico / datos.

## Veredicto

Stage 2.5 **resuelve el bug del optimizador** (entregable mecánico verificado
por closing-loop) pero **no valida la capa de masas** contra literatura: el
sesgo real-data persiste, ahora aislado como el problema dominante. La
conclusión del deepwork ("capa no publicable") se mantiene, pero el diagnóstico
es radicalmente más preciso: de los tres problemas que estaban confundidos
(degeneración M↔deltas, falta de leverage, bias estructural 1/M), el real es
una **sistemática de datos/modelo en el régimen con señal**, separable ahora
que el optimizador funciona y Pallas está identificado como leverage-limited.

## Limitaciones / notas

- La variante `--noise realistic` del closing-loop **no es válida
  cuantitativamente**: inyecta ruido por cuadratura RA/Dec que incluye el error
  sistemático across-scan (~433 mas), inconsistente con el whitening along-scan
  (~2 mas) del likelihood (χ²_red sale ~4000). Para atribuir ruido vs
  sistemática hay que muestrear de la covarianza 2D real. Pendiente.
- σ_M por curvatura no es robusta en el régimen plano/degenerado (Pallas).
  Usar un cruce Δχ²=1 sobre el perfil sería más fiable.

## Próximos pasos sugeridos

- **Investigar el sesgo real-data de Hygiea** (9× alto): per-observation
  residuals post-encuentro, estabilidad target-a-target del mínimo perfilado,
  scan de χ²(masa) sobre datos reales (¿unimodal?), perturbador secundario.
  Se solapa con Track B1 (outliers Stage 2).
- Arreglar el muestreo de ruido del closing-loop (covarianza 2D) para separar
  ruido de sistemática.
- σ_M robusta vía perfil Δχ²=1.
