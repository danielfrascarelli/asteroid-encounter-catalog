# Track A Stage 2 — Multi-target joint fit: diagnóstico

> Estado: **gate FAIL**. Fecha: 2026-05-29.
> Entregables: [src/mass/forward_model_joint_multitarget.py](../src/mass/forward_model_joint_multitarget.py),
> [scripts/mass/fit_mass_gaia_multitarget.py](../scripts/mass/fit_mass_gaia_multitarget.py),
> [tests/test_forward_model_joint_multitarget.py](../tests/test_forward_model_joint_multitarget.py),
> [data/output/stage_a2_multitarget_validation.csv](../data/output/stage_a2_multitarget_validation.csv).

## Objetivo

Romper la degeneración hipotética M ↔ Δ-elementos del joint fit single-target
(que sesgaba las masas en Stage 4, ver [mass_layer_validation.md](mass_layer_validation.md))
compartiendo una sola `log10_M_perturber` entre los N targets de un mismo
cuerpo grande, con 6 deltas orbitales libres por target. La consistencia
inter-target debería restringir la masa más que cualquier target aislado:
si los deltas son ortogonales a la masa, σ_M baja por ~√N y la masa fitted
debería converger al valor literatura.

## Setup

- Forward model: `residuals_joint_multitarget` — vector de parámetros
  `(log10_M, da_1..dM_1, da_2..dM_2, ...)` de dimensión `1 + 6N`. Misma masa
  compartida, deltas por-target independientes.
- Likelihood: Mahalanobis 2D (AL+AC), priors `default`.
- Calibradores con N ≥ 2 targets `ok` en
  [stage4_validation_summary.csv](../data/output/stage4_validation_summary.csv):
  - **Pallas** (perturber 2): 5 targets (28036, 47563, 59882, 60093, 73243).
  - **Hygiea** (perturber 10): 5 targets (4803, 16772, 45989, 47605, 58775).
  - **Ceres** (perturber 1): solo 1 target `ok` → multi-target no aplica.
  - **Vesta** (perturber 4): 0 fits exitosos → fuera.

Tests sintéticos: 6/6 pass, incluido el que recupera M compartida + da_rel
per-target dentro de tolerancia con señal de masa time-varying (Heaviside) y
offset orbital constante por target.

## Resultado

| Perturber | N | M_fit multi (kg) | σ_M (kg) | ratio multi/lit | ratio single/lit | z | χ²_red | gate |
|-----------|---|------------------|----------|-----------------|------------------|------|--------|------|
| Pallas    | 5 | 1.171×10²⁰       | 1.27×10¹⁸| **0.5711**      | 0.5711           | −17.0| 0.57   | FAIL |
| Hygiea    | 5 | 1.919×10¹⁹       | 1.06×10¹⁷| **0.2313**      | 0.2340           | −16.0| 74.0   | FAIL |

Criterio: ≥ 2/3 calibradores con |z| < 3. Resultado: **0/2**.

## Lectura

El hallazgo decisivo no es que falle el gate, sino **cuánto** falla:

1. **La masa multi-target es idéntica a la single-target.** Para Pallas el
   ratio coincide a 4 decimales (0.5711 = 0.5711); para Hygiea cae dentro del
   1 % (0.2313 vs 0.2340). Compartir M entre 5 targets **no movió la masa**.
   Esto **refuta la hipótesis de la degeneración M ↔ deltas**: si los deltas
   estuvieran absorbiendo señal de masa de forma idiosincrática por target,
   el fit conjunto la habría recuperado. No lo hace.

2. **Los deltas ajustados son diminutos** (da_rel ~ 10⁻⁷–10⁻⁶, es decir
   metros sobre una órbita de ~3 AU). No están absorbiendo deflección: ya
   eran ~0 en el single-target. El espacio de los deltas no es el problema.

3. **El bias es coherente y sistemático por target**, no ruido que promedia.
   Cada target del mismo perturber arrastra el mismo factor de subestimación.
   Combinado con el escalado con 1/M_real observado en Stage 4
   (0.77 Ceres / 0.57 Pallas / 0.23 Hygiea), apunta a un sesgo **en el
   modelo de deflección / los datos**, no en la parametrización orbital.

4. **σ_M colapsa irrealmente** (Pallas 1 %, Hygiea 0.5 %) y `jtj_condition`
   ~10¹³. El fit conjunto está sobre-confiado: la consistencia inter-target
   aprieta el error formal pero alrededor del valor **equivocado**. El χ²_red
   de Hygiea (74) está dominado por el target 45989 (χ²_red single = 160),
   que el fit conjunto no aísla.

## Veredicto

**Stage A2 falla el gate.** El multi-target no corrige el bias porque la
degeneración que pretendía romper no es el mecanismo causal. Pasar a
**Stage A3 (OU drift)** con una advertencia importante: la premisa de A3
—que reducir los grados de libertad de los deltas de 6 a 2 evita que
absorban señal— queda **debilitada** por este resultado, dado que los deltas
aquí ya son ~0 y no absorben nada. Antes de invertir 2–3 semanas en A3,
conviene reorientar el diagnóstico hacia el **modelo de deflección mismo**:

- ¿El escalado 1/M_real sugiere saturación del efecto de deflección modelado,
  o una geometría de encuentro mal capturada para targets pequeños?
- ¿Hay una normalización / signo / light-time en el cálculo de la deflección
  que subestima sistemáticamente la amplitud?
- Validar el forward model contra una deflección N-body sintética inyectada
  con masa conocida sobre los mismos transits (closing-the-loop test): si el
  forward model no recupera la masa inyectada, el bug está ahí y ni A2 ni A3
  lo arreglan.

Recomendación: **no arrancar A3 a ciegas**; primero el closing-the-loop test
del modelo de deflección sobre datos sintéticos (es ~1 día y discrimina entre
"bug del forward model" y "limitación del dataset DR3").
