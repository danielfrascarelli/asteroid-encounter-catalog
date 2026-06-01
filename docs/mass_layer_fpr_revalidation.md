# Re-validación de la capa de masas con Gaia FPR (FPR_INGEST_PLAN Stage 4)

> El gate científico del plan. Pregunta: ¿el baseline extendido de FPR (~66 vs
> ~34 meses) reabre la determinación de masas que se cerró con DR3
> ([docs/mass_layer_track_a_closure.md](mass_layer_track_a_closure.md),
> `project_mass_optimizer_bug`)?
>
> **Respuesta corta: no con la metodología actual.** FPR aporta los datos
> post-encuentro que a DR3 le faltaban, pero el χ²(masa) sigue siendo **no
> unimodal** con el mínimo lejos de la literatura — la misma no-identificabilidad
> que cerró DR3. El entregable real son las Stages 0–3 (ingesta FPR reusable);
> Stage 4 es un resultado negativo honesto, no una reapertura.

Corridas: 2026-06-01, branch `feat/fpr-stage4-revalidation`, sobre el archivo
Gaia en vivo. Artefactos en `data/output/{dr3,fpr}_scan/` y
`data/output/fits_{dr3,fpr}/` (gitignored).

---

## 1. El dato sí mejora: FPR llena el arco post-encuentro

`fit_mass_gaia_loo` (LOO orbit + ventana de masa) sobre encuentros de Big-4
cerca del **final** de la ventana DR3 (2017-05-28):

| perturber→target | fecha | release | obs ventana (pre+post) | masa fiteada |
|------------------|-------|---------|------------------------|--------------|
| (2) Pallas → 28036 | 2017-03-01 | DR3 | 62 + **0** | — (sin post, no fittea) |
| (2) Pallas → 28036 | 2017-03-01 | FPR | 61 + **148** | 1.07e20 ± 0.02 dex |
| (2) Pallas → 73243 | 2017-05-08 | DR3 | 62 + **0** | — (sin post, no fittea) |
| (2) Pallas → 73243 | 2017-05-08 | FPR | 62 + **312** | 2.7e15 ± **9.3 dex** |

La premisa del plan se confirma **a nivel de datos**: encuentros que DR3 no
podía fittear (0 observaciones post-encuentro, el encuentro cae al borde del fin
de misión) ahora tienen arcos post de cientos de tránsitos en FPR.

## 2. Pero la masa no se vuelve defendible

Dos modos de falla, ambos heredados del cierre de DR3:

**(a) Restringido a la ventana de deflexión (±90 d) → no identificable.**
Window-scan del χ²(masa) profiled (mismo diagnóstico que cerró DR3):

| caso | release | M_min/lit | χ²ᵣ@min | χ²ᵣ@lit | Δχ²(lit−min) | σ | unimodal | mín. locales |
|------|---------|-----------|---------|---------|--------------|---|----------|--------------|
| Ceres → 18937 | DR3 | 0.42× | 0.61 | 0.94 | 10.3 | 3.2 | **No** | 4 |
| Pallas → 28036 | FPR | **8.66×** | 0.62 | 0.64 | 1.1 | 1.0 | **No** | — |

En FPR (Pallas 28036) el mínimo barre hasta **8.7× la masa de literatura**, la
curva es **multimodal**, y la masa de literatura queda a sólo **1σ** del mínimo:
estadísticamente indistinguible → **no identificable**. Es la misma firma
(χ²(masa) dentado/multimodal, mínimo errante) que cerró DR3.

**(b) Abierto al arco post completo → drift orbital domina.**
Pallas 73243 (FPR, arco post de 2.7 años): la órbita LOO (anclada en datos
pre-encuentro) predice el arco post a sólo ~418 mas RMS — el **drift orbital
acumulado** sobre el arco largo (gotcha de Kepler/arco-corto ya documentado en
CLAUDE.md) swampa la señal de deflexión (~40 mas) y deja la masa sin constreñir
(±9.3 dex). Cuando el arco se usa entero, las deltas orbitales libres absorben
el drift a una masa arbitraria — exactamente la no-identificabilidad de
[docs/mass_layer_stage_a2_6_realdata_bias.md](mass_layer_stage_a2_6_realdata_bias.md).

Y en los casos donde sí da un número con σ chica (Pallas 28036 LOO: 1.07e20),
el ratio fit/lit ≈ **0.52×** — **el mismo sesgo** que DR3 reportó para ese
perturber (0.57×, ver `stage4_validation_summary.csv`). El baseline extra no
remueve el sistemático.

## 3. Conclusión y camino

FPR **no reabre** masas defendibles con la maquinaria LOO/profiled actual. El
factor limitante no es cantidad de datos sino la **no-identificabilidad
metodológica**: la degeneración masa ↔ drift orbital sobre el arco. Reabrir
masas requeriría una **solución global simultánea** (órbitas + masas ajustadas
conjuntamente sobre todos los encuentros, como Fuentes-Muñoz 2024/2025), no el
ajuste por-encuentro de este pipeline — fuera del scope de la ingesta FPR.

**Lo que sí queda entregado y es reusable:**
- Ingesta FPR completa e intercambiable con DR3 vía flag `release` (Stages 0–3).
- Toda la capa de masas (`fit_mass_gaia_loo`, `realdata_mass_scan`,
  `_build_bundle`) acepta `--release` y corre sobre FPR sin más cambios.
- Evidencia cuantitativa de que FPR-solo (con esta metodología) no alcanza —
  sin sobre-afirmar ni retractar (cf. (111) Ate, `project_ate_mass_result`).

**Caveats de esta evaluación:**
- Muestra chica (4 encuentros LOO + 2 window-scans), limitada por compute. La
  señal es consistente y reproduce el cierre DR3, pero no es un barrido
  exhaustivo de los 231 pares de Fuentes-Muñoz.
- No se intentó la solución global simultánea (el método que sí funciona en la
  literatura); la conclusión es sobre *esta* metodología, no sobre FPR en
  abstracto.
- Targets con baseline pre-encuentro corto (p.ej. 18937, encuentro 2015-06 con
  misión arrancando 2014-08) son marginales para LOO en **ambos** releases.
