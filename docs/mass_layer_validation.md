# Stage 4 — Validación contra masas conocidas (DAWN / Goffin / Vernazza)

> Gate de publicación: re-fittear (1) Ceres, (4) Vesta, (2) Pallas,
> (10) Hygiea con el pipeline completo (joint + Mahalanobis 2D) y
> verificar que las masas literatura quedan dentro de **|z| < 3** del
> fit, donde `z = (M_fit − M_lit) / sqrt(σ_fit² + σ_lit²)`.

**Branch**: `track2/stage4-literature-validation`
**Fecha**: 2026-05-29
**Outputs**:
- `data/output/stage4_validation_summary.csv` (todos los fits)
- `data/output/stage4_validation/fit_*.json` (uno por par)

---

## TL;DR

**El gate de Stage 4 falla**: 0 / 11 fits exitosos caen dentro de
`|z| < 3` de la masa literatura. El pipeline **subestima
sistemáticamente** las masas conocidas en un factor que crece a
medida que la masa real disminuye:

| Calibrator | M_lit (kg) | M_fit típico (kg) | ratio | bias factor |
|---|---|---|---|---|
| (1) Ceres | 4.71×10²⁰ | 3.64×10²⁰ | 0.77 | 1.30× under |
| (4) Vesta | 2.59×10²⁰ | (no exitoso) | — | — |
| (2) Pallas | 2.05×10²⁰ | 1.17×10²⁰ | 0.57 | 1.75× under |
| (10) Hygiea | 8.30×10¹⁹ | 1.96×10¹⁹ | 0.24 | **4.2× under** |

Las incertidumbres reportadas (`σ_fit`) son **artificialmente
estrechas**: `jtj_condition ~ 10¹²` indica que la columna de masa del
Jacobiano está casi totalmente correlada con los 6 deltas orbitales,
lo que infla la confianza aparente. Los `z`-scores reportados son
todos `|z| > 6` y hasta `|z| ≈ 26` para Ceres — la asimetría no es
estadística, es **sistemática del modelo**.

**Conclusión científica honesta**: el pipeline en su estado actual
**no puede reportar masas absolutas con calidad publicable**. La
señal de masa que recupera es real pero biased — los 6 deltas
orbitales absorben parte de la deflección inducida por la masa,
especialmente para perturbers pequeños.

---

## Mecanismo del bias

El forward model joint tiene 7 parámetros: 1 masa + 6 Δ-elementos
orbitales del target. Esos 6 deltas absorben drift orbital (la razón
original del Stage 1), pero **también absorben parcialmente la
deflección gravitatoria** del encuentro porque ambas señales
afectan la trayectoria del target.

Esto se confirma en los condicionamientos del Jacobiano:

- Ceres / 18937: `jtj_condition` no reportado directamente pero
  `log10_mass_sigma = 0.0012` con n_joint = 234 → uncertainty
  artificialmente baja.
- Pallas / 28036: `jtj_condition = 8.2×10¹²` → cerca de singular.

Cuando dos columnas del Jacobiano están muy correladas, el optimizer
encuentra una dirección de mínimo donde `(M, deltas)` se compensan,
y mueve la masa solo lo necesario para minimizar χ² conjunto. Si la
señal de masa puede absorberse en deltas con menor "costo" en priors
que en la masa, el optimizer prefiere absorberla — y resulta una
masa underestimated.

**Magnitud del bias y dependencia con la masa**: el ratio M_fit/M_lit
es 0.77 (Ceres) → 0.57 (Pallas) → 0.24 (Hygiea). Es decir, **a menor
masa real, más fracción se absorbe en los deltas**. Esto tiene
sentido: el signal-to-(absorbed-by-deltas) ratio escala con la masa.

---

## Tabla detallada (11 fits exitosos)

```
perturber  target  n_joint  chi2_red    mass_fit            sigma_fit       z         ratio
1 Ceres    18937   234      1.13        3.64e20             1.04e18         -25.83    0.77
2 Pallas   28036    62      0.36        1.17e20             9.13e17         -17.30    0.57
2 Pallas   47563    86      0.70        1.17e20             1.70e18         -16.71    0.57
2 Pallas   59882    42      0.45        1.17e20             1.22e19         -6.66     0.57
2 Pallas   60093    42      0.45        1.17e20             4.21e18         -13.45    0.57
2 Pallas   73243    62      0.75        1.17e20             9.53e18         -8.17     0.57
10 Hygiea  4803     48      0.20        1.94e19             2.29e17         -15.87    0.23
10 Hygiea  16772    69      0.95        1.94e19             8.93e17         -15.51    0.23
10 Hygiea  45989   157    160.3*        2.04e19             6.35e17         -15.45    0.25
10 Hygiea  47605   111     10.9*        1.96e19             1.32e18         -15.04    0.24
10 Hygiea  58775    36      0.50        1.97e19             4.86e18         -10.06    0.24
```

*Los dos Hygiea con χ²_red alto siguen subestimando — el bias no se
explica por mal ajuste sino por absorción estructural.*

**Vesta**: ningún fit exitoso. De 8 targets seleccionados, 4 cayeron
por encuentros muy tempranos en Gaia DR3 (no hay LOO baseline de 180
días) y 4 con cero transitos Gaia para ese target. Se necesita una
selección más cuidadosa de targets Vesta o más datos (FPR/DR4).

---

## Consistencia inter-target

Para Pallas los 5 fits exitosos dan **exactamente la misma masa**
(1.17e20 kg, ±0.05% entre fits). Para Hygiea los 5 dan 1.94–2.04×10¹⁹
(±5%). Esto significa que **el fit tiene un attractor estable** —
la combinación H-magnitude initial + dinámica + drift converge al
mismo punto independientemente del target. No es ruido estadístico:
es un **bias estructural reproducible**.

Eso es bueno en el sentido de que el bias se puede **calibrar y
corregir**, pero implica también que las "detecciones" de masa con
chi²_red < 1 (Stage 2/3) son consistentes con esa misma absorción
estructural y no necesariamente con detecciones reales.

---

## Implicaciones para Stage 0/1/2/3 (auto-citas)

1. **(111) Ate, (206) Hersilia, (124) Alkeste / 3294** previamente
   "detectados": las masas reportadas (5.4e17, 2.6e17, 5.0e17 kg)
   son ahora **doblemente cuestionables**:
   - Specificity test (Stage 3): nulls reproducen el mismo orden de
     magnitud de masa.
   - Calibración Stage 4: el pipeline subestima por 1.3–4× — las
     masas verdaderas (si existen) podrían ser 2–4× mayores; pero
     dado que el bias depende de la masa real (que no conocemos),
     no se puede corregir por escala simple.

2. **(19) Fortuna / 13346 y (49) Pales / 94474** (los 2 que pasaron
   specificity en χ²): la magnitud de masa reportada está sujeta al
   mismo bias. Sí hay señal **estadística** de deflección, pero el
   valor absoluto está subestimado por factor desconocido.

---

## ¿Se puede arreglar?

Posibles direcciones (todas requieren trabajo adicional fuera del
alcance Stage 4):

1. **Tightening de priors orbitales**: los priors gaussianos sobre
   los 6 deltas son anchos (`σ_da_rel = 2e-4`, etc.) porque
   permitían absorber drift Gaia. Apretándolos al nivel de
   uncertainty real de MPCORB (~10× más estrechos), los deltas
   tendrían menos espacio para absorber masa. Riesgo: los fits ya
   bien ajustados de Stage 2 se desestabilizan.

2. **Forward model físico**: en lugar de deltas absolutos, modelar
   drift como un proceso de Ornstein-Uhlenbeck con varianza
   conocida. Reduce los grados de libertad efectivos.

3. **Multi-target fit conjunto**: re-fittear todos los targets de
   Pallas (o Ceres) simultáneamente, compartiendo masa y dejando
   solo los deltas por-target libres. La masa estaría restringida
   por la consistencia entre targets.

4. **Cambiar dataset**: Gaia FPR / DR4 ofrecerá ~10× más obs por
   target. Con más datos, el optimizer tendrá menos espacio para
   confundir mass con deltas — pero esto no llega antes de 2026.

5. **Reportar solo masas relativas**: si el pipeline tiene un bias
   reproducible (Pallas siempre 0.57×, Hygiea siempre 0.24×), se
   puede calibrar y reportar ratios. Útil para perturbers en la
   misma familia dinámica pero no resuelve la magnitud absoluta.

---

## Decisión sobre Track 2

**Stage 4 falla el gate de publicación.** La capa de masas del
pipeline en su estado actual:

- Detecta señal real de deflección en algunos casos (Fortuna, Pales).
- **No** puede reportar masas absolutas con confianza publicable.
- **No** reproduce masas conocidas dentro de 3σ — bias sistemático
  de factor 1.3× a 4×.

Próximos pasos honestos:

1. **Pausar Track 2** hasta nuevo trabajo en el modelo (priors
   tighter, multi-target joint, o DR4).
2. **Documentar limitaciones** en `FROZEN_RUN.md` / `ROADMAP.md`:
   el catálogo de encuentros (Track 1) es publicable; la
   determinación de masas no.
3. **Retractar** las afirmaciones cuantitativas Stage 0/1/2 sobre
   (111) Ate (5.43e17 kg) y otros — siguen siendo registros del
   desarrollo del pipeline pero no detecciones publicables.

---

## Entregables Stage 4

- [x] `scripts/mass/run_stage4_validation.py`
- [x] `data/output/stage4_validation_summary.csv` (31 filas, 11 ok)
- [x] `data/output/stage4_validation/` (JSON per-par)
- [x] `docs/mass_layer_validation.md` (este documento)
