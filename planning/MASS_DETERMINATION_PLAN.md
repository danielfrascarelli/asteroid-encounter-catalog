# Plan: Determinación de masas de asteroides estilo Fuentes-Muñoz

> **Estado:** 🟡 ACTIVO
> **Última actualización:** 2026-06-01
> Plan para convertir el catálogo de encuentros en **determinaciones de masa
> publicables**, replicando la metodología de **solución global simultánea**
> (órbitas + masa por mínimos cuadrados sobre el arco completo, con la
> covarianza along-scan de Gaia). Criterio de éxito: reproducir las masas de los
> 4 calibradores (Ceres/Vesta/Pallas/Hygiea) dentro de σ y producir ≥1 masa
> nueva defendible.
>
> Sustituye a `ROADMAP_TO_MASSES.md` (borrado): su enfoque por-encuentro fue
> refutado por el cierre Track A.

## Por qué hace falta esto (no es más datos, es arquitectura)

La capa de masas está cerrada para DR3 ([docs/mass_layer_track_a_closure.md](../docs/mass_layer_track_a_closure.md))
y FPR-solo tampoco la reabre con el método actual
([docs/mass_layer_fpr_revalidation.md](../docs/mass_layer_fpr_revalidation.md)).
La raíz es la **degeneración masa↔drift orbital**: nuestro `fit_mass_gaia_loo`
ajusta órbita y masa en **pasos secuenciales** (fija la órbita con datos
pre-encuentro, luego lee la masa del post), y al separarlos tira la información
que distingue una deflexión real de un error de órbita que crece con el tiempo →
χ²(masa) multimodal, mínimo errante, no identificable.

Fuentes-Muñoz (y todo el campo: OrbFit, JPL) resuelve **órbitas + masa juntas**,
en un único ajuste de mínimos cuadrados sobre el arco completo, con el Jacobiano
de las **ecuaciones variacionales** (∂obs/∂elementos y ∂obs/∂GM) y stackeando
**muchos** test-asteroids por perturber. La degeneración se maneja dentro de la
covarianza, no se descarta. Eso es lo que hay que construir.

---

## Tabla de estado

| # | tarea | fase | estado | entregable / gate |
|---|-------|------|--------|-------------------|
| T1 | Esqueleto `orbdet` + primitivas matemáticas (constants, kepler, frames, time) | 0 | ✅ | `src/orbdet/` · #70 (65 tests) |
| T2 | Modelo dinámico N-cuerpos (rebound: Sol+planetas+asteroides grandes+perturber) | 0 | ✅* | `src/orbdet/dynamics.py` · validado vs límite dos-cuerpos a 1e-8 AU; *cross-check Horizons marcado, pendiente de entorno con acceso JPL |
| T3 | Ecuaciones variacionales (∂estado/∂elementos y ∂estado/∂GM) | 0 | ✅ | `src/orbdet/variational.py`: STM analítica (rebound add_variation) + ∂x/∂elementos (Φ·J_elem) coincide con FD a <1e-6; ∂x/∂GM por DF central con meseta de Richardson |
| T4 | Modelo de observación + covarianza along-scan anisotrópica | 0 | ✅ | `src/orbdet/observation.py`: estado→ICRS→RA/Dec + light-time iterativa + covarianza AL/AC anisotrópica; gate verde (chain N-cuerpos vs oráculo kepleriano <0.1 mas; ruido AL → χ²/obs≈1) |
| T5 | Corrector diferencial (OD por mínimos cuadrados, arco completo) | 0 | ✅ | `src/orbdet/least_squares.py` (LM genérico) + `orbit_determination.py`; gate verde: recupera órbita sintética sin ruido (χ²<1e-6) y con ruido AL (χ²_red≈1, <5σ) |
| — | *(alternativa)* evaluar integrar OrbFit de terceros | 0 | ⬜ | spike de decisión (ver abajo) |
| T6 | Ajuste conjunto órbita+masa de un perturber | 1 | ✅ | `src/orbdet/mass_determination.py`; closing-loop verde: masa sintética inyectada recuperada ratio≈1.0 (sin ruido <2e-3; con ruido AL dentro de 3σ, σ informativa) |
| T7 | Stacking multi-asteroide (GM compartido, N targets) | 1 | ✅ | `mass_determination.determine_shared_mass` (sistema en flecha 1+6N); gate verde: σ(GM)∝1/√N (s2/s1≈1/√2, s4/s1≈0.5 a <5%) |
| T8 | Modelo de fuerzas + pesos completo (efemérides, debiasing, outliers) | 1 | ⬜ | χ²_red ≈ 1 en datos reales |
| T9 | Adaptador FPR → motor (obs + covarianza por tránsito) | 2 | ⬜ | corre Big-4 end-to-end sobre FPR |
| T10 | Validación contra literatura (4 calibradores + Fuentes-Muñoz + Goffin/Galád) | 2 | ⬜ | \|z\| < 3 en los 4 calibradores |
| T11 | Corrida de producción + catálogo de masas + writeup | 2 | ⬜ | ≥1 masa nueva defendible |

---

## Fase 0 — Motor de determinación de órbitas propio ("OrbFit de cero")

> **Esta fase ES "implementar OrbFit".** Una librería de astrodinámica
> autocontenida en `src/orbdet/`, testeable de forma aislada, que NO importa
> otros módulos del pipeline (contrato de aislamiento verificado por test).
> Depende solo de numpy/scipy/astropy/rebound. Es además el hogar canónico de la
> matemática del proyecto (consolida, con el tiempo, `time_utils`,
> `astrometry/transforms`, `propagate/kepler`+`nbody`, `mass/forward_model_*`).

### T1 — Esqueleto + primitivas matemáticas
**Entregable:** `src/orbdet/{constants,kepler,frames,time_scales}.py` + `tests/orbdet/`.
**Gate:** round-trip elementos↔estado; invariantes (energía, momento angular,
período) dentro de tolerancia; conversiones de tiempo coinciden con astropy; un
test de aislamiento falla si `orbdet` importa cualquier `src.*` no-orbdet.

### T2 — Modelo dinámico N-cuerpos
**Entregable:** `src/orbdet/dynamics.py` — ecuaciones de movimiento vía rebound,
conjunto de perturbers configurable (Sol + 8 planetas + asteroides grandes + el
perturber bajo estudio).
**Gate:** propagar un asteroide numerado y reproducir JPL Horizons a ≲ pocas mas
sobre la ventana Gaia. **Depende de:** T1.

### T3 — Ecuaciones variacionales (el corazón del método)
**Entregable:** `src/orbdet/variational.py` — parciales ∂x(t)/∂x₀ (matriz de
transición de estado) y ∂x(t)/∂GM_perturber, junto con la dinámica de T2.

**Decisión de diseño (2026-06-01): enfoque (a) — diferencias finitas para
∂x/∂GM.**
- **∂x/∂elementos (∂x/∂x₀):** analítico vía `sim.add_variation()` de rebound
  (partículas variacionales de 1er orden respecto a las condiciones iniciales).
- **∂x/∂GM_perturber:** **diferencias finitas** — propagar a GM±δ y diferencia
  central, con δ relativo elegido (~1e-3·GM) y verificación de convergencia
  (Richardson / barrido de δ).
- **Por qué (a) y no (b) (partícula variacional analítica respecto a la masa):**
  desbloquea T4/T5/T6 end-to-end cuanto antes con código simple y robusto; el
  costo extra de 2 propagaciones por evaluación de la parcial de masa es
  aceptable en esta fase. La variante analítica (b) queda como **optimización
  posterior** (más exacta y rápida en producción, bastante más código y sutil),
  a reconsiderar si el costo o el ruido numérico de las DF se vuelven limitantes
  en T7 (stacking multi-asteroide).

**Gate:** (i) ∂x/∂elementos analítico coincide con su diferencia finita a < 1e-6
relativo; (ii) ∂x/∂GM por DF es estable bajo refinamiento de δ (meseta de
Richardson). **Depende de:** T2.

### T4 — Modelo de observación + covarianza anisotrópica
**Entregable:** `src/orbdet/observation.py` — estado heliocéntrico → baricéntrico
ICRS → RA/Dec desde la posición de Gaia, corrección de light-time iterativa, y la
**covarianza along-scan/across-scan anisotrópica** de Gaia (AL ~0.2–2 mas, AC
~cientos). Reutiliza la lógica de proyección AL que ya funciona
(`al_residuals_and_weights`).
**Gate:** los residuos de una órbita conocida quedan al nivel del ruido AL.
**Depende de:** T1. (Esta es la mitad Gaia-específica que ya tenemos resuelta en
parte — el lever que OrbFit de stock probablemente no modela.)

### T5 — Corrector diferencial (OD por mínimos cuadrados)
**Entregable:** `src/orbdet/least_squares.py` + `orbit_determination.py` —
Gauss-Newton/Levenberg-Marquardt para los 6 elementos sobre el **arco completo**
(sin split LOO), con Jacobiano analítico de T3 y pesos de T4.
**Gate:** recupera una órbita conocida a partir de observaciones sintéticas con
ruido; converge robustamente. **Depende de:** T3, T4.

### Alternativa: integrar OrbFit de terceros (spike de decisión)
Antes de comprometerse a T1–T5, un spike barato (~1–2 días):
1. **Showstopper:** ¿OrbFit soporta covarianza astrométrica **anisotrópica**
   (pesos AL/AC por obs)? Sin eso, Gaia pierde el lever → opción muerta.
2. Si pasa: reproducir una masa de literatura con datos **terrestres** que OrbFit
   come nativo (Goffin 2014 / Galád 2002), aislando "¿funciona la determinación
   de masa?" de "¿podemos adaptar Gaia?".
**Resultado:** documentar build-from-scratch (Fase 0) vs integrar OrbFit. La
intuición actual (ver historial) es que construir es más limpio porque ya
tenemos la mitad Gaia-específica y OrbFit nos haría rehacerla contra Fortran.

---

## Fase 1 — Determinación de masa

### T6 — Ajuste conjunto órbita + masa (un perturber)
**Entregable:** GM del perturber como parámetro libre **en el mismo sistema de
mínimos cuadrados** que los 6 elementos del target; covarianza masa↔elementos
propagada. `src/orbdet/mass_determination.py`.
**Gate:** closing-loop — inyectar una masa en datos sintéticos y recuperarla con
ratio ≈ 1.0 y σ realista (lo que el LOO secuencial nunca logró sobre datos
reales). **Depende de:** T5.

### T7 — Stacking multi-asteroide
**Entregable:** GM compartido entre muchos test-asteroids, parámetros orbitales
por target, resueltos en un solo sistema de ecuaciones normales.
**Gate:** σ(GM) baja como ~1/√N al agregar targets; la multimodalidad de
1-encuentro (que vimos en FPR) se rompe. **Depende de:** T6.

### T8 — Modelo de fuerzas y pesos completo
**Entregable:** efemérides planetarias JPL de alta precisión, perturbers
asteroidales grandes en el modelo, rechazo de outliers, modelo de error/debiasing
de Gaia.
**Gate:** χ²_red ≈ 1 en ajustes reales (no inflado ni deflactado). **Paralelo a**
T6/T7.

---

## Fase 2 — Datos y validación

### T9 — Adaptador FPR → motor
**Entregable:** alimentar el motor con obs FPR + covarianza por tránsito
(ingesta ya existe, flag `release=fpr`); filtrar `is_rejected`.
**Gate:** corre un Big-4 end-to-end sobre FPR. **Depende de:** T5.

### T10 — Validación contra literatura
**Entregable:** reproducir masas de Ceres/Vesta/Pallas/Hygiea; cruzar
Fuentes-Muñoz 2024/25 (231 masas) y Goffin/Galád.
**Gate:** |z| < 3 para los 4 calibradores; acuerdo con un número significativo de
masas de Fuentes-Muñoz. **Depende de:** T7, T9.

### T11 — Producción + catálogo + writeup
**Entregable:** corrida sobre los perturbers viables, catálogo de masas nuevas
con incertidumbres, y writeup publicable.
**Gate:** ≥1 masa nueva defendible (pasa specificity test, consistente con
literatura donde solapa). **Depende de:** T10.

---

## Fuera de scope
- Detección/caracterización de encuentros (ya hecho y congelado, ver `FROZEN_RUN.md`).
- DR4 (cuando salga, el motor lo soportará vía el flag `release`).
- Reescribir el pipeline de detección o el dashboard.

## Riesgos principales
- **Parcial ∂obs/∂GM** en rebound (T3): puede requerir partícula variacional
  custom; diferencias finitas como fallback de arranque.
- **Efemérides planetarias** de precisión suficiente para no contaminar el ajuste.
- **Cómputo** del solve multi-asteroide (T7) con muchos targets × arco largo.
- **Decisión OrbFit** (spike): si el showstopper anisotrópico no se resuelve,
  build-from-scratch es el único camino.
