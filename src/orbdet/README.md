# orbdet — motor de determinación de órbitas y masas

Implementación propia ("OrbFit de cero") de un motor de determinación de órbitas
por mínimos cuadrados, y sobre él la determinación de masas de asteroides
perturbadores por el método de encuentros cercanos (estrategia tipo
Fuentes-Muñoz). Implementa el plan completo [`planning/MASS_DETERMINATION_PLAN.md`](../../planning/MASS_DETERMINATION_PLAN.md) (T1–T11).

## Por qué existe (y por qué aislado)

La capa de masas vieja (`scripts/mass/fit_mass_gaia_loo.py`) ajusta órbita y masa
en **pasos secuenciales** y no resuelve la degeneración masa↔drift orbital → χ²
multimodal, no identificable (ver `docs/mass_layer_track_a_closure.md`). El
método correcto es un **ajuste simultáneo** de órbitas + masa sobre el arco
completo, con el Jacobiano de las **ecuaciones variacionales**. Eso requiere un
motor de astrodinámica limpio y testeable, no más parches sobre el pipeline.

## Contrato de aislamiento (verificado por test)

`orbdet` **no importa ningún otro módulo `src.*`** del proyecto. Solo depende de:

- librería estándar de Python,
- `numpy`, `scipy`, `astropy`, `rebound`, `assist` (efeméride DE440 + 16 perturbadores).

Internamente usa **imports relativos** (`from .constants import ...`).
`tests/orbdet/test_isolation.py` falla si algún archivo de `orbdet` viola esto.
Esto garantiza que el motor se pueda testear y razonar de forma independiente, y
que sea el **hogar canónico** de la matemática del proyecto (con el tiempo
absorbe `time_utils`, `astrometry/transforms`, `propagate/kepler`, etc.).

## Convenciones de unidades (internas)

| magnitud | unidad |
|----------|--------|
| distancia | AU |
| tiempo | días (JD TDB salvo indicación) |
| ángulo | radianes |
| velocidad | AU/día |
| GM | AU³/día² |

Conversiones a km/grados/ISO solo en los bordes (entrada/salida).

## Módulos

| módulo | estado | rol |
|--------|--------|-----|
| `constants.py` | ✅ T1 | constantes físicas (GM_sun, AU, c, oblicuidad) — única fuente de verdad |
| `time_scales.py` | ✅ T1 | conversiones JD entre TDB/TCB/UTC (vía astropy) |
| `frames.py` | ✅ T1 | rotaciones y eclíptica↔ecuatorial |
| `kepler.py` | ✅ T1 | dos cuerpos: elementos↔estado, ecuación de Kepler, propagación |
| `dynamics.py` | ✅ T2 | fuerzas N-cuerpos (rebound) |
| `variational.py` | ✅ T3 | ∂estado/∂elementos (Φ·J_elem analítico) y ∂estado/∂GM (FD + Richardson) |
| `observation.py` | ✅ T4 | estado→RA/Dec + light-time + covarianza along-scan anisotrópica |
| `least_squares.py` / `orbit_determination.py` | ✅ T5 | corrector diferencial (Levenberg-Marquardt) |
| `mass_determination.py` | ✅ T6/T7/T8 | ajuste conjunto órbitas+masa, stacking 1+6N, covarianza en bloques por FOV + calibración del piso |
| `dynamics_assist.py` | ✅ T8 | modelo de fuerzas ASSIST (DE440 + GR EIH + 16 perturbadores) |
| `gaia_adapter.py` | ✅ T9 | adaptador datos reales: σ_AL, MPCORB→elementos, grupos FOV |

**Estado:** T1–T11 completas (plan cerrado, PR #80). Validado sobre Gaia FPR: 4/4
calibradores |z|<3 + masa nueva (16) Psyche = 2.43×10¹⁹ kg ±3.3%. Ver
[`docs/orbdet_engine_status.md`](../../docs/orbdet_engine_status.md) y
[`docs/mass_determination_results.md`](../../docs/mass_determination_results.md).
