# orbdet — motor de determinación de órbitas y masas

Implementación propia ("OrbFit de cero") de un motor de determinación de órbitas
por mínimos cuadrados, y sobre él la determinación de masas de asteroides
perturbadores por el método de encuentros cercanos (estrategia tipo
Fuentes-Muñoz). Es la **Fase 0** de [`planning/MASS_DETERMINATION_PLAN.md`](../../planning/MASS_DETERMINATION_PLAN.md).

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
- `numpy`, `scipy`, `astropy`, `rebound`.

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
| `dynamics.py` | ⬜ T2 | fuerzas N-cuerpos (rebound) |
| `variational.py` | ⬜ T3 | ∂estado/∂elementos y ∂estado/∂GM |
| `observation.py` | ⬜ T4 | estado→RA/Dec + covarianza along-scan |
| `least_squares.py` / `orbit_determination.py` | ⬜ T5 | corrector diferencial |
| `mass_determination.py` | ⬜ T6/T7 | ajuste conjunto órbitas+masa, stacking |
