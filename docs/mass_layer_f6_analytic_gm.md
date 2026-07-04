# F6 — Parcial ∂x/∂GM analítica (partícula variacional de masa)

> **Estado:** ✅ HECHO para el backend `rebound` (2026-07-04).
> 🔒 BLOQUEADO por la API de la efeméride para el backend `assist` (producción),
> donde la parcial sigue siendo por diferencias finitas.
> Toca `src/orbdet/variational.py` y `src/orbdet/mass_determination.py`.

## Objetivo

Reemplazar la parcial ∂x/∂GM del perturbador —hasta ahora por **diferencias
finitas centrales** (dos propagaciones completas del arco por evaluación de
Jacobiano)— por una **partícula variacional analítica** integrada junto a la
trayectoria en una sola propagación por sentido.

**Gate** (de `planning/MASS_FUTURE_WORK.md` §F6):
- parcial analítica vs FD **< 1e-6 relativo** ✅
- reduce el nº de propagaciones ✅ (dos props FD → cero props extra en `rebound`)

## Formulación

El objetivo (no masivo) obedece `d²r/dt² = a(r, t; GM_p)`, con la fuerza del
perturbador estudiado `a_p = GM_p (r_p − r)/|r_p − r|³` entre los términos. La
sensibilidad `s(t) = ∂[r, v]/∂GM_p` satisface la ecuación variacional de primer
orden

```
ds_r/dt = s_v
ds_v/dt = (∂a/∂r) · s_r + ∂a/∂GM_p
```

donde `∂a/∂GM_p` es el forzamiento directo (analítico) y `∂a/∂r` es el Jacobiano
de **toda** la fuerza (Sol + planetas + GR + perturbadores + su acoplamiento
mutuo). Integrar esta ecuación con el mismo integrador que la órbita da `s(t)`
exacta (a orden de máquina del integrador), sin el sesgo de paso O(δ²) de la FD
ni sus dos propagaciones.

## Implementación (backend `rebound`)

En el backend `orbdet.dynamics` el Sol y los planetas son **partículas masivas de
REBOUND**, así que el Jacobiano completo `∂a/∂r` lo conoce el integrador
variacional de REBOUND. Se usa la partícula variacional de masa nativa:

```python
var = sim.add_variation(order=1)
var.vary(perturber_particle_index, "m")   # ∂/∂ mass_msun del perturbador
sim.integrate(dt)
# var.particles[test_idx] = ∂[r, v]_test / ∂ mass_msun
```

Como `sim.G = GM_SUN` y las masas van en M_sun, `GM = GM_SUN · mass_msun`, de modo
que `∂x/∂GM = (∂x/∂mass_msun) / GM_SUN` — misma convención y unidades que la FD.
La variacional captura además el acoplamiento perturbador↔resto del sistema (el
perturbador desvía a los planetas y estos al objetivo), idéntico a lo que la FD
recupera al re-propagar todo.

Función: `orbdet.variational.partial_wrt_gm_variational`. Cableada en
`mass_determination._forward_al` (rama `rebound`) tras el flag
`_ModelConfig.gm_variational` (default `True`); la FD queda disponible con
`gm_variational=False` para validación.

## Validación

`tests/orbdet/test_variational.py::test_dgm_variational_matches_fd`: la parcial
variacional coincide con la FD central **extrapolada por Richardson** (δ y δ/2, para
quitar la truncación O(δ²)) a **< 1e-6 relativo** por época. La suite completa
`tests/orbdet/` (124 tests) y los gates sintéticos T6/T7 (`test_mass_determination`,
`test_stacking`) pasan con la variacional activada en el backend `rebound`.

## Por qué ASSIST queda por fuera (y qué haría falta)

Bajo `orbdet.dynamics_assist` el Sol, los planetas y la corrección GR (EIH) son
**fuerzas adicionales de la efeméride** (`assist.Extras`), no partículas de
REBOUND. El integrador variacional de REBOUND **no propaga** las partículas
variacionales a través de esas fuerzas, así que `var.vary(..., "m")` allí omitiría
el término dominante `∂a/∂r` (el tirón solar) y daría una sensibilidad errónea. Por
eso el backend `assist` —el de producción— mantiene la FD.

Camino futuro para ASSIST (no bloqueante; sigue siendo baja prioridad): integrar la
misma ecuación variacional pero obteniendo el producto Jacobiano-vector
`(∂a/∂r)·s_r` por **diferencia finita direccional de la aceleración** de ASSIST
(dos evaluaciones de fuerza por paso, no dos propagaciones completas), o
reimplementar analíticamente el Jacobiano de Sol/planetas/GR. Ambas exceden el
alcance de F6 y su ganancia en producción es marginal (la parcial de masa es 2 de
~15 propagaciones por objetivo; el grueso son las 12 de ∂x/∂elementos).
