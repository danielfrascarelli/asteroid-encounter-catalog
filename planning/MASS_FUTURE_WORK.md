# Determinación de masas — ítems abiertos y mejoras

> **Estado:** 🟡 CASI CERRADO — F1–F7 hechos; F6 cerrado para el backend `rebound`
> (2026-07-04), bloqueado en `assist`; solo queda F8 (espera DR4).
> **Última actualización:** 2026-07-04.
> Trabajo futuro sobre el motor de masas `src/orbdet/`. El plan original
> (T1–T11) está completo; sus resultados están en
> [`docs/mass_determination_results.md`](../docs/mass_determination_results.md) y la
> arquitectura en [`docs/orbdet_engine_status.md`](../docs/orbdet_engine_status.md).
> Este documento lista únicamente lo que queda por hacer.

## Tabla de estado

| # | ítem | prioridad | estado | gate de aceptación |
|---|------|-----------|--------|--------------------|
| F1 | σ externa por-perturbador (jackknife/bootstrap) | alta | ✅ | σ_jackknife reportada para los 16 (PR #83); reemplaza la σ formal donde la excede |
| F2 | Regularización del par masa↔órbita en perturbadores débiles | alta | ✅ | criterio de identificabilidad (measured/not_identifiable/cota) por curvatura de χ² (PR #83) |
| F3 | Acotar el sesgo medio −4 % de los calibradores | media | ✅ (hipótesis refutada) | fondo 16→35 mueve masas <0.25%, f_sys 4.16%→4.26%: el −4% NO es incompletitud del fondo. Ver [`docs/mass_f3_background_extension.md`](../docs/mass_f3_background_extension.md) (PR #89) |
| F4 | Catálogo de perturbadores más allá de los 16 de `sb441-n16.bsp` | media | ✅ | rama custom (flag `--perturber-orbit-source`); (19) Fortuna χ²_red=0.977, (9) Metis 0.981. Ver [`docs/mass_layer_f4_design.md`](../docs/mass_layer_f4_design.md) (PR #89) |
| F5 | Extender el cruce Fuentes-Muñoz a perturbadores débiles | media | ✅ | 10/10 medidas en \|z\|<3 con σ jackknife (vs 5/10 formal). Ver [`docs/mass_crosscheck_jack.md`](../docs/mass_crosscheck_jack.md) (PR #90) |
| F6 | ∂x/∂GM analítico (partícula variacional de masa) | baja | ✅ backend `rebound` / 🔒 `assist` bloqueado | parcial analítica vs FD < 1e-6 ✅; ahorra 2 props/Jacobiano. `rebound` usa `partial_wrt_gm_variational`; `assist` sigue FD (fuerzas de efeméride no propagan variacionales). Ver [`docs/mass_layer_f6_analytic_gm.md`](../docs/mass_layer_f6_analytic_gm.md) |
| F7 | Cruce con masas terrestres (Goffin 2014, Galád 2002) donde solape | baja | ✅ | ratio y z reportados para los perturbadores en común (PR #84) |
| F8 | Perturbadores target-limited (Pallas) con Gaia DR4 | baja | ⬜ (espera datos) | N(Pallas) ≥ 20 con baseline DR4 |

## Detalle

### F1 — σ externa por-perturbador
**Problema medido.** Con N ≥ 20 la σ de Fisher cae por debajo de 0.2 %
([resultados §Dependencia con N](../docs/mass_determination_results.md)), pero el
sesgo real frente a literatura es 1–6 % en los calibradores y mayor en perturbadores
débiles. La σ formal no captura el error de regresión masa↔órbita.
**Acción.** Estimar σ por jackknife (dejar-un-objetivo-fuera) o bootstrap sobre el
conjunto de objetivos por perturbador; reportar la mayor entre σ_formal y σ_externa.
**Gate.** σ_jackknife disponible para los 16 perturbadores; los z del cruce
Fuentes-Muñoz recomputados con ella.

### F2 — Regularización masa↔órbita en perturbadores débiles
**Problema medido.** 6 perturbadores con χ²_red ≈ 1 y N ≥ 20 dan ratio fit/DE441 en
[0.39, 0.72]; el ajuste explica la astrometría con menor masa y órbita reajustada
(regresión hacia masa nula cuando la deflexión queda bajo el ruido por-encuentro).
**Acción.** Evaluar prior/penalización sobre los Δ-elementos del objetivo, o un test
de identificabilidad (curvatura de χ² respecto a la masa) que decida si el perturbador
admite medida o sólo cota inferior.
**Gate.** Los 6 perturbadores con ratio < 0.8 o bien recuperan la masa de referencia
dentro de la σ de F1, o se reportan como cota inferior con criterio explícito.

### F3 — Sesgo medio −4 % de los calibradores
**Problema medido.** Ratio medio de Ceres/Vesta/Hygiea (N ≥ 20) = 0.96; el déficit de
4 % es sistemático respecto a DAWN/Vernazza.
**Hipótesis.** Completitud del fondo: el modelo incluye 16 perturbadores; los menores
no modelados desvían las órbitas de los objetivos.
**Acción.** Ampliar el conjunto de perturbadores de fondo y medir el cambio en f_sys.
**Gate.** f_sys (RMS de ratio − 1 sobre calibradores) por debajo de 4.2 %.

### F4 — Perturbadores fuera de los 16 de la efeméride
**Restricción actual.** El motor toma la órbita del perturbador de `sb441-n16.bsp`
(su masa es el parámetro libre). Sólo hay 16 cuerpos con órbita en esa efeméride.
**Acción.** Integrar la órbita del perturbador en el propio modelo N-cuerpos (en vez de
leerla de la efeméride) para admitir cuerpos sin entrada en sb441-n16.
**Gate.** ≥ 1 perturbador no presente en sb441-n16 ajustado con χ²_red ∈ [0.95, 1.05].

### F5 — Extensión del cruce Fuentes-Muñoz
Una vez disponible la σ de F1, recomputar z para los 12 perturbadores no-calibradores
y reportar cuántos quedan en |z| < 3. Hoy 4/8 fiables están en |z| < 3 con σ formal,
limitado por la subestimación de σ en perturbadores débiles.

### F6 — Parcial de masa analítica ✅ (backend `rebound`)
**Hecho (2026-07-04).** `orbdet.variational.partial_wrt_gm_variational` integra
∂x/∂GM con la partícula variacional de masa de REBOUND (`Variation.vary(i, "m")`)
en una sola propagación por sentido; cableada en `mass_determination._forward_al`
(rama `rebound`) tras `gm_variational=True`. Gate cumplido: coincide con la FD
(extrapolada Richardson) a **< 1e-6** (`test_dgm_variational_matches_fd`) y ahorra
las dos propagaciones de la FD. Detalle: [`docs/mass_layer_f6_analytic_gm.md`](../docs/mass_layer_f6_analytic_gm.md).

**Pendiente / bloqueado en `assist` (producción).** Bajo ASSIST el Sol/planetas/GR
son fuerzas adicionales de la efeméride que no propagan partículas variacionales de
REBOUND, así que allí ∂x/∂GM (y ∂x/∂elementos) siguen por diferencias finitas
centrales sobre `propagate_assist`. Camino futuro (no bloqueante): integrar la
variacional con el producto Jacobiano-vector `(∂a/∂r)·s` por diferencia finita
direccional de la *aceleración* de ASSIST (2 evals de fuerza/paso, no 2
propagaciones), o Jacobiano analítico de Sol/planetas/GR. Ganancia marginal en
producción (la parcial de masa es 2 de ~15 propagaciones; el grueso son las 12 de
∂x/∂elementos).

### F7 — Cruce con masas terrestres
Goffin (2014) y Galád (2002) reportan masas por astrometría terrestre. Comparar donde
solape con los 16 perturbadores aporta una referencia independiente de Gaia/DE441.

### F8 — Perturbadores target-limited
Pallas tiene 6–7 encuentros < 0.05 AU en el catálogo DR3 completo. Gaia DR4 (baseline y
densidad mayores) aumentará el número de objetivos. El motor soporta DR4 por el flag
`release`.

## Fuera de scope

- Reabrir el enfoque por-encuentro LOO (`src/mass/`): cerrado y superado por el ajuste
  simultáneo (ver [`docs/mass_layer_track_a_closure.md`](../docs/mass_layer_track_a_closure.md)).
- Detección/caracterización de encuentros: congelada (ver [`FROZEN_RUN.md`](../FROZEN_RUN.md)).
