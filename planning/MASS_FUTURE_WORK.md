# Determinación de masas — ítems abiertos y mejoras

> **Estado:** ⬜ PENDIENTE (ninguno bloquea los resultados ya medidos).
> **Última actualización:** 2026-06-30.
> Trabajo futuro sobre el motor de masas `src/orbdet/`. El plan original
> (T1–T11) está completo; sus resultados están en
> [`docs/mass_determination_results.md`](../docs/mass_determination_results.md) y la
> arquitectura en [`docs/orbdet_engine_status.md`](../docs/orbdet_engine_status.md).
> Este documento lista únicamente lo que queda por hacer.

## Tabla de estado

| # | ítem | prioridad | estado | gate de aceptación |
|---|------|-----------|--------|--------------------|
| F1 | σ externa por-perturbador (jackknife/bootstrap) | alta | ⬜ | σ_jackknife reportada para los 16; reemplaza la σ formal donde la excede |
| F2 | Regularización del par masa↔órbita en perturbadores débiles | alta | ⬜ | los 6 perturbadores con ratio < 0.8 dejan de regresar a la baja, o se declaran cota inferior |
| F3 | Acotar el sesgo medio −4 % de los calibradores | media | ⬜ | f_sys medido < 4.2 % al extender el fondo de perturbadores |
| F4 | Catálogo de perturbadores más allá de los 16 de `sb441-n16.bsp` | media | ⬜ | ≥ 1 perturbador fuera de los 16 ajustado con χ²_red ∈ [0.95, 1.05] |
| F5 | Extender el cruce Fuentes-Muñoz a perturbadores débiles | media | ⬜ | reportar z con σ de F1 para los 12 no-calibradores |
| F6 | ∂x/∂GM analítico (partícula variacional sobre ASSIST) | baja | ⬜ | parcial analítica coincide con FD a < 1e-6 relativo; reduce nº de propagaciones |
| F7 | Cruce con masas terrestres (Goffin 2014, Galád 2002) donde solape | baja | ⬜ | ratio y z reportados para los perturbadores en común |
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

### F6 — Parcial de masa analítica
**Estado actual.** Bajo ASSIST, las parciales (∂x/∂elementos y ∂x/∂GM) se calculan por
diferencias finitas centrales sobre `propagate_assist` (las partículas variacionales de
rebound no propagan a través de las fuerzas de la efeméride).
**Acción.** Implementar la partícula variacional respecto a la masa o un esquema
adjunto, reduciendo el número de propagaciones por evaluación de Jacobiano.
**Gate.** Parcial analítica vs FD < 1e-6 relativo; reducción medible de tiempo de cómputo.

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
