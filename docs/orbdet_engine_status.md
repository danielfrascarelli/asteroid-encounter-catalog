# Motor `orbdet` — estado y arquitectura

> **Estado:** T1–T11 completas (PR #80 en `main`). T10: las 4 masas calibradoras se
> recuperan con |z| < 3 sobre Gaia FPR (N ≥ 20), χ²_red ≈ 1, tras modelar la
> correlación intra-tránsito de los CCDs de Gaia (covarianza en bloques por FOV, piso
> autocalibrado). Con N = 7 el ratio fit/lit de los calibradores cae en [1.12, 1.29];
> con N ≥ 20 en [0.943, 0.990]. T11: (16) Psyche = 2.43×10¹⁹ kg, σ_stat 3.3 %
> (ratio 1.020 frente a DE441; 1.014 frente a Fuentes-Muñoz 2025).
> **Última actualización:** 2026-06-30.
> Resultados: [`mass_determination_results.md`](mass_determination_results.md).
> Trabajo futuro: [`planning/MASS_FUTURE_WORK.md`](../planning/MASS_FUTURE_WORK.md).

## Por qué existe

La capa de masas vieja (`src/mass`, LOO secuencial) quedó **cerrada**: no produce
masas identificables (χ²(masa) multimodal, |z| > 6 en calibradores) en DR3 ni en FPR
([`docs/mass_layer_track_a_closure.md`](mass_layer_track_a_closure.md),
[`docs/mass_layer_fpr_revalidation.md`](mass_layer_fpr_revalidation.md)). La raíz
es la **degeneración masa↔drift orbital**: el LOO ajusta órbita y masa en pasos
secuenciales y, al separarlos, descarta la información que distingue una deflexión
real de un error de órbita que crece con el tiempo → χ²(masa) multimodal, mínimo
errante, no identificable.

`src/orbdet` es el motor que resuelve **órbitas + masa juntas** en un único ajuste
de mínimos cuadrados sobre el arco completo (estrategia Fuentes-Muñoz / OrbFit /
JPL), con el Jacobiano de las ecuaciones variacionales y stacking de muchos
test-asteroids por perturbador. La degeneración se maneja **dentro de la
covarianza**, no se descarta.

## Contrato de aislamiento

`orbdet` **no importa ningún otro módulo del proyecto** (`src.detect`, `src.mass`,
`src.propagate`, …). Depende solo de numpy/scipy/astropy/rebound/assist. Verificado por
[`tests/orbdet/test_isolation.py`](../tests/orbdet/test_isolation.py). Esto permite
testear el motor de forma totalmente independiente del pipeline y será el hogar
canónico de la matemática del proyecto.

Convenciones internas: AU, día, radián, M_sun; GM en AU³/día²; JD TDB; marco
eclíptico J2000 baricéntrico para la dinámica; ICRS para la observación.

## Módulos (`src/orbdet/`)

| Módulo | Rol | Tarea |
|--------|-----|-------|
| `constants.py` | Constantes físicas (GM_sun, c, oblicuidad, conversiones masa↔GM) | T1 |
| `kepler.py` | Mecánica de dos cuerpos: elementos↔estado, propagación, `dstate_delements` (Jacobiano analítico ∂[r,v]/∂elementos) | T1, T3 |
| `frames.py` | Rotaciones y derivadas (`rotation_*`, `drotation_*`), eclíptica↔ecuatorial | T1, T3 |
| `time_scales.py` | Conversiones TCB/TDB/UTC vía astropy; época Gaia → JD TDB | T1 |
| `dynamics.py` | Dinámica N-cuerpos (rebound): Sol + planetas + asteroides perturbadores | T2 |
| `variational.py` | **Ecuaciones variacionales**: matriz de transición de estado Φ(t) analítica (`add_variation`), ∂x/∂elementos (Φ·J_elem), ∂x/∂GM por diferencias finitas + Richardson | T3 |
| `observation.py` | **Modelo de observación**: estado→ICRS→RA/Dec, light-time iterativa, covarianza along-scan/across-scan anisotrópica, Jacobiano observacional | T4 |
| `least_squares.py` | Corrector **Levenberg-Marquardt** genérico (residuos blanqueados, covarianza = (JᵀJ)⁻¹) | T5 |
| `orbit_determination.py` | **OD por mínimos cuadrados** de los 6 elementos sobre el arco completo | T5 |
| `mass_determination.py` | **Ajuste conjunto órbita+masa** (T6) y **stacking multi-objetivo** (T7, sistema en flecha 1+6N); `backend="assist"` con parciales por FD sobre ASSIST; **covarianza diagonal en bloques por FOV** (`_block_whiten`) + **calibración del piso sistemático** (`calibrate_sys_floor`) para χ²_red≈1 | T6, T7, T8 |
| `dynamics_assist.py` | **Modelo de fuerzas ASSIST**: efeméride JPL DE440 + 8 planetas/Luna/Plutón + GR (EIH) + 16 perturbadores asteroidales masivos con masa variable | T8 |
| `gaia_adapter.py` | **Adaptador de datos reales**: σ_AL por proyección de la covarianza (RA,Dec) de Gaia, MPCORB→`KeplerElements`, armonización de épocas N-cuerpos, **agrupación por cruce FOV** (`fov_groups_from_epochs`, para la covarianza en bloques), ensamblado de `TargetObservations` | T9 |

## Estado por tarea y gates verificados

| # | Tarea | Estado | Gate (verde) |
|---|-------|--------|--------------|
| T1 | Esqueleto + primitivas | ✅ | round-trips, invariantes, aislamiento |
| T2 | Dinámica N-cuerpos | ✅ | límite dos-cuerpos vs Kepler a 1e-8 AU; cross-check Horizons cerrado por T8 (ASSIST vs Horizons 0.17 mas/900 d) |
| T3 | Ecuaciones variacionales | ✅ | ∂x/∂elem analítico vs FD <1e-6; meseta de Richardson para ∂x/∂GM |
| T4 | Observación + covarianza AL | ✅ | chain N-cuerpos vs oráculo kepleriano <0.1 mas; ruido AL → χ²/obs≈1 |
| T5 | Corrector diferencial | ✅ | recupera órbita sintética: sin ruido χ²<1e-6, con ruido χ²_red≈1, <5σ |
| T6 | Ajuste conjunto órbita+masa | ✅ | **closing-loop**: masa inyectada ratio≈1.0 (sin ruido <2e-3; con ruido <3σ, σ informativa) |
| T7 | Stacking multi-asteroide | ✅ | **σ(GM)∝1/√N** (s2/s1≈1/√2, s4/s1≈0.5 a <5%) |
| T8 | Fuerzas + pesos completos | ✅ | ASSIST vs Horizons 0.17 mas; **χ²_red≈1 sobre datos reales** vía covarianza en bloques por FOV (piso `s_c` autocalibrado) + clip 4σ. Big-4 FPR: χ²_red∈[0.97,1.00] |
| T9 | Adaptador FPR/DR3 → motor | ✅ | adaptador + `fov_groups_from_epochs`; `scripts/mass/orbdet_fit_realdata.py` corre **Big-4 end-to-end sobre FPR real** |
| T10 | Validación literatura | ✅ | 4/4 con |z| < 3 (N ≥ 20): Ceres z = −1.01, Vesta −1.30, Hygiea −0.13, Pallas +2.67 (N = 6). Ratio fit/lit en [0.943, 0.990] para los 3 con N ≥ 20. Cruce Fuentes-Muñoz 2025: Psyche ratio 1.014 (z = +0.25) |
| T11 | Producción + catálogo | ✅ | barrido de 16 perturbadores (`build_mass_catalog.py`, `mass_determination_results.md`). (16) Psyche = 2.43×10¹⁹ kg, σ_stat 3.3 % (ratio 1.020 frente a DE441). 6 perturbadores con deflexión débil: ratio en [0.39, 0.72] (ver trabajo futuro) |

PRs de la sesión 2026-06-28: T3 #73, T4 #74, T5 #75, T6 #76, T7 #77.
Sesión 2026-06-29: maquinaria T8/T9 (`dynamics_assist`, `gaia_adapter`); run Big-4
FPR + covarianza en bloques por FOV (`_block_whiten`, `calibrate_sys_floor`) → T10.

## Determinabilidad de masas: qué muestran los resultados

- **Closing-loop sintético (T6/T7).** El ajuste simultáneo recupera una masa inyectada
  con ratio 1.0 ± 2×10⁻³ (sin ruido) y dentro de 3σ con ruido AL; la incertidumbre
  escala como 1/√N al apilar objetivos (s2/s1 ≈ 1/√2, s4/s1 ≈ 0.5, error < 5 %).
- **Datos reales (FPR).** El run Big-4 con N ≥ 20 recupera las 4 masas calibradoras
  con |z| < 3 y χ²_red ≈ 1. La causa del cierre de Track A queda acotada al método
  secuencial (LOO), no al leverage de la astrometría.
- **Modelo de error.** Con covarianza diagonal, σ(masa) se subestima en un factor 1.66
  (los 7 CCDs correlacionados por cruce FOV, ICC = 0.32, se tratan como
  independientes). La covarianza en bloques por FOV lo corrige. El sesgo del estimador
  sobre la geometría real es compatible con cero (closing-loop, ratio medio 0.997 sobre
  3 semillas).
- **Independencia del catálogo de encuentros.** El motor `orbdet` es ortogonal a la
  detección/caracterización (72.236.904 encuentros, congelados en `FROZEN_RUN.md`).

## Trabajo futuro

T1–T11 están cerradas (tabla arriba). Los ítems abiertos (σ externa por-perturbador,
regularización masa↔órbita, acotar el sesgo de −4 %, perturbadores fuera de los 16,
DR4) están en [`planning/MASS_FUTURE_WORK.md`](../planning/MASS_FUTURE_WORK.md). El
cruce con Fuentes-Muñoz 2025 está hecho (resultados en
[`mass_determination_results.md`](mass_determination_results.md)).

## Maquinaria de datos reales (sesión 2026-06-29)

- **Modelo de fuerzas ASSIST (T8):** reemplaza los planetas integrados libremente
  por la efeméride JPL DE440 leída en cada paso + GR (EIH) + 16 perturbadores
  asteroidales masivos (masa variable para el perturbador bajo estudio). Bajo
  ASSIST las partículas variacionales de rebound no propagan a través de las
  fuerzas de la efeméride, así que **todas** las parciales (∂x/∂elementos y
  ∂x/∂masa) salen por diferencias finitas centrales sobre `propagate_assist`.
  Gate verde: vs Horizons **0.17 mas** sobre 900 d (vs 239.69 mas de planetas
  libres, **1404×**); closing-loop ASSIST recupera la masa inyectada a ratio≈1.
- **Adaptador Gaia (T9):** σ_AL proyectando la elipse de covarianza (RA,Dec)
  completa de Gaia (sistemática + aleatoria) sobre la dirección de barrido;
  MPCORB (grados)→`KeplerElements`; armonización de épocas propagando con el
  propio N-cuerpos (no Kepler de dos cuerpos); ensamblado de `TargetObservations`.
  Unit-tested contra la forma cuadrática explícita y un round-trip sintético.
- **Covarianza en bloques por FOV (T8/T10):** los CCDs de un mismo cruce de plano
  focal (~7, separados ~5 s) comparten error sistemático → residuos correlacionados
  (ICC≈0.32 medido en datos reales). `_block_whiten` aplica `C_bloque = diag(σ_AL²)
  + s_c²·11ᵀ` por grupo FOV (Cholesky por bloque); `calibrate_sys_floor` fija el piso
  correlacionado `s_c` por bisección para χ²_red ≈ 1. Corrige la subestimación de
  σ(masa) (factor 1.66 con covarianza diagonal). `fov_groups_from_epochs` arma los
  grupos. Verificado contra `C⁻¹` explícito (`tests/orbdet/test_block_covariance.py`).
- **Dependencia nueva:** `assist>=1.1` (requiere `rebound` 4.x). Efemérides
  (`linux_p1550p2650.440`, `sb441-n16.bsp`, ~750 MB) en `$ORBDET_EPHEM_DIR`
  (default `data/raw/ephem`), no versionadas.

## Notas operativas

- Lint/format/mypy: `docker run --rm -v "$PWD":/app -w /app gaia-asteroid-encounters <cmd>`
  (el servicio `pipeline` de compose corre contra el source baked en la imagen).
- La suite `tests/orbdet/` completa puede colgar **localmente** por la descarga
  IERS/efemérides de astropy sin red; correr archivos sueltos con
  `-p no:cacheprovider`. En CI (con red) corre normal.
- Tests pesados (rebound) marcados `slow`; no se deseleccionan por defecto (solo
  `horizons`, que requiere red, lo está).
