# Motor `orbdet` — estado y arquitectura

> **Estado:** Fase 0 completa (T1–T5) + núcleo de Fase 1 (T6, T7) + **maquinaria
> de datos reales (T8 modelo de fuerzas ASSIST + T9 adaptador Gaia) construida y
> validada** (vs JPL Horizons a sub-mas). Falta el run end-to-end sobre FPR real y
> la validación contra literatura (T10–T11), donde se decide el veredicto científico.
> **Última actualización:** 2026-06-29.
> Roadmap detallado: [`planning/MASS_DETERMINATION_PLAN.md`](../planning/MASS_DETERMINATION_PLAN.md).

## Por qué existe

La capa de masas vieja (`src/mass`, LOO secuencial) quedó **cerrada**: no determina
ninguna masa defendible en DR3 ni en FPR
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
| `mass_determination.py` | **Ajuste conjunto órbita+masa** (T6) y **stacking multi-objetivo** (T7, sistema en flecha 1+6N); `backend="assist"` con parciales por FD sobre ASSIST | T6, T7, T8 |
| `dynamics_assist.py` | **Modelo de fuerzas state-of-the-art** (ASSIST): efeméride JPL DE440 + 8 planetas/Luna/Plutón + GR (EIH) + 16 perturbadores asteroidales masivos con masa variable | T8 |
| `gaia_adapter.py` | **Adaptador de datos reales**: σ_AL por proyección de la covarianza (RA,Dec) de Gaia, MPCORB→`KeplerElements`, armonización de épocas N-cuerpos, ensamblado de `TargetObservations` | T9 |

## Estado por tarea y gates verificados

| # | Tarea | Estado | Gate (verde) |
|---|-------|--------|--------------|
| T1 | Esqueleto + primitivas | ✅ | round-trips, invariantes, aislamiento |
| T2 | Dinámica N-cuerpos | ✅ | límite dos-cuerpos vs Kepler a 1e-8 AU (cross-check Horizons pendiente de red) |
| T3 | Ecuaciones variacionales | ✅ | ∂x/∂elem analítico vs FD <1e-6; meseta de Richardson para ∂x/∂GM |
| T4 | Observación + covarianza AL | ✅ | chain N-cuerpos vs oráculo kepleriano <0.1 mas; ruido AL → χ²/obs≈1 |
| T5 | Corrector diferencial | ✅ | recupera órbita sintética: sin ruido χ²<1e-6, con ruido χ²_red≈1, <5σ |
| T6 | Ajuste conjunto órbita+masa | ✅ | **closing-loop**: masa inyectada ratio≈1.0 (sin ruido <2e-3; con ruido <3σ, σ informativa) |
| T7 | Stacking multi-asteroide | ✅ | **σ(GM)∝1/√N** (s2/s1≈1/√2, s4/s1≈0.5 a <5%) |
| T8 | Fuerzas + pesos completos | 🟡 | **maquinaria + validación verde**: ASSIST vs Horizons 0.17 mas (1404× mejor que planetas libres), closing-loop ASSIST ratio≈1. Falta χ²_red≈1 + outliers/debiasing sobre datos reales |
| T9 | Adaptador FPR/DR3 → motor | 🟡 | **adaptador construido + unit-tested sintético**. Falta el script IO que corre un Big-4 end-to-end sobre FPR real |
| T10 | Validación literatura | ⬜ | \|z\|<3 en los 4 calibradores |
| T11 | Producción + catálogo | ⬜ | ≥1 masa nueva defendible |

PRs de la sesión 2026-06-28: T3 #73, T4 #74, T5 #75, T6 #76, T7 #77.
Sesión 2026-06-29: maquinaria T8 (`dynamics_assist`) + T9 (`gaia_adapter`).

## Qué significa para la determinabilidad de masas

- **Sintéticamente está demostrado** que el ajuste **conjunto** rompe la
  degeneración que hundió al LOO: recupera una masa inyectada a ratio≈1.0 (T6) y
  la incertidumbre baja como 1/√N al apilar objetivos (T7). Es el mecanismo
  Fuentes-Muñoz funcionando en principio.
- **El veredicto sobre datos reales sigue pendiente.** El closing-loop usa un
  encuentro cercano *construido* con leverage garantizado. Si DR3/FPR real tiene
  leverage suficiente lo decide **T10** (reproducir los 4 calibradores con
  |z|<3), no este test sintético. El cierre Track A advirtió que el leverage de
  DR3 puede ser intrínsecamente bajo — eso no lo cambia la metodología, solo
  más/mejores datos (FPR + stacking masivo).
- **El catálogo geométrico del README (72.236.904 encuentros) no se toca** — el
  motor `orbdet` es ortogonal a la detección/caracterización (congeladas en
  `FROZEN_RUN.md`).

## Próximos pasos (Fase 1 restante + Fase 2)

- **Run end-to-end (cuello de botella):** script IO en `scripts/mass` que carga
  FPR/DR3 real (parquet, `src/ingest`, flag `release`), filtra `is_rejected`, y
  alimenta `gaia_adapter.build_target_observations` → `determine_shared_mass`
  (`backend="assist"`, fondo de `big_asteroid_perturbers`) para un Big-4. Cierra
  los gates literales de T8 (χ²_red≈1 + rechazo de outliers/debiasing) y T9
  (Big-4 end-to-end). El IO vive fuera de `orbdet` por el contrato de aislamiento.
- **T10:** reproducir Ceres/Vesta/Pallas/Hygiea; cruzar Fuentes-Muñoz 2024/25.
- **T11:** corrida de producción, catálogo de masas nuevas, writeup.

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
