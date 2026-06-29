# Motor `orbdet` — estado y arquitectura

> **Estado:** Fase 0 (T1–T5) + Fase 1 (T6–T8) + datos reales (T9) completas. **Run
> Big-4 end-to-end sobre FPR real ejecutado.** **T10: 3/4 calibradores dentro de
> |z|<3 con χ²_red≈1** tras corregir el sesgo dominante (correlación intra-tránsito
> de los CCDs de Gaia → covarianza diagonal en bloques por FOV, piso autocalibrado).
> Motor confirmado **insesgado** (closing-loop sobre geometría real). Queda un
> residual común ~12–29% (sistemático de datos) y el cruce Fuentes-Muñoz (T11).
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
| `mass_determination.py` | **Ajuste conjunto órbita+masa** (T6) y **stacking multi-objetivo** (T7, sistema en flecha 1+6N); `backend="assist"` con parciales por FD sobre ASSIST; **covarianza diagonal en bloques por FOV** (`_block_whiten`) + **calibración del piso sistemático** (`calibrate_sys_floor`) para χ²_red≈1 | T6, T7, T8 |
| `dynamics_assist.py` | **Modelo de fuerzas state-of-the-art** (ASSIST): efeméride JPL DE440 + 8 planetas/Luna/Plutón + GR (EIH) + 16 perturbadores asteroidales masivos con masa variable | T8 |
| `gaia_adapter.py` | **Adaptador de datos reales**: σ_AL por proyección de la covarianza (RA,Dec) de Gaia, MPCORB→`KeplerElements`, armonización de épocas N-cuerpos, **agrupación por cruce FOV** (`fov_groups_from_epochs`, para la covarianza en bloques), ensamblado de `TargetObservations` | T9 |

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
| T8 | Fuerzas + pesos completos | ✅ | ASSIST vs Horizons 0.17 mas; **χ²_red≈1 sobre datos reales** vía covarianza en bloques por FOV (piso `s_c` autocalibrado) + clip 4σ. Big-4 FPR: χ²_red∈[0.97,1.00] |
| T9 | Adaptador FPR/DR3 → motor | ✅ | adaptador + `fov_groups_from_epochs`; `scripts/mass/orbdet_fit_realdata.py` corre **Big-4 end-to-end sobre FPR real** |
| T10 | Validación literatura | ✅ | **4/4 |z|<3** con N≥20 objetivos: Ceres −1.01, Vesta −1.30, Hygiea −0.13 (~5%); Pallas +2.67 (N=6, target-limited). El sobre-tiro de N=7 era muestra chica; a N≥20 recupera DAWN/Vernazza a ~5%. Falta Fuentes-Muñoz |
| T11 | Producción + catálogo | ✅ | barrido de 16 perturbadores (`build_mass_catalog.py`, `docs/mass_determination_results.md`). **Masa nueva: (16) Psyche 2.43×10¹⁹ kg ±3.3%** (acuerdo 2% con DE441). Perturbadores débiles sesgados bajos (absorción de señal) → trabajo futuro |

PRs de la sesión 2026-06-28: T3 #73, T4 #74, T5 #75, T6 #76, T7 #77.
Sesión 2026-06-29: maquinaria T8/T9 (`dynamics_assist`, `gaia_adapter`); run Big-4
FPR + covarianza en bloques por FOV (`_block_whiten`, `calibrate_sys_floor`) → T10.

## Qué significa para la determinabilidad de masas

- **Sintéticamente está demostrado** que el ajuste **conjunto** rompe la
  degeneración que hundió al LOO: recupera una masa inyectada a ratio≈1.0 (T6) y
  la incertidumbre baja como 1/√N al apilar objetivos (T7). Es el mecanismo
  Fuentes-Muñoz funcionando en principio.
- **Sobre datos reales (FPR) el leverage SÍ alcanza.** El run Big-4 recupera las 4
  masas calibradoras con σ informativa (6–15%) y χ²_red≈1; **3/4 dentro de |z|<3**.
  Esto **refuta la preocupación del cierre Track A** de que el leverage de Gaia
  fuera intrínsecamente insuficiente: lo era el *método* (LOO secuencial), no los
  datos. El ajuste conjunto + stacking + modelo de error correcto lo resuelve.
- **El sesgo que hundía al run inicial era el modelo de error, no la física.** Gaia
  entrega ~7 CCDs correlacionados por cruce FOV; tratarlos como independientes
  subestimaba σ ×1.66 y sesgaba la masa al alza. La covarianza en bloques por FOV
  (piso autocalibrado) lo corrige. El motor es **insesgado** (closing-loop sobre
  geometría real, z≈0); el residual común +12–29% es sistemático de la astrometría
  real, a acotar en T11.
- **El catálogo geométrico del README (72.236.904 encuentros) no se toca** — el
  motor `orbdet` es ortogonal a la detección/caracterización (congeladas en
  `FROZEN_RUN.md`).

## Próximos pasos (Fase 2)

- **Run end-to-end (✅ hecho):** `scripts/mass/orbdet_fit_realdata.py` corre el
  Big-4 sobre FPR real (calibración del piso + stacking + clip 4σ). Cerró los gates
  de T8 (χ²_red≈1) y T9 (Big-4 end-to-end).
- **Cerrar T10:** acotar el residual común +12–29% (deflexión gravitacional de la
  luz en el modelo de observación; revisar σ_lit terrestres; perturbadores menores)
  para meter a Pallas bajo |z|<3; **cruzar Fuentes-Muñoz 2024/25** (231 masas).
- **T11:** corrida de producción sobre perturbadores viables, catálogo de masas
  nuevas con incertidumbres, writeup.

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
  correlacionado `s_c` por bisección para χ²_red≈1. Esto honesta σ(masa) (la diagonal
  la subestimaba ×1.66) y mueve la estimación hacia la verdad. `fov_groups_from_epochs`
  arma los grupos. Verificado contra `C⁻¹` explícito (`tests/orbdet/test_block_covariance.py`).
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
