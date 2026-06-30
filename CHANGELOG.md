# Changelog

## 2026-06-30 — Determinación de masas con motor propio (`orbdet`, T1–T11 ✅, PR #80)

### Features
- **Motor `src/orbdet/` — "OrbFit de cero"** (aislado del resto de `src.*`, verificado por test). Determinación de órbitas + masa por mínimos cuadrados sobre el arco completo:
  - `kepler/frames/time_scales/constants` (T1), `dynamics` N-cuerpos rebound (T2), `variational` (T3, ∂x/∂elem analítico + ∂x/∂GM por FD/Richardson), `observation` (T4, RA/Dec + light-time + covarianza along-scan anisotrópica), `least_squares`/`orbit_determination` (T5, Levenberg-Marquardt).
  - `mass_determination` (T6/T7/T8): ajuste **conjunto** órbita+masa, stacking multi-objetivo (sistema en flecha 1+6N), **covarianza diagonal en bloques por FOV** (`_block_whiten`) + autocalibración del piso sistemático (`calibrate_sys_floor`) para χ²_red≈1.
  - `dynamics_assist` (T8): modelo de fuerzas ASSIST (DE440 + GR EIH + 16 perturbadores asteroidales); vs Horizons 0.17 mas sobre 900 d.
  - `gaia_adapter` (T9): σ_AL por proyección de la covarianza Gaia, MPCORB→elementos, grupos por cruce FOV (`fov_groups_from_epochs`).
- **Scripts de masas** (`scripts/mass/`): `orbdet_fit_realdata.py` (Big-4 end-to-end sobre FPR real, stacking + clip 4σ + paralelización ~6×), `build_mass_catalog.py` (catálogo con `σ_total=√(σ_stat²+(f_sys·M)²)`).

### Resultados científicos
- **Validación (T10): 4/4 calibradores dentro de |z|<3** sobre Gaia FPR con N≥20 objetivos (Ceres/Vesta/Hygiea a ~5%; Pallas +2.67, target-limited a N=6). **Refuta el cierre Track A**: el problema era el método (LOO secuencial), no el leverage de Gaia.
- **Producción (T11): masa nueva defendible — (16) Psyche = 2.43×10¹⁹ kg ±3.3%** (acuerdo 2% con DE441). Hallazgo: perturbadores con deflexión débil se sesgan bajos (absorción de señal masa↔órbita).
- Detalle: [docs/mass_determination_results.md](docs/mass_determination_results.md), [docs/orbdet_engine_status.md](docs/orbdet_engine_status.md).

### Dependencias
- `assist>=1.1` (requiere `rebound` 4.x); efemérides (`linux_p1550p2650.440`, `sb441-n16.bsp`, ~750 MB) en `$ORBDET_EPHEM_DIR` (no versionadas).

## 2026-05-24 — Corrida congelada + catálogo híbrido + cierre capa de masas (LOO)

### Features
- **Catálogo congelado** ([FROZEN_RUN.md](FROZEN_RUN.md)): 72.236.904 candidatos geométricos a 0.05 AU sobre 98.775 numerados, con hashes de inputs/outputs.
- **Caracterización por streaming**: `encounters_characterized_full.parquet` (72M filas, 18.9% Gaia-observables) sin OOM.
- **Catálogo híbrido N-body** (`encounters_catalog_hybrid_stageb.parquet`): re-refinamiento N-body del subset crítico (`q_min<1.8 ∨ e_max>0.3`, 12%); `refinement_method` por fila.
- **Recall del prefiltro medido**: 76.4% en cola adversa; fix radial-overlap → 100% ([docs/prefilter_recall.md](docs/prefilter_recall.md)).
- **Validación consolidada** ([docs/literature_validation.md](docs/literature_validation.md)): gate 4/4, Fienga 3/4, Galád 4/4, Fuentes-Muñoz 2025 11.804/40.004 confirmaciones.
- **Reorganización de `scripts/`** en `{ingest,pipeline,mass,validate,bench,dev}/`.

### Capa de masas (enfoque viejo LOO)
- **Track A CERRADO**: el LOO secuencial no determina masas en DR3 (no-identificabilidad masa↔drift); FPR-solo tampoco lo reabre con ese método. Posteriormente **superado** por el motor `orbdet` (ver entrada 2026-06-30).

## 2026-05-18 — Validación contra literatura + propagación N-body

### Features
- **Multi-snapshot MPCORB** ([src/ingest/mpcorb_archive.py](src/ingest/mpcorb_archive.py)): descubrir y seleccionar el snapshot MPCORB de época más cercana a la ventana temporal. Reduce error de Kepler de ~0.03 AU a < 0.001 AU para la ventana Gaia.
- **download_mpcorb_historical.py**: descarga snapshots históricos de MPCORB.DAT vía Wayback Machine.
- **Propagación N-body con `rebound`** ([src/propagate/nbody.py](src/propagate/nbody.py)): WHFast (default) o IAS15, Sol + planetas + Ceres/Vesta/Pallas/Hygiea opcionales.
- **Cache de trayectorias** ([src/propagate/cache.py](src/propagate/cache.py)): `np.memmap` para trayectorias (T, N, 3); streaming a disco evita OOM en runs grandes.
- **Validador Fienga 2003** ([scripts/validate_fienga_2003.py](scripts/validate_fienga_2003.py)): cross-match contra VizieR J/A+A/406/751.
- **Validador Galád 2002** ([scripts/validate_galad_2002.py](scripts/validate_galad_2002.py)): cross-match contra HTML del paper.
- **Validador JPL Horizons** ([scripts/validate_jpl_horizons.py](scripts/validate_jpl_horizons.py)): cross-check 3-way (nuestro vs literatura vs JPL).
- **Comparador Kepler vs rebound** ([scripts/compare_kepler_vs_rebound.py](scripts/compare_kepler_vs_rebound.py)): diff de catálogos + deltas estadísticos.
- **Makefile** con atajos para los comandos más comunes.

### Fixes
- **Bug: refinement Kepler sobre catálogos rebound** ([src/detect/refine.py](src/detect/refine.py)): `refine_candidates` ahora acepta `positions` + `time_grid` y, cuando los recibe, hace interpolación cuadrática sobre el cache N-body en lugar de re-propagar con Kepler. Antes el catálogo rebound era *efectivamente* un catálogo Kepler.
- **Bug: parallel scan con `positions=memmap` se colgaba en 100k+** ([src/detect/parallel.py](src/detect/parallel.py)): cuando positions es un `np.memmap`, pasar solo el filename a los workers (cada worker re-abre la memmap, OS comparte páginas via page cache).
- **Bug: parallel scan con prefilter activo + N pequeño se colgaba** ([src/detect/parallel.py](src/detect/parallel.py)): cuando `pairs.nbytes > 1 MB`, escribir a tempfile y pasar path en lugar de array pickleado por initargs.
- **Bug: OOM en integración rebound 100k** ([src/propagate/cache.py](src/propagate/cache.py), [src/propagate/nbody.py](src/propagate/nbody.py)): la trayectoria (29.5 GB float32) se escribía dos veces (in-memory + disco). Ahora streamea directo al memmap mediante un parámetro `out=` en `propagate_grid_nbody`.

### Tests
- 3 archivos nuevos en `tests/`: `test_mpcorb_archive.py` (5 tests), `test_cache.py` (6 tests), `test_parallel_memmap.py` (1 integration), `test_refine_cache.py` (2 tests), `test_parallel_pairs_spill.py` (1 test). Total: 198 tests passing (era 183).

### CI
- ruff + black + mypy todos verdes.

### Docs
- [ROADMAP.md](ROADMAP.md): estado actual + bugs conocidos.
- [VALIDATION_SUMMARY.md](VALIDATION_SUMMARY.md): tabla de validación contra Fienga + Galád + JPL.
- [CONTRIBUTING.md](CONTRIBUTING.md): patrones para extender el pipeline.
- [README.md](README.md): nueva sección de validación con 4 niveles.

### Resultados
- Catálogo a 0.05 AU: 4 035 700 encuentros (rebound) / 4 036 495 (Kepler).
- Validación: 4/4 Fienga + 4/4 Galád = 8/8 matched, MAE vs JPL = 0.00020 AU (ours) / 0.00004 AU (literatura).
- Outlier de Hygiea-Birkle (1.2 mAU vs literatura con Kepler): se cierra con rebound + cache-aware refinement → ~0.012 AU = JPL.
