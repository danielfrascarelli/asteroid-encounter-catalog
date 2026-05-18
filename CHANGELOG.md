# Changelog

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
