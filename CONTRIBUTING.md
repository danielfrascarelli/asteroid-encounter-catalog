# Contributing — Catálogo de Encuentros Cercanos (Gaia DR3)

> Patrones para extender el pipeline. Si vas a agregar código nuevo, ojalá esto te ahorre tiempo.

---

## Reglas que el proyecto sigue

- Todo corre dentro de **Docker** (`docker compose run --rm pipeline ...`). No se asume Python local.
- Internamente los **tiempos son JD TDB** y las **distancias son AU**. Las conversiones a/desde UTC, TCB, km, deg sólo ocurren en interfaces de entrada/salida.
- **Polars** se prefiere sobre pandas. **Logging** (no `print`). **Type hints** en funciones públicas.
- El stack está fijado en `pyproject.toml`; no agregar dependencias nuevas sin justificación.
- `ruff` + `black` + `mypy` deben pasar en CI (`docker compose run --rm pipeline sh -c "ruff check . && black --check . && mypy src scripts"`).
- Cualquier cambio debe pasar los tests: `docker compose run --rm test pytest tests/ -v`.

## Cómo extender el pipeline

### Agregar una fuente de datos externa nueva

Caso de uso típico: tenés un catálogo (VizieR, HTML, FTP) que querés cruzar contra el catálogo del pipeline.

1. Agregá una entrada bajo `sources:` en [config.yaml](config.yaml):

   ```yaml
   sources:
     mi_fuente:
       url: "https://ejemplo.org/data.csv"
       output_filename: "mi_fuente.parquet"
   ```

2. Agregá un dataclass `MiFuenteSourceConfig` en [src/utils/config.py](src/utils/config.py) y registrá en `SourcesConfig` + `_require` + `_build`. Seguí el patrón de `Fienga2003SourceConfig` o `Galad2002SourceConfig`.

3. Creá `scripts/download_mi_fuente.py` que lea `cfg.sources.mi_fuente`, baje, parsee y escriba a `data/raw/<output_filename>` + un sidecar `*_metadata.json`.

4. (Opcional) Si la fuente tiene encuentros para cruzar, creá `scripts/validate_mi_fuente.py` siguiendo el patrón de `scripts/validate_fienga_2003.py`:

   - Filtrá a la ventana Gaia (`cfg.time_window`).
   - Particioná por `Impact <= cfg.detection.threshold_au` (esperados vs flybys más anchos).
   - Cross-match `(perturber, target)` ± tolerancia temporal.
   - Escribí `data/output/<name>_matches.csv` y `_misses.csv`.

### Agregar un nuevo método de propagación

1. Implementá un módulo en [src/propagate/](src/propagate/) que exponga `propagate_grid(elements, time_grid) -> ndarray (T, N, 3)` o un iterator `(t, pos)`.
2. Registralo en [src/propagate/grid.py:propagate_full_grid](src/propagate/grid.py) con un branch sobre `cfg.propagation.method`.
3. Si la propagación es cara (e.g. N-body), aceptá un buffer de salida `out: np.ndarray | None` y enrutalo desde [src/propagate/cache.py](src/propagate/cache.py) para streaming a disco (evita OOM en runs grandes).
4. Cache key: agregá un campo a `build_cache_key` en `cache.py` con los parámetros distintivos del método.

### Agregar un snapshot MPCORB histórico

```bash
docker compose run --rm pipeline python -m scripts.download_mpcorb_historical --year 2018 --month 6
```

Esto descarga del Wayback Machine y deja el archivo en `data/raw/mpcorb_archive/MPCORB_<YYYYMMDD>.DAT` con su sidecar JSON. El pipeline auto-selecciona el snapshot con época más cercana al centro de la ventana `cfg.time_window`.

### Agregar un nuevo umbral / configuración temporal

Usá `config.local.yaml` (gitignored) para overrides locales sin tocar `config.yaml`:

```yaml
# config.local.yaml
detection:
  threshold_au: 0.05
output:
  filename: "encounters_catalog_005au"
```

Se mergea automáticamente sobre `config.yaml`.

### Agregar tests

Los tests viven en [tests/](tests/) y usan pytest. Mock las llamadas a red. Para tests que requieran datos reales chicos, fabricá los archivos en `tmp_path` (ver `tests/test_mpcorb_archive.py` como ejemplo).

Correr:

```bash
docker compose run --rm test pytest tests/test_mi_modulo.py -v
```

## Convenciones de git

- **Branch descriptiva** (`feat/...`, `fix/...`, `docs/...`). Nunca en `main`.
- **Autor**: `Daniel Frascarelli <dsanfra@gmail.com>` (todos los commits). No agregar `Co-Authored-By: Claude`.
- **No commitear automáticamente**: presentar diff al usuario y esperar OK explícito.
- **No saltar hooks** (`--no-verify`).
- Una vez confirmado, abrir PR contra `main` en lugar de mergear directo.

## Pitfalls específicos del dominio

- **Escalas de tiempo**: TCB (Gaia) ≠ TDB (JPL) ≠ UTC. Diferencia ~segundos, pero importante a 0.001 AU. Convertir explícitamente con `astropy.time.Time(..., scale=...)`.
- **Frame**: MPCORB usa heliocéntrico; Gaia reporta baricéntrico. Diferencia ~0.005 AU.
- **Light-time correction**: el proyecto detecta encuentros 3D geométricos, NO aparentes desde Gaia. No aplicar light-time aquí.
- **Época de elementos orbitales**: Kepler 2-cuerpos puro propagado >5 años acumula error >0.01 AU. Usar el snapshot histórico apropiado (`scripts/download_mpcorb_historical.py`) o cambiar a propagación N-body (`config.propagation.method: rebound`).
- **Memoria**: matriz `(T, N, 3)` para 100k asteroides × 25k pasos × float32 = 30 GB. Usar `np.memmap` (lo hace `src/propagate/cache.py` automáticamente).

## Layout reciente del repo

```
src/
├── detect/          # KD-tree scan + refinement + parallel
├── ingest/          # MPCORB parser + multi-snapshot archive + Gaia
├── propagate/       # Kepler 2-body + rebound N-body + cache
├── characterize/    # Velocidades, magnitudes, observabilidad Gaia
├── catalog/         # Schema, writer, query API
├── dashboard/       # Streamlit (pendiente)
└── utils/           # Config tipado, tiempos
scripts/
├── download_*.py    # MPC, Gaia, Fienga, Galád, MPCORB histórico
├── validate_*.py    # Cross-match contra catálogos + JPL Horizons
├── run_pipeline.py  # Entry point
└── compare_kepler_vs_rebound.py
tests/               # 194 tests, todos pasan
```

## Para más contexto

- [CLAUDE.md](CLAUDE.md) — contexto para asistentes IA.
- [ROADMAP.md](ROADMAP.md) — fases del plan (1–7) + mejoras post-MVP.
- [VALIDATION_SUMMARY.md](VALIDATION_SUMMARY.md) — resultados de validación al día.
