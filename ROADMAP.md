# Roadmap de implementación — Catálogo de Encuentros Cercanos (Gaia DR3)

> Planificación por hitos/fases. Cada fase tiene un objetivo concreto, entregables verificables y criterios de aceptación explícitos. El orden es estrictamente secuencial salvo donde se indica.

---

## Estado actual (2026-05-18)

**Fases 1–7 originales completas.** Pipeline end-to-end operativo: ingesta MPCORB + Gaia DR3 → propagación → detección KD-tree → caracterización → catálogo Parquet → dashboard Streamlit. Última corrida sobre **98.775 asteroides numerados (a∈[1.5, 4.0] AU)** en la ventana Gaia DR3 (2014-07-25 → 2017-05-28) a umbral **0.05 AU** produce **4.036.495 encuentros** en ~23 min (Kepler, 28 workers).

**Mejoras post-MVP implementadas:**

- **Multi-snapshot MPCORB indexado por época** ([src/ingest/mpcorb_archive.py](src/ingest/mpcorb_archive.py), [scripts/download_mpcorb_historical.py](scripts/download_mpcorb_historical.py)): descarga snapshots históricos vía Wayback Machine, auto-selecciona el de época más cercana al centro de la ventana temporal. Reduce error de Kepler de ~0.03 AU (época 2026 → 2017) a <0.001 AU (época 2015 → 2017).
- **Validación contra literatura**:
  - [scripts/validate_fienga_2003.py](scripts/validate_fienga_2003.py) — Fienga et al. (2003), VizieR `J/A+A/406/751`. 100% match (1/1 a 0.01 AU, 4/4 a 0.05 AU).
  - [scripts/validate_galad_2002.py](scripts/validate_galad_2002.py) — Galád & Gray (2002), parseado del HTML del paper (no está en VizieR). 100% match (4/4 a 0.05 AU).
  - [scripts/validate_jpl_horizons.py](scripts/validate_jpl_horizons.py) — Cross-check 3-way (nuestro vs literatura vs JPL Horizons). MAE(ours−JPL)=0.0002 AU, MAE(lit−JPL)=0.00004 AU.
- **Propagación N-body** ([src/propagate/nbody.py](src/propagate/nbody.py)): REBOUND con WHFast, Sol+Júpiter+Saturno como cuerpos masivos, asteroides como test particles. Opcional: Ceres/Vesta/Pallas/Hygiea como perturbadores (`config.propagation.rebound.include_major_asteroids`). Integra 100k asteroides × 25k pasos en ~9 min.
- **Cache de trayectorias** ([src/propagate/cache.py](src/propagate/cache.py)): persiste `(T, N, 3)` float32 en `np.memmap` (29.5 GB para 100k×25k), cache hit en <1 s. La integración streamea directo al disco para evitar OOM.

**Bugs conocidos:**
- `src/detect/parallel.py` se cuelga cuando: (a) `pairs` del prefilter es grande (subset chico ⇒ pickle de pares × 28 workers); (b) `positions=memmap` con N=100k+ a umbral 0.05 AU (sospecha: page-fault thrashing al leer el cache de 30 GB en 28 procesos). Workaround: `n_workers=1`. Fix pendiente.

**Para detalle completo de cada fase, ver secciones siguientes.**

---

## Fase 1 — Infraestructura y datos

**Objetivo**: Tener el entorno operativo, los datos en disco y las conversiones de tiempo correctas.

### Tareas

#### 1.1 Setup del entorno (Docker)

Todo el proyecto se ejecuta dentro de Docker. No se requiere Python local ni entorno virtual en el host.

**Archivos de infraestructura:**
- [ ] `Dockerfile` — imagen basada en `python:3.11-slim`, usa `uv` con `UV_SYSTEM_PYTHON=1` para instalar dependencias sin venv. Incluye gcc/g++ para extensiones compiladas (scipy, rebound).
- [ ] `docker-compose.yml` — tres servicios:
  - `pipeline`: runner principal, monta `./data` y `./logs` como volúmenes.
  - `dashboard`: corre Streamlit en `0.0.0.0:8501`, depende del servicio `pipeline`.
  - `test`: ejecuta `pytest tests/ -v`.
- [ ] `pyproject.toml` — declarar todas las dependencias del stack bajo `[project.dependencies]` y extras `[dev]` para herramientas de calidad. La imagen se reconstruye automáticamente cuando este archivo cambia.
- [ ] `ruff` + `black` + `mypy` configurados en `pyproject.toml` bajo `[tool.*]`.
- [ ] `.gitignore` que excluya `data/raw/`, `data/cache/`, `data/output/`, `logs/`, `*.local.yaml`, `.venv/`.
- [ ] `.dockerignore` que excluya `data/`, `logs/`, `*.local.yaml`, `.git/`, `__pycache__/`, `*.pyc`.

**Comandos de referencia (todos desde el host):**
```bash
docker compose build                                    # Construir la imagen
docker compose run --rm pipeline python -m scripts.download_mpcorb
docker compose run --rm pipeline python -m scripts.run_pipeline --config config.yaml
docker compose run --rm test                            # pytest
docker compose up dashboard                             # Streamlit en localhost:8501
```

**Verificación:**
- [ ] `docker compose build` termina sin errores.
- [ ] `docker compose run --rm pipeline python -c "import astropy, polars, scipy; print('OK')"` imprime `OK`.
- [ ] `docker compose run --rm test pytest tests/ --collect-only` lista los tests sin errores de importación.

#### 1.2 Módulo `src/utils/config.py`
- [ ] Leer `config.yaml` con `PyYAML` y exponer un objeto tipado (dataclass o Pydantic).
- [ ] Fallar rápido si faltan claves obligatorias.
- [ ] Soportar override con `config.local.yaml` (gitignored).

#### 1.3 Módulo `src/utils/time_utils.py`
- [ ] Función `utc_to_tdb(t: str | Time) -> Time` con escala explícita.
- [ ] Función `tcb_to_tdb(jd_tcb: float) -> float`.
- [ ] Tests unitarios en `tests/test_time_utils.py` con valores de referencia de `astropy`.
- [ ] Documentar en docstring la diferencia TCB vs TDB (~1.6 ms/día, ~1.8 s/año).

#### 1.4 Descarga de MPCORB (`src/ingest/mpcorb.py`)
- [ ] Script `scripts/download_mpcorb.py` que descarga desde `config.sources.mpcorb.url`.
- [ ] Registrar hash MD5 y fecha de descarga en `data/raw/mpcorb_metadata.json`.
- [ ] Parser que lee el formato de ancho fijo MPCORB.DAT y devuelve un `polars.DataFrame`.
- [ ] Columnas mínimas: `number`, `designation`, `H`, `G`, `a`, `e`, `i`, `Omega`, `omega`, `M`, `epoch_jd`.
- [ ] Filtrar por `subset.only_numbered` y `subset.semimajor_axis_au` según config.
- [ ] Test: verificar que (1) Ceres, (4) Vesta, (2) Pallas están presentes y sus `a` son correctos (±0.001 AU).

#### 1.5 Descarga de Gaia SSO (`src/ingest/gaia_sso.py`)
- [ ] Script `scripts/download_gaia_sso.py` que querea `gaiadr3.sso_observation` vía TAP.
- [ ] Guardar en `data/raw/gaia_sso.parquet`.
- [ ] Columnas mínimas según `config.sources.gaia_sso.columns`.
- [ ] Verificar que `epoch` está en TCB y documentarlo en el schema del archivo.
- [ ] Test básico: contar que hay observaciones para Ceres (`number_mp = 1`).

### Entregables
- `Dockerfile` + `docker-compose.yml` funcionales.
- `pyproject.toml` con dependencias completas.
- `data/raw/MPCORB.DAT` + metadatos de descarga.
- `data/raw/gaia_sso.parquet` con observaciones Gaia.
- `src/utils/time_utils.py` testeado.
- `src/utils/config.py` funcional.

### Criterios de aceptación
- `docker compose build` resuelve sin errores.
- `docker compose run --rm test pytest tests/test_time_utils.py` pasa.
- `docker compose run --rm test pytest tests/test_ingest.py` pasa: Ceres, Vesta, Pallas en MPCORB con `a` correctos.
- `docker compose run --rm pipeline ruff check .` sin errores.

---

## Fase 2 — Propagación orbital

**Objetivo**: Poder propagar cualquier asteroide del catálogo a una grilla temporal densa y verificar la precisión contra JPL Horizons.

### Tareas

#### 2.1 Propagador Kepler (`src/propagate/kepler.py`)
- [ ] Función `mean_to_eccentric_anomaly(M, e, tol=1e-12)` (iteración de Newton).
- [ ] Función `elements_to_state_vector(a, e, i, Omega, omega, E) -> (pos, vel)` en AU y AU/día.
- [ ] Función `propagate_kepler(elements, epoch_jd, target_jd_tdb_array) -> ndarray[N_times, 3]`.
- [ ] Frame de referencia: heliocéntrico eclíptico J2000 (lo que usa MPCORB); documentar la convención.
- [ ] Conversión heliocéntrico → baricéntrico: sumar posición del Sol desde efemérides (`astropy.coordinates.solar_system_ephemeris`).

#### 2.2 Grilla temporal (`src/propagate/grid.py`)
- [ ] Generar array de JD TDB desde `config.time_window.start` a `config.time_window.end` con paso `config.propagation.time_step_hours`.
- [ ] Función `propagate_all(df_elements, jd_grid) -> ndarray[N_asteroids, N_times, 3]`.
- [ ] Estrategia de memoria: no materializar el array completo si N_asteroids × N_times × 3 × 8 bytes > umbral configurable. Fallback a procesamiento en chunks.
- [ ] Cachear el resultado en `data/cache/positions_{hash}.npy` donde `hash` = hash de los parámetros de propagación + hash del subset de elementos.

#### 2.3 Validación contra JPL Horizons (`src/propagate/validate.py`)
- [ ] Para un conjunto de asteroides de prueba (Ceres, Vesta, Pallas + 2-3 MBAs aleatorios), queréar JPL Horizons con `astroquery.jplhorizons` en 5 épocas distribuidas por el período Gaia.
- [ ] Comparar posiciones baricéntricas propagadas vs Horizons.
- [ ] Aceptar error ≤ 0.001 AU (1000 km) para el Kepler de 2 cuerpos en el período 2014–2017.
- [ ] Loguear los errores por asteroide y época.

#### 2.4 Propagador rebound (opcional, post-MVP)
- [ ] Integración con `rebound` usando integrador `ias15`.
- [ ] Incluir planetas mayores como perturbadores (`config.propagation.rebound.include_planets`).
- [ ] Misma interfaz que el propagador Kepler para que sean intercambiables.
- [ ] Comparar resultados con Kepler: cuantificar mejora en error residual.

### Entregables
- `src/propagate/kepler.py` con funciones puras testeadas.
- `src/propagate/grid.py` con cache en disco.
- `src/propagate/validate.py` con reporte de errores contra Horizons.
- `tests/test_propagation.py` con casos para Ceres y al menos 2 MBAs.

### Criterios de aceptación
- Error de propagación Kepler vs Horizons ≤ 0.001 AU para los asteroides de prueba en 2014–2017.
- Cache funcional: segunda corrida no requiere recomputar si los parámetros no cambiaron.
- `pytest tests/test_propagation.py -v` pasa.

---

## Fase 3 — Detección de encuentros (MVP)

**Objetivo**: Detectar encuentros en un subset pequeño (~1000 asteroides) y verificar que el algoritmo es correcto antes de escalar.

### Tareas

#### 3.1 Filtro orbital previo (`src/detect/prefilter.py`)
- [ ] Función `compatible_pairs(df_elements) -> list[tuple[int, int]]` que devuelve pares con `|Δa| ≤ config.detection.prefilter.semimajor_diff_max_au` y `Δi ≤ config.detection.prefilter.inclination_diff_max_deg`.
- [ ] Implementar con operaciones vectorizadas en polars (no loops Python).
- [ ] Loguear cuántos pares surviven el filtro (fracción del total N²/2).

#### 3.2 KD-tree por step temporal (`src/detect/kdtree_scan.py`)
- [ ] Para cada step de la grilla, construir `scipy.spatial.cKDTree` con posiciones de los asteroides en ese step.
- [ ] Queréar vecinos dentro de `config.detection.threshold_au` para cada asteroide.
- [ ] Filtrar por los pares compatibles del prefilter.
- [ ] Acumular candidatos: lista de `(id_1, id_2, jd_step, dist_approx)`.

#### 3.3 Refinamiento temporal (`src/detect/refine.py`)
- [ ] Para cada candidato grueso, propagar los dos cuerpos en una sub-grilla densa centrada en `jd_step`, con paso `config.detection.refinement.fine_time_step_seconds` y ventana `config.detection.refinement.window_hours`.
- [ ] Encontrar el mínimo real de distancia por interpolación cuadrática sobre los 3 puntos más cercanos al mínimo de la sub-grilla.
- [ ] Confirmar o descartar: solo conservar si distancia mínima ≤ `threshold_au`.
- [ ] Registrar: `(id_1, id_2, jd_min_tdb, dist_min_au, relative_velocity_au_per_day)`.

#### 3.4 Pipeline de detección completo (`src/detect/pipeline.py`)
- [ ] Orquestar prefilter → KD-tree scan → refinamiento.
- [ ] Input: `df_elements`, `positions_array`, `jd_grid`.
- [ ] Output: `polars.DataFrame` con columnas del encuentro.
- [ ] Evitar duplicados: cada par `(id_1, id_2)` con `id_1 < id_2`.

#### 3.5 Test sobre subset de 1000 asteroides
- [ ] Correr pipeline con `config.subset.max_asteroids: 1000` (incluir Ceres y Vesta forzosamente).
- [ ] Verificar que el encuentro Ceres–Parthenope del config de validación aparece.
- [ ] Comparar tiempos de ejecución y uso de memoria.

### Entregables
- `src/detect/prefilter.py`, `kdtree_scan.py`, `refine.py`, `pipeline.py`.
- `tests/test_detection.py` con test de regresión para encuentros conocidos.
- Resultado en `data/output/encounters_test_1000.parquet`.

### Criterios de aceptación
- Encuentros de validación del config aparecen en el output del subset de 1000.
- No hay pares duplicados en el output.
- `pytest tests/test_detection.py::test_known_encounter` pasa.
- Pipeline de 1000 asteroides corre en < 10 minutos en hardware modesto (no optimizado aún).

### Cómo testear y validar la Fase 3

**Tests unitarios (sin datos reales, rápidos):**
```bash
docker compose run --rm -v ./src:/app/src -v ./tests:/app/tests \
  pipeline pytest tests/test_detection.py -v
```
Los 28 tests cubren prefilter, KD-tree scan, refinamiento cuadrático y el pipeline
completo usando asteroides sintéticos construidos analíticamente.  No requieren
`data/raw/` ni red.

**Validación sobre datos reales (subset de 1000 asteroides del cinturón principal):**
```bash
# Asume que MPCORB.DAT y gaia_sso.parquet ya están en data/raw/
docker compose run --rm pipeline python - <<'EOF'
import polars as pl
from src.ingest.mpcorb import parse_mpcorb
from src.propagate.grid import make_time_grid
from src.detect.pipeline import detect_encounters
from astropy.time import Time

elements = (
    parse_mpcorb("data/raw/MPCORB.DAT", semimajor_min_au=2.0, semimajor_max_au=3.5)
    .head(1000)
)
# Gaia DR3 window (TDB approximation)
t_start = Time("2014-07-25", scale="tdb").jd
t_end   = Time("2016-07-25", scale="tdb").jd
grid    = make_time_grid(t_start, t_end, step_hours=1.0)

encounters = detect_encounters(elements, grid, threshold_au=0.01)
print(encounters.sort("dist_au").head(20))
encounters.write_parquet("data/output/encounters_test_1000.parquet", compression="zstd")
EOF
```

**Qué verificar en el output:**
- Al menos un encuentro encontrado dentro del periodo Gaia (era de ~3 años es suficiente).
- `number_1 < number_2` para cada fila (garantizado por el prefilter de índices `triu`).
- `dist_au` ≤ 0.01 para todas las filas.
- `rel_vel_au_day` > 0 (velocidades relativas físicamente plausibles: 0.001–0.1 AU/día).
- Sin filas duplicadas: `len(encounters) == encounters.unique(["number_1", "number_2"]).len()`.

---

## Fase 4 — Paralelización y corrida completa

**Objetivo**: Escalar el pipeline al dataset completo (~100k asteroides numerados) con uso eficiente de CPU y memoria.

### Tareas

#### 4.1 Paralelización temporal (`src/detect/parallel.py`)
- [ ] Dividir la grilla temporal en chunks de `config.parallel.chunk_size_days` días.
- [ ] Distribuir chunks entre workers con `multiprocessing.Pool` (o `dask` según config).
- [ ] Cada worker procesa su chunk de forma independiente y retorna una lista de candidatos.
- [ ] Merge de candidatos al final, desduplicar, luego refinamiento.
- [ ] Manejar `n_workers: "auto"` usando `os.cpu_count()`.

#### 4.2 Gestión de memoria
- [ ] Benchmarkar uso de memoria para N = 10k, 50k, 100k asteroides × 3 años de grilla.
- [ ] Si el array de posiciones no cabe en RAM, implementar procesamiento por bloques de asteroides (chunked propagation).
- [ ] Documentar el compromiso memoria/velocidad en el README.

#### 4.3 Corrida completa
- [ ] Correr pipeline con todos los asteroides numerados del cinturón principal (`subset.only_numbered: true`, `semimajor_axis_au: [1.5, 4.0]`).
- [ ] Verificar que Ceres, Vesta, Pallas, Hygiea aparecen en el catálogo de encuentros (obligatorio).
- [ ] Registrar tiempo de ejecución total y peak de memoria.

#### 4.4 Logging y progreso
- [ ] `tqdm` en el loop principal de steps temporales.
- [ ] Log INFO al inicio de cada chunk (worker, rango de fechas, N candidatos encontrados).
- [ ] Guardar log en `logs/run_{name}.log`.

### Entregables
- `src/detect/parallel.py` funcional.
- Corrida completa con output en `data/output/encounters_full.parquet`.
- Reporte de performance (tiempo, memoria) en el log.

### Criterios de aceptación
- Ceres, Vesta, Pallas, Hygiea aparecen en el catálogo de encuentros.
- No hay condición de carrera entre workers (verificado con seed fijo y corridas repetidas).
- La corrida completa termina en tiempo razonable (< 4h en 8 cores).

---

## Fase 5 — Caracterización de encuentros

**Objetivo**: Enriquecer cada encuentro con propiedades físicas y observacionales.

### Tareas

#### 5.1 Geometría del encuentro (`src/characterize/geometry.py`)
- [ ] Velocidad relativa en el mínimo: `|v1 - v2|` en AU/día, convertir a m/s para output.
- [ ] Dirección del encuentro: ángulo entre el vector de posición relativa y el de velocidad relativa.
- [ ] Parámetro de impacto mínimo (si el encuentro fuera gravitacionalmente focalizado).

#### 5.2 Visibilidad desde Gaia (`src/characterize/observability.py`)
- [ ] Calcular la posición de Gaia en el momento del encuentro (efemérides Gaia de `astropy` o JPL).
- [ ] Ángulo de elongación solar del punto del encuentro desde Gaia.
- [ ] ¿Estaba el par dentro de la zona de exclusión solar de Gaia (< 45°)? Flag booleano.
- [ ] ¿Estaba el par en el rango de magnitud observable por Gaia (G < 21)?

#### 5.3 Propiedades físicas de los asteroides (`src/characterize/physical.py`)
- [ ] Diámetro estimado desde magnitud absoluta H y albedo: `D = (1329 / sqrt(p)) * 10^(-H/5)` km.
- [ ] Usar `config.characterize.default_albedo` como albedo por defecto.
- [ ] Clasificación aproximada: MBA, NEA, Troyana, Centauro (desde `a` y `e`).
- [ ] Masa estimada si el asteroide está en la lista de masas conocidas de JPL SBDB.

#### 5.4 Ensamblado del registro de encuentro (`src/characterize/encounter.py`)
- [ ] Schema final del registro: `(id_1, id_2, designation_1, designation_2, jd_tdb_min, date_utc_min, dist_min_au, dist_min_km, v_rel_au_per_day, v_rel_m_per_s, H_1, H_2, D1_km, D2_km, gaia_observable, solar_elongation_deg, class_1, class_2)`.
- [ ] Función `characterize_encounter(row, df_elements) -> dict`.
- [ ] Procesar todos los encuentros del DataFrame de detección.

### Entregables
- `src/characterize/` con los módulos anteriores.
- `tests/test_characterize.py` con casos unitarios para las fórmulas.
- DataFrame de encuentros caracterizados.

### Criterios de aceptación
- Diámetros de Ceres (940 km) y Vesta (525 km) calculados con < 5% de error.
- Velocidades relativas en el rango físico plausible para MBAs (< 10 km/s).
- `pytest tests/test_characterize.py` pasa.

---

## Fase 6 — Catálogo final

**Objetivo**: Construir la base de datos consultable con el catálogo completo y garantizar su reproducibilidad.

### Tareas

#### 6.1 Schema del catálogo (`src/catalog/schema.py`)
- [ ] Definir schema Parquet tipado (columnas, tipos, descripciones).
- [ ] Incluir columna `run_id` para identificar la corrida.

#### 6.2 Escritura del catálogo (`src/catalog/writer.py`)
- [ ] Escribir `data/output/encounters_catalog.parquet` con compresión `zstd`.
- [ ] Escribir sidecar `data/output/encounters_catalog_metadata.json` con:
  - Versión de MPCORB (hash + fecha de descarga).
  - Config completo usado.
  - Versiones de las dependencias clave (astropy, scipy, polars).
  - Timestamp de la corrida.
  - Número de encuentros totales.

#### 6.3 API de consulta (`src/catalog/query.py`)
- [ ] Función `load_catalog(path) -> polars.DataFrame`.
- [ ] Función `filter_encounters(df, min_dist_au, max_dist_au, date_start, date_end, body_ids)`.
- [ ] Función `top_encounters(df, n, by) -> polars.DataFrame` donde `by` puede ser `dist_min_au`, `v_rel`, etc.

#### 6.4 Tests de regresión finales (`tests/test_validation.py`)
- [ ] Test que Ceres–Parthenope, Vesta–Thetis, Pallas–Thalia aparecen en el catálogo (pares del `config.validation.known_pairs`).
- [ ] Test que el catálogo no tiene pares duplicados.
- [ ] Test que el sidecar de metadatos existe y es JSON válido.

### Entregables
- `data/output/encounters_catalog.parquet` (catálogo final).
- `data/output/encounters_catalog_metadata.json` (metadatos de reproducibilidad).
- `src/catalog/` completo.
- `pytest tests/test_validation.py` pasa.

### Criterios de aceptación
- Pares de validación presentes.
- Catálogo legible con `polars.read_parquet()` sin dependencias adicionales.
- Metadatos completos que permiten reproducir la corrida.

---

## Fase 7 — Dashboard y validación cruzada

**Objetivo**: Interfaz de exploración interactiva y validación del catálogo contra la literatura.

### Tareas

#### 7.1 Dashboard Streamlit (`src/dashboard/app.py`)
- [ ] Carga del catálogo al inicio; caché con `@st.cache_data`.
- [ ] **Filtros** en sidebar: rango de distancia mínima, rango de fechas, tipo de cuerpo, observable/no observable desde Gaia.
- [ ] **Tabla de encuentros** filtrada con columnas clave.
- [ ] **Histograma** de distribución de distancias mínimas.
- [ ] **Scatter** velocidad relativa vs distancia mínima.
- [ ] **Mapa del cinturón** (proyección a/e o a/i): puntos coloreados por distancia del encuentro.
- [ ] **Vista detalle** de un encuentro individual: trayectorias propagadas en ±24h, geometría.
- [ ] Botón de exportar CSV para el subset filtrado.

#### 7.2 Validación cruzada con literatura
- [ ] Cruzar catálogo generado con la lista de pares de Goffin (2014) para el período Gaia.
- [ ] Cruzar con Fuentes-Muñoz et al. (2024): verificar que los encuentros usados para determinación de masas aparecen.
- [ ] Documentar discrepancias y posibles causas (diferencia en threshold, método de propagación, etc.).
- [ ] Escribir breve análisis en `notebooks/validation_literature.ipynb`.

#### 7.3 Análisis exploratorio (notebooks)
- [ ] `notebooks/encounter_statistics.ipynb`: distribución de distancias, velocidades, épocas, clasificación de asteroides.
- [ ] `notebooks/big_encounters.ipynb`: encuentros más cercanos, con los cuerpos más grandes, entre NEAs.

#### 7.4 README final
- [ ] Actualizar README.md con resultados reales (número de encuentros encontrados, encuentro más cercano, etc.).
- [ ] Agregar sección de resultados preliminares con una figura.
- [ ] Verificar que todos los comandos del README funcionan desde cero.

### Entregables
- `src/dashboard/app.py` funcional.
- `notebooks/validation_literature.ipynb` con análisis.
- `notebooks/encounter_statistics.ipynb` y `notebooks/big_encounters.ipynb`.
- README actualizado con resultados reales.

### Criterios de aceptación
- `streamlit run src/dashboard/app.py` levanta sin errores.
- Encuentros de Fuentes-Muñoz et al. (2024) o Goffin (2014) identificados en el catálogo.
- README ejecutable de principio a fin desde un entorno limpio.

---

## Dependencias entre fases

```
Fase 1 (Infraestructura)
    └── Fase 2 (Propagación)
            └── Fase 3 (Detección MVP)
                    └── Fase 4 (Paralelización + escala)
                            └── Fase 5 (Caracterización)
                                    └── Fase 6 (Catálogo final)
                                                └── Fase 7 (Dashboard + validación)
```

---

## Gates de calidad entre fases

| Gate | Condición para avanzar a la siguiente fase |
|------|---------------------------------------------|
| 1 → 2 | Tests de ingest pasan; Ceres y Vesta en MPCORB con `a` correctos |
| 2 → 3 | Error de propagación Kepler ≤ 0.001 AU para asteroides de prueba |
| 3 → 4 | Encuentros de validación detectados en subset de 1000; sin duplicados |
| 4 → 5 | Ceres, Vesta, Pallas, Hygiea en catálogo; corrida completa terminada |
| 5 → 6 | Diámetros y velocidades en rangos físicos plausibles |
| 6 → 7 | `pytest tests/test_validation.py` pasa; catálogo legible sin dependencias |

---

## Parámetros clave de config.yaml por fase

| Fase | Parámetros relevantes |
|------|-----------------------|
| 1 | `sources.mpcorb`, `sources.gaia_sso`, `subset.only_numbered`, `subset.semimajor_axis_au` |
| 2 | `propagation.method`, `propagation.time_step_hours`, `time_window`, `propagation.cache_results` |
| 3 | `detection.threshold_au`, `detection.prefilter`, `detection.kdtree`, `detection.refinement` |
| 4 | `parallel.enabled`, `parallel.n_workers`, `parallel.backend`, `parallel.chunk_size_days` |
| 5 | `characterize.*` |
| 6 | `output.format`, `output.compression`, `output.include_metadata`, `validation.known_pairs` |
| 7 | — (solo lectura del catálogo) |
