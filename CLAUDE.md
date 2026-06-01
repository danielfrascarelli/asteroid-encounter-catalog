# CLAUDE.md

> Contexto del proyecto para asistentes de IA (Claude Code, Cursor, etc.).
> Este archivo orienta a cualquier agente que trabaje sobre este repositorio.

---

## Resumen del proyecto

**Nombre**: Catálogo de Encuentros Cercanos entre Asteroides (Gaia DR3)

**Objetivo**: Construir un pipeline que detecta sistemáticamente pares de asteroides que pasaron a menos de un umbral de distancia uno del otro durante el período de observación de la misión Gaia (julio 2014 – mayo 2017 en DR3), generando un catálogo nuevo de encuentros cercanos con su geometría y metadatos.

**Autor**: Estudiante de astronomía, Ingeniero de sistemas, perfil técnico Python avanzado.

**Duración estimada**: 5–7 semanas.

**Foco**: Programación de alto nivel (data engineering, algoritmos espaciales, paralelización) **y física tan rigurosa como haga falta para obtener resultados científicamente publicables**. No hay límite autoimpuesto de simplicidad física: cuando un resultado lo exija (p. ej. determinación de masas de asteroides), se implementa el modelado completo en vez de atajos — determinación de órbitas por mínimos cuadrados sobre el arco completo, ecuaciones variacionales para las parciales (∂obs/∂elementos y ∂obs/∂GM), ajuste simultáneo órbitas+masa, y modelo de errores astrométricos de Gaia (covarianza along-scan anisotrópica). Se parte de elementos orbitales públicos y librerías existentes (rebound, astropy), pero se las extiende —o se integran herramientas de astrodinámica de terceros— donde la ciencia lo requiera. El criterio de éxito es la **defendibilidad científica del resultado**, no la simplicidad de la implementación.

---

## Pregunta científica que responde

> Dados los ~150.000 asteroides observados por Gaia, ¿qué pares pasaron a menos de X distancia uno del otro durante el período de observación? ¿Cómo se distribuyen estos encuentros estadísticamente? ¿Hay encuentros llamativos entre asteroides grandes, NEAs, o miembros de la misma familia dinámica?

---

## Stack tecnológico

### Core
- **Python** ≥ 3.11
- **astropy**: tiempo, coordenadas, unidades
- **astroquery**: queries a JPL Horizons y MPC
- **numpy / scipy**: cálculo numérico, KD-trees (`scipy.spatial`)
- **polars** / **duckdb**: manejo de datos a escala (preferidos sobre pandas)
- **rebound** o **poliastro**: propagación orbital (ver sección de algoritmos)

### Visualización
- **plotly**: gráficos interactivos
- **streamlit**: dashboard final de exploración

### Infraestructura
- **multiprocessing** / **dask**: paralelización
- **rich** / **tqdm**: logging y progreso
- **pytest**: testing
- **uv** o **poetry**: gestión de entorno y dependencias

---

## Fuentes de datos

| Dato | Fuente | Formato | Notas |
|------|--------|---------|-------|
| Observaciones Gaia SSO | Gaia Archive (`gea.esac.esa.int`) | CSV / parquet | Tabla `gaiadr3.sso_observation` |
| Elementos orbitales (bulk) | MPC (`minorplanetcenter.net`) | MPCORB.DAT | Descarga única, offline |
| Validación de órbitas puntuales | JPL Horizons (`astroquery.jplhorizons`) | API online | Solo casos específicos, rate-limited |
| Masas conocidas (validación) | JPL SBDB | API online | Para identificar perturbers conocidos |

Los datos crudos viven en `data/raw/` y **no** se versionan (ver `.gitignore`).

---

## Decisiones de diseño clave

### 1. Umbral de encuentro: 0.05 AU por defecto
Compromiso entre permisividad estadística (capturar suficientes eventos) y relevancia física (encuentros usables para futuras determinaciones de masa). Configurable vía `config.yaml`.

### 2. Propagación orbital: Kepler analítico → rebound
- **MVP**: Kepler de dos cuerpos puro (Sol + asteroide). Suficiente para detección a 0.01 AU en escalas de pocos años.
- **Refinamiento opcional**: `rebound` con planetas mayores para reducir falsos positivos cerca de resonancias.

### 3. Subset inicial: asteroides numerados
Empezar solo con asteroides con número MPC asignado (mejor calidad orbital). Expandir a designaciones provisionales después.

### 4. Algoritmo de detección: KD-tree por step temporal
- Discretizar tiempo en pasos de Δt (típicamente 1 hora).
- En cada step, construir KD-tree 3D con posiciones de todos los asteroides.
- Query de vecinos para cada asteroide; filtrar por umbral.
- Refinar candidatos con sub-grid temporal más fino alrededor del mínimo aparente.

### 5. Tiempos: siempre JD TDB internamente
Toda la lógica interna usa Julian Date en escala TDB (Barycentric Dynamical Time). Conversiones a/desde TCB (lo que Gaia reporta) y UTC se hacen **solo** en interfaces de entrada/salida. **Nunca mezclar floats de JD sin saber su escala**.

---

## Estructura del repositorio

```
.
├── CLAUDE.md                    # Este archivo
├── README.md                    # Descripción para humanos
├── pyproject.toml
├── config.yaml                  # Parámetros del pipeline
├── data/
│   ├── raw/                     # Datos descargados (gitignored)
│   ├── cache/                   # Computaciones intermedias (gitignored)
│   └── output/                  # Catálogo final de encuentros
├── src/
│   ├── ingest/                  # Descarga y parseo de fuentes externas
│   ├── propagate/               # Propagación orbital
│   ├── detect/                  # Algoritmo de detección de encuentros
│   ├── characterize/            # Cálculo de propiedades de encuentros
│   ├── catalog/                 # Construcción y consulta del catálogo final
│   ├── dashboard/               # Streamlit app
│   └── utils/                   # Tiempo, geometría, IO
├── notebooks/                   # Exploración interactiva (NO pipeline)
├── tests/
└── scripts/                     # Entry points CLI
```

---

## Convenciones

### Unidades
- **Distancias**: AU en todo el código interno. Convertir a km solo para output user-facing.
- **Tiempos**: JD TDB internamente. Strings ISO solo para input/output.
- **Ángulos**: radianes internamente. Grados solo para input/output.
- **Velocidades**: AU/día internamente.
- Usar `astropy.units` cuando haya ambigüedad o en código no-hot-path. Evitar en bucles internos por overhead.

### Estilo
- Type hints obligatorios en funciones públicas.
- Docstrings estilo NumPy.
- Funciones puras donde sea posible.
- Side effects (IO, logging) explícitos.
- `ruff` + `black` para formato.

### Logging
- `logging` estándar (no `print`).
- Nivel INFO para hitos del pipeline, DEBUG para detalle.
- Progreso de loops largos con `tqdm`.

### Reproducibilidad
- Seed fijo (`config.yaml`) para cualquier randomness.
- Datos intermedios cacheados con hash del input + parámetros.
- Versión / fecha de descarga de MPCORB registrada en el output.

---

## Algoritmos: descripción de alto nivel

### Detección eficiente de encuentros

**Problema**: 150.000 asteroides × N pasos temporales ⇒ ~22.500 millones de pares por step → infactible naive.

**Solución en capas**:
1. **Filtro orbital previo**: descartar pares con semiejes mayores muy distintos (no pueden encontrarse físicamente). Reduce el espacio de búsqueda en órdenes de magnitud.
2. **KD-tree espacial por step**: O(N log N) construcción, O(log N) por query.
3. **Refinamiento temporal**: para cada candidato grueso, sub-grid temporal denso alrededor del mínimo aparente.
4. **Paralelización**: bloques temporales son independientes → `multiprocessing.Pool` o `dask`.

### Caracterización de encuentros

Para cada par detectado, calcular:
- Distancia mínima y época exacta (interpolación cuadrática en el mínimo).
- Velocidad relativa en el encuentro.
- Geometría respecto a Gaia (¿observable? ¿ángulo de fase favorable?).
- Magnitudes y diámetros estimados de ambos cuerpos.

---

## Pitfalls / gotchas específicos del dominio

⚠️ **Escalas de tiempo**: TCB (lo que reporta Gaia) ≠ TDB (lo que usan efemérides JPL) ≠ UTC. Diferencias del orden de segundos, pero importantes a precisión sub-AU. Siempre convertir explícitamente con `astropy.time.Time`.

⚠️ **Frame de referencia**: posiciones baricéntricas (centro de masas del sistema solar) vs heliocéntricas (centro del Sol). Diferencia ~0.005 AU. MPCORB usa heliocéntrico, Gaia reporta baricéntrico.

⚠️ **Light-time correction**: la posición "observada" por Gaia es la del asteroide cuando emitió la luz, no cuando Gaia la recibió. Para detección de encuentros entre asteroides (problema geométrico puro) **no aplica** esta corrección, pero es fácil olvidarlo si se compara con posiciones observacionales.

⚠️ **Época de los elementos orbitales**: MPCORB da elementos en una época fija. Propagar muchos años desde esa época con Kepler puro acumula error. Para Gaia (período corto, ~3 años) está bien, pero verificar contra JPL Horizons.

⚠️ **Memoria**: una matriz N×N de distancias para 150k asteroides son ~90 GB en float32. **No materializarla**. Trabajar con KD-trees y queries puntuales.

⚠️ **MPCORB se actualiza**: registrar el hash/fecha de descarga en cada corrida. Resultados no son bit-exactly reproducibles entre versiones del catálogo.

⚠️ **Encuentros vs co-localizaciones aparentes**: dos asteroides pueden aparecer cercanos en RA/Dec desde la Tierra pero estar a AU de distancia en 3D. Este proyecto detecta encuentros **3D reales**, no apariencias en el plano del cielo.

---

## Glosario de términos del dominio

| Término | Significado |
|---------|-------------|
| **SSO** | Solar System Object |
| **NEA** | Near-Earth Asteroid |
| **MBA** | Main Belt Asteroid |
| **MPC** | Minor Planet Center |
| **MPCORB** | Archivo bulk de elementos orbitales del MPC |
| **Elementos keplerianos** | a (semieje mayor), e (excentricidad), i (inclinación), Ω (nodo ascendente), ω (argumento del perihelio), M (anomalía media) en una época dada |
| **AU** | Unidad Astronómica (~1.5×10⁸ km) |
| **TCB / TDB / UTC** | Escalas de tiempo (Barycentric Coordinate Time / Barycentric Dynamical Time / Coordinated Universal Time) |
| **Tránsito (transit)** | Paso de un objeto por el plano focal de Gaia (~40 s) |
| **Perturber / test asteroid** | En estudios de masas, el cuerpo grande que perturba y el pequeño que es perturbado |

---

## Comandos comunes

Todo se ejecuta dentro de Docker. No se requiere Python local.

```bash
# Build de la imagen (requerido tras cambios en pyproject.toml)
docker compose build

# Descarga inicial de datos
docker compose run --rm pipeline python -m scripts.ingest.download_mpcorb
docker compose run --rm pipeline python -m scripts.ingest.download_gaia_sso

# Pipeline completo
docker compose run --rm pipeline python -m scripts.pipeline.run_pipeline --config config.yaml

# Solo detección (asume propagación ya hecha)
docker compose run --rm pipeline python -m scripts.pipeline.detect_deflections

# Fitting de masas
docker compose run --rm pipeline python -m scripts.mass.fit_mass_gaia_loo
docker compose run --rm pipeline python -m scripts.mass.summarize_loo_fits

# Validación contra literatura
docker compose run --rm pipeline python -m scripts.validate.validate_goffin_2014
docker compose run --rm pipeline python -m scripts.validate.validate_literature

# Dashboard
docker compose up dashboard      # abre http://localhost:8501

# Tests
docker compose run --rm test
docker compose run --rm test pytest tests/test_detection.py::test_known_encounter

# Formato y lint
docker compose run --rm pipeline ruff check . --fix
docker compose run --rm pipeline black .

# Shell interactivo dentro del contenedor
docker compose run --rm pipeline bash
```

### Estructura de scripts

```
scripts/
├── ingest/      # descarga de datos externos (MPCORB, Gaia SSO, literatura)
├── pipeline/    # ejecución del pipeline (detección, caracterización, masas)
├── mass/        # fitting y análisis de masas (LOO, linear, perturbers)
├── validate/    # validación contra literatura y JPL Horizons
├── bench/       # benchmarks de rendimiento y experimentos
└── dev/         # sanity checks, demos y tests de desarrollo
```

---

## Hitos / roadmap

Ver [ROADMAP.md](ROADMAP.md) para el detalle completo de cada fase (tareas, entregables, criterios de aceptación).

| # | Fase | Descripción | Gate de salida |
|---|------|-------------|----------------|
| 1 | **Infraestructura y datos** | Setup del entorno, descarga de MPCORB y Gaia SSO, módulo de conversión de tiempos | Ceres y Vesta en MPCORB con `a` correctos |
| 2 | **Propagación orbital** | Propagador Kepler 2-cuerpos, grilla temporal con cache, validación contra JPL Horizons | Error ≤ 0.001 AU para asteroides de prueba |
| 3 | **Detección MVP** | Filtro orbital previo, KD-tree por step temporal, refinamiento sobre subset de ~1000 asteroides | Encuentros de validación detectados; sin duplicados |
| 4 | **Escala y paralelización** | Paralelización por bloques temporales, corrida completa sobre todos los numerados | Ceres, Vesta, Pallas, Hygiea en el catálogo |
| 5 | **Caracterización** | Geometría del encuentro, visibilidad desde Gaia, propiedades físicas estimadas | Diámetros y velocidades en rangos físicos plausibles |
| 6 | **Catálogo final** | Schema Parquet tipado, sidecar de metadatos, API de consulta, tests de regresión | `pytest tests/test_validation.py` pasa |
| 7 | **Dashboard y validación** | Streamlit app, cruce con Goffin (2014) y Fuentes-Muñoz (2024), README final | Dashboard levanta; encuentros de literatura identificados |

---

## Validación

Encuentros conocidos en la literatura (para tests de regresión):
- Pares reportados en Goffin (2014) — usados para determinación de masas.
- Encuentros con (1) Ceres, (4) Vesta, (2) Pallas, (10) Hygiea en el período Gaia (deben aparecer **obligatoriamente** en el catálogo, dado el tamaño de su esfera de Hill).

---

## Referencias clave

- **Tanga et al. (2023)**: "Gaia Data Release 3. The Solar System survey", A&A 674, A12.
- **Goffin (2014)**: "Asteroid mass determinations", A&A 565, A56.
- **Fuentes-Muñoz et al. (2024)**: "231 asteroid masses from Gaia FPR", LPSC #2388.
- **Documentación Gaia DR3**: `gea.esac.esa.int/archive/documentation`.
- **MPCORB format**: `minorplanetcenter.net/iau/info/MPOrbitFormat.html`.

---

## Notas para el agente de IA

- Cuando se discutan tiempos, **siempre preguntar/declarar la escala** (TCB/TDB/UTC/JD/MJD). No asumir.
- Frente a operaciones N² obvias, **sugerir indexación o filtrado previo** antes de proponer paralelización bruta.
- Si se sugiere usar pandas, **proponer polars/duckdb** como alternativa más performante en este contexto.
- No introducir nuevas dependencias sin justificación; el stack está pensado para ser estable.
- Los notebooks son para exploración, **no** para pipeline. Lógica reutilizable debe migrar a `src/`.
- Los datos crudos de MPCORB y Gaia pueden ser grandes (GB). Antes de cargarlos enteros, considerar streaming o queries SQL.
- Cuando se ajusten parámetros del pipeline, hacerlo en `config.yaml`, **no** hardcodeado.

---

## Convenciones de colaboración con el agente

### Git
- **Autor de todos los commits**: `Daniel Frascarelli <dsanfra@gmail.com>` — cuenta GitHub: `danielfrascarelli`.
- **No agregar `Co-Authored-By`** en ningún commit. El agente no debe aparecer como co-autor.
- Antes de commitear, verificar que author y committer tengan el email correcto (`dsanfra@gmail.com`). Si el git config local usa otro email, sobreescribir con `GIT_COMMITTER_EMAIL` al hacer el commit.
- No usar `--no-verify` ni saltar hooks salvo indicación explícita.

### Flujo de trabajo con Git
- **Nunca trabajar en `main`**. Ante cualquier cambio, crear una branch descriptiva (ej. `feat/propagation`, `fix/kdtree-duplicates`, `docs/roadmap-update`).
- **No commitear ni pushear automáticamente** tras editar archivos. Presentar los cambios al usuario y esperar confirmación explícita antes de hacer `git commit`.
- Una vez confirmado el commit, abrir un PR hacia `main` en lugar de mergear directamente.

### Estilo de trabajo
- Respuestas cortas y directas. Sin resumen al final de cada respuesta ("lo que hice fue...").
- Operar sin pausas para preguntas de clarificación; tomar la decisión razonable y continuar.
- Todo el proyecto corre dentro de Docker. No sugerir ni usar entornos Python locales.
