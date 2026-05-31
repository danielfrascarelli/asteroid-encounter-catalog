# 🌌 Catálogo de Encuentros Cercanos entre Asteroides

> Detección sistemática de pares de asteroides que pasaron cerca uno del otro durante el período de observación de la misión Gaia (DR3), construyendo un catálogo nuevo de encuentros cercanos con su geometría y propiedades físicas.

![Python](https://img.shields.io/badge/python-3.11+-blue.svg)
![Status](https://img.shields.io/badge/status-completado-brightgreen.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)

---

## 📖 Sobre el proyecto

La misión **Gaia** de la ESA observó cientos de miles de asteroides con precisión astrométrica sub-miliarcosegundo entre 2014 y 2017. Esos datos, publicados en **Gaia DR3**, abren una ventana única para estudiar la dinámica del cinturón principal: por primera vez podemos rastrear con altísima precisión el movimiento de ~150.000 asteroides simultáneamente.

Este proyecto pregunta algo simple pero técnicamente desafiante:

> *De todos los asteroides observados por Gaia, ¿cuáles pasaron cerca unos de otros (en 3D, no aparentemente) durante el período de la misión?*

La respuesta es valiosa por varias razones:

- 🔭 **Construye un catálogo nuevo** de encuentros cercanos validable contra la literatura.
- 📊 **Permite análisis estadísticos** sobre la frecuencia y distribución de estos eventos.
- 🪨 **Identifica candidatos** para futuros estudios de determinación de masas (los encuentros son la herramienta clásica para "pesar" asteroides grandes).
- 💻 **Es un excelente ejercicio** de ingeniería de software aplicada a un dataset astronómico real a escala.

## ✨ Qué hace el pipeline

1. **Descarga** elementos orbitales del MPC y observaciones de Gaia DR3.
2. **Propaga** las órbitas de ~150.000 asteroides en una grilla temporal densa.
3. **Detecta** eficientemente todos los pares cuya separación 3D bajó de un umbral configurable.
4. **Refina** cada candidato con sub-muestreo temporal alrededor del mínimo de distancia.
5. **Caracteriza** cada encuentro: distancia mínima, época exacta, velocidad relativa, geometría observacional, magnitudes.
6. **Cataloga** los resultados en una base de datos consultable.
7. **Visualiza** todo en un dashboard interactivo donde puedes filtrar, explorar y exportar.

## 🚀 Inicio rápido

### Requisitos previos

- Docker + Docker Compose
- ~10 GB de espacio en disco (datos crudos + cache)
- 8 GB de RAM mínimo (16 GB recomendado para subset completo)

No se requiere Python local — todo corre dentro del contenedor.

### Instalación

```bash
git clone https://github.com/tu-usuario/asteroid-encounters.git
cd asteroid-encounters
docker compose build
```

### Descarga de datos

```bash
# 1. MPCORB — elementos orbitales (actuales, ~90 MB comprimidos)
docker compose run --rm pipeline python -m scripts.ingest.download_mpcorb

# 2. Snapshot histórico 2015 (recomendado: reduce error Kepler de ~30 mAU a <1 mAU)
# El pipeline auto-selecciona el snapshot con época más cercana al centro de la ventana.
docker compose run --rm pipeline python -m scripts.ingest.download_mpcorb_historical --year 2015 --month 6
```

> **Nota**: `download_gaia_sso` no es requerido para el pipeline de detección. El pipeline lee únicamente MPCORB.

### Pipeline completo

```bash
# Detección + caracterización completa (todos los asteroides numerados)
docker compose run --rm pipeline python -m scripts.pipeline.run_pipeline
docker compose run --rm pipeline python -m scripts.pipeline.characterize_catalog
```

Para caracterizar el catálogo congelado completo (72 M filas) sin OOM, usar el
modo streaming (chunked, RAM acotada por chunk):

```bash
docker compose run --rm pipeline python -m scripts.pipeline.characterize_catalog \
    --input data/output/encounters_catalog_hybrid_stageb.parquet --streaming on
# → data/output/encounters_characterized_full.parquet (72.236.904 filas, 18,9% Gaia-observables)
```

Produce el catálogo caracterizado y su sidecar de metadatos. Para los counts y
hashes de la corrida congelada, ver [`FROZEN_RUN.md`](FROZEN_RUN.md).

### Subset rápido para pruebas

```bash
# Configurar para corrida rápida (~5000 asteroides, ~10 minutos)
cp config.yaml config.local.yaml
# editar config.local.yaml → subset.max_asteroids: 5000

docker compose run --rm pipeline python -m scripts.pipeline.run_pipeline --config config.local.yaml
docker compose run --rm pipeline python -m scripts.pipeline.characterize_catalog --config config.local.yaml
```

### Explorar resultados

```bash
docker compose up dashboard
```

Y abrir `http://localhost:8501` en el navegador.

## 🔧 Configuración

Todos los parámetros del pipeline están en `config.yaml`. Los más importantes:

| Parámetro | Default | Descripción |
|-----------|---------|-------------|
| `detection.threshold_au` | `0.05` | Umbral de distancia para considerar "encuentro" |
| `detection.time_step_hours` | `1.0` | Resolución temporal de la grilla de propagación |
| `subset.only_numbered` | `true` | Limitarse a asteroides con número MPC asignado |
| `subset.max_asteroids` | `null` | Tope opcional para tests rápidos |
| `parallel.n_workers` | `auto` | Procesos paralelos (auto = número de cores) |
| `propagation.method` | `kepler` | `kepler` (rápido) o `rebound` (preciso) |

Ver `config.yaml` para la lista completa con comentarios.

## 🏗️ Arquitectura del pipeline

```
┌────────────────┐    ┌────────────────┐    ┌─────────────────┐
│   MPC / JPL    │    │  Gaia Archive  │    │   JPL Horizons  │
│   (orbital     │    │  (epochs +     │    │   (validación   │
│   elements)    │    │   asteroid     │    │   puntual)      │
│                │    │   list)        │    │                 │
└────────┬───────┘    └────────┬───────┘    └────────┬────────┘
         │                     │                     │
         ▼                     ▼                     ▼
   ┌─────────────────────────────────────────────────────┐
   │              src/ingest (descarga y parseo)         │
   └────────────────────────┬────────────────────────────┘
                            ▼
   ┌─────────────────────────────────────────────────────┐
   │       src/propagate (órbitas → grilla 3D temporal)  │
   └────────────────────────┬────────────────────────────┘
                            ▼
   ┌─────────────────────────────────────────────────────┐
   │     src/detect (KD-tree → candidatos → refinamiento)│
   └────────────────────────┬────────────────────────────┘
                            ▼
   ┌─────────────────────────────────────────────────────┐
   │   src/characterize (geometría, velocidad, contexto) │
   └────────────────────────┬────────────────────────────┘
                            ▼
   ┌─────────────────────────────────────────────────────┐
   │   src/catalog (DuckDB) ──► src/dashboard (Streamlit)│
   └─────────────────────────────────────────────────────┘
```

## 📁 Estructura del repositorio

```
.
├── CLAUDE.md                # Contexto para asistentes de IA
├── README.md                # Este archivo
├── config.yaml              # Parámetros del pipeline
├── data/                    # Datos (gitignored)
├── src/                     # Código fuente modular
├── notebooks/               # Exploración interactiva
├── tests/                   # Suite de tests
└── scripts/                 # Entry points CLI
```

Ver `CLAUDE.md` para una descripción detallada de cada módulo.

## 📊 Resultados (corrida congelada)

> **Alcance**: la corrida congelada produce un **catálogo de candidatos
> geométricos bajo supuestos congelados** (rebound coarse-scan + refinamiento
> Kepler 2-cuerpos, prefilter heurístico). NO es un catálogo completo, NO es
> un catálogo de masas, y NO está validado a precisión sub-cadencia.  Leer
> [`FROZEN_RUN.md` § "Scope and limits"](FROZEN_RUN.md) antes de citar
> cualquier número de abajo.

| Métrica | Valor |
|---------|-------|
| MPCORB snapshot | `MPCORB_20160217.DAT` |
| Ventana temporal | 2014-07-25 → 2017-05-28 (Gaia DR3) |
| Umbral de detección | **0.05 AU** |
| Scan / refine | rebound (whfast, Sun+Jupiter+Saturn) / **Kepler 2-cuerpos** |
| **Candidatos geométricos (Kepler-refined)** | **72.236.904** |
| Separación mínima Kepler-refined | **6.6 × 10⁻⁶ AU** (≈ 988 km, bajo el modelo Kepler) |
| Candidatos < 0.001 AU | 26.038 |
| Candidatos < 0.01 AU | 2.833.425 |

Hashes de inputs/outputs y la tabla completa de claims en
[`FROZEN_RUN.md`](FROZEN_RUN.md).

### Top-3 separaciones mínimas (bajo el refinamiento Kepler 2-cuerpos)

| Rank | Cuerpo 1 | Cuerpo 2 | Distancia Kepler (AU) |
|------|----------|----------|-----------------------|
| 1 | (153222) 2000 YD43 | (238587) 2004 YX3 | 6.6 × 10⁻⁶ |
| 2 | (15072) Landolt | (387599) 2001 XF180 | 1.2 × 10⁻⁵ |
| 3 | (270730) 2002 QE130 | (366918) 2005 UC211 | 1.5 × 10⁻⁵ |

Estos son los pares con menor separación bajo el modelo de refinamiento
Kepler. La distancia mínima real bajo dinámica N-body completa puede
diferir, especialmente para órbitas con alta excentricidad o cerca de
resonancias con Júpiter.

### Notas honestas sobre alcance

- **Catálogo geométrico de candidatos**: defensible como tal. Los 72 M pares
  son los que pasan el filtro bajo la pipeline congelada. NO se puede afirmar
  completitud ni precisión μAU global (validación limitada a 8 pares + grilla
  1 h de JPL).
- **Completitud del prefiltro (audit #2): cuantificada.** El prefiltro orbital
  (|Δa|≤0.5 AU ∧ |Δi|≤30°) tiene **76,4 % de recall en la cola adversa**
  (alta-e/alta-i): pierde ~143 k encuentros reales, casi todos por el corte
  |Δa|≤0.5 (ciego a la excentricidad). Fix recomendado y verificado: prefiltro
  de solapamiento radial → 100 % recall. Detalle en
  [`docs/prefilter_recall.md`](docs/prefilter_recall.md). No usar la palabra
  "completo" sobre el catálogo congelado.
- **Capa de masas**: NO publicable. El archivo
  `mass_followup_candidates.csv` (41 filas) son *targets de seguimiento*,
  no medidas. El test de specificidad da 0/41
  ([`encounter_analysis/DETECTIONS.md`](encounter_analysis/DETECTIONS.md))
  y los chi²_red del LOO batch dan mediana 425 / max 7.2 × 10⁵ — el
  modelo está mal especificado para esto. La capa requiere joint orbit +
  mass con covarianza AL real (audit #6, semanas de trabajo).
- **Catálogos caracterizados**: `encounters_characterized_full.parquet`
  (72.236.904 filas, con observabilidad Gaia + magnitudes/diámetros, generado
  vía streaming) cubre el catálogo congelado completo. El antiguo
  `encounters_characterized.parquet` (158.672 filas) es una referencia de
  desarrollo. La corrida congelada (72M) es la única canónica.
- (2) Pallas (i = 34.9°) aparece con 47 encuentros en el catálogo (el más
  cercano a 6.3 × 10⁻³ AU) — su alta inclinación orbital la mantiene
  parcialmente separada del plano del cinturón.

## ✅ Validación

El pipeline se valida en tres niveles:

### 1. Tests de regresión (CI)

```bash
docker compose run --rm test pytest tests/ -v   # ~300 tests, todos pasan (regresión, no validación científica completa)
```

### 2. Cross-match contra catálogos publicados

Resumen consolidado en [`docs/literature_validation.md`](docs/literature_validation.md).
Catálogos independientes de encuentros conocidos:

- **Gate de cuerpos grandes** (criterio del CLAUDE.md): (1) Ceres, (2) Pallas,
  (4) Vesta, (10) Hygiea presentes con conteos 352/47/458/162 — test de
  regresión `tests/test_validation.py::TestFrozenMajorBodyGate`.
- **Fienga et al. (2003)** [A&A 406, 751] — **3/4** en ventana (residual mediano
  52 µAU; el par (804,733) es un detection gap near-threshold).
- **Galád & Gray (2002)** [A&A 391, 1115] — **4/4** encuentros de Hygiea (µAU–20 µAU).
- **Fuentes-Muñoz et al. (2025)** [AJ 170, 353] — **11.804 / 40.004 (29,5 %)**
  de los pares perturbador→objetivo de su Tabla 5 (Gaia FPR) están presentes en
  el catálogo DR3. Es un **lower bound** (FPR ajusta sobre baseline completo →
  la mayoría de los encuentros caen fuera de la ventana DR3); los presentes son
  confirmaciones independientes.

```bash
docker compose run --rm pipeline python -m scripts.ingest.download_fienga_2003
docker compose run --rm pipeline python -m scripts.ingest.download_galad_2002
docker compose run --rm pipeline python -m scripts.ingest.download_fuentes_munoz
# Después del pipeline:
docker compose run --rm pipeline python -m scripts.validate.validate_fienga_2003
docker compose run --rm pipeline python -m scripts.validate.validate_galad_2002
docker compose run --rm pipeline python -m scripts.validate.validate_fuentes_munoz_2025
```

> Goffin (2014) no se puede cruzar par-a-par: su catálogo VizieR `J/A+A/565/A56`
> contiene sólo tablas de masas, no la lista de encuentros (ver el doc consolidado).

### 3. Spot-check contra JPL Horizons (ground truth)

Cada par matcheado se vuelve a consultar contra JPL Horizons (DE440 + N-body completo) en una ventana fina alrededor de la fecha:

```bash
docker compose run --rm pipeline python -m scripts.validate.validate_jpl_horizons
```

Resultados (8 pares de literatura): |nuestro − JPL| ≤ ~5 × 10⁻⁶ AU **a la cadencia
de muestreo de JPL** (1 h / 30 min) — no es una prueba de precisión sub-cadencia
global (ver [`FROZEN_RUN.md`](FROZEN_RUN.md) límite 3).

### 4. Multi-snapshot MPCORB (precisión histórica)

La propagación Kepler 2-cuerpos acumula error si la época de los elementos orbitales está lejos de la ventana temporal de interés. Para mitigar:

```bash
# Descarga un snapshot histórico del Wayback Machine al año pedido
docker compose run --rm pipeline python -m scripts.ingest.download_mpcorb_historical --year 2015 --month 6
```

El pipeline auto-detecta los snapshots disponibles en `data/raw/mpcorb_archive/` y selecciona el de época más cercana al centro de la ventana. Para la ventana Gaia (2014–2017), un snapshot 2015 reduce el error de ~0.03 AU (MPCORB actual ≈ 2026) a < 0.001 AU.

## 🔁 Reproducibilidad

Todo corre dentro de Docker; no se requiere Python local. Para reproducir el
catálogo congelado desde cero:

```bash
# 0. Build de la imagen
docker compose build

# 1. Datos: MPCORB actual + snapshot histórico cercano a la ventana Gaia
docker compose run --rm pipeline python -m scripts.ingest.download_mpcorb
docker compose run --rm pipeline python -m scripts.ingest.download_mpcorb_historical --year 2016 --month 2
#    (la corrida congelada usó MPCORB_20160217; el pipeline auto-selecciona el
#     snapshot de época más cercana al centro de la ventana — ver FROZEN_RUN.md)

# 2. Pipeline de detección + caracterización (streaming para el catálogo completo)
docker compose run --rm pipeline python -m scripts.pipeline.run_pipeline --config config.local.yaml
docker compose run --rm pipeline python -m scripts.pipeline.characterize_catalog \
    --input data/output/encounters_catalog_hybrid_stageb.parquet --streaming on

# 3. Validación contra literatura
docker compose run --rm pipeline python -m scripts.ingest.download_fienga_2003
docker compose run --rm pipeline python -m scripts.ingest.download_galad_2002
docker compose run --rm pipeline python -m scripts.ingest.download_fuentes_munoz
docker compose run --rm pipeline python -m scripts.validate.validate_fienga_2003
docker compose run --rm pipeline python -m scripts.validate.validate_galad_2002
docker compose run --rm pipeline python -m scripts.validate.validate_fuentes_munoz_2025

# 4. Gate de regresión de los 4 cuerpos grandes (requiere el catálogo congelado local)
RUN_REAL_CATALOG_TESTS=1 docker compose run --rm -e RUN_REAL_CATALOG_TESTS=1 \
    pipeline pytest tests/test_validation.py::TestFrozenMajorBodyGate -q

# 5. Dashboard
docker compose up dashboard      # http://localhost:8501
```

**Notas de reproducibilidad** (ver [`FROZEN_RUN.md`](FROZEN_RUN.md) para hashes):

- La réplica bit-a-bit requiere el mismo commit de código y las mismas versiones
  de dependencias listadas en el sidecar de provenance.
- MPCORB se actualiza: cada corrida registra el hash/fecha del snapshot usado.
- Los datos crudos y los outputs (`data/raw`, `data/output`) están gitignored;
  se regeneran con los comandos de arriba.

## 🛠️ Desarrollo

```bash
# Tests
docker compose run --rm test

# Tests específicos
docker compose run --rm test pytest tests/test_ingest_gaia_sso.py -v

# Formato y lint
docker compose run --rm pipeline ruff check . --fix
docker compose run --rm pipeline black .

# Type checking
docker compose run --rm pipeline mypy src/
```

### Debugear con VS Code

Para adjuntar el debugger de VS Code a un proceso corriendo en Docker:

```bash
# 1. Levantar el contenedor en modo debug (queda esperando la conexión)
docker compose -f docker-compose.yml -f docker-compose.debug.yml \
  run --rm --service-ports pipeline -m scripts.ingest.download_gaia_sso --config config.yaml

# 2. En VS Code: Run & Debug → "Docker: Attach to pipeline" → ▶
```

El proceso no arranca hasta que VS Code se conecte. Podés poner breakpoints en cualquier archivo de `src/` o `scripts/`.

## 📚 Referencias

- Tanga, P., et al. (2023). *Gaia Data Release 3. The Solar System survey*. A&A 674, A12. [[link]](https://www.aanda.org/articles/aa/full_html/2023/06/aa43796-22/aa43796-22.html)
- Goffin, E. (2014). *Astrometric asteroid masses: a simultaneous determination*. A&A 565, A56.
- Fuentes-Muñoz, O., Farnocchia, D., Giorgini, J. D., & Park, R. S. (2025). *Asteroid Mass Estimation by Mutual Perturbations during Close Encounters after Gaia FPR*. AJ 170, 353. [DOI:10.3847/1538-3881/ae0cc9](https://doi.org/10.3847/1538-3881/ae0cc9)
- Documentación oficial de Gaia DR3: [gea.esac.esa.int](https://gea.esac.esa.int/archive/documentation/)
- Formato MPCORB: [minorplanetcenter.net](https://www.minorplanetcenter.net/iau/info/MPOrbitFormat.html)

## 📝 Notas

Este proyecto es un trabajo **académico/educativo**. No pretende reemplazar pipelines profesionales de servicios de efemérides (JPL, IMCCE), sino servir como ejercicio práctico de astronomía computacional y data engineering aplicado.

Si te resulta útil o tenés sugerencias, sentite libre de abrir un issue o un PR.

## 📄 Licencia

MIT — ver `LICENSE`.

## 🙏 Agradecimientos

- Misión Gaia (ESA) por hacer públicos los datos.
- Minor Planet Center y JPL por mantener efemérides accesibles.
- La comunidad de `astropy`, `rebound`, `poliastro` y `duckdb`.
