# 🌌 Catálogo de Encuentros Cercanos entre Asteroides

> Detección sistemática de pares de asteroides que pasaron cerca uno del otro durante el período de observación de la misión Gaia (DR3), construyendo un catálogo nuevo de encuentros cercanos con su geometría y propiedades físicas.

![Python](https://img.shields.io/badge/python-3.11+-blue.svg)
![Status](https://img.shields.io/badge/status-en%20desarrollo-yellow.svg)
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

- Python ≥ 3.11
- ~10 GB de espacio en disco (datos crudos + cache)
- 8 GB de RAM mínimo (16 GB recomendado para subset completo)

### Instalación

```bash
# Clonar el repo
git clone https://github.com/tu-usuario/asteroid-encounters.git
cd asteroid-encounters

# Crear entorno e instalar dependencias
uv sync
# o alternativamente:
# poetry install
```

### Descarga de datos

```bash
# MPCORB (~50 MB, varios minutos)
python -m scripts.download_mpcorb

# Observaciones Gaia SSO (~varios GB, dependiendo del subset)
python -m scripts.download_gaia_sso --subset numbered
```

### Primer test sobre subset pequeño

```bash
# Configurar para corrida rápida (1000 asteroides)
cp config.yaml config.local.yaml
# editar config.local.yaml → subset.max_asteroids: 1000

# Correr pipeline
python -m scripts.run_pipeline --config config.local.yaml
```

Deberías ver salida estilo:
```
[INFO] Cargando 1000 asteroides desde MPCORB...
[INFO] Propagando órbitas (Δt = 1h, 3 años)...
[INFO] Construyendo KD-trees por step temporal...
[INFO] Detectando encuentros (umbral = 0.01 AU)...
[INFO] Refinando 47 candidatos...
[INFO] 38 encuentros confirmados → data/output/encounters.parquet
```

### Explorar resultados

```bash
streamlit run src/dashboard/app.py
```

Y abrir `http://localhost:8501` en el navegador.

## 🔧 Configuración

Todos los parámetros del pipeline están en `config.yaml`. Los más importantes:

| Parámetro | Default | Descripción |
|-----------|---------|-------------|
| `detection.threshold_au` | `0.01` | Umbral de distancia para considerar "encuentro" |
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

## ✅ Validación

El pipeline se valida contra encuentros conocidos en la literatura:

- Encuentros con asteroides masivos (Ceres, Vesta, Pallas, Hygiea) en el período Gaia.
- Pares reportados por Goffin (2014) y Fuentes-Muñoz et al. (2024).

Los tests de regresión verifican que estos casos aparezcan en cada corrida:

```bash
pytest tests/test_validation.py -v
```

## 🛠️ Desarrollo

```bash
# Tests
pytest tests/ -v

# Formato y lint
ruff check . --fix
black .

# Type checking
mypy src/
```

Antes de hacer commit, asegurate de que pasen los tests y el lint:

```bash
make check  # tests + lint + format check
```

## 📚 Referencias

- Tanga, P., et al. (2023). *Gaia Data Release 3. The Solar System survey*. A&A 674, A12. [[link]](https://www.aanda.org/articles/aa/full_html/2023/06/aa43796-22/aa43796-22.html)
- Goffin, E. (2014). *Astrometric asteroid masses: a simultaneous determination*. A&A 565, A56.
- Fuentes-Muñoz, O., et al. (2024). *Asteroid masses from Gaia FPR*. LPSC #2388.
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
