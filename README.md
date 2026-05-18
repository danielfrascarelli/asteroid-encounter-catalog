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
# MPCORB (~50 MB, varios minutos)
docker compose run --rm pipeline python -m scripts.download_mpcorb

# Observaciones Gaia SSO (~varios GB, puede tardar horas)
# Los chunks se guardan en data/cache/gaia_sso_chunks/ — si se interrumpe,
# el siguiente run retoma desde el último chunk completado.
docker compose run --rm pipeline python -m scripts.download_gaia_sso --config config.yaml
```

### Pipeline completo

```bash
# Detección + caracterización completa (todos los asteroides numerados)
docker compose run --rm pipeline python -m scripts.run_pipeline
docker compose run --rm pipeline python -m scripts.characterize_catalog
```

Produce `data/output/encounters_characterized.parquet` (~119k filas) y el sidecar de metadatos.

### Subset rápido para pruebas

```bash
# Configurar para corrida rápida (~5000 asteroides, ~10 minutos)
cp config.yaml config.local.yaml
# editar config.local.yaml → subset.max_asteroids: 5000

docker compose run --rm pipeline python -m scripts.run_pipeline --config config.local.yaml
docker compose run --rm pipeline python -m scripts.characterize_catalog --config config.local.yaml
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

## 📊 Resultados (corrida completa sobre todos los asteroides numerados)

| Métrica | Valor |
|---------|-------|
| Asteroides procesados | ~99.999 numerados (MPCORB) |
| Ventana temporal | 2014-07-25 → 2017-05-28 (Gaia DR3) |
| Umbral de detección | 0.01 AU |
| **Encuentros detectados** | **119.546** |
| Gaia-observables (elong > 45°, mag < 21) | 50.473 |
| Encuentro más cercano | **0.000043 AU** (≈ 6.434 km) |
| Velocidad relativa (rango) | 0.032 – 25.23 km/s |
| Diámetro cuerpo 1 (rango) | 1 – 795 km |
| Cuerpos grandes confirmados | Ceres (2 enc.), Vesta (8 enc.), Hygiea (1 enc.) |

### Encuentros destacados

| Rank | Cuerpo 1 | Cuerpo 2 | Distancia (AU) | Fecha | Vel (km/s) |
|------|----------|----------|----------------|-------|------------|
| 1 | 193507 | 343572 | 0.000043 | 2016-07-23 | 0.24 |
| 2 | 63313 | 197297 | 0.000056 | 2015-02-10 | 0.65 |
| 3 | 78160 | 176588 | 0.000061 | 2016-11-03 | 0.52 |

### Nota sobre (2) Pallas

Pallas tiene inclinación i = 34.9° (la mayor entre los asteroides masivos), lo que mantiene su órbita bien separada del plano del cinturón principal durante la ventana Gaia. No se detectaron encuentros < 0.01 AU con Pallas — este es un resultado físicamente correcto, no un bug del pipeline.

## ✅ Validación

El pipeline se valida contra encuentros conocidos en la literatura:

- Encuentros con asteroides masivos (Ceres, Vesta, Pallas, Hygiea) en el período Gaia.
- Pares reportados por Goffin (2014) y Fuentes-Muñoz et al. (2024).
- Validación contra JPL Horizons: encuentro más cercano confirmado con posición baricéntrica.

Los tests de regresión verifican que estos casos aparezcan en cada corrida:

```bash
docker compose run --rm test pytest tests/test_validation.py -v
```

Validación cruzada con literatura:

```bash
docker compose run --rm pipeline python -m scripts.validate_literature
```

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
  run --rm --service-ports pipeline -m scripts.download_gaia_sso --config config.yaml

# 2. En VS Code: Run & Debug → "Docker: Attach to pipeline" → ▶
```

El proceso no arranca hasta que VS Code se conecte. Podés poner breakpoints en cualquier archivo de `src/` o `scripts/`.

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
