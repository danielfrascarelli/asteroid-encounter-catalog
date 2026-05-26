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
docker compose run --rm pipeline python -m scripts.download_mpcorb

# 2. Snapshot histórico 2015 (recomendado: reduce error Kepler de ~30 mAU a <1 mAU)
# El pipeline auto-selecciona el snapshot con época más cercana al centro de la ventana.
docker compose run --rm pipeline python -m scripts.download_mpcorb_historical --year 2015 --month 6
```

> **Nota**: `download_gaia_sso` no es requerido para el pipeline de detección. El pipeline lee únicamente MPCORB.

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

La corrida canónica está documentada en [`FROZEN_RUN.md`](FROZEN_RUN.md) — éste
es el único documento que debe citarse para reproducir o referenciar el
catálogo geométrico.

| Métrica | Valor |
|---------|-------|
| MPCORB snapshot | `MPCORB_20160217.DAT` |
| Ventana temporal | 2014-07-25 → 2017-05-28 (Gaia DR3) |
| Umbral de detección | **0.05 AU** |
| Scan / refine | rebound (whfast, Sun+Jupiter+Saturn) / Kepler |
| **Encuentros detectados** | **72.236.904** |
| Encuentro más cercano | **6.6 × 10⁻⁶ AU** (≈ 988 km) |
| Encuentros < 0.001 AU | 26.038 |
| Encuentros < 0.01 AU | 2.833.425 |

Hashes de inputs/outputs y la tabla completa de claims en
[`FROZEN_RUN.md`](FROZEN_RUN.md).

### Top-3 encuentros más cercanos

| Rank | Cuerpo 1 | Cuerpo 2 | Distancia (AU) |
|------|----------|----------|----------------|
| 1 | (153222) 2000 YD43 | (238587) 2004 YX3 | 6.6 × 10⁻⁶ |
| 2 | (15072) Landolt | (387599) 2001 XF180 | 1.2 × 10⁻⁵ |
| 3 | (270730) 2002 QE130 | (366918) 2005 UC211 | 1.5 × 10⁻⁵ |

### Notas

- El catálogo geométrico es el único output considerado publicable. El
  archivo `publishable_mass_candidates.csv` (41 filas) son **candidatos**
  para fitting de masas, no masas medidas; la capa de mass-fitting es
  exploratoria (ver caveats en `FROZEN_RUN.md`).
- El catálogo viejo `encounters_catalog.parquet` (158.672 filas, umbral
  0.05 AU sobre un subset menor de asteroides) y su versión caracterizada
  `encounters_characterized.parquet` quedan como referencia histórica;
  la corrida congelada (72M) los reemplaza.
- (2) Pallas (i = 34.9°) aparece con 47 encuentros en el catálogo (el más
  cercano a 6.3 × 10⁻³ AU) — su alta inclinación orbital la mantiene
  parcialmente separada del plano del cinturón.

## ✅ Validación

El pipeline se valida en tres niveles:

### 1. Tests de regresión (CI)

```bash
docker compose run --rm test pytest tests/ -v   # 198 tests, todos pasan
```

### 2. Cross-match contra catálogos publicados

Dos catálogos independientes de encuentros conocidos se descargan automáticamente:

- **Fienga et al. (2003)** [A&A 406, 751] — predicciones N-body 2003–2022, VizieR `J/A+A/406/751`.
- **Galád & Gray (2002)** [A&A 391, 1115] — candidatos para determinación de masas, parseados del HTML del artículo.

Para correr:

```bash
docker compose run --rm pipeline python -m scripts.download_fienga_2003
docker compose run --rm pipeline python -m scripts.download_galad_2002

# Después del pipeline:
docker compose run --rm pipeline python -m scripts.validate_fienga_2003
docker compose run --rm pipeline python -m scripts.validate_galad_2002
```

Resultados a 0.05 AU sobre el catálogo Kepler de 4M encuentros: **100% match** (4/4 Fienga + 4/4 Galád en la ventana Gaia).

### 3. Spot-check contra JPL Horizons (ground truth)

Cada par matcheado se vuelve a consultar contra JPL Horizons (DE440 + N-body completo) en una ventana fina alrededor de la fecha:

```bash
docker compose run --rm pipeline python -m scripts.validate_jpl_horizons
```

Resultados: MAE(nuestro − JPL) = 0.0002 AU; MAE(literatura − JPL) = 0.00004 AU. JPL agrees con la literatura, confirmando ambos catálogos como referencias confiables.

### 4. Multi-snapshot MPCORB (precisión histórica)

La propagación Kepler 2-cuerpos acumula error si la época de los elementos orbitales está lejos de la ventana temporal de interés. Para mitigar:

```bash
# Descarga un snapshot histórico del Wayback Machine al año pedido
docker compose run --rm pipeline python -m scripts.download_mpcorb_historical --year 2015 --month 6
```

El pipeline auto-detecta los snapshots disponibles en `data/raw/mpcorb_archive/` y selecciona el de época más cercana al centro de la ventana. Para la ventana Gaia (2014–2017), un snapshot 2015 reduce el error de ~0.03 AU (MPCORB actual ≈ 2026) a < 0.001 AU.

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
