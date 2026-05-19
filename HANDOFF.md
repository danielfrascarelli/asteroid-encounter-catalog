# Handoff — Estado del análisis de encuentros novedosos (2026-05-19)

> Documento de continuidad. Si retomás el trabajo más tarde o lo retoma otra persona,
> este archivo describe **qué se hizo**, **qué quedó pendiente** y **cómo seguir**.

---

## Contexto rápido

El proyecto detecta encuentros cercanos entre asteroides durante la ventana Gaia DR3
(2014-07-25 → 2017-05-28). El pipeline produce ~4 millones de encuentros a 0.05 AU.

Después de cruzar contra toda la literatura publicada (Fienga 2003, Galád 2002,
Goffin 2014, Fuentes-Muñoz 2024), **119.545 encuentros son novedosos** (no aparecen
en ningún paper). De esos identificamos **41 candidatos viables a determinación de
masa asteroidal** usando astrometría Gaia DR3.

El objetivo científico: pesar asteroides cuya masa no se conoce todavía.

---

## Lo que se hizo (PRs)

### PRs mergeados a main

- **PR #6** ([feat/encounter-analysis-filter](https://github.com/danielfrascarelli/asteroid-encounter-catalog/pull/6)) — `encounter_analysis/`
  - `filter_candidates.py` — filtra los 119.545 encuentros a 379 relevantes
  - `README.md` — documenta los criterios de filtrado

### PRs abiertos (esperando merge)

- **PR #1** ([feat/threshold-0.05-default](https://github.com/danielfrascarelli/asteroid-encounter-catalog/pull/1)) — cambia threshold default 0.01 → 0.05 AU
- **PR #2** ([refactor/no-pipeline-defaults](https://github.com/danielfrascarelli/asteroid-encounter-catalog/pull/2)) — quita defaults de `detect_encounters()`
- **PR #3** ([docs/fix-data-download-instructions](https://github.com/danielfrascarelli/asteroid-encounter-catalog/pull/3)) — fixes README/ROADMAP
- **PR #4** ([fix/wayback-429-retry](https://github.com/danielfrascarelli/asteroid-encounter-catalog/pull/4)) — fix HTTP 429 en download_mpcorb_historical
- **PR #5** ([docs/analysis-plan-and-papers](https://github.com/danielfrascarelli/asteroid-encounter-catalog/pull/5)) — ANALYSIS_PLAN.md + 6 PDFs de referencia
- **PR #7** ([feat/novel-encounter-analysis-scripts](https://github.com/danielfrascarelli/asteroid-encounter-catalog/pull/7)) — pipeline completo de análisis Cat A/B + verificación Gaia + demo de deflexión

### Scripts agregados (en PR #7)

- `scripts/validate_novel_a.py` — valida 7 encuentros con Ceres/Vesta contra JPL Horizons. **MAE = 3.21e-4 AU** ✅ (umbral 5e-3)
- `scripts/analyze_mass_candidates.py` — rankea candidatos Cat B por deflexión esperada
- `scripts/check_gaia_observations.py` — verifica que Gaia haya observado el target ±180 días del encuentro (TAP async, sin truncar)
- `scripts/demo_ate_deflection.py` — demo end-to-end intentando detectar la perturbación de (111) Ate sobre 2000_nt3

### Outputs generados en `data/output/`

- `cat_a_jpl_validation.csv` — 7 filas, validación Cat A contra JPL
- `mass_candidates.csv` — 100 candidatos Cat B rankeados por deflexión
- `gaia_observations_check.csv` — 100 candidatos con conteo de transits Gaia antes/después
- `publishable_mass_candidates.csv` — **41 candidatos viables**, ordenados por δ
- `ate_2000nt3_residuals.csv` — residuales de 76 observaciones de 2000_nt3 alrededor del encuentro de Ate

### Reportes

- `ANALYSIS_PLAN.md` — plan original + resultados de Cat A y Cat B
- `encounter_analysis/RESULTS.md` — reporte científico final de los 41 candidatos

---

## El resultado científico

**41 encuentros novedosos potencialmente publicables como nuevas determinaciones
de masa asteroidal.** Todos cumplen:

1. Señal de deflexión esperada ≥ 100 μas (la precisión por tránsito de Gaia)
2. Al menos 3 observaciones Gaia ANTES del encuentro
3. Al menos 3 observaciones Gaia DESPUÉS del encuentro
4. Validado contra JPL Horizons que la geometría del pipeline es correcta (MAE < 5e-4 AU)

Top 3 destacados:

| # | Perturber (a pesar) | Target | Fecha | δ (μas) | obs antes/después | mag |
|---|---------------------|--------|-------|---------|-------------------|-----|
| 1 | **(57) Mnemosyne** | 2008_ef40 | 2016-08-26 | 22.712 | 24 / 49 | 20.0 |
| 2 | **(111) Ate** | 2000_nt3 | 2016-06-08 | 4.906 | 39 / 37 | **18.6** ← brillante |
| 3 | **(241) Germania** | 2000_jc23 | 2016-06-27 | 3.495 | 29 / **204** | 19.5 |

Ver lista completa en `data/output/publishable_mass_candidates.csv`.

---

## Lo que quedó pendiente

### 1. Detección real de la perturbación (CRÍTICO para publicación)

El script `demo_ate_deflection.py` intenta detectar la perturbación de (111) Ate
sobre 2000_nt3 propagando la órbita MPCORB con Kepler 2-cuerpos y comparando contra
las observaciones Gaia. **No funciona** por una razón esperada:

- MPCORB.DAT snapshot disponible: 2012-09-18 (epoch_jd = 2456200.5)
- Encuentro: 2016-06-08
- Propagación Kepler 2-cuerpos sobre ~4 años acumula error de ~875 arcsec
- Señal de Ate: ~5 mas
- **Ratio señal/ruido: 5 mas / 875.000 mas ≈ 6e-6** → indetectable así

La técnica clásica de mass determination resuelve esto con un **fit de órbita a las
observaciones pre-encuentro**:

1. Tomar las primeras N observaciones de 2000_nt3 antes del encuentro (~30-90 días)
2. Fitear una órbita Kepler a esas observaciones con `scipy.optimize.least_squares`
   (6 parámetros: a, e, i, Ω, ω, M₀)
3. Propagar esa órbita forward al instante de cada observación post-encuentro
4. Los residuales (observado - predicho) deberían ser ~0 antes y ~deflexión-de-Ate después

**Esfuerzo estimado**: 2-3 horas. Requiere coordinar:
- `src/propagate/kepler.py` (ya existe)
- `scipy.optimize.least_squares`
- Las transformaciones de coordenadas en `demo_ate_deflection.py` (heliocéntrico
  eclíptico → barycéntrico ICRS → línea de vista desde Gaia → RA/Dec)
- Manejo de incertidumbre por observación (usar el error formal de Gaia si está)

**Alternativa más barata**: bajar un snapshot histórico MPCORB de 2015 (con
`scripts.download_mpcorb_historical --year 2015 --month 6`) y reintentar el demo.
La propagación sería de ~1 año en vez de 4 → error de Kepler ~30-50 mas, todavía
arriba de la señal pero mucho mejor. Hoy intentamos ese download y dio HTTP 429
(rate limit). El fix está implementado (PR #4) pero no probado.

### 2. Cross-check contra masas anunciadas recientemente

Algunos de los 41 candidatos podrían ya tener masas publicadas en preprints
recientes que no están en los catálogos viejos. Vale la pena hacer una búsqueda
en ADS / arXiv por cada perturber del top 10:

```
(57) Mnemosyne mass
(111) Ate mass
(241) Germania mass
...
```

Y descartar los que ya tengan masa medida.

### 3. Propagación N-body para los top candidatos

`config.propagation.method: rebound` activa N-body con Sol + Júpiter + Saturno
y opcionalmente Ceres/Vesta/Pallas/Hygiea. Re-correr el detector con esto sobre
los pares Cat B confirmaría que las geometrías son robustas (no artefactos de
Kepler).

### 4. Estimación de incertidumbre

Cada masa potencial necesita una barra de error. Eso requiere:
- Densidad asumida (usamos 1.5 g/cm³ - varía entre 1.0 y 3.5 según taxonomía)
- Incertidumbre en la geometría del encuentro
- Incertidumbre astrométrica de Gaia (depende de magnitud y N transits)

### 5. Publicación / paper

Si los 3-5 top candidatos se detectan con un fit real (item 1), hay material
para un paper técnico:
- "A catalog of N novel asteroid mass-determination candidates from Gaia DR3"
- Estructura: pipeline, validación (Cat A), candidatos (Cat B), detección
  ejemplo (uno o dos casos resueltos), tabla completa.

---

## Cómo reanudar

### Setup

Ya está todo configurado. El proyecto corre dentro de Docker:

```bash
cd /home/daniel/Documents/gaia-project/gaia
git checkout feat/novel-encounter-analysis-scripts   # o main si PR #7 ya mergeó
docker compose build   # solo si cambiaste pyproject.toml
```

### Estado actual de los datos

Lo que está en `data/raw/`:
- `MPCORB.DAT` — snapshot actual de MPCORB (varios MB)
- `mpcorb_archive/MPCORB_20120918.DAT` — snapshot histórico 2012 (existente)
- ❌ `mpcorb_archive/MPCORB_2015*.DAT` — **falta**. Download falló por 429.
- `gaia_sso.parquet` — descarga incompleta (1103 asteroides solamente).
  No es necesario para el pipeline, sólo para `gaia_has_target` en analyze_mass_candidates.

Lo que está en `data/output/`:
- Todos los archivos listados arriba ya están generados.

### Comandos para regenerar todo

```bash
# 1. (opcional) Descargar snapshot 2015 para mejor propagación
docker compose run --rm pipeline python -m scripts.download_mpcorb_historical --year 2015 --month 6

# 2. Filtrar candidatos relevantes (rápido, sin red)
docker compose run --rm pipeline python encounter_analysis/filter_candidates.py

# 3. Rankear Cat B por señal de deflexión (rápido, sin red)
docker compose run --rm pipeline python -m scripts.analyze_mass_candidates --top-n 100

# 4. Verificar observabilidad contra Gaia DR3 (~10 s, requiere red)
docker compose run --rm pipeline python -m scripts.check_gaia_observations

# 5. (opcional) Validar Cat A contra JPL Horizons (~25 s)
docker compose run --rm pipeline python -m scripts.validate_novel_a

# 6. Demo del intento de detección Ate→2000_nt3
docker compose run --rm pipeline python -m scripts.demo_ate_deflection
```

### Próximo paso concreto a implementar

**Crear `scripts/fit_pre_encounter_orbit.py`** que haga el fit Kepler a las
observaciones pre-encuentro del target. Esqueleto:

```python
"""Fit pre-encounter Kepler orbit to Gaia observations of a target asteroid.

For a (perturber, target, date) tuple from the candidate catalog:
1. Download Gaia transits in [encounter - half_window, encounter - blackout].
2. Initial guess: MPCORB osculating elements.
3. scipy.optimize.least_squares over (a, e, i, Ω, ω, M₀) minimising
   Σ residuals_ra² + residuals_dec² (mas).
4. Propagate fitted orbit forward to all post-encounter observations.
5. Output: post-encounter residuals vs time.

If residuals show systematic offset > Gaia precision (~100 μas), the
encounter signal is detected.
"""

import scipy.optimize
from src.propagate.kepler import kepler_to_cartesian
# ... reuse coordinate transforms from demo_ate_deflection.py
```

Tiempo estimado: 2-3 horas implementación + 1 hora pruebas y debugging.

---

## Branches y archivos relevantes

| Archivo | Branch | Estado |
|---------|--------|--------|
| `encounter_analysis/filter_candidates.py` | main | mergeado (PR #6) |
| `encounter_analysis/README.md` | main | mergeado (PR #6) |
| `encounter_analysis/RESULTS.md` | feat/novel-encounter-analysis-scripts | PR #7 abierto |
| `scripts/validate_novel_a.py` | feat/novel-encounter-analysis-scripts | PR #7 abierto |
| `scripts/analyze_mass_candidates.py` | feat/novel-encounter-analysis-scripts | PR #7 abierto |
| `scripts/check_gaia_observations.py` | feat/novel-encounter-analysis-scripts | PR #7 abierto |
| `scripts/demo_ate_deflection.py` | feat/novel-encounter-analysis-scripts | PR #7 abierto |
| `ANALYSIS_PLAN.md` | feat/novel-encounter-analysis-scripts | PR #7 abierto |
| `HANDOFF.md` (este archivo) | feat/novel-encounter-analysis-scripts | PR #7 abierto |

Para empezar limpio: mergear todos los PRs abiertos en orden (#1, #2, #3, #4, #5, #7).

---

## Decisiones tomadas en esta sesión que vale la pena recordar

- **Threshold default 0.05 AU** (antes 0.01). Esto da más encuentros pero menos precisión por encuentro. Para mass determination, el filtro `dist_au < 0.02` en filter_candidates.py recupera la precisión.
- **No defaults en `detect_encounters()`**: todos los parámetros vienen de config. Si querés correr con valores distintos, editá `config.yaml` o pasá `config.local.yaml`.
- **`download_gaia_sso` NO es prerequisito del pipeline.** Solo lee MPCORB.
- **Gaia `epoch` column en `sso_observation` es days since J2010.0 TCB**, no JD como dice el metadata. Tener cuidado al trabajar con esa tabla.
- **TAP sync limita a 2000 rows.** Usar `launch_job_async` para queries grandes.
- **Densidad asumida ρ = 1.5 g/cm³** para mass-from-diameter. Es un C-type/MBA típico. Puede variar factor ~2.

---

## Contacto / referencias

- Repo: <https://github.com/danielfrascarelli/asteroid-encounter-catalog>
- Papers de validación: ver `papers/` (6 PDFs incluidos)
- Gaia DR3 archive: <https://gea.esac.esa.int/archive/>
- ADQL docs: <https://gea.esac.esa.int/archive/documentation/GDR3/Gaia_archive/chap_datamodel/>
