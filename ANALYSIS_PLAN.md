# Plan: Análisis de Encuentros Novedosos — Categoría A y B

> ✅ Implementado. Scripts en `scripts/validate_novel_a.py` y `scripts/analyze_mass_candidates.py`.

---

## Contexto

El pipeline detectó 119.545 encuentros no documentados en ninguna literatura científica. Después de filtrar por `gaia_observable`, `dist_au < 0.02`, `diameter_1_km > 30` y `rel_vel_km_s < 8.0`, quedaron 379 encuentros ordenados por `deflection_score`.

Input: `data/output/relevant_novel_encounters.csv`

---

## Categoría A — Encuentros con Ceres/Vesta (calibración de precisión)

**7 encuentros** · Perturbers: (1) Ceres, (4) Vesta · Masas conocidas con alta precisión (misión Dawn)

**Objetivo**: Usar como benchmark de precisión del pipeline. Las masas conocidas permiten calcular la deflexión astrométrica predicha y comparar nuestra distancia/época detectada contra JPL Horizons como ground truth.

### Encuentros

| Perturber | Target | Fecha | dist (AU) | vel (km/s) |
|-----------|--------|-------|-----------|------------|
| (4) Vesta | 2001_uy75 | 2015-11-02 | 0.00695 | 2.46 |
| (4) Vesta | 1997_ce20 | 2016-10-19 | 0.00930 | 3.37 |
| (1) Ceres | 2005_uh343 | 2017-01-08 | 0.00460 | 4.64 |
| (4) Vesta | 2000_sa336 | 2017-02-28 | 0.00412 | 5.62 |
| (4) Vesta | 2001_su160 | 2016-12-14 | 0.00610 | 4.91 |
| (4) Vesta | 2000_am232 | 2016-12-19 | 0.00717 | 4.84 |
| (1) Ceres | 2005_qf142 | 2017-01-22 | 0.00869 | 6.62 |

### Tareas

1. Crear `scripts/validate_novel_a.py`:
   - Para cada par, consultar JPL Horizons en ventana ±2 días, pasos de 30 min
   - Reutilizar `_horizons_vectors()` y `_jpl_min_distance()` de `scripts/validate_jpl_horizons.py`
   - Calcular `our_dist` vs `jpl_dist`, diferencia y MAE global
   - Calcular deflexión astrométrica esperada: δ ≈ 2GM/(v²·b) [rad → μas]
     - Masas: Ceres = 9.384e20 kg, Vesta = 2.591e20 kg (Dawn)
   - Output: `data/output/cat_a_jpl_validation.csv`

2. Criterio de éxito: MAE(nuestro − JPL) < 0.005 AU para los 7 encuentros

### Resultados ✅ (2026-05-19)

| Perturber | Target | Nuestro (AU) | JPL (AU) | Δ (AU) | δ esperada (μas) |
|-----------|--------|-------------|----------|--------|-----------------|
| Vesta | 2001_uy75 | 0.006954 | 0.006952 | +1.8e-6 | 1,129,425 |
| Vesta | 1997_ce20 | 0.009304 | 0.009310 | −6.6e-6 | 451,426 |
| Ceres | 2005_uh343 | 0.004601 | 0.004008 | +5.9e-4 | 1,742,683 |
| Vesta | 2000_sa336 | 0.004118 | 0.004011 | +1.1e-4 | 367,039 |
| Vesta | 2001_su160 | 0.006096 | 0.006240 | −1.4e-4 | 324,803 |
| Vesta | 2000_am232 | 0.007174 | 0.006806 | +3.7e-4 | 283,348 |
| Ceres | 2005_qf142 | 0.008688 | 0.007662 | +1.0e-3 | 453,078 |

**MAE(nuestro − JPL) = 3.21e-4 AU** ← bien por debajo del umbral 0.005 AU ✅

Todos los 7 encuentros tienen señales de deflexión >> 100 μas (rango: 283k–1.7M μas).
Estos encuentros con Ceres/Vesta son detectables trivialmente — sirven como benchmark de precisión.

---

## Categoría B — Encuentros con masas desconocidas (candidatos a determinación de masa)

**8 encuentros** · Perturbers sin masa publicada

**Objetivo**: Identificar qué encuentros podrían producir una NUEVA medición de masa asteroidal a partir de la astrometría de Gaia. Filtro clave: señal de deflexión esperada > precisión de Gaia por tránsito (~100 μas).

### Encuentros

| Perturber | Target | Fecha | dist (AU) | vel (km/s) | Score |
|-----------|--------|-------|-----------|------------|-------|
| (57) Mnemosyne | 2008_ef40 | 2016-08-26 | 0.00244 | 2.85 | 1.57e8 |
| (68) Leto | 2002_vy95 | 2015-04-06 | 0.00522 | 2.74 | 8.30e7 |
| (165) Loreley | 1996_tf50 | 2014-12-08 | 0.00254 | 2.82 | 4.59e7 |
| (26) Proserpina | 2010_xj25 | 2015-10-13 | 0.00505 | 2.54 | 4.48e7 |
| (511) Davida | 2002_tf | 2015-01-13 | 0.00599 | 5.03 | 4.21e7 |
| **(111) Ate** | **2000_nt3** | **2016-06-08** | **0.000472** | **5.76** | 3.39e7 ← más cercano |
| (113) Amalthea | 2001_vr121 | 2016-11-24 | 0.00910 | 1.12 | 3.00e7 |
| (46) Hestia | Sitensky | 2016-01-14 | 0.00977 | 1.15 | 2.69e7 |

### Tareas

1. Crear `scripts/analyze_mass_candidates.py`:
   - Para cada perturber: estimar masa a partir del diámetro con densidad media ρ ≈ 1.5 g/cm³
     - M = ρ · (4/3)π · (D/2)³
   - Calcular ángulo de deflexión: δ ≈ 2GM/(v²·b) [rad → μas]
   - Comparar δ contra umbral Gaia (100 μas por tránsito)
   - Verificar si el target aparece en `data/raw/gaia_sso.parquet` (si existe)
   - Output ranked: `data/output/mass_candidates.csv`

2. Misma validación JPL Horizons que Categoría A para los 8 pares

3. Nota especial (111) Ate + 2000_nt3: con b = 0.000472 AU = 70.600 km, la deflexión
   es probablemente >> 100 μas. Calcular explícitamente y destacar en el output.

4. Criterio de éxito: tabla ordenada con columnas
   `perturber, target, date, dist_au, mass_est_kg, deflection_muas, gaia_has_target, viable`

### Resultados ✅ (2026-05-19, top 20, sin JPL)

| Rank | Perturber | Target | dist (AU) | δ (μas) | viable |
|------|-----------|--------|-----------|---------|--------|
| 1 | (57) Mnemosyne | 2008_ef40 | 0.00244 | 22,712 | ✅ |
| 2 | (68) Leto | 2002_vy95 | 0.00522 | 11,999 | ✅ |
| 3 | (18) Melpomene | 2002_tq275 | 0.00535 | 7,139 | ✅ |
| 4 | (165) Loreley | 1996_tf50 | 0.00254 | 6,633 | ✅ |
| 5 | (26) Proserpina | 2010_xj25 | 0.00505 | 6,482 | ✅ |
| 6 | (511) Davida | 2002_tf | 0.00599 | 6,082 | ✅ |
| 7 | (42) Isis | 1998_mf36 | 0.00501 | 5,863 | ✅ |
| 8 | (68) Leto | 2000_sr27 | 0.00739 | 5,583 | ✅ |
| 9 | (110) Lydia | 2000_aj60 | 0.00430 | 5,357 | ✅ |
| **10** | **(111) Ate** | **2000_nt3** | **0.000472** | **4,906** | **✅ ← más cercano** |
| 11–20 | ... | ... | ... | 2,728–4,505 | ✅ |

**Todos los 20 candidatos tienen δ > 2700 μas (27× el umbral Gaia de 100 μas).**

Nota: `gaia_has_target = False` para todos porque el archivo `data/raw/gaia_sso.parquet` local
solo contiene 1103 asteroides (subconjunto del catálogo DR3 completo de ~156k). Este resultado
no es definitivo — requiere verificar contra el catálogo DR3 completo en el Gaia Archive.

---

## Infraestructura a reutilizar

| Necesidad | Fuente |
|-----------|--------|
| Consulta JPL Horizons | `scripts/validate_jpl_horizons.py` → `_horizons_vectors()`, `_jpl_min_distance()` |
| Diámetro → km | `src/characterize/physical.py` → `diameter_km()` |
| Observabilidad Gaia | `src/characterize/observability.py` → `is_gaia_observable()` |
| Datos de entrada | `data/output/relevant_novel_encounters.csv` |

## Archivos creados ✅

- `scripts/validate_novel_a.py`
- `scripts/analyze_mass_candidates.py`
- `data/output/cat_a_jpl_validation.csv` (generado)
- `data/output/mass_candidates.csv` (generado)

## Verificación final

```bash
docker compose run --rm pipeline python -m scripts.validate_novel_a
docker compose run --rm pipeline python -m scripts.analyze_mass_candidates
```
