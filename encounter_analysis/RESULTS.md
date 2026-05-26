# Resultados — Candidatos a Determinación de Masa Asteroidal

> Resultados del análisis de encuentros novedosos detectados en este pipeline.
> Última corrida: 2026-05-19.

---

## Resumen ejecutivo

De los **119.545 encuentros novedosos** en el catálogo (ningún paper publicado los menciona),
identificamos **41 candidatos prioritarios** para determinación de masa asteroidal usando
astrometría de Gaia DR3.

> ⚠️ **Nota metodológica**: los 41 candidatos están priorizados por el score de deflexión
> impulsiva (δ = 2GM/v²b ≥ 100 μas) y cobertura Gaia bracketing. Esto es una condición
> necesaria pero **no suficiente** para una detección publicable. Un resultado de masa real
> requiere un fit LOO/N-body AL-weighted con mejora estadística robusta sobre el modelo
> masa=0. Hasta ahora solo (111) Ate + (18105) 2000 NT3 cuenta con ese análisis completo.

| Etapa | Encuentros | Threshold |
|-------|-----------|-----------|
| Corrida exploratoria (MPCORB ~99k obj.) | 4.036.495 | 0.05 AU |
| Corrida de producción actual (MPCORB 433k obj.) | 119.546 | 0.01 AU |
| No documentados en literatura (producción) | 119.545 | — |
| Físicamente relevantes (filtros astrofísicos) | 379 | — |
| Top 100 por score de deflexión | 100 | — |
| **Prioritarios (δ ≥ 100 μas + Gaia obs en ambos lados)** | **41** | — |
| **Con fit LOO/N-body completo (detección real)** | **1** ← (111) Ate | — |

---

## Pipeline de análisis

### Paso 1 — Filtrado físico ([filter_candidates.py](filter_candidates.py))

Lee `data/output/novel_encounters_not_in_literature.csv` y aplica:
- `gaia_observable == true` (geometría observable por Gaia)
- `dist_au < 0.02` (encuentro genuinamente cercano)
- `diameter_1_km > 30` (perturber con masa apreciable)
- `rel_vel_km_s < 8.0` (deflexión significativa)

Score: `D₁³ / (dist_au × v²)` — proporcional a la deflexión gravitacional esperada.

**Output**: `relevant_novel_encounters.csv` (379 filas).

### Paso 2 — Estimación de señal ([scripts/analyze_mass_candidates.py](../scripts/analyze_mass_candidates.py))

Para cada candidato:
- Masa estimada del perturber: M = ρ · (4/3)π(D/2)³ con ρ = 1.5 g/cm³
- Deflexión: δ = 2GM/(v²b) → conversión a μas
- Flag `viable` si δ ≥ 100 μas (precisión de Gaia por tránsito)

**Output**: `mass_candidates.csv` (top N por score, default 20).

### Paso 3 — Verificación de observabilidad ([scripts/check_gaia_observations.py](../scripts/check_gaia_observations.py))

Para cada candidato, query directo al Gaia Archive vía TAP (async, sin truncar):
- ¿Cuántas observaciones tiene el target en ±180 días alrededor del encuentro?
- ¿Cuántas son ANTES del encuentro (con 7 días de blackout)?
- ¿Cuántas DESPUÉS?
- Flag `viable_obs` si ≥ 3 antes Y ≥ 3 después (mínimo clásico para fit astrométrico)

**Output**:
- `gaia_observations_check.csv` (todos los candidatos analizados)
- `mass_followup_candidates.csv` (solo los viables)

### Paso 4 — Validación pipeline ([scripts/validate_novel_a.py](../scripts/validate_novel_a.py))

Los 7 encuentros con Ceres/Vesta (masas conocidas) sirven como benchmark:
- MAE(nuestro − JPL) = 3.21e-4 AU << umbral 0.005 AU ✅
- Confirma que la geometría detectada del pipeline es precisa

---

## Top 10 candidatos publicables

| Rank | Perturber (a pesar) | Target (observado) | Fecha | dist (AU) | δ (μas) | obs antes | obs después | mag |
|------|---------------------|---------------------|-------|-----------|---------|-----------|-------------|-----|
| 1 | **(57) Mnemosyne** | 2008_ef40 | 2016-08-26 | 0.00244 | **22,712** | 24 | 49 | 20.0 |
| 4 | (165) Loreley | 1996_tf50 | 2014-12-08 | 0.00254 | 6,633 | 9 | 23 | 20.4 |
| 10 | **(111) Ate** | 2000_nt3 | 2016-06-08 | **0.000472** | 4,906 | 39 | 37 | **18.6** |
| 13 | (124) Alkeste | 2000_qs165 | 2016-04-22 | 0.00329 | 4,411 | 68 | 76 | 19.4 |
| 14 | (113) Amalthea | 2001_vr121 | 2016-11-24 | 0.00910 | 4,341 | 22 | 8 | 19.8 |
| 15 | (46) Hestia | Sitensky | 2016-01-14 | 0.00977 | 3,893 | 17 | 37 | 18.5 |
| 17 | **(241) Germania** | 2000_jc23 | 2016-06-27 | 0.00712 | 3,495 | 29 | **204** | 19.5 |
| 18 | (511) Davida | 2003_sm90 | 2014-11-19 | 0.00903 | 3,250 | 10 | 65 | 19.8 |
| 23 | (206) Hersilia | 1999_vf5 | 2017-02-12 | 0.00818 | 2,341 | 24 | 26 | 19.7 |
| 28 | (19) Fortuna | 2000_ad1 | 2016-07-27 | 0.00406 | 1,921 | 27 | 43 | 19.6 |

### Standouts

- **(57) Mnemosyne**: Señal masiva de 22,712 μas (227× la precisión de Gaia). 73 observaciones bracketing. **Candidato #1**.
- **(111) Ate + 2000_nt3**: El encuentro más cercano del catálogo (0.000472 AU = 70.600 km). Target brillante (mag 18.6) con 76 observaciones. **Mejor SNR esperado.**
- **(241) Germania**: 204 observaciones después del encuentro — cobertura astrométrica excepcional.
- **(83) Beatrix + Cunitza** (rank #41 por δ, pero target mag 16.3 — el más brillante de la lista).

---

## Ranking completo

Los 41 candidatos viables están en `data/output/mass_followup_candidates.csv`,
ordenados por deflexión esperada (descendente).

Para reproducir:

```bash
# Paso 1: filtrar candidatos relevantes
docker compose run --rm pipeline python encounter_analysis/filter_candidates.py

# Paso 2: rankear por señal de deflexión
docker compose run --rm pipeline python -m scripts.analyze_mass_candidates --top-n 100

# Paso 3: verificar contra Gaia DR3
docker compose run --rm pipeline python -m scripts.check_gaia_observations
```

---

## Notas y caveats

1. **Masas estimadas** del perturber se basan en densidad asumida ρ = 1.5 g/cm³.
   Para los asteroides reales la densidad varía entre ~1.0 (C-type carbonáceo) y ~3.5
   (S-type, M-type metálico) g/cm³. La masa real puede variar en factor ~2.

2. **Deflexión** calculada con la fórmula lineal δ = 2GM/(v²b), válida para
   ángulos pequeños. Para los encuentros más cercanos puede subestimar ligeramente
   la deflexión real, pero el orden de magnitud es correcto.

3. **Precisión Gaia 100 μas** es el valor típico per-transit para un MBA. Targets más
   brillantes alcanzan mejor precisión; los más oscuros pueden estar limitados a
   ~500 μas por tránsito. Múltiples tránsitos promediados reducen el ruido como √N.

4. **Sesgo de catálogo**: nuestros encuentros usan propagación Kepler con MPCORB
   actual. Verificamos contra JPL (MAE = 3.2e-4 AU) para Ceres/Vesta. Para los Cat B
   con encuentros muy cercanos (ej. (111) Ate a 70.600 km), una validación N-body
   con `--propagation rebound` daría confirmación adicional.

5. **Validación de masa real**: ninguno de estos 41 candidatos tiene masa publicada,
   pero algunos podrían estar siendo trabajados en pipelines profesionales no
   públicos (ej. Fuentes-Muñoz et al. están analizando Gaia FPR). Cruzar con
   anuncios de masas recientes sería el siguiente paso antes de cualquier publicación.

---

## Próximos pasos posibles

- [ ] Cross-check N-body (rebound) para los top 10 candidatos
- [ ] Estimación de incertidumbre orbital del target durante el encuentro
- [ ] Generar fit astrométrico real con las observaciones Gaia descargadas
- [ ] Comparar con masas publicadas en Fuentes-Muñoz (2024) y otros
- [ ] Análisis específico por subgrupo dinámico (familias colisionales)
