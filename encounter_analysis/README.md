# Análisis de Encuentros Novedosos

Scripts para filtrar y analizar los encuentros detectados por el pipeline que no aparecen en ninguna literatura científica publicada.

## Contexto

El pipeline produce ~4M encuentros a 0.05 AU. De esos, ~119.545 no están documentados en Fienga (2003), Galád & Gray (2002) ni Goffin (2014). La mayoría son irrelevantes para astronomía de precisión: cuerpos pequeños, encuentros lejanos, velocidades altas. Este directorio contiene el código para destilar ese conjunto en candidatos útiles.

## Scripts

### `filter_candidates.py`

Lee `data/output/novel_encounters_not_in_literature.csv` y aplica criterios astrofísicos para producir `data/output/relevant_novel_encounters.csv`.

```bash
# Corrida estándar (parámetros por defecto)
docker compose run --rm pipeline python encounter_analysis/filter_candidates.py

# Ajustar umbrales
docker compose run --rm pipeline python encounter_analysis/filter_candidates.py \
    --max-dist-au 0.01 \
    --min-diameter-km 50 \
    --max-vel-km-s 6.0

# Incluir encuentros no observables por Gaia
docker compose run --rm pipeline python encounter_analysis/filter_candidates.py --all-observability
```

## Criterios de filtrado

### 1. Observabilidad Gaia (`gaia_observable == true`)

**Por qué**: Si Gaia no observó al asteroide target cerca de la fecha del encuentro, el encuentro no produce señal astrométrica medible. La condición de observabilidad es:
- Elongación solar > 45° (el objeto no está demasiado cerca del Sol para Gaia)
- Magnitud aparente < 21 (dentro del límite de detección de Gaia)

### 2. Distancia máxima (`dist_au < 0.02`)

**Por qué**: La deflexión gravitacional escala como 1/b (distancia de máximo acercamiento). A 0.05 AU la mayoría de encuentros producen deflexiones de < 1 μas, indetectables. El corte en 0.02 AU (≈ 3 millones de km, ~8× la distancia Tierra–Luna) garantiza que la perturbación gravitacional sea al menos del orden de la precisión de Gaia (~0.1 mas por tránsito).

### 3. Diámetro mínimo del perturber (`diameter_1_km > 30`)

**Por qué**: La deflexión escala como M ∝ D³. Un asteroide de 30 km tiene ~14.000 veces más masa que uno de 3 km. Por debajo de 30 km la masa es tan pequeña que la deflexión es indetectable con cualquier instrumento actual, incluso en un encuentro muy cercano.

### 4. Velocidad relativa máxima (`rel_vel_km_s < 8.0`)

**Por qué**: La deflexión escala como 1/v². A velocidades altas el tiempo de influencia gravitacional es muy corto. Un encuentro a 0.01 AU con v = 15 km/s produce ~4× menos deflexión que el mismo encuentro con v = 7.5 km/s. Encuentros rápidos son poco útiles para determinación de masa.

## Score de deflexión

```
deflection_score = diameter_1_km³ / (dist_au × rel_vel_km_s²)
```

Derivado de la fórmula de deflexión gravitacional:

```
δ ≈ 2GM / (v² · b)
```

donde `M ∝ D³` (densidad constante), `b = dist_au`, `v = rel_vel_km_s`. El score es proporcional a la deflexión esperada hasta una constante que depende de la densidad. Se usa para ordenar los candidatos de mayor a menor utilidad para determinación de masa.

## Categorías de output

El CSV de salida incluye una columna `mass_unknown`:

- **`mass_unknown = false` (Categoría A)**: el perturber es un cuerpo con masa publicada (Ceres, Vesta, etc.). Útiles como **benchmark de precisión del pipeline** — la deflexión predicha con la masa conocida puede compararse contra la detectada por Gaia.

- **`mass_unknown = true` (Categoría B)**: el perturber no tiene masa publicada. Son los **candidatos genuinos a nueva determinación de masa** — si la deflexión del target es detectable en los datos Gaia, podría derivarse la masa del perturber.

## Columnas del output

| Columna | Descripción |
|---------|-------------|
| `number_1`, `designation_1` | Perturber (cuerpo más grande / más masivo) |
| `diameter_1_km` | Diámetro estimado del perturber (km) |
| `number_2`, `designation_2` | Target (cuerpo perturbado) |
| `diameter_2_km` | Diámetro estimado del target (km) |
| `date_utc` | Fecha del máximo acercamiento (UTC) |
| `dist_au`, `dist_km` | Distancia mínima |
| `rel_vel_km_s` | Velocidad relativa en el encuentro |
| `class_1`, `class_2` | Clasificación orbital (MBA, NEA, etc.) |
| `solar_elongation_deg` | Elongación solar en la fecha del encuentro |
| `gaia_observable` | Si cumple criterios de observabilidad Gaia |
| `deflection_score` | Score de priorización (ver fórmula arriba) |
| `mass_unknown` | `true` = candidato a nueva determinación de masa |
