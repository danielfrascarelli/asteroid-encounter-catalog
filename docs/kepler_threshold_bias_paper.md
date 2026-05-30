# El sesgo del threshold Kepler cerca de 0.05 AU — nota técnica (borrador)

> Borrador de **Track B Stage 3** del [FOLLOWUP_PLAN.md](../FOLLOWUP_PLAN.md).
> Caracteriza cuantitativamente los 25,283 cruces del threshold Kepler→N-body
> reportados en [FROZEN_RUN.md](../FROZEN_RUN.md) y **corrige su interpretación**.
> Análisis reproducible:
> [scripts/validate/analyze_kepler_threshold_bias.py](../scripts/validate/analyze_kepler_threshold_bias.py),
> figuras en [notebooks/kepler_threshold_bias_analysis.ipynb](../notebooks/kepler_threshold_bias_analysis.ipynb).

## Resumen

El catálogo de detección se construye con un refinamiento Kepler de dos cuerpos
y un umbral de 0.05 AU. Un subconjunto (8.73 M pares, seleccionado por criterio
orbital `q_min < 1.8 AU OR e_max > 0.3`) se re-refinó con N-body. De esos,
**25,283 (0.29 %) cruzan el umbral** al pasar a N-body: Kepler los puso <0.05 AU,
N-body los recomputa ≥0.05 AU. **Cero** cruzan en la dirección opuesta.

`FROZEN_RUN.md` leyó esto como *"un pequeño sesgo sistemático de
sobre-detección de Kepler cerca del threshold, sin falsos negativos"*. **Esa
lectura es engañosa en dos puntos**, y corregirla es la contribución de esta nota:

1. **La asimetría "0 cruces hacia abajo" es censura, no medición.** El catálogo
   contiene *únicamente* pares con Kepler <0.05 AU (es el umbral de detección;
   ver tabla de FROZEN_RUN: 100 % <0.05). Los pares que Kepler colocó *por
   encima* de 0.05 AU nunca se escribieron, así que es **imposible** observarlos
   cruzando hacia abajo. La tasa de falsos negativos **no es 0** — es
   **no medible** con este catálogo.
2. **La corrección N-body no es sistemáticamente hacia arriba.** El Δdist =
   d_Nbody − d_Kepler tiene **media −1.3×10⁻⁵ AU y mediana −1.0×10⁻⁶ AU**
   (ligeramente *negativa*: N-body tiende a *acercar*), con **std 4.3×10⁻⁴ AU**.
   Es decir: la corrección está **dominada por scatter simétrico**, no por un
   desplazamiento unidireccional. Los cruces hacia arriba ocurren porque el
   scatter empuja por encima de 0.05 a los pares que estaban *justo por debajo*;
   sus espejos (que cruzarían hacia abajo) caen fuera del catálogo.

## Datos y método

- **Catálogo**: `data/output/encounters_catalog_hybrid_stageb.parquet`, columnas
  `dist_au_kepler`, `dist_au_nbody`, `rel_vel_au_day`, `refinement_method`.
- **Filas refinadas N-body**: 8,728,509 (todas con Kepler <0.05 AU).
- **Elementos orbitales** para el corte por banda: cache MPCORB
  `data/cache/nbody_validation/mpcorb_stageb_elements.parquet`.
- Todo se agrega con DuckDB (no se materializa el catálogo en memoria).

## Resultados

### 1. Tasa de cruce concentrada en el último bin

Casi todos los cruces ocurren en los 0.005 AU justo debajo del umbral:

| banda d_Kepler (AU) | n pares | media Δdist (AU) | tasa de cruce↑ |
|---|---:|---:|---:|
| [0.000, 0.040) | 5.51 M | −3 a −15 ×10⁻⁶ | 0.000 % |
| [0.040, 0.045) | 1.49 M | −1.7×10⁻⁵ | 0.0039 % |
| [0.045, 0.050) | 1.67 M | −1.8×10⁻⁵ | **1.51 %** |

Físicamente esperado: sólo un par cuya distancia Kepler ya está dentro de
~1 std (4×10⁻⁴ AU) del umbral puede ser empujado al otro lado por la corrección.

### 2. Δdist: scatter simétrico, no sesgo unidireccional

Globalmente Δdist tiene media ≈ mediana ≈ 0 a nivel de 10⁻⁵–10⁻⁶ AU, frente a
una std de 4.3×10⁻⁴ AU (≈ 65,000 km). La media *ligeramente negativa* implica
que, si hay un sesgo neto, Kepler **sub-estima** levemente la proximidad — lo
contrario de "sobre-detección". La sobre-detección observada es un efecto de
**borde + scatter**, no un sesgo físico direccional.

### 3. Dependencia con velocidad relativa

La tasa de cruce crece monótonamente con v_rel (quintiles de Kepler<0.05):

| quintil v_rel | rango (AU/d) | tasa de cruce↑ |
|---|---|---:|
| q1 (lento) | [0.000, 0.0020] | 0.150 % |
| q3 | [0.0029, 0.0037] | 0.314 % |
| q5 (rápido) | [0.0047, 0.039] | 0.403 % |

Encuentros rápidos son más sensibles a la perturbación planetaria que Kepler
omite: un pequeño cambio de geometría/timing desplaza más el mínimo aparente.

### 4. Banda orbital de los cruces

De los 25,283 cruces, **22,979 (90.9 %)** tienen `q_min < 1.8 AU` (mediana
q_min = 1.75 AU); **8,220 (32.5 %)** tienen `e_max > 0.3` (mediana e_max =
0.256). Los cruces se concentran en la población de **perihelio bajo / e alta**,
exactamente donde la aproximación de dos cuerpos es menos fiable (perturbación
solar+planetaria fuerte cerca del perihelio). Esto **valida a posteriori** el
criterio de selección del subset N-body (`q_min<1.8 OR e_max>0.3`).

## Implicaciones para usuarios del catálogo

1. **Cerca del threshold (0.045–0.050 AU), ~1.5 % de las detecciones Kepler son
   espurias** (N-body las saca de 0.05 AU). Para ciencia sensible a la
   completitud/pureza cerca del corte, usar siempre el catálogo híbrido N-body.
2. **La tasa de falsos negativos es comparable, no nula.** Por simetría del
   scatter, debería haber un número similar de pares que Kepler puso *justo por
   encima* de 0.05 y que N-body bajaría — pero el catálogo los excluye. Reclamar
   "sin falsos negativos" es incorrecto.
3. **El sesgo es peor en órbitas de perihelio bajo / alta excentricidad y en
   encuentros rápidos.** Un usuario interesado en NEAs / Marte-crossers debe
   tratar las detecciones Kepler cerca del threshold con especial cautela.

## Experimento para medir la tasa de falsos negativos (trabajo futuro)

La cantidad faltante —cuántos encuentros reales (<0.05 AU N-body) pierde el
prefiltro Kepler— requiere re-refinar con N-body una muestra de pares con
**d_Kepler ∈ [0.05, ~0.06] AU**, que el catálogo actual descarta. Plan:

1. Re-correr la detección con threshold 0.06 AU (sólo el prefiltro Kepler) sobre
   una muestra orbitalmente representativa.
2. N-body-refinar los pares en [0.05, 0.06] y contar cuántos bajan de 0.05.
3. Combinar con la tasa de sobre-detección de esta nota para una matriz de
   confusión completa del prefiltro cerca del threshold.

## ¿Publicar?

Como **nota técnica / apéndice metodológico** del catálogo: sí, aporta una
caracterización honesta del prefiltro. Como paper independiente: marginal —
el resultado central (scatter de la corrección de dos cuerpos cerca de un corte
de distancia produce sobre-detección de borde) es esperable y probablemente ya
discutido en la literatura de catálogos de close-approach (revisar Fienga,
JPL CNEOS). **Decisión pendiente** (humana): apéndice del catálogo vs. nota
técnica standalone. Recomendación: apéndice, con el experimento de falsos
negativos como motivación para DR4.

## Reproducir

```bash
docker compose run --rm pipeline python -m scripts.validate.analyze_kepler_threshold_bias \
  --out-prefix data/output/kepler_bias/threshold
# -> threshold_summary.json, threshold_bands.csv
```
