# El sesgo del threshold Kepler cerca de 0.05 AU — apéndice de métodos del catálogo

> **Track B Stage 3** del follow-up post-deepwork (disuelto; ver [ROADMAP.md](../ROADMAP.md) § "Estado actual"). **Decisión
> (2026-05-30): este material va como apéndice de métodos/caveats del paper del
> catálogo, NO como nota standalone** (ver "Decisión de publicación" al final).
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

## Tasa de falsos negativos — MEDIDA (2026-05-31, Track C2)

La cantidad que faltaba —cuántos encuentros reales (<0.05 AU N-body) censura el
catálogo Kepler— ahora está **medida**
([scripts/validate/measure_threshold_false_negatives.py](../scripts/validate/measure_threshold_false_negatives.py),
artefactos en `data/output/kepler_false_negatives/`). Método: detección Kepler a
threshold **0.06 AU** sobre **10.000** cuerpos numerados (a∈[1.5,4.0], muestra
seeded del snapshot congelado), aislando los pares con `d_Kepler ∈ [0.05, 0.06)`
(los que el catálogo de 0.05 descarta), y re-refinando cada uno bajo N-body
completo (±12 h, IAS15-grade, Sol+Júpiter+Saturno+4 mayores).

Resultado sobre **17.469 pares de banda** (0 fallidos):

- **122 cruzan hacia abajo** (N-body < 0.05 AU) → **tasa de falsos negativos en
  la banda = 0.70 %** [IC95 0.59 %, 0.83 %].
- 48 de esos 122 cruces son `near_boundary` (mínimo posiblemente fuera de la
  ventana ±12 h); excluyéndolos, un piso conservador es 74/17.469 = **0.42 %**.
  La tasa real cae en **~0.4–0.7 %** de la banda.
- `Δdist = d_Nbody − d_Kepler`: mediana **−3×10⁻⁷ AU**, media −2×10⁻⁵, **σ =
  3.7×10⁻⁴ AU** — scatter simétrico dominante, consistente con la §2 (el efecto
  es ruido cerca del corte, NO un sesgo unidireccional, en **ambos** lados).

**Matriz de confusión cerca de 0.05 AU (cerrada):** los cruces son simétricos y
scatter-dominados. Hacia arriba: ~1.5 % en el bin de borde [0.045,0.05) (§1).
Hacia abajo: ~0.4–0.7 % de [0.05,0.06). El catálogo, al contener sólo
Kepler<0.05, **observa los cruces hacia arriba pero censura los de abajo** — y
ahora sabemos que esos no observados pesan ~0.4–0.7 % de la banda adyacente.

**Extrapolación catalog-wide** (orden de magnitud): la banda [0.05,0.06) escala
∝N²; 17.469 pares en 10k cuerpos ⇒ ~3.5×10⁷ en los 449k numerados ⇒ **~1.5–2.5×10⁵
encuentros reales <0.05 AU censurados** por el corte Kepler (independiente del
déficit de recall del prefiltro, [prefilter_recall.md](prefilter_recall.md)).
Caveat: sujeto al near_boundary y a que estos pares además deben pasar el prefiltro.

## Decisión de publicación

### Actualización 2026-05-31 (tras medir los falsos negativos, Track C2)

El **gatillo** que la decisión original dejó abierto —medir la tasa de falsos
negativos y cerrar la matriz de confusión— **se cumplió** (sección anterior). El
número genuinamente novel ya no está sin medir: la matriz cerca de 0.05 AU está
completa (cruces simétricos ~1.5 % arriba / ~0.4–0.7 % abajo, scatter-dominados,
σ=3.7×10⁻⁴). Eso **vuelve viable una nota técnica standalone** (p. ej. RNAAS):
"crossing rates simétricos Kepler↔N-body en un umbral de distancia, medidos sobre
un catálogo Gaia DR3 de 72 M encuentros". Sigue siendo física de bajo perfil
(sesgo de selección genérico), así que la recomendación es **opcional**: si se
quiere un output de publicación de bajo costo, el material ya está; si no, queda
como apéndice de métodos del paper del catálogo. **Decisión de escribirla o no:
del autor** (es una llamada de alcance de publicación, no técnica).

### Decisión original (2026-05-30)

**Decidido: apéndice de métodos/caveats del paper del catálogo, NO nota
standalone.** Razones científicas:

1. **No es física novel.** El esquema Kepler→N-body de dos etapas es práctica
   estándar, y el efecto —errores de dos cuerpos cerca de un corte de distancia
   producen cruces censurados/asimétricos— es un caso genérico de sesgo de
   selección en un umbral (familia Eddington/Malmquist). Una búsqueda de
   literatura no encontró una caracterización prominente *específica* del efecto,
   pero el mecanismo es esperable.
2. **El valor es específico de este catálogo**: cuantifica la pureza cerca del
   corte (~1.5 % espurio en el bin de borde) y dónde (q bajo / e alta / rápidos).
   Eso es exactamente lo que documenta una sección de métodos/caveats — sirve
   directamente a los usuarios del catálogo.
3. **Como standalone sería débil**: el número genuinamente novel —la tasa de
   falsos negativos— está **sin medir** (requiere re-refinar pares con Kepler ∈
   [0.05, 0.06] AU, ausentes del catálogo). Un referee lo vería incompleto.

**Acción tomada**: corregida la afirmación errónea "no false negatives" en
[FROZEN_RUN.md](../FROZEN_RUN.md) (era censura, no medición). Este doc + el
notebook quedan como fuente del apéndice cuando se escriba el paper del catálogo.

**Gatillo para reconsiderar standalone**: si en el futuro se corre el experimento
de falsos negativos (re-refinar [0.05, 0.06] AU) y se obtiene una matriz de
confusión completa del prefiltro, el conjunto sí justificaría una nota técnica
(p. ej. RNAAS).

## Reproducir

```bash
docker compose run --rm pipeline python -m scripts.validate.analyze_kepler_threshold_bias \
  --out-prefix data/output/kepler_bias/threshold
# -> threshold_summary.json, threshold_bands.csv
```
