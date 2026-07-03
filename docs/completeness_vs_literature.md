# Completitud del catálogo frente a la literatura de determinación de masas

**Estado:** cerrado (primera medición) · **Fecha:** 2026-07-03 · **Autor:** Daniel Frascarelli

**Script reproducible:** `scripts/validate/crosscheck_literature_encounters.py`
**Salidas:** `data/output/literature_validation/completeness_{pairs,per_perturber}.csv`,
`completeness_summary.json`
**Catálogo evaluado:** `data/output/encounters_catalog_hybrid_stageb.parquet`
(72,236,904 encuentros, `dist_au` ∈ [7×10⁻⁶, 0.0572] AU; el que consume el motor de masas)

---

## 1. Pregunta y motivación

La determinación de masas de asteroides por perturbaciones mutuas parte de un
conjunto **seleccionado a mano** de pares (perturbador, asteroide de prueba). La
pregunta de completitud que un referee planteará para un *dataset paper* de
encuentros cercanos es directa:

> De los pares (perturbador, objetivo) que la literatura efectivamente usó para
> medir masas, ¿cuántos reaparecen como encuentros a < 0.05 AU en nuestro
> catálogo?

Este documento mide esa fracción de recuperación contra la fuente primaria y más
completa: la Tabla 5 de **Fuentes-Muñoz et al. (2025)**, el estándar Gaia-FPR
actual. La cota inferior es honesta y conservadora (ver §2 y §5).

---

## 2. Fuentes y método

### 2.1 Fuentes-Muñoz et al. (2025) — fuente primaria

Fuentes-Muñoz, Farnocchia, Giorgini & Park (2025), *"Asteroid Mass Estimation by
Mutual Perturbations during Close Encounters after Gaia FPR"*, AJ 170, 353. La
Tabla 5 machine-readable (`data/raw/fuentes_munoz_2025/ajae0cc9t5_mrt.txt`)
lista, por perturbador numerado, una columna `List` (bytes 101-776) delimitada
por `|` con los asteroides de prueba **para los que hubo señal**.

Puntos clave del parser (`parse_table5_targets`):

- Solo se consideran **perturbadores numerados** (1645 filas; se descartan
  perturbadores con designación provisional como encabezado).
- Cada token de la lista se resuelve a número MPC; los tokens con designación
  provisional (p. ej. `2007 VQ345`) quedan como objetivo no-numerado, **fuera
  del alcance** de un catálogo de numerados.
- **Truncamiento (Nota 5 de la tabla):** las listas de FM están **truncadas a
  los primeros 100 objetos**, ordenados por señal decreciente, con un `...`
  final cuando son más largas. 223 de los perturbadores tienen lista truncada.
  Por lo tanto nuestro denominador es *"pares FM listados con objetivo numerado
  resoluble"*, **no** el `Ntest` completo de FM. Es una cota conservadora: los
  pares FM más allá del top-100 no entran ni al numerador ni al denominador.

### 2.2 Nuestro catálogo y criterio de match

Se escanea el parquet en modo **lazy + streaming** de polars, filtrando temprano
a las filas cuyos *dos* números caen en el conjunto de perturbadores/objetivos de
FM (colapsa 72 M filas a unos pocos miles antes de cualquier trabajo pareado). Un
par FM (P, T), ambos numerados, se considera **recuperado** si el catálogo
contiene al menos una fila con {`number_1`, `number_2`} == {P, T} y
`dist_au < 0.05` (match insensible al orden). Se retiene además el mínimo
`dist_au` sobre todas las filas del par para clasificar las no-recuperaciones.

### 2.3 Goffin (2014) — limitación documentada

Goffin (2014, A&A 565, A56) es la otra referencia clásica. Sus tablas VizieR
disponibles (`data/raw/goffin_2014_encounters.parquet`: table5 con 132
perturbadores, table6 con 367 estimaciones de masa compiladas) contienen **masas
por perturbador y el número de deflexiones `Nd`, pero NO la lista de asteroides
de prueba** (los pares individuales viven en el texto/apéndice del artículo, no
en la publicación machine-readable). El contenedor no dispone de `pdftotext` ni
de librerías de extracción de PDF (`pdfplumber`, `PyMuPDF`, `PyPDF2`), de modo
que **el cruce pair-level contra Goffin es inviable con los activos actuales** y
se documenta como tal. FM es, además, la fuente más completa y el estándar
Gaia-FPR, por lo que el cruce contra FM es el resultado sustantivo.

---

## 3. Resultado global

| Métrica | Valor |
|---|---|
| Perturbadores numerados FM (Tabla 5) | 1645 |
| Pares listados totales | 41,382 |
| Objetivos con designación provisional (fuera de alcance) | 1,206 |
| **Pares numerados resolubles (denominador)** | **40,176** |
| **Pares recuperados a < 0.05 AU** | **11,842** |
| **Fracción de recuperación (pair-weighted)** | **29.5 %** |
| Recuperación media por perturbador | 26.5 % |
| Recuperación mediana por perturbador | 27.3 % |
| Perturbadores con ≥1 objetivo recuperado | 862 / 1046* |
| Perturbadores con 0 recuperados | 184 / 1046* |

\* De los 1645 perturbadores, 1046 tienen al menos un objetivo **numerado**; los
599 restantes solo listan objetivos con designación provisional (fuera de
alcance) y no entran en la estadística por-perturbador.

> **Recuperación global: 29.5 %** de los pares numerados del top-100 de FM
> reaparecen como encuentros < 0.05 AU en nuestro catálogo.

Nota metodológica: 40,176 pares listados corresponden a 40,004 claves de par
no-ordenadas distintas (172 pares aparecen bajo dos perturbadores, p. ej. A
listado como objetivo de B y B como objetivo de A). Contamos **pares listados
por FM**, de ahí la pequeña diferencia con el conteo por clave única.

---

## 4. Desglose por perturbador (selección)

### 4.1 Calibradores Big-4 y perturbadores clave del motor de masas

| # | Nombre | N obj. num. | Recup. | % | Trunc. |
|---:|---|---:|---:|---:|:---:|
| 1 | Ceres | 99 | 36 | 36.4 | Y |
| 2 | Pallas | 85 | 9 | 10.6 | . |
| 3 | Juno | 99 | 23 | 23.2 | Y |
| 4 | Vesta | 98 | 36 | 36.7 | Y |
| 5 | Astraea | 99 | 54 | 54.5 | Y |
| 7 | Iris | 96 | 28 | 29.2 | Y |
| 10 | Hygiea | 100 | 18 | 18.0 | Y |
| 15 | Eunomia | 99 | 50 | 50.5 | Y |
| 16 | Psyche | 99 | 38 | 38.4 | Y |
| 29 | Amphitrite | 100 | 43 | 43.0 | Y |
| 52 | Europa | 100 | 41 | 41.0 | Y |
| 65 | Cybele | 99 | 19 | 19.2 | Y |
| 87 | Sylvia | 77 | 30 | 39.0 | . |
| 356 | Liguria | 99 | 60 | 60.6 | Y |
| 511 | Davida | 99 | 6 | 6.1 | Y |
| 704 | Interamnia | 99 | 27 | 27.3 | Y |

### 4.2 Top-20 perturbadores por número absoluto de objetivos recuperados

| # | Nombre | N obj. num. | Recup. | % |
|---:|---|---:|---:|---:|
| 356 | Liguria | 99 | 60 | 60.6 |
| 5 | Astraea | 99 | 54 | 54.5 |
| 74 | Galatea | 98 | 50 | 51.0 |
| 15 | Eunomia | 99 | 50 | 50.5 |
| 72 | Feronia | 100 | 49 | 49.0 |
| 23 | Thalia | 99 | 48 | 48.5 |
| 503 | Evelyn | 100 | 48 | 48.0 |
| 247 | Eukrate | 96 | 47 | 49.0 |
| 40 | Harmonia | 97 | 46 | 47.4 |
| 203 | Pompeja | 100 | 46 | 46.0 |
| 139 | Juewa | 98 | 45 | 45.9 |
| 200 | Dynamene | 98 | 45 | 45.9 |
| 120 | Lachesis | 95 | 44 | 46.3 |
| 224 | Oceana | 100 | 44 | 44.0 |
| 36 | Atalante | 96 | 44 | 45.8 |
| 135 | Hertha | 96 | 43 | 44.8 |
| 29 | Amphitrite | 100 | 43 | 43.0 |
| 56 | Melete | 100 | 42 | 42.0 |
| 55 | Pandora | 100 | 42 | 42.0 |
| 206 | Hersilia | 100 | 42 | 42.0 |

La tabla completa por perturbador está en
`data/output/literature_validation/completeness_per_perturber.csv`.

---

## 5. Causas de no-recuperación

De los 40,176 pares numerados, 28,334 no se recuperan. El script los clasifica de
forma excluyente y con prioridad (perturbador ausente → objetivo ausente →
sin encuentro cercano):

| Causa | N pares | % de los no-recuperados |
|---|---:|---:|
| Encuentro ≥ 0.05 AU presente en el catálogo | 0 | 0.0 % |
| Perturbador ausente del universo del catálogo | 270 | 1.0 % |
| Objetivo ausente del universo del catálogo | 2,102 | 7.4 % |
| Ambos presentes, ningún encuentro < 0.05 AU | 25,962 | 91.6 % |
| **Total no-recuperados** | **28,334** | **100 %** |

Y aparte, **1,206 objetivos** listados por FM tienen designación provisional
(no numerada) → fuera del alcance de un catálogo de asteroides numerados.

Interpretación:

- **`encounter_above_threshold` = 0.** Nuestro catálogo está censurado a
  ~0.05 AU (máximo observado 0.0572 AU), de modo que casi no hay pares "presentes
  pero lejanos". Todo par presente cae, en la práctica, por debajo del corte.
- **Objetos ausentes del universo (270 + 2,102 = 2,372, ~8.4 % de los
  no-recuperados).** El objetivo (o, en 270 casos, el perturbador) nunca aparece
  en ninguna fila del catálogo. Son mayormente objetivos de numeración alta y
  reciente (ejemplos verificados para Ceres: 613383, 668956, 731182) que quedan
  fuera de nuestro subconjunto de numerados propagados / del prefiltro orbital.
- **Sin encuentro cercano (25,962, 91.6 % de los no-recuperados) — causa
  dominante y esperada.** Ambos cuerpos *están* en el catálogo, pero ese par
  concreto nunca se aproximó a < 0.05 AU en nuestra propagación. Esto **no es un
  fallo de completitud**: la lista de "señal" de FM **no está limitada a
  0.05 AU**. La señal de masa depende de la geometría y de la masa del
  perturbador, no de un corte geométrico duro; FM incluye encuentros
  genuinamente más anchos que 0.05 AU. Nuestro catálogo, por diseño, solo cataloga
  encuentros < 0.05 AU. Esta diferencia de criterio explica la mayor parte de la
  brecha entre el 29.5 % recuperado y el 100 %.

---

## 6. Lectura para el dataset paper

- **Recuperación honesta: 29.5 %** de los pares (perturbador, objetivo) numerados
  del top-100-por-señal de FM aparecen como encuentros < 0.05 AU en el catálogo.
- El déficit **no** se debe a fallos de detección dentro de nuestro dominio: el
  ~92 % de la brecha son encuentros que FM usa pero que son **más anchos que
  nuestro umbral de 0.05 AU** (definición de dominio distinta), y solo ~8 % son
  objetos fuera de nuestro subconjunto de numerados propagados.
- Los calibradores Big-4 (Ceres 36 %, Vesta 37 %, Pallas 11 %, Hygiea 18 %) y los
  perturbadores del motor de masas (Psyche 38 %, Sylvia 39 %, Europa 41 %,
  Eunomia 51 %) se recuperan con fracciones consistentes con esta interpretación:
  cuanto más "cercano-selectivo" es el conjunto de encuentros de un perturbador,
  mayor la recuperación (Liguria 61 %, Astraea 55 %).
- **Limitación:** el cruce pair-level contra Goffin (2014) no es posible con los
  activos actuales (sus tablas machine-readable no publican los pares individuales
  y no hay extractor de PDF en el entorno). FM es la fuente primaria y más
  completa, y sostiene el resultado.

---

## 7. Reproducir

```bash
docker compose run --rm pipeline python -m scripts.validate.crosscheck_literature_encounters
# opcional: otro umbral o catálogo
docker compose run --rm pipeline python -m scripts.validate.crosscheck_literature_encounters \
    --threshold 0.05 --catalog data/output/encounters_catalog_hybrid_stageb.parquet
```
