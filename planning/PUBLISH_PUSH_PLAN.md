# Empuje a publicación — plan maestro multi-frente

> **Estado:** 🟡 ACTIVO
> **Última actualización:** 2026-07-02
> Plan para llevar el proyecto de "resultados validados" a "material publicable".
> Criterio de éxito: identificar y documentar lo genuinamente nuevo (dataset +
> methods), y producir al menos una masa que no exista en la literatura (F4).

## Diagnóstico de partida (por qué este plan)

Las masas medidas hasta hoy son **validación** (calibradores conocidos, semillas
DE441, Psyche ya cubierto por Fuentes-Muñoz 2025), no descubrimiento. Lo
genuinamente nuevo que ya tenemos:

1. **Catálogo de encuentros 3D con presupuesto de completitud medido** (censura
   ~0.70 %, recall prefiltro 76 % en el tail adverso). Nadie publica un catálogo
   de encuentros con su propio budget de falsos negativos cuantificado.
2. **Marco metodológico de masas**: covarianza en bloques por FOV (ICC=0.32
   medido, factor 1.66 en σ), criterio de identificabilidad (measured /
   not_identifiable / cota) + σ jackknife.

Lo que falta para un "titular de descubrimiento": (a) una masa nueva de verdad
(F4, perturbador fuera de los 16), y (b) rastrear el catálogo por eventos únicos.

## Tabla de estado

| # | Frente | tipo | estado | gate |
|---|--------|------|--------|------|
| P1 | Minería de eventos únicos en el catálogo | análisis | ✅ | `docs/notable_encounters.md`. Sin grande-grande D≳100km nuevo (2 eventos, ambos conocidos). Mejor candidato masa fuera de los 16: (9) Metis |
| P2 | F4 — masa nueva (perturbador fuera de los 16) | código + run | ✅ | rama custom OK. (19) Fortuna: 1.13e19±2.2e18 kg (z=+1.25 vs FM), χ²_red=0.977. (9) Metis: 4.74e18±1.53e18 kg (ratio 0.73, z=−1.14 vs FM), χ²_red=0.981. Gate metodológico VERDE — pero AMBAS ya están en Goffin+FM, NO son masas nuevas; σ del 20-32%, medidas coarse |
| P3 | F5 — cruce Fuentes-Muñoz con σ jackknife | análisis | ✅ | `docs/mass_crosscheck_jack.md`. 10/10 medidas en \|z\|<3 con σ jack (vs 5/10 formal) |
| P4 | Dataset paper — draft + consolidación | escritura | 🟡 | scaffold en `docs/dataset_paper_draft.md`; §4 pendiente de llenar con P1 |
| P5 | F3 — sesgo −4 % (fondo de perturbadores) | run | ✅ | HIPÓTESIS REFUTADA: fondo 16→35 (20 cuerpos FM 2025) mueve masas <0.25%; f_sys 4.16%→4.26%. El −4% NO es incompletitud del fondo. Negativo limpio |
| P6 | Merge PRs #83/#84 + higiene de branches | git | ⬜ | (requiere confirmación del usuario) |

## Conclusión sobre novedad (2026-07-03)

El gate F4 quedó verde **metodológicamente** (el motor ajusta un cuerpo fuera de
los 16 de la efeméride, χ²_red=0.977). Pero Fortuna ya tiene masa en Goffin (2014)
y Fuentes-Muñoz (2025) → **no es una masa nueva**. FM 2025 publicó 231 masas, que
cubren esencialmente todos los MBA grandes clásicos, así que "una masa que nadie
tenía" es improbable por la vía de perturbadores grandes. La minería (P1) tampoco
halló un encuentro grande-grande nuevo llamativo. **La novedad real y defendible
es el dataset (catálogo con completitud medida) + el marco metodológico** (cov en
bloques por FOV, σ jackknife, identificabilidad), NO un descubrimiento puntual.
Nota técnica: la masa de Fortuna sale ~32 % alta vs literatura, consistente dentro
de la σ jackknife (z=+1.25) pero con σ del 20 % — Fortuna no es una medida fina.

## Ejecución en paralelo

Frentes con archivos disjuntos → se implementan en paralelo sobre el working tree:

- **P1** escribe `scripts/bench/mine_notable_encounters.py` + `docs/notable_encounters.md`.
- **P2** edita `scripts/mass/orbdet_fit_realdata.py` + `scripts/mass/find_extra_perturber_candidates.py`.
- **P3** escribe `scripts/mass/crosscheck_fuentes_munoz_jack.py` + `docs/mass_crosscheck_jack.md`.

Frentes que requieren runs largos (Docker + Gaia TAP, horas) se secuencian tras
revisar el código: **P2 run**, **P5**.

## Detalle por frente

### P1 — Minería de eventos únicos
Fuente: `data/output/encounters_characterized_full.parquet` (diámetros, clases,
H, distancias, v_rel). Buscar:
- Encuentros grande-grande (ambos D ≳ 50 km) — raros, alto interés.
- Pares en la misma región dinámica (proximidad en a/e/i como proxy de familia).
- Extremos: mínima distancia absoluta, máxima/mínima v_rel, encuentros lentos
  (candidatos naturales a determinación de masa futura).
- Cruce contra pares ya reportados (Goffin 2014, Fuentes-Muñoz) para separar
  lo nuevo de lo conocido.

### P2 — F4 masa nueva
Diseño completo en [`docs/mass_layer_f4_design.md`](../docs/mass_layer_f4_design.md).
Rama custom en la capa IO (no toca `src/orbdet/`): órbita del 17º desde Horizons,
masa-semilla por H, fondo = los 16 completos, N_active=17. Candidatos: (24) Themis,
(532) Herculina, (29) Amphitrite, (354) Eleonora.

### P3 — F5 cruce con σ jackknife
`mass_catalog_jack.csv` ya trae σ_jack y z_total. Producir la tabla de cruce
Fuentes-Muñoz 2025 para los 12 no-calibradores con la σ externa; contar cuántos
quedan en |z|<3 vs los 4/8 con σ formal.

### P4 — Dataset paper
Outline + secciones: datos (Gaia DR3 window), método de detección (KD-tree +
refinamiento N-cuerpos), presupuesto de completitud (censura + recall medidos),
tabla de eventos notables (de P1). Es lo más maduro para publicar ya.

### P5 — F3 sesgo −4 %
Extender el fondo de perturbadores más allá de los 16 y medir el cambio en f_sys
sobre los calibradores. Run largo; se secuencia tras P2.

### P6 — Git
2 PRs abiertas: #83 (jackknife+flag), #84 (Goffin). Mergear y ordenar branches
**requiere confirmación explícita del usuario** (convención del proyecto).
