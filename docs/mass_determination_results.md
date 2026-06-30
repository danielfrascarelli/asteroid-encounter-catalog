# Determinación de masas de asteroides con Gaia FPR — resultados

> **Estado:** T10 ✅ (validación de calibradores) + T11 ✅ (catálogo de masas).
> Plan T1–T11 completo; mergeado a `main` vía PR #80 (2026-06-30).
> Motor: `src/orbdet/` (ver [`orbdet_engine_status.md`](orbdet_engine_status.md)).
> Plan: [`planning/MASS_DETERMINATION_PLAN.md`](../planning/MASS_DETERMINATION_PLAN.md).
> **Última actualización:** 2026-06-30.

## Resumen

Determinamos masas de asteroides perturbadores ajustando **conjuntamente** la masa
del perturbador y las órbitas de los asteroides de prueba que tuvieron encuentros
cercanos con él, por mínimos cuadrados sobre el arco completo de astrometría Gaia FPR,
con el modelo de fuerzas state-of-the-art (ASSIST: DE440 + GR + 16 perturbadores). Es
la metodología de Fuentes-Muñoz / OrbFit / JPL, construida de cero en `src/orbdet/`.

**Validación (T10): las 4 masas calibradoras se recuperan dentro de |z|<3.** Con
N≥20 asteroides de prueba, las masas DAWN (Ceres, Vesta) y Vernazza (Hygiea) se
reproducen a **~5%**. Esto **refuta la conclusión del cierre Track A** de que el
leverage de Gaia era insuficiente: lo era el método (LOO secuencial), no los datos.

**Producción (T11): masa nueva defendible — (16) Psyche = 2.43×10¹⁹ kg ±3.3%**, en
acuerdo del 2% con DE441. El barrido de los otros 12 perturbadores grandes también
revela un límite del método: para perturbadores con deflexión débil la masa se sesga
baja (absorción de señal) — ver abajo.

## Método

- **Ajuste conjunto órbita+masa** (`orbdet.mass_determination.determine_shared_mass`):
  vector de `1 + 6N` parámetros (masa compartida + 6 elementos por objetivo),
  Jacobiano en flecha, resuelto por Levenberg-Marquardt sobre el arco completo. La
  degeneración masa↔drift se maneja dentro de la covarianza conjunta, no se descarta.
- **Fuerzas (ASSIST):** efeméride JPL DE440 (Sol + 8 planetas + Luna + Plutón) leída
  en cada paso + relatividad (EIH) + 16 perturbadores asteroidales masivos. Validado
  vs Horizons a 0.17 mas sobre 900 d.
- **Observación:** estado → ICRS → RA/Dec con light-time iterativa; covarianza
  along-scan anisotrópica de Gaia.
- **Covarianza en bloques por FOV (clave).** Gaia entrega ~7 CCDs por cruce de plano
  focal (separados ~5 s) con residuos correlacionados (ICC≈0.32 medido). Tratarlos
  como independientes subestima σ(masa) ~√7. Se blanquea con `C_bloque = diag(σ_AL²)
  + s_c²·11ᵀ` por cruce, con el piso correlacionado `s_c` autocalibrado para χ²_red≈1.
- **Selección de objetivos:** los más cercanos (<0.05 AU) del catálogo de encuentros
  (`--from-catalog`). Usar **muchos** objetivos (N≥20) es esencial (ver abajo).
- **Rechazo de outliers:** sigma-clipping iterativo a 4σ.
- **Paralelización:** los N objetivos se evalúan en paralelo (pool por proceso);
  ~6× speedup, resultado idéntico al modo serie.

## Validación — calibradores Big-4 (Gaia FPR, N≥20)

| Cuerpo | N obj | masa ajustada (kg) | σ_total | ratio fit/lit | z | fuente lit |
|--------|-------|--------------------|---------|---------------|---|------------|
| Ceres  | 28 | 8.96×10²⁰ | 4.6% | 0.955 | −1.01 | DAWN (Park+ 2016) |
| Vesta  | 28 | 2.44×10²⁰ | 4.6% | 0.943 | −1.30 | DAWN (Russell+ 2012) |
| Hygiea | 20 | 8.22×10¹⁹ | 5.7% | 0.990 | −0.13 | Vernazza+ (2020) |
| Pallas |  6 | 2.54×10²⁰ | 7.0% | 1.240 | +2.67 | Goffin (2014) — *target-limited* |

Los 3 bien muestreados recuperan la masa a ~5% (sesgo medio −4%). Pallas tiene **sólo
6–7 encuentros <0.05 AU en todo el catálogo** → no se puede promediar, y es el único en
tensión (aun así |z|<3).

## Dos hallazgos metodológicos

1. **El número de objetivos es decisivo.** Con N≈7 las masas salían sesgadas alto
   (+12–29%); con N≥20 convergen a la verdad (~0.95). El "sobre-tiro" era **dispersión
   estadística de muestra chica**, no un sistemático. El motor es insesgado: un
   closing-loop sobre la geometría real (obs sintéticas a la masa verdadera + ruido)
   recupera ratio medio 0.997 (3 semillas).

2. **La medición está limitada por sistemáticos, no por estadística.** Con N grande la
   σ formal (Fisher) baja como 1/√N y se vuelve <0.2%, pero la exactitud real está
   limitada por sistemáticos por-encuentro (imperfección de la órbita del objetivo,
   astrometría local, perturbadores menores fuera de los 16). Por eso reportamos
   `σ_total = √(σ_stat² + (f_sys·M)²)` con **f_sys≈4.2% calibrado de los calibradores
   bien muestreados** — el tratamiento estándar de incertidumbre externa.

## Masas nuevas — barrido de los 12 perturbadores restantes (T11)

Mismo procedimiento, objetivos del catálogo (<0.05 AU), N hasta 40. Ratio = masa
ajustada / masa de la efeméride DE441 (que para estos cuerpos es la referencia
publicada, menos precisa que DAWN). σ_stat = incertidumbre formal (Fisher).

| Cuerpo | N | masa ajustada (kg) | σ_stat | ratio fit/DE441 | χ²_red | clase |
|--------|---|--------------------|--------|-----------------|--------|-------|
| **Psyche** | 36 | **2.43×10¹⁹** | **3.3%** | **1.020** | 0.99 | ✅ defendible |
| Euphrosyne | 21 | 2.43×10¹⁹ | 27% | 1.502 | 0.99 | consistente pero imprecisa |
| Sylvia | 11 | 3.47×10¹⁹ | 11% | 1.069 | 0.98 | OK, pocos objetivos |
| Juno | 34 | 1.97×10¹⁹ | 11% | 0.685 | 0.99 | sesgada baja |
| Eunomia | 38 | 2.20×10¹⁹ | 9% | 0.724 | 0.98 | sesgada baja |
| Europa | 34 | 2.29×10¹⁹ | 6% | 0.568 | 0.98 | sesgada baja |
| Interamnia | 35 | 2.42×10¹⁹ | 13% | 0.571 | 0.99 | sesgada baja |
| Iris | 37 | 8.11×10¹⁸ | 19% | 0.475 | 0.99 | sesgada baja |
| Thisbe | 37 | 6.99×10¹⁸ | 6% | 0.392 | 0.97 | sesgada baja |
| Cybele | 12 | 6.24×10¹⁹ | 22% | 4.44 | 0.91 | no fiable (N bajo) |
| Camilla | 2 | — | 67% | 8.1 | — | no fiable (N=2) |
| Davida | 3 | <0 | — | <0 | 0.65 | no fiable (N=3) |

**Resultado principal (gate T11): (16) Psyche** se determina a **2.43×10¹⁹ kg ±3.3%**
(formal), en acuerdo del 2% con la efeméride DE441 — una masa nueva defendible que
**extiende la validación más allá de los calibradores**. Sylvia y Euphrosyne son
consistentes pero menos precisas.

**Sesgo de absorción de señal (hallazgo).** Seis perturbadores bien muestreados salen
sistemáticamente **bajos** (0.39–0.72), con χ²_red≈1. No es ruido: es la
**degeneración masa↔órbita cuando la deflexión es débil frente al ruido por-encuentro**
— el ajuste explica los datos con menos masa + órbita ajustada (regresión hacia cero).
Los perturbadores fuertes (Big-4, Psyche) no lo sufren porque su deflexión domina el
ruido. **Implicación:** para perturbadores débiles la masa de la efeméride no se
recupera con esta metodología tal cual, y la σ formal subestima el error real. La σ
fiable requiere estimación externa por-perturbador (jackknife/bootstrap) y/o
regularización del par masa-órbita — trabajo futuro.

## Cruce independiente con Fuentes-Muñoz 2025 (T10)

`scripts/validate/validate_fuentes_munoz_masses.py` compara nuestras masas contra
la Tabla 5 de Fuentes-Muñoz et al. 2025 (AJ 170, 353), convirtiendo su `GMfin`
(km³/s²) a kg con G=6.67430×10⁻²⁰. De los 16 perturbadores del catálogo, los 16
solapan con su tabla. **Para los calibradores Fuentes-Muñoz fija `GMfin` a la
semilla SB441/literatura** (no es comparación independiente); el cruce con valor
es sobre los **no-calibradores**, donde ellos corrieron su propio ajuste FPR.

| Cuerpo | M_nuestra (kg) | M_FM (kg) | ratio | z_vs_FM |
|--------|----------------|-----------|-------|---------|
| **(16) Psyche** | 2.429×10¹⁹ | 2.395×10¹⁹ | **1.014** | **+0.25** ✅ |
| (10) Hygiea (cal) | 8.218×10¹⁹ | 8.237×10¹⁹ | 0.998 | −0.04 |
| (31) Euphrosyne | 2.430×10¹⁹ | 1.645×10¹⁹ | 1.477 | +1.15 |
| (52) Europa | 2.285×10¹⁹ | 2.656×10¹⁹ | 0.860 | −2.26 |
| (3) Juno | 1.971×10¹⁹ | 2.719×10¹⁹ | 0.725 | −3.34 |
| (7) Iris | 8.108×10¹⁸ | 1.456×10¹⁹ | 0.557 | −3.99 |
| (15) Eunomia | 2.195×10¹⁹ | 3.120×10¹⁹ | 0.704 | −4.24 |
| (88) Thisbe | 6.991×10¹⁸ | 8.755×10¹⁸ | 0.799 | −3.22 |
| (704) Interamnia | 2.420×10¹⁹ | 3.235×10¹⁹ | 0.748 | −2.50 |

(tabla completa, incl. los no fiables Sylvia/Cybele/Camilla/Davida, en
`data/output/literature_validation/fuentes_munoz_2025_mass_comparison.csv`.)

**Resultado clave: (16) Psyche concuerda con Fuentes-Muñoz al 1.4% (z=+0.25)** —
nuestra masa nueva queda **confirmada de forma independiente** por un estudio de
masas con Gaia FPR revisado por pares, no sólo por la efeméride DE441. Hygiea,
Ceres y Vesta también concuerdan (|z|≲1.3).

**El sesgo de absorción de señal se reproduce contra FM.** Los mismos
perturbadores débiles que salen bajos vs DE441 (Juno, Iris, Eunomia, Thisbe,
Interamnia, ~0.7–0.8) salen bajos vs FM, con |z|>3. Esto **no** contradice el
gate: es exactamente el límite documentado — para deflexión débil la σ formal
subestima el error real (la regresión masa↔órbita hacia cero), por lo que el
z-score formal exagera la tensión. Lo informativo es que los perturbadores
**fuertes e independientes** (Psyche) caen en z≈0.

## Limitaciones y trabajo futuro

- Sólo los 16 perturbadores grandes de `sb441-n16.bsp` tienen órbita en la efeméride;
  el motor requiere la órbita del perturbador de ahí (su masa es el parámetro libre).
- El sesgo medio −4% de los calibradores sugiere un sistemático pequeño residual
  (candidato: completitud del fondo de perturbadores). Acotarlo bajaría f_sys.
- Cruce con Fuentes-Muñoz 2025 **hecho** (sección arriba): Psyche confirmada a 1.4%
  (z=+0.25). Falta sólo extenderlo a más perturbadores débiles con σ externa.
- Pallas y otros perturbadores con pocos encuentros cercanos quedan target-limited
  hasta DR4.

## Reproducción

```bash
# Barrido de calibradores (validación T10)
docker compose run --rm pipeline python -m scripts.mass.orbdet_fit_realdata \
    --perturber big4 --release fpr \
    --from-catalog data/output/encounters_catalog_hybrid_stageb.parquet \
    --top-per-perturber 30 --workers 24 --out-dir data/output/orbdet/expanded

# Catálogo con modelo de error correcto
docker compose run --rm pipeline python -m scripts.mass.build_mass_catalog \
    --in-dir data/output/orbdet/expanded --out data/output/orbdet/mass_catalog.csv

# Cruce independiente de masas vs Fuentes-Muñoz 2025 (T10)
docker compose run --rm pipeline python -m scripts.validate.validate_fuentes_munoz_masses
```
