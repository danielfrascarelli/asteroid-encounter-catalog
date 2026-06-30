# Determinación de masas de asteroides con Gaia FPR — resultados

> **Estado:** completo. Motor `src/orbdet/`; arquitectura y gates de aceptación en
> [`orbdet_engine_status.md`](orbdet_engine_status.md). Ítems abiertos y mejoras en
> [`planning/MASS_FUTURE_WORK.md`](../planning/MASS_FUTURE_WORK.md).
> **Última actualización:** 2026-06-30.

## Resumen

Se determinan masas de asteroides perturbadores por ajuste simultáneo de la masa del
perturbador y las órbitas de los asteroides de prueba con encuentros < 0.05 AU, por
mínimos cuadrados sobre el arco completo de astrometría Gaia FPR. El modelo de fuerzas
es ASSIST (efeméride DE440; correcciones relativistas EIH; 16 perturbadores
asteroidales masivos). La implementación es propia (`src/orbdet/`), sin dependencia de
software de determinación de órbitas de terceros.

Resultados medidos:

- Las 4 masas calibradoras (Ceres, Vesta, Hygiea, Pallas) se recuperan con |z| < 3
  respecto a la literatura.
- Con N ≥ 20 objetivos, el ratio ajuste/literatura de Ceres, Vesta (DAWN) e Hygiea
  (Vernazza 2020) cae en [0.943, 0.990]; la diferencia con la literatura es 1.0–5.7 %.
- (16) Psyche: M = 2.43 × 10¹⁹ kg, σ_stat = 3.3 %. Ratio respecto a DE441 = 1.020;
  ratio respecto a Fuentes-Muñoz et al. (2025) = 1.014, z = +0.25.
- 6 perturbadores cuya deflexión queda bajo el ruido astrométrico por-encuentro
  arrojan ratio en [0.39, 0.72] (sesgo a la baja; mecanismo en §Sesgo a la baja).

El cierre previo de la capa de masas (Track A) atribuyó el fracaso al leverage de los
datos DR3. Estos resultados, obtenidos sobre FPR con ajuste simultáneo, acotan la causa
al método secuencial (LOO orbit→mass), no al leverage de la astrometría.

## Método

- **Ajuste simultáneo órbita+masa** (`orbdet.mass_determination.determine_shared_mass`):
  vector de 1 + 6N parámetros (masa compartida más 6 elementos por objetivo), Jacobiano
  en flecha, resuelto por Levenberg-Marquardt sobre el arco completo. La correlación
  masa↔drift orbital queda en la covarianza conjunta.
- **Fuerzas (ASSIST):** efeméride JPL DE440 (Sol, 8 planetas, Luna, Plutón) evaluada en
  cada paso; relatividad EIH; 16 perturbadores asteroidales masivos. Error de
  propagación frente a JPL Horizons: 0.17 mas sobre 900 d.
- **Observación:** estado → ICRS → (RA, Dec) con corrección de light-time iterativa;
  covarianza along-scan/across-scan anisotrópica de Gaia.
- **Covarianza en bloques por FOV.** Cada cruce de plano focal de Gaia produce 7 CCDs
  (separación ≈ 5 s) con residuos correlacionados (ICC = 0.32 medido sobre datos
  reales). Tratarlos como independientes subestima σ(masa) en un factor 1.66 (N_efectivo
  ≈ 0.36 N). Se blanquea con `C_bloque = diag(σ_AL²) + s_c²·11ᵀ` por cruce; el piso
  correlacionado s_c se calibra por bisección hasta χ²_red = 1.
- **Selección de objetivos:** encuentros < 0.05 AU del catálogo congelado
  (`--from-catalog`). El número de objetivos condiciona el sesgo (§Dependencia con N).
- **Rechazo de outliers:** sigma-clipping iterativo a 4σ.
- **Paralelización:** los N objetivos se evalúan en un pool por proceso (contexto
  `fork`, efeméride compartida por copy-on-write); factor 6 de aceleración respecto al
  modo serie, con resultado idéntico (test de equivalencia).

## Validación — calibradores Big-4 (Gaia FPR, N ≥ 20)

| Cuerpo | N obj | masa ajustada (kg) | σ_total | ratio fit/lit | z | fuente lit |
|--------|-------|--------------------|---------|---------------|---|------------|
| Ceres  | 28 | 8.96×10²⁰ | 4.6 % | 0.955 | −1.01 | DAWN (Park+ 2016) |
| Vesta  | 28 | 2.44×10²⁰ | 4.6 % | 0.943 | −1.30 | DAWN (Russell+ 2012) |
| Hygiea | 20 | 8.22×10¹⁹ | 5.7 % | 0.990 | −0.13 | Vernazza+ (2020) |
| Pallas |  6 | 2.54×10²⁰ | 7.0 % | 1.240 | +2.67 | Goffin (2014) |

Ceres, Vesta e Hygiea tienen N ≥ 20 y ratio en [0.943, 0.990]; el sesgo medio de los
tres es −4 %. Pallas tiene 6–7 encuentros < 0.05 AU en el catálogo completo (limitado
por objetivos), ratio 1.240, z = +2.67.

## Dependencia con N y origen del error

1. **Número de objetivos.** Con N ≈ 7, el ratio de los calibradores cae en [1.12, 1.29];
   con N ≥ 20 cae en [0.943, 0.990]. La diferencia es consistente con dispersión
   estadística de muestra pequeña. Un closing-loop sobre la geometría real (observaciones
   sintéticas generadas a la masa de referencia, con ruido del modelo, procesadas por el
   mismo pipeline) recupera ratio medio 0.997 sobre 3 semillas: el sesgo del estimador es
   compatible con cero.

2. **Límite por sistemáticos.** Con N grande, la σ formal (Fisher) escala como 1/√N y cae
   por debajo de 0.2 %, pero la exactitud queda acotada por sistemáticos por-encuentro
   (imperfección de la órbita del objetivo, astrometría local, perturbadores fuera de los
   16 modelados). Se reporta σ_total = √(σ_stat² + (f_sys·M)²) con f_sys = 4.2 %, igual a
   la RMS de (ratio − 1) de los calibradores con N ≥ 20.

## Barrido de los 12 perturbadores restantes

Mismo procedimiento; objetivos del catálogo (< 0.05 AU), N hasta 40. Ratio = masa
ajustada / masa de la efeméride DE441 (referencia publicada para estos cuerpos).
σ_stat = incertidumbre formal de Fisher.

| Cuerpo | N | masa ajustada (kg) | σ_stat | ratio fit/DE441 | χ²_red |
|--------|---|--------------------|--------|-----------------|--------|
| Psyche | 36 | 2.43×10¹⁹ | 3.3 % | 1.020 | 0.99 |
| Euphrosyne | 21 | 2.43×10¹⁹ | 27 % | 1.502 | 0.99 |
| Sylvia | 11 | 3.47×10¹⁹ | 11 % | 1.069 | 0.98 |
| Juno | 34 | 1.97×10¹⁹ | 11 % | 0.685 | 0.99 |
| Eunomia | 38 | 2.20×10¹⁹ | 9 % | 0.724 | 0.98 |
| Europa | 34 | 2.29×10¹⁹ | 6 % | 0.568 | 0.98 |
| Interamnia | 35 | 2.42×10¹⁹ | 13 % | 0.571 | 0.99 |
| Iris | 37 | 8.11×10¹⁸ | 19 % | 0.475 | 0.99 |
| Thisbe | 37 | 6.99×10¹⁸ | 6 % | 0.392 | 0.97 |
| Cybele | 12 | 6.24×10¹⁹ | 22 % | 4.44 | 0.91 |
| Camilla | 2 | — | 67 % | 8.1 | — |
| Davida | 3 | < 0 | — | < 0 | 0.65 |

**(16) Psyche:** M = 2.43 × 10¹⁹ kg, σ_stat = 3.3 %, N = 36, χ²_red = 0.99. Ratio
respecto a DE441 = 1.020. Es la única masa fuera de los calibradores con σ_stat < 5 % y
χ²_red en [0.97, 1.00]. Sylvia (ratio 1.069, σ 11 %) y Euphrosyne (ratio 1.502, σ 27 %)
son consistentes con DE441 dentro de su σ pero con menor precisión.

### Sesgo a la baja en perturbadores con deflexión débil

6 perturbadores con N ≥ 20 y χ²_red ≈ 1 arrojan ratio en [0.39, 0.72] (Juno, Eunomia,
Europa, Interamnia, Iris, Thisbe). El mecanismo es la degeneración masa↔órbita cuando la
señal de deflexión es comparable o inferior al ruido por-encuentro: el ajuste reproduce
la astrometría con menor masa y órbita reajustada (regresión hacia masa nula). Los
perturbadores cuya deflexión supera el ruido (Big-4, Psyche) no presentan este efecto.
Consecuencia operativa: para perturbadores débiles, (i) la masa de la efeméride no se
recupera con esta metodología sin regularización adicional, y (ii) la σ formal subestima
el error. La cuantificación de σ por-perturbador (jackknife/bootstrap) está en
[`planning/MASS_FUTURE_WORK.md`](../planning/MASS_FUTURE_WORK.md).

## Cruce con Fuentes-Muñoz (2025)

`scripts/validate/validate_fuentes_munoz_masses.py` compara las masas del catálogo con
la Tabla 5 de Fuentes-Muñoz et al. (2025, AJ 170, 353), convirtiendo su GMfin (km³/s²) a
masa con G = 6.67430 × 10⁻²⁰ km³ kg⁻¹ s⁻². Los 16 perturbadores del catálogo solapan con
su tabla. Para los calibradores, Fuentes-Muñoz fija GMfin a la semilla SB441/literatura;
esas filas no constituyen comparación independiente. El cruce con valor independiente es
sobre los no-calibradores, donde ese trabajo ejecutó su propio ajuste FPR.

| Cuerpo | M_orbdet (kg) | M_FM (kg) | ratio | z_vs_FM |
|--------|---------------|-----------|-------|---------|
| (16) Psyche | 2.429×10¹⁹ | 2.395×10¹⁹ | 1.014 | +0.25 |
| (10) Hygiea (cal) | 8.218×10¹⁹ | 8.237×10¹⁹ | 0.998 | −0.04 |
| (31) Euphrosyne | 2.430×10¹⁹ | 1.645×10¹⁹ | 1.477 | +1.15 |
| (52) Europa | 2.285×10¹⁹ | 2.656×10¹⁹ | 0.860 | −2.26 |
| (3) Juno | 1.971×10¹⁹ | 2.719×10¹⁹ | 0.725 | −3.34 |
| (7) Iris | 8.108×10¹⁸ | 1.456×10¹⁹ | 0.557 | −3.99 |
| (15) Eunomia | 2.195×10¹⁹ | 3.120×10¹⁹ | 0.704 | −4.24 |
| (88) Thisbe | 6.991×10¹⁸ | 8.755×10¹⁸ | 0.799 | −3.22 |
| (704) Interamnia | 2.420×10¹⁹ | 3.235×10¹⁹ | 0.748 | −2.50 |

(Tabla completa, incluidos los casos no fiables Sylvia/Cybele/Camilla/Davida, en
`data/output/literature_validation/fuentes_munoz_2025_mass_comparison.csv`.)

(16) Psyche concuerda con Fuentes-Muñoz con ratio 1.014, z = +0.25: la masa coincide con
un ajuste FPR independiente dentro de 1.4 %. Hygiea, Ceres y Vesta concuerdan con
|z| ≤ 1.3. Los 6 perturbadores con sesgo a la baja frente a DE441 muestran el mismo signo
frente a FM (ratio 0.56–0.80), con |z| > 3; esto es consistente con la subestimación de σ
descrita en §Sesgo a la baja: el z formal sobreestima la tensión porque σ_FM y σ_orbdet
no incorporan el error de regresión masa↔órbita.

## Reproducción

```bash
# Barrido de calibradores (validación)
docker compose run --rm pipeline python -m scripts.mass.orbdet_fit_realdata \
    --perturber big4 --release fpr \
    --from-catalog data/output/encounters_catalog_hybrid_stageb.parquet \
    --top-per-perturber 30 --workers 24 --out-dir data/output/orbdet/expanded

# Catálogo de masas con modelo de error completo
docker compose run --rm pipeline python -m scripts.mass.build_mass_catalog \
    --in-dir data/output/orbdet/expanded --out data/output/orbdet/mass_catalog.csv

# Cruce de masas vs Fuentes-Muñoz (2025)
docker compose run --rm pipeline python -m scripts.validate.validate_fuentes_munoz_masses
```
