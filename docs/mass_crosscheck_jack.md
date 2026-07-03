# Cruce de masas vs Fuentes-Muñoz (2025) con σ externa por jackknife (F5)

**Estado**: cerrado — cruce ejecutado sobre el catálogo con jackknife.
**Fecha**: 2026-07-02
**Frente**: P3 / ítem F5
**Script**: `scripts/mass/crosscheck_fuentes_munoz_jack.py`
**Insumo**: `data/output/orbdet/mass_catalog_jack.csv` (16 perturbadores)
**Referencia**: Fuentes-Muñoz et al. (2025), AJ 170, 353, Tabla 5 (`GMfin`),
parseada de `data/raw/fuentes_munoz_2025/ajae0cc9t5_mrt.txt`
(GM → masa con G = 6.67430e-20 km³ kg⁻¹ s⁻²).

---

## Qué se hizo y por qué

El cruce previo (`validate_fuentes_munoz_masses.py`) mide la consistencia con
Fuentes-Muñoz usando `sigma_total_kg` del catálogo **no-jack**, cuya parte
estadística es la σ **formal** de Fisher (diagonal de la covarianza). Esa σ baja
como 1/√N y **subestima** el error real: ignora la regresión masa↔órbita que sólo
se ve dejando encuentros afuera. La σ por **jackknife** (F1) reincorpora esa
regresión, así que el z construido sobre ella es el defendible.

Este ítem (F5) recomputa el cruce con la σ jackknife y pone los dos z lado a lado:

- `z_formal = (M − M_FM) / √(σ_formal² + σ_FM²)`
- `z_jack   = (M − M_FM) / √(σ_total²  + σ_FM²)`, donde `σ_total` ya combina
  σ_jack con el piso sistemático (`sigma_total_kg` del catálogo jack).

### Referencia usada: disponible vs faltante

- **Fuentes-Muñoz Tabla 5 SÍ está en el repo** y cubre los **16** perturbadores
  del catálogo (parseada del MRT descargado). Es la referencia usada aquí, y es
  independiente para los no-calibradores.
- La columna `ref_mass_kg`/`ref_source` del propio CSV **no** es Fuentes-Muñoz
  para los no-calibradores: trae la semilla `DE441 ephemeris (seed)` (y
  DAWN/Goffin/Vernazza para los 4 calibradores). No se usó como referencia del
  cruce; se prefirió la Tabla 5 de FM, que es un ajuste FPR independiente.
- **Advertencia**: para los Big-4 (Ceres, Vesta, Pallas, Hygiea) FM ancla `GMfin`
  a la semilla literatura/SB441, así que esas filas **no** son comparación
  independiente — recuperan los mismos valores que ya calibramos. El cruce real
  son los no-calibradores.

---

## Tabla de resultados

`st`: estado de identificabilidad (F2). `M` = measured (snr_jack ≥ 3, medida
genuina); `bnd` = not_identifiable (la deflexión queda bajo el ruido
por-encuentro → **cota, no medida**); `?` = unknown (ajuste sin jackknife).
`cal` marca calibrador. Masas en kg.

| # | nombre | cal | st | M_ours | σ_formal | σ_jack | M_FM | ratio | z_formal | z_jack | status |
|---:|:---|:--:|:--:|---:|---:|---:|---:|---:|---:|---:|:---|
| 1 | Ceres | Y | M | 9.014e+20 | 1.632e+19 | 3.783e+19 | 9.385e+20 | 0.960 | −2.27 | −0.70 | measured |
| 2 | Pallas | Y | M | 2.542e+20 | 1.416e+19 | 1.052e+19 | 2.049e+20 | 1.240 | +3.46 | +2.78 | measured |
| 3 | Juno | . | M | 1.971e+19 | 2.073e+18 | 5.759e+18 | 2.719e+19 | 0.725 | −3.60 | −1.29 | measured |
| 4 | Vesta | Y | M | 2.436e+20 | 4.538e+18 | 8.477e+18 | 2.590e+20 | 0.941 | −3.38 | −1.16 | measured |
| 7 | Iris | . | bnd | 8.108e+18 | 1.572e+18 | 8.294e+18 | 1.456e+19 | 0.557 | −4.08 | −0.78 | not_identifiable |
| 10 | Hygiea | Y | M | 8.218e+19 | 3.097e+18 | 6.864e+18 | 8.237e+19 | 0.998 | −0.06 | −0.02 | measured |
| 15 | Eunomia | . | M | 2.195e+19 | 1.962e+18 | 5.705e+18 | 3.120e+19 | 0.704 | −4.68 | −1.60 | measured |
| 16 | Psyche | . | M | 2.429e+19 | 8.130e+17 | 2.056e+18 | 2.395e+19 | 1.014 | +0.40 | +0.15 | measured |
| 31 | Euphrosyne | . | bnd | 2.430e+19 | 6.636e+18 | 1.674e+19 | 1.645e+19 | 1.477 | +1.17 | +0.47 | not_identifiable |
| 52 | Europa | . | M | 2.285e+19 | 1.311e+18 | 2.754e+18 | 2.656e+19 | 0.860 | −2.78 | −1.27 | measured |
| 65 | Cybele | . | bnd | 4.111e+19 | 1.388e+19 | 1.407e+19 | 1.452e+19 | 2.831 | +1.91 | +1.88 | not_identifiable |
| 87 | Sylvia | . | bnd | 3.472e+19 | 3.894e+18 | 1.332e+19 | 1.424e+19 | 2.438 | +5.20 | +1.53 | not_identifiable |
| 88 | Thisbe | . | M | 6.991e+18 | 4.444e+17 | 1.358e+18 | 8.755e+18 | 0.799 | −3.81 | −1.27 | measured |
| 107 | Camilla | . | ? | 1.748e+20 | 1.166e+20 | — | 1.156e+19 | 15.122 | +1.40 | +1.40 | unknown |
| 511 | Davida | . | bnd | −1.058e+20 | 1.574e+19 | 4.696e+19 | 2.927e+19 | −3.616 | −8.58 | −2.86 | not_identifiable |
| 704 | Interamnia | . | M | 2.420e+19 | 3.017e+18 | 7.888e+18 | 3.235e+19 | 0.748 | −2.64 | −1.02 | measured |

Salidas de datos:
`data/output/literature_validation/fuentes_munoz_jack_mass_comparison.csv` y
`.../fuentes_munoz_jack_mass_summary.json`.

---

## Conteo final: |z| < 3

Sólo cuentan como test de consistencia las masas `measured` (10 de 16). Las
`not_identifiable` son cotas y se **excluyen** del conteo; `unknown` (Camilla, sin
jackknife) también.

| población | con σ formal | con σ jackknife |
|:---|:--:|:--:|
| **Todas las measured (overlap)** | **5 / 10** | **10 / 10** |
| **No-calibradores measured** (independiente) | **3 / 6** | **6 / 6** |

Nota de continuidad con el reporte previo: el `mass_summary.json` no-jack reportaba
"4/8 dentro de |z|<3" sobre el bloque `independent_reliable` (8 objetos, umbral de
fiabilidad por nº de objetivos). Aquí el bloque independiente son las 6
no-calibradoras `measured` (el criterio de inclusión pasó de "reliable" a
"measured" según F2). En ambos marcos el resultado es el mismo: **con la σ formal
alrededor de la mitad de las masas quedan en tensión aparente (|z|>3), y con la σ
jackknife todas las medidas se vuelven consistentes con Fuentes-Muñoz**.

---

## Lectura honesta

- La σ jackknife es sistemáticamente **mayor** que la formal (típicamente ×2–3;
  p. ej. Ceres 3.78e19 vs 1.63e19, Vesta 8.48e18 vs 4.54e18). Absorber la
  regresión masa↔órbita ensancha la barra de error y elimina las tensiones
  espurias de la σ formal. Esto no es "hacer pasar" resultados: es reconocer que
  la σ formal subestima el error real, cosa que ya se veía en que los propios
  calibradores Big-4 (masa verdadera conocida) quedaban a |z|>3 con σ formal
  (Pallas +3.46, Vesta −3.38).
- Las cotas siguen siendo cotas. Iris, Euphrosyne, Cybele, Sylvia y Davida son
  `not_identifiable`: su deflexión no supera el ruido por-encuentro (snr_jack < 3).
  Que su `z_jack` caiga dentro de ±3 **no** las convierte en medidas — la σ jack es
  tan grande que casi cualquier valor sería "consistente". Davida además tiene masa
  ajustada negativa (no física) con sólo 3 objetivos: es ruido, no señal.
- Camilla (107, `unknown`) tiene 2 objetivos y no se le corrió jackknife; su ratio
  ×15 vs FM es puramente un artefacto de muestra mínima, no un resultado.
- Conclusión defendible: de los 12 no-calibradores, **6 son masas medidas y las 6
  son consistentes con Fuentes-Muñoz bajo la σ externa por jackknife** (Juno,
  Vesta† , Eunomia, Psyche, Europa, Thisbe, Interamnia — descontando calibradores);
  las restantes son cotas o muestras insuficientes. La σ jackknife es la barra de
  error que hay que reportar.

† Vesta es calibrador; las 6 measured no-calibradoras son Juno, Eunomia, Psyche,
Europa, Thisbe e Interamnia.
