# F3 — Extensión del fondo de perturbadores: ¿explica el sesgo medio −4 %?

> **Estado:** ✅ COMPLETO (Ceres/Vesta/Hygiea corridos, fondo 16 vs 35). 2026-07-03.
> **Resultado:** hipótesis refutada — extender el fondo no reduce el sesgo −4 %
> (f_sys 4.158 % → 4.257 %). Gate no cumplido por esta vía.
> Ítem F3 de [`planning/MASS_FUTURE_WORK.md`](../planning/MASS_FUTURE_WORK.md).
> **Gate:** f_sys (RMS de ratio−1 sobre calibradores con N≥20) por debajo de 4.2 %.

## Motivación

El ajuste conjunto órbita+masa (motor `orbdet`, backend ASSIST) recupera los
calibradores Ceres/Vesta/Hygiea con un ratio medio fit/referencia ≈ 0.96 con
N ≥ 20 objetivos — un déficit sistemático de **−4 %** frente a las masas de
DAWN/Vernazza. El χ²_red ≈ 1 en todos, así que no es un problema de escala de
error: la astrometría se explica con **menos masa** de la esperada.

**Hipótesis (completitud del fondo).** El modelo dinámico incluye solo los 16
perturbadores asteroidales masivos de la efeméride `sb441-n16.bsp`. Los asteroides
masivos **no modelados** (17º en adelante) también desvían las órbitas de los
objetivos; al no estar en el modelo, su tirón se absorbe re-ajustando la órbita del
objetivo y **restando** señal a la masa del perturbador bajo estudio → sesgo a la
baja. Si la hipótesis es correcta, extender el fondo debería subir los ratios hacia
1.0 y bajar f_sys.

## Método

Se extiende el fondo de 16 a 16+N cuerpos agregando los **N asteroides más masivos
fuera de los 16**, con:

- **masa** de Fuentes-Muñoz et al. (2025), Tabla 5, columna `GMfin`
  (M = GM/G, G = 6.67430×10⁻²⁰ km³ kg⁻¹ s⁻²), y
- **órbita** (elementos osculadores heliocéntricos eclípticos en la época común
  del ajuste) desde **JPL Horizons** — mismo tratamiento que los 16 (estado inicial
  desde efeméride/Horizons, integrado luego bajo Sol+planetas+GR vía ASSIST).

Cada cuerpo extra entra como una partícula masiva más de rebound. La fuerza
`ASTEROIDS` de la efeméride está excluida en el motor, así que **no hay doble
conteo**: ningún asteroide (ni los 16 ni los extra) aparece en las fuerzas de la
efeméride. Implementado en la capa IO
([`scripts/mass/orbdet_fit_realdata.py`](../scripts/mass/orbdet_fit_realdata.py),
helpers `_fm_extra_perturbers` / `_extended_background`, flag `--extra-background N`),
reusando la maquinaria de perturbador custom de F4. El motor `src/orbdet/` **no se
tocó**.

### Cuerpos agregados al fondo (N=20, los más masivos de FM 2025 fuera de los 16)

| # | nombre | M_FM (kg) |
|---|--------|-----------|
| 29 | Amphitrite | 1.358×10¹⁹ |
| 6 | Hebe | 1.237×10¹⁹ |
| 451 | Patientia | 1.197×10¹⁹ |
| 532 | Herculina | 9.959×10¹⁸ |
| 324 | Bamberga | 9.808×10¹⁸ |
| 19 | Fortuna | 8.524×10¹⁸ |
| 22 | Kalliope | 7.764×10¹⁸ |
| 130 | Elektra | 7.437×10¹⁸ |
| 13 | Egeria | 7.230×10¹⁸ |
| 48 | Doris | 7.150×10¹⁸ |
| 24 | Themis | 6.842×10¹⁸ |
| 702 | Alauda | 6.827×10¹⁸ |
| 423 | Diotima | 6.783×10¹⁸ |
| 354 | Eleonora | 6.604×10¹⁸ |
| 9 | Metis | 6.482×10¹⁸ |
| 39 | Laetitia | 6.139×10¹⁸ |
| 45 | Eugenia | 6.032×10¹⁸ |
| 92 | Undina | 5.939×10¹⁸ |
| 372 | Palma | 5.512×10¹⁸ |
| 375 | Ursula | 5.497×10¹⁸ |

Para referencia: el menor de los 16 de fondo (Iris) tiene M ≈ 1.3×10¹⁹ kg y el
menor calibrador (Hygiea) ≈ 8.3×10¹⁹ kg. Los 20 extra suman ≈ 1.6×10²⁰ kg,
repartidos por todo el cinturón principal.

## Resultados

Los 3 calibradores con N ≥ 20 objetivos, ajustados con el fondo de 16 y con el
fondo extendido de 35 (idénticos objetivos, corte de rechazo y semilla; sin
jackknife, que no interviene en f_sys):

| # | cuerpo | N | masa fondo-16 (kg) | masa fondo-35 (kg) | Δmasa | ratio-16 | ratio-35 | χ²_red-35 |
|---|--------|---|--------------------|--------------------|-------|----------|----------|-----------|
| 1 | Ceres  | 35 | 9.0144×10²⁰ | 9.0122×10²⁰ | −0.02 % | 0.9606 | 0.9604 | 0.989 |
| 4 | Vesta  | 35 | 2.4360×10²⁰ | 2.4303×10²⁰ | −0.23 % | 0.9405 | 0.9384 | 0.964 |
| 10 | Hygiea | 20 | 8.2178×10¹⁹ | 8.2326×10¹⁹ | +0.18 % | 0.9901 | 0.9919 | 0.980 |

**f_sys** = RMS(ratio − 1) sobre los 3 calibradores:

| fondo | f_sys | gate < 4.2 % |
|-------|-------|--------------|
| 16 (base) | **4.158 %** | ✅ (marginal) |
| 35 (extendido) | **4.257 %** | ❌ |

## Conclusión

**La hipótesis de completitud del fondo queda refutada.** Agregar los 20
asteroides más masivos de FM 2025 al modelo dinámico (fondo 16 → 35, +1.6×10²⁰ kg
repartidos por el cinturón) desplaza las masas de los calibradores en **< 0.25 %** —
tres órdenes de magnitud por debajo del déficit del −4 % que se buscaba explicar.
El ratio fit/referencia se queda clavado en ~0.96 (Ceres, Vesta) y ~0.99 (Hygiea),
y f_sys **no baja: sube marginalmente** (4.158 % → 4.257 %), cruzando el gate hacia
el lado equivocado.

Lectura física: los perturbadores no modelados tiran de los objetivos, pero su
efecto sobre la masa reconstruida del calibrador es **despreciable** frente al
sesgo observado. La deflexión de un fondo lejano y disperso se promedia a casi cero
sobre el conjunto de objetivos; lo que sesga la masa a la baja **no** es el fondo
faltante. Candidatos remanentes para el −4 % (fuera del alcance de F3): imperfección
de la órbita de cada objetivo de prueba (su propio arco Gaia es corto), sistemáticos
astrométricos locales por-encuentro, o un piso de correlación intra-tránsito no
capturado del todo por el modelo de bloques por FOV. El gate de F3 (f_sys < 4.2 %
por esta vía) **no se cumple**, y la evidencia indica que esta vía no es la
correcta. Nota: el fondo-16 base ya está en 4.158 %, apenas por debajo del 4.2 %,
de modo que el −4 % de sesgo medio no es reducible extendiendo el catálogo de
perturbadores.

## Comandos exactos

```bash
# Ceres (perturbador 1) con fondo extendido de 20 cuerpos, N≥20 objetivos:
docker compose run --rm pipeline python -m scripts.mass.orbdet_fit_realdata \
    --perturber 1 --release fpr \
    --from-catalog data/output/encounters_catalog_hybrid_stageb.parquet \
    --top-per-perturber 40 --extra-background 20 --jackknife --workers 6 \
    --out data/output/orbdet/f3_extbg/ceres_extbg20_fpr.json

# Vesta (4):
docker compose run --rm pipeline python -m scripts.mass.orbdet_fit_realdata \
    --perturber 4 --release fpr \
    --from-catalog data/output/encounters_catalog_hybrid_stageb.parquet \
    --top-per-perturber 40 --extra-background 20 --jackknife --workers 6 \
    --out data/output/orbdet/f3_extbg/vesta_extbg20_fpr.json

# Hygiea (10):
docker compose run --rm pipeline python -m scripts.mass.orbdet_fit_realdata \
    --perturber 10 --release fpr \
    --from-catalog data/output/encounters_catalog_hybrid_stageb.parquet \
    --top-per-perturber 40 --extra-background 20 --jackknife --workers 6 \
    --out data/output/orbdet/f3_extbg/hygiea_extbg20_fpr.json

# f_sys (fondo 16 vs extendido):
docker compose run --rm pipeline python -m scripts.mass.f3_fsys \
    --jack-catalog data/output/orbdet/mass_catalog_jack.csv \
    --ext-json data/output/orbdet/f3_extbg/ceres_extbg20_fpr.json \
    --ext-json data/output/orbdet/f3_extbg/vesta_extbg20_fpr.json \
    --ext-json data/output/orbdet/f3_extbg/hygiea_extbg20_fpr.json
```
