# Roadmap: De Candidatos a Masas Publicables

> Plan detallado del trabajo que falta para convertir los 24 perturbers
> novedosos en mediciones de masa publicables.

---

## Resumen ejecutivo

**Estado actual**: 41 candidatos viables identificados, 24 perturbers genuinamente
novedosos (no en literatura), pero las detecciones de perturbación NO son
encounter-specific con los métodos actuales.

**Lo que falta**: un fit conjunto orbit-mass usando N-body, que separe la
perturbación del encuentro del drift orbital. **Estimación total: 4-6 semanas
de trabajo enfocado**.

**Output final**: un catálogo de ~10-20 nuevas masas asteroidales con
incertidumbres, listo para publicar.

---

## Por qué los métodos actuales no alcanzan

El problema con la cadena actual:

```
Gaia observa → comparamos contra Horizons → vemos un shift
```

**Limitación**: Horizons usa su propio fit orbital del target, que se hizo SIN
incluir la masa del perturber. Ese fit puede tener residuales de tamaño similar
a la perturbación que queremos medir. Los shifts coinciden con el encuentro
pero también con muchas otras fechas (null test: 0/41 strongly specific).

**Solución**: hacer NUESTRO PROPIO fit orbital, donde la masa del perturber sea
un parámetro libre. Si la masa que mejor ajusta es no-cero, eso ES la detección.

```
Gaia observa → fit conjunto (orbit + perturber mass) → masa con incertidumbre
```

---

## FASE 1 — Core fitter (1-2 semanas)

### Objetivo
Implementar `scripts/fit_perturber_mass.py` que para un par (perturber, target,
fecha) devuelva la masa estimada del perturber con su incertidumbre.

### Tareas

#### 1.1 Coordinate transformations (1-2 días)
- [ ] `src/astrometry/transforms.py`: heliocentric_ecliptic → barycentric_ICRS
- [ ] Light-time correction iterativa:
  ```
  t_retarded = t_obs - |r_target(t_retarded) - r_gaia(t_obs)| / c
  Iterar 2-3 veces hasta converger.
  ```
- [ ] Stellar aberration (Gaia velocity, ~30 km/s)
- [ ] Test: con MPCORB de un asteroide conocido, reproducir RA/Dec de
      Horizons a ~1 mas

#### 1.2 N-body propagator con perturber configurable (2-3 días)
- [ ] Wrapper sobre `src/propagate/nbody.py` que acepte:
      - target asteroide (test particle)
      - planetas mayores (siempre)
      - big-4 asteroides (Ceres/Pallas/Vesta/Hygiea, masas conocidas)
      - **perturber** (con masa como input, valor a optimizar)
- [ ] Devolver heliocentric position at array of epochs
- [ ] Test: integrar Ceres como perturber, recuperar el efecto sobre un test particle

#### 1.3 Forward model: epochs → RA/Dec (1 día)
- [ ] Función que toma (orbital elements, perturber mass, epochs, gaia positions)
      y devuelve (predicted RA, predicted Dec)
- [ ] Internamente: propagate → light-time → aberration → transforms

#### 1.4 Residual function + least-squares fit (2-3 días)
- [ ] Free parameters: 6 elementos orbitales del target + masa del perturber
- [ ] Residuales: (obs_RA - pred_RA) × cos(dec), (obs_Dec - pred_Dec) en mas
- [ ] Weights: 1/σ²_per_transit (necesita σ_Gaia, ver Fase 2)
- [ ] `scipy.optimize.least_squares` con method='lm' o 'trf'
- [ ] Initial guess: MPCORB elements + mass from diameter (ρ=1.5)
- [ ] Covariance matrix → uncertainties via Jacobian

#### 1.5 CLI script (1 día)
- [ ] `scripts/fit_perturber_mass.py --perturber 511 --target 115180 --date 2014-11-19`
- [ ] Output: `data/output/mass_fit_<pert>_<target>.csv` con mass, sigma_mass,
      reduced_chi2, ra_residuals, dec_residuals
- [ ] Logging detallado para debugging

### Entregable Fase 1
Un fit que para (511) Davida + 2003_sm90 devuelva la masa de Davida con
incertidumbre. Debe corresponder con Goffin 2014 (~3.5e19 kg) dentro de
factor ~2.

### Riesgos
- **Convergencia**: el espacio de parámetros tiene 7 dimensiones. Puede haber
  múltiples mínimos locales. Mitigar con buena initial guess y posibles
  re-starts.
- **Performance**: cada evaluación del forward model implica una integración
  N-body (~1-10 segundos). Un fit típico necesita ~100-1000 evaluaciones.
  Optimizar: caché de trayectorias planetarias, paralelización.

---

## FASE 2 — Gaia observation uncertainties (3-5 días)

### Objetivo
Bajar las observaciones Gaia con su error formal per-transit, ya que afecta
directamente los pesos del fit.

### Tareas
- [ ] Actualizar `scripts/check_gaia_observations.py` (y derivados) para
      incluir las columnas de error de la tabla `gaiadr3.sso_observation`.
      Las columnas relevantes son `ra_error_systematic`, `dec_error_systematic`,
      `ra_error_random`, `dec_error_random` y la matriz de covarianza si está
      disponible.
- [ ] Test: para targets brillantes (mag<18), σ ~ 0.1-0.3 mas; para débiles
      (mag>20), σ ~ 1-3 mas
- [ ] Integrar errores como weights en el fit
- [ ] (Opcional) Sumar systematic en cuadratura

### Entregable Fase 2
Un fit que pondera correctamente las observaciones según su precisión.

---

## FASE 3 — Validación con calibration set (3-5 días)

### Objetivo
Aplicar Fase 1 a los 5 perturbers con masa conocida y verificar reproducción.

### Tareas
- [ ] Fit de cada uno de los 5:
      - (19) Fortuna (3 encounters en el catálogo)
      - (46) Hestia
      - (165) Loreley
      - (241) Germania
      - (511) Davida
- [ ] Comparar con literatura (Fienga 2003, Galad 2002, Goffin 2014, Fuentes-Muñoz 2024)
- [ ] Para cada uno: tabla con (mass_fit ± σ_fit) vs (mass_pub ± σ_pub)
- [ ] Calibration metric: bias = mean(fit/pub - 1), σ_calib = std

### Criterios de éxito
- |bias| < 0.20 (≤ 20% systematic)
- σ_calib < 0.30 (≤ 30% scatter)

### Si NO se cumplen
- **Bias grande**: revisar coordinate transforms (probable: aberración o
  light-time mal calculados)
- **Scatter grande**: revisar weights (probable: Gaia errors mal asignados)
  o presencia de perturbers no modelados

### Entregable Fase 3
Documento `data/output/calibration_results.csv` y conclusión "el método está
calibrado dentro de X%".

---

## FASE 4 — Aplicación a candidatos novedosos (5-7 días)

### Objetivo
Aplicar el fit a los 24 perturbers novedosos y producir el catálogo de masas.

### Tareas
- [ ] Loop sobre los 24 perturbers, cada uno con el encuentro de mayor δ esperado
- [ ] Para cada uno: fit + uncertainty + diagnostics (chi², residual plot)
- [ ] **MCMC para los top candidatos** (5-10 mejores) usando `emcee`:
      - mejor caracterización de la posterior
      - detecta degeneracies entre masa y elementos
- [ ] Outlier rejection: σ-clipping iterativo en los residuales
- [ ] Output: `data/output/mass_catalog_v1.csv` con
      `perturber, target, encounter_date, mass_kg, mass_sigma_kg,
       chi2_red, n_transits, fit_quality_flag`

### Sub-tareas paralelas
- [ ] Per-candidate diagnostic plots (residuales pre/post-fit, posterior MCMC)
- [ ] Identificar masas robustas vs marginal vs no-detection

### Entregable Fase 4
Catálogo de masas: probablemente **10-20 masas robustas** + algunas marginales
+ no-detections. Algunas de las 24 candidatos pueden quedar sin masa medible
si el fit no converge o la incertidumbre supera el valor central.

---

## FASE 5 — Cross-check + writeup (1-2 semanas)

### Objetivo
Confirmar novedad de cada masa medida y preparar para publicación.

### Tareas

#### 5.1 Cross-check literatura completo (3 días)
- [ ] Descargar y cargar:
      - Goffin (2014) full catalog
      - Fuentes-Muñoz et al. (2024) — masas de Gaia FPR  
      - Park et al. (2021) DE440 — masas de planetas/asteroides
- [ ] Para cada masa nueva, búsqueda directa en ADS / arXiv del perturber
- [ ] Distinguir entre "verdaderamente novel" y "preprint reciente"

#### 5.2 Catálogo final + paper (5-7 días)
- [ ] Tabla principal: perturber number, name, D_km, mass_new ± σ,
      mass_pub ± σ (si existe), encounter date, target, n_transits, χ²
- [ ] Figuras:
      - Histograma de masas
      - Mass vs diameter (para inferir densidad)
      - Calibration: fit_mass vs pub_mass para los 5 conocidos
      - Residual plots de 2-3 ejemplos
- [ ] Methods section: pipeline, fit, validation
- [ ] Discussion: comparison with similar Gaia work, implications

#### 5.3 Reproducibilidad (2 días)
- [ ] README ejecutable: cómo correr todo desde cero
- [ ] Docker setup verificado
- [ ] Tests unitarios para los componentes clave del fit
- [ ] Citation file (CITATION.cff)

### Entregable Fase 5
Manuscript draft listo para revisión + código publicable.

---

## Cronograma estimado

```
Semana 1-2:  Fase 1 (core fitter)
Semana 2:    Fase 2 (Gaia errors)         | paralelo con final de Fase 1
Semana 3:    Fase 3 (calibración)
Semana 4:    Fase 4 (24 candidatos)
Semana 5-6:  Fase 5 (validation + paper)
```

**Total: 5-6 semanas** asumiendo trabajo de tiempo completo. Para trabajo
de medio tiempo (~20 h/sem), duplica: ~10-12 semanas.

---

## Riesgos principales

### Técnicos

| Riesgo | Impacto | Mitigación |
|--------|---------|------------|
| El fit no converge para targets débiles | Alto | Mejor initial guess vía MPCORB + outlier rejection |
| Residuales sistemáticos > masa | Alto | Cross-check con Cat A (Ceres/Vesta) |
| Performance: fit muy lento | Medio | Caché de planetary motion, paralelización |
| Light-time / aberration con bugs | Alto | Test contra Horizons (intent #3 ya probó esto) |

### Científicos

| Riesgo | Impacto | Mitigación |
|--------|---------|------------|
| Otros perturbers contaminan signal | Medio | Incluirlos en el N-body, o fitting iterativo |
| Errors Gaia mal estimados | Medio | Boost factor empírico desde Cat A |
| Masas verdaderamente novedosas tienen baja SNR | Medio | Reportar como upper limits |
| Algún candidato ya tiene masa anunciada | Bajo | Cross-check periódico con ADS/arXiv |

### De proyecto

- **Scope creep**: agregar más candidatos / tests → mantenerse focusado
- **Burnout en Fase 1**: la complejidad de las transformaciones de coordenadas
  puede ser tediosa → dividir en sub-tareas con tests intermedios

---

## Dependencies externas

- Gaia DR3 archive (TAP) — disponible, sin auth
- JPL Horizons (astroquery) — disponible
- REBOUND (ya en pyproject.toml) — instalado en Docker
- emcee (para MCMC, agregarlo si Fase 4 lo requiere) — pip install simple

**No requiere acceso a datos privados ni colaboraciones externas**.

---

## Cómo arrancar (próxima sesión)

```bash
# 1. Mergear todos los PRs abiertos
gh pr merge 1 --merge && gh pr merge 2 --merge && ...

# 2. Crear branch para Fase 1
git checkout -b feat/fit-perturber-mass

# 3. Empezar con coordinate transforms (Fase 1.1)
#    Archivo: src/astrometry/transforms.py
#    Test: tests/test_astrometry_transforms.py

# 4. Iterativo: implementar cada sub-tarea con tests
```

---

## Versión mínima viable (si querés un demo más chico)

Si querés un proof-of-concept en 1 semana en vez del catálogo completo:

**Mini-roadmap (5 días)**:
1. Día 1-2: coordinate transforms + N-body propagator
2. Día 3: forward model + fit script básico (sin MCMC)
3. Día 4: aplicarlo a (511) Davida (calibration target con masa conocida)
4. Día 5: validar la masa recuperada vs Goffin 2014, escribir notas

Output: una sola masa fitted (Davida) + demostración del método funcionando.
Si funciona, escalar a las otras 28 ya es ingeniería pura.

Este mini-roadmap es la prueba de concepto antes de invertir en el catálogo
completo.

---

## Lo que NO requiere este roadmap

(para evitar confusión sobre el scope)

- ❌ Mejorar el catálogo de encuentros (ya está al 100%)
- ❌ Filtrar más candidatos (los 24 novel son suficientes)
- ❌ Recorrer los 119k novedosos (ya validamos el subset relevante)
- ❌ Mejor propagación Kepler (suficiente para detección, no para masas)
- ❌ Cross-check de literatura masivo (basta con los catálogos cargados)

El cuello de botella es **únicamente el fit conjunto**.

---

## Conclusión

Este pipeline está al ~70% del camino a publicación. Los componentes faltantes
están bien identificados y técnicamente accesibles. **No hay incertidumbre
científica fundamental — es trabajo de ingeniería astrométrica pura**.

Una iteración focused de 5-6 semanas convierte los 24 candidatos en un
catálogo de masas publicable. El mini-roadmap de 5 días resuelve la cuestión
"¿funciona?" antes de comprometerse al scope completo.
