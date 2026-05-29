# Stage 2 — Covarianza Gaia AL/AC: diseño

> Plan técnico para reemplazar el likelihood AL-1D actual por una
> formulación Mahalanobis 2D que use ambas componentes (along-scan y
> across-scan) con la matriz de covarianza Gaia 2×2 por observación.

**Branch**: `track2/stage2-gaia-covariance`
**Inputs**: `docs/mass_layer_joint_diagnostic.md` (los 3 outliers χ²_red ≥ 10).
**Outputs esperados**: `src/mass/likelihood_al.py` + tests + diagnóstico
re-corrida del batch.

---

## Diagnóstico del estado actual

`src/mass/forward_model_joint.py::al_residuals_and_weights` proyecta el
residuo tangencial `(dRA*, dDec)` (mas) sobre el versor along-scan
`e_AL = (sin PA, cos PA)`:

```
r_AL = dRA* sin(PA) + dDec cos(PA)
σ²_AL = e_AL^T Σ_RA,Dec e_AL   (sistematic + random sumados en cuadratura)
```

y descarta la proyección perpendicular `r_AC`. Justificación implícita:
`σ_AC ≫ σ_AL` (~10×) en Gaia, así que AC aporta poca información de
parámetro. Pero **el χ² test sí se evalúa contra el modelo en las dos
componentes**: si el modelo no predice bien la parte AC, descartarla
oculta misfit, **no** lo elimina del fit (porque el optimizer minimiza
solo AL, no ve el residuo AC), pero distorsiona el χ²_red reportado.

Hipótesis: para los 3 outliers Stage 1 (n_obs ≫ 50, χ²_red 26–192) el
residuo restante AL contiene sistemática 2D que el modelo joint no
puede absorber con sus 7 parámetros, y AL-only **infla** el χ²_red al
no permitir que la sistemática cancele parcialmente entre AL y AC.

## Diseño del nuevo likelihood

### Vector residuo por observación

Para cada observación `i`, formar el vector 2D en el plano tangente
del cielo:

```
δ_i = (dRA*_i, dDec_i)^T  [mas]
```

donde `dRA* = (RA_obs − RA_pred) cos(Dec_pred)` (componente tangencial,
ya manejada por `residuals_mas`).

### Matriz de covarianza por observación

La covarianza total Σ_i suma sistematic + random:

```
Σ_i = Σ_sys_i + Σ_rand_i
```

donde, para cada componente,

```
Σ_X = | σ_RA²        ρ σ_RA σ_Dec |
      | ρ σ_RA σ_Dec  σ_Dec²       |
```

(elementos: `ra_error_*_mas`, `dec_error_*_mas`, `ra_dec_correlation_*`,
ya disponibles vía la query TAP existente en
`scripts/mass/fit_mass_gaia_loo.py:115-117`).

### Residual Mahalanobis 2D

```
χ²_i = δ_i^T Σ_i⁻¹ δ_i
```

Para inyectarlo en `scipy.optimize.least_squares`, factorizar
`Σ_i⁻¹ = L_i^T L_i` (Cholesky 2×2) y devolver el vector "blanqueado":

```
r_i = L_i δ_i  ∈ R²
```

de modo que `||r_i||² = χ²_i`. El vector completo es la concatenación de
los `r_i` por observación más el vector de priors (sin cambio).

`least_squares` minimiza `(1/2) Σ ||r_i||²`, que es el negativo del
log-likelihood Mahalanobis salvo constantes. Equivalente a fit por
máxima verosimilitud bajo errores Gaussianos 2D.

### Inversión robusta de Σ 2×2

Para una matriz `[[a, b], [b, c]]` con det `d = a c − b²`:

```
Σ⁻¹ = (1/d) [[c, −b], [−b, a]]
```

Numéricamente: si `d < ε` por alguna sistemática mal reportada (ρ → ±1),
caer a un fallback diagonal `diag(1/a, 1/c)` y emitir un warning.
Tests deben cubrir este caso.

### Cholesky 2×2 explícito

Para `Σ⁻¹ = L^T L` con L triangular superior:

```
L11 = sqrt(1/d * c)
L12 = (1/d * −b) / L11
L22 = sqrt(1/d * a − L12²)
```

Implementación vectorizada con numpy: cada observación produce un par
de números `(r_AL', r_AC')` y se apilan en un vector de longitud `2 N`.

## Cambios al forward model

### Nueva función `mahalanobis_residuals_2d`

```python
def mahalanobis_residuals_2d(
    dra_mas: np.ndarray,        # shape (N,) = (RA_obs - RA_pred) * cos(Dec_pred)
    ddec_mas: np.ndarray,       # shape (N,)
    ra_err_sys: np.ndarray,     # shape (N,)
    dec_err_sys: np.ndarray,
    corr_sys: np.ndarray,       # shape (N,), valores en [-1, 1]
    ra_err_rand: np.ndarray,
    dec_err_rand: np.ndarray,
    corr_rand: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return whitened residuals (shape (2N,)) and per-obs chi² (shape (N,))."""
```

Devuelve `(r_whitened, chi2_per_obs)`. Equivale a `r_whitened = L · δ`
con L Cholesky de Σ⁻¹.

### Wrapper de compatibilidad

Conservar `al_residuals_and_weights` sin cambio (otros call sites lo
usan, incluido `fit_mass_gaia_loo.py` simple-fit). Agregar parámetro
opcional `likelihood: Literal["al", "mahalanobis2d"]` en
`residuals_joint` y `fit_joint`; default `"al"` para no romper
reproducibilidad. Tests Stage 1 deben seguir verdes.

## Tests sintéticos (criterio de aceptación)

`tests/test_likelihood_al.py`:

1. **Diagonal sin correlación**: σ_RA = σ_Dec = 5 mas, ρ = 0, residuos
   gaussianos N(0, σ²) → `mean(chi²_per_obs) ≈ 2` (2 grados de libertad
   por observación), `std` consistente.

2. **Recuperación de chi² Mahalanobis conocido**: residuo `(3, 4)` con
   σ_RA = 1, σ_Dec = 2, ρ = 0 → χ² = 9 + 4 = 13.

3. **Correlación ρ = 0.5**: residuo `(1, 1)` con σ_RA = σ_Dec = 1 →
   χ² calculado por la fórmula cerrada `(δᵀ Σ⁻¹ δ)` debe coincidir
   con el output de `mahalanobis_residuals_2d` dentro de 1e-12.

4. **Suma sistematic+random**: dos contribuciones diagonales se suman
   correctamente en cuadratura (σ_total² = σ_sys² + σ_rand²).

5. **Fallback degenerate**: ρ = 0.9999 (det ≈ 0) → no NaN, warning
   emitido, vector `r_whitened` finito.

6. **Equivalencia con AL en límite σ_AC → ∞**: si σ_RA = σ_Dec con
   `corr` ajustado para que `σ_AC = 1000 σ_AL`, χ² Mahalanobis 2D
   debe coincidir con `(r_AL/σ_AL)²` dentro de 1e-3.

   *Nota*: este límite no es trivial de construir porque σ_AL y σ_AC
   se derivan de (σ_RA, σ_Dec, ρ, PA). Construir el caso con
   PA = 45° y σ_RA = σ_Dec ajustando ρ para que el eje propio
   coincida con el eje AC.

## Plan operativo

| Sub-paso | Output | ~Tiempo |
|---|---|---|
| 2.1 Implementar `src/mass/likelihood_al.py` | módulo Mahalanobis 2D | 2 h |
| 2.2 Tests `tests/test_likelihood_al.py` | 6 tests pasando | 2 h |
| 2.3 Integrar en `forward_model_joint.py` con flag | residuals_joint con `likelihood=` | 1 h |
| 2.4 Tests existentes `tests/test_forward_model_joint.py` | 4/4 siguen verdes | 30 min |
| 2.5 Re-correr 3 outliers Stage 1 con Mahalanobis 2D | tabla comparativa | 30 min cómputo |
| 2.6 Re-correr batch 41 candidatos completo | `loo_batch_results_joint_mahal.csv` | ~5 min cómputo |
| 2.7 Diagnóstico `docs/mass_layer_stage2_diagnostic.md` | reporte χ²_red, Δmasas, conclusión | 2 h |
| 2.8 PR + merge | merge | 30 min |

**Criterios de aceptación finales**:

- Tests sintéticos (1)-(6) verdes.
- Re-corrida 41-candidatos: χ²_red mediano joint+Mahalanobis ≤ AL-1D
  mediano (1.52). Si **mejor**, conclusión clara. Si **igual**,
  conclusión "AL-1D era suficiente" — también informativo.
- Para los 3 outliers: documentar si bajan a < 10 (ideal) o si siguen
  altos (entonces el problema no es covarianza sino sistemática
  no modelada — gate hacia Stage 3 honesto).
- PR mergeada a main.

## Riesgos y mitigaciones

- **Riesgo**: el factor Cholesky por observación añade latency. Mitigación:
  vectorizar todo con numpy (det / inv 2×2 cerrado, sin loops). Profile
  contra el AL-1D actual; aceptar hasta 2× slowdown porque el batch
  total es 5 min.
- **Riesgo**: ρ del catálogo Gaia no está acotado en [−1, 1] por bugs
  de export. Mitigación: clamp con warning antes de usarlo en det.
- **Riesgo**: el `least_squares` con vector más largo (2N + 6 vs N + 6)
  reasigna grados de libertad → χ²_red se computa con `n_obs_effective`
  distinto. Mitigación: definir `n_obs_effective = 2N` para Mahalanobis
  2D y documentarlo en el output.
