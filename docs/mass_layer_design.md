# Mass layer design — joint orbit + mass fit

**Fecha**: 2026-05-26
**Stage**: Track 2 / Stage 1.2

## Objetivo

Reemplazar el ajuste secuencial actual:

1. ajustar 6 elementos orbitales del target,
2. congelarlos,
3. ajustar solo `M_perturber`,

por un ajuste conjunto de 7 parametros:

```text
θ = (log10_M_perturber, da_rel, de, di_deg, dOmega_deg, domega_deg, dM_deg)
```

El cambio necesario no es un nuevo propagador. El propagador N-body y el
forward astrometrico existentes son reutilizables; el cambio es la
parametrizacion del optimizador y la forma de los residuales.

## Contrato propuesto

Nuevo modulo:

```text
src/mass/forward_model_joint.py
```

Tipos/funciones iniciales:

```python
@dataclass(frozen=True)
class GaiaObservationBundle:
    jd_tdb: np.ndarray
    gaia_xyz_bary: np.ndarray
    ra_deg: np.ndarray
    dec_deg: np.ndarray
    position_angle_scan_deg: np.ndarray
    ra_error_systematic_mas: np.ndarray
    dec_error_systematic_mas: np.ndarray
    ra_dec_correlation_systematic: np.ndarray
    ra_error_random_mas: np.ndarray
    dec_error_random_mas: np.ndarray
    ra_dec_correlation_random: np.ndarray


@dataclass(frozen=True)
class JointFitPriors:
    sigma_da_rel: float
    sigma_de: float
    sigma_di_deg: float
    sigma_dOmega_deg: float
    sigma_domega_deg: float
    sigma_dM_deg: float
    log10_mass_bounds: tuple[float, float]


def apply_target_deltas(target_elements: dict, params: np.ndarray) -> dict:
    ...


def residuals_joint(
    params: np.ndarray,
    target_elements: dict,
    perturber_elements: dict,
    obs: GaiaObservationBundle,
    priors: JointFitPriors,
    background_elements: dict[str, dict] | None = None,
    dt_days: float = 1.0,
    integrator: str = "whfast",
) -> np.ndarray:
    ...
```

`residuals_joint` debe devolver un vector concatenado:

```text
[r_AL / sigma_AL, parameter_priors]
```

Esto permite usar `scipy.optimize.least_squares(method="trf")` igual que el
codigo actual, pero la masa y el drift orbital compiten dentro de una misma
likelihood.

## Parametrizacion

Se recomienda:

- `log10_M_perturber` en vez de masa lineal, porque las masas plausibles cubren
  muchos ordenes de magnitud.
- `da_rel = Δa / a`, no `Δa` absoluto, para estabilizar escalas entre objetos.
- Deltas angulares en grados, normalizados a `[-180, 180)` despues de aplicar.
- `e + de` clamped por bounds del optimizador a `[0, 0.999)`, no por clipping
  silencioso dentro del forward model.

Aplicacion:

```text
a_au       = a0 * (1 + da_rel)
e          = e0 + de
i_deg      = i0 + di_deg
Omega_deg  = wrap(Omega0 + dOmega_deg)
omega_deg  = wrap(omega0 + domega_deg)
M_deg      = wrap(M0 + dM_deg)
```

## Priors y bounds iniciales

Hasta disponer de covarianzas reales por objeto, usar priors conservadores y
explicitos:

| parametro | sigma inicial | bounds iniciales |
|---|---:|---:|
| `log10_M_perturber` | sin prior | `[14, 23]` |
| `da_rel` | `2e-4` | `±1e-3` |
| `de` | `5e-4` | `±2e-3` |
| `di_deg` | `0.05` | `±0.25` |
| `dOmega_deg` | `0.2` | `±1.0` |
| `domega_deg` | `0.2` | `±1.0` |
| `dM_deg` | `0.5` | `±2.5` |

Estos valores deben quedar en `JointFitPriors`, no hard-codeados dentro del
residual. Si el fit empuja repetidamente contra bounds, eso debe aparecer en
diagnosticos y bloquear claims de masa.

## Ventana de datos

Para la primera implementacion:

- Mantener la logica LOO para construir una orbita inicial razonable.
- Usar el fit conjunto sobre la ventana que hoy usa Fase B, pero incluir
  observaciones pre y post fuera del blackout.
- Excluir `abs(days_from_encounter) < blackout_days` para evitar que el
  minimo geometrico exacto domine con modelado temporal imperfecto.
- No restar un offset post-pre dentro del residual conjunto; si hace falta
  un offset instrumental, debe ser parametro explicito de nuisance en una
  etapa posterior.

## Diagnosticos obligatorios

El wrapper debe persistir al menos:

- `mass_kg`, `log10_mass`, incertidumbre formal e incertidumbre inflada.
- Los seis deltas orbitales ajustados.
- `chi2_red_joint`.
- `n_obs_fit`, `n_pre`, `n_post`.
- `active_bounds`: lista de parametros a menos de 1% del bound.
- `corr_mass_orbit_proxy`: condicion o rango singular de `JᵀJ` como proxy de
  degeneracion.
- `rms_al_mas` global, pre y post.

## Tests minimos

1. `apply_target_deltas` aplica `da_rel` y wrap angular correctamente.
2. Con `log10_M` fijo y deltas cero, `residuals_joint` reproduce el residual
   AL del forward model actual.
3. Dataset sintetico sin perturber masivo: el fit debe devolver masa cerca del
   bound bajo o no significativa, y recuperar deltas orbitales inyectados.
4. Dataset sintetico con masa conocida y drift orbital conocido: recupera masa
   y deltas dentro de tolerancia numerica.

Los tests 1-2 pueden ser unitarios puros. Los tests 3-4 pueden usar un sistema
reducido con `include_planets=("sun",)` para mantenerse rapidos y offline.

## Riesgos

- El forward model actual integra una simulacion nueva por evaluacion; un fit
  7D puede ser bastante mas caro que el fit 1D.
- Sin covarianzas orbitales reales, los priors siguen siendo heuristicas.
- Si `JᵀJ` queda mal condicionado, el resultado cientificamente correcto es
  reportar degeneracion masa-orbita, no forzar una masa.
- La falta de cache local para Gaia TAP limita reproducibilidad de corridas
  batch; debe resolverse antes de publicar resultados de masa.
