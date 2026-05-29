# Mass layer audit — Track 2 Stage 1.1

**Fecha**: 2026-05-26
**Estado**: inventario inicial del modelo actual

## Scope

Este documento audita la capa actual de ajuste de masas antes de introducir el
fit conjunto `M + Δorbita`. El foco es el flujo activo en
[`scripts/mass/fit_mass_gaia_loo.py`](../scripts/mass/fit_mass_gaia_loo.py),
que es el pipeline más avanzado frente a los prototipos previos
`fit_perturber_mass.py` y `fit_mass_linear.py`.

## Flujo actual

1. Descarga todas las observaciones Gaia DR3 SSO del target desde TAP
   (`fetch_gaia_full`).
2. Convierte epochs Gaia TCB a JD TDB y usa `x_gaia/y_gaia/z_gaia` como
   posición baricéntrica del observador.
3. Selecciona el snapshot MPCORB más cercano al encuentro si no se pasa
   `--mpcorb`.
4. Fase A (`fit_orbit_loo`): ajusta 6 elementos orbitales del target usando
   solo observaciones pre-encuentro fuera de la ventana LOO.
5. Fase B (`fit_mass_from_window`): congela esos 6 elementos ajustados y
   optimiza solo `log10(M_perturber)` sobre observaciones post-encuentro.
6. Escribe un JSON por par `data/output/fit_<perturber>_<target>_loo.json`.
   `scripts/mass/summarize_loo_fits.py` consolida esos JSON en
   `data/output/loo_batch_results.csv`.

## Modelo fisico actual

La prediccion astrometrica esta encapsulada en
[`src/astrometry/forward_model.py`](../src/astrometry/forward_model.py):

- Entrada: elementos MPCORB del target, elementos del perturber, masa del
  perturber, epochs Gaia y posicion baricentrica de Gaia.
- Propagacion: `propagate_target_with_perturber` en
  [`src/propagate/nbody_perturber.py`](../src/propagate/nbody_perturber.py).
- Cuerpos: planetas por defecto, perturber masivo con masa ajustable, target
  como particula sin masa, y asteroides de fondo opcionales con masas fijas.
- Salida: RA/Dec predichos en el frame astrometrico baricentrico usado por
  Gaia DR3 SSO, con correccion iterativa de light-time.

## Parametros ajustados

Fase A ajusta 6 parametros absolutos del target:

- `a_au`, `e`, `i_deg`, `Omega_deg`, `omega_deg`, `M_deg`.
- Bounds hard-coded alrededor de MPCORB: `a ± 0.05 AU`, `e ± 0.05`,
  `i ± 2 deg`, `Omega/omega ± 10 deg`, `M ± 15 deg`.
- Regularizacion gaussiana ad hoc via `reg_sigma =
  [5e-4, 5e-4, 0.1, 0.3, 0.3, 0.5]`.
- La masa del perturber se fija en cero durante este ajuste orbital.

Fase B ajusta 1 parametro:

- `log10(M_perturber)`, con bounds `[14, 23]` en kg.
- Los 6 elementos del target quedan congelados en el valor de Fase A.
- La incertidumbre de masa sale de `(JᵀJ)^-1 * chi2_red`; no propaga la
  incertidumbre/covarianza del ajuste orbital de Fase A.

## Likelihood y pesos

El codigo ya usa una mejora importante frente a un chi2 2D ingenuo:

- `al_residuals_and_weights` proyecta `(ΔRA, ΔDec)` sobre la direccion
  along-scan usando `position_angle_scan`.
- `sigma_AL` se calcula proyectando la covarianza Gaia RA/Dec sistematica y
  random sobre esa misma direccion.
- El residuo optimizado es escalar: `r_AL / sigma_AL`.

Limitacion: esto no es todavia un likelihood Mahalanobis 2D completo. La
componente across-scan se descarta, lo cual es pragmatico para SSO, pero Stage
2 debe decidir si conviene modelar la matriz 2x2 completa por observacion.

## Supuestos criticos

- El target y el perturber comparten epoch MPCORB dentro de `1e-3` dias.
- El target es siempre una particula test; su masa no retroalimenta la
  dinamica.
- Las masas de asteroides de fondo son fijas y vienen de `_MAJOR_ASTEROIDS`.
- La Fase A usa solo pre-encuentro para no absorber la senial post-encuentro.
- La Fase B resta el promedio pre-encuentro de `r_AL` como baseline residual
  antes de evaluar el post-encuentro.
- La seleccion de snapshot MPCORB reduce drift por epoch, pero no usa
  covarianzas orbitales reales de MPCORB/Gaia.

## Debilidades que explican el fallo actual

1. El ajuste no es conjunto. La masa se estima despues de congelar la orbita,
   por lo que drift orbital residual y masa siguen siendo confundibles.
2. La covarianza entre los 6 elementos orbitales y `M_perturber` no existe en
   el resultado; la incertidumbre formal de masa es demasiado optimista.
3. Los priors de elementos son heuristicas, no incertidumbres observacionales
   del objeto.
4. El baseline post-pre en Fase B solo corrige un offset AL medio; no captura
   drift temporal dentro de la ventana post-encuentro.
5. El batch depende de consultas TAP live para observaciones Gaia, salvo que se
   introduzca un cache explicito.

## Punto de insercion para Stage 1.2/1.3

El lugar correcto para el rediseño es un modulo nuevo
`src/mass/forward_model_joint.py`, no mas logica dentro del script CLI.

Interfaz recomendada:

```python
def predict_residuals_joint(
    params: np.ndarray,
    target_elements: dict,
    perturber_elements: dict,
    obs: GaiaObservationBundle,
    background_elements: dict[str, dict] | None = None,
) -> np.ndarray:
    ...
```

Parametros iniciales:

- `log10_M_perturber`
- `da_rel` o `da_au`
- `de`
- `di_deg`
- `dOmega_deg`
- `domega_deg`
- `dM_deg`

El wrapper debe aplicar los deltas al target, llamar al `forward_model`
existente y devolver residuales AL normalizados sobre una ventana comun. Esto
reutiliza la propagacion N-body ya validada y cambia solo la parametrizacion y
la likelihood del ajuste.

## Implicacion para criterios de aceptacion

Stage 1 no debe considerarse cerrado solo por producir masas. Debe demostrar
que el fit conjunto baja `chi2_red` de ventana de forma sustancial y que los
parametros orbitales absorben drift sin colapsar `M_perturber` contra los
bounds. Si `chi2_red` sigue en cientos, el siguiente blocker probable es la
covarianza Gaia/AL o sistematica observacional no modelada.
