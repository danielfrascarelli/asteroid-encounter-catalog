# F4 — Perturbadores fuera de los 16 de la efeméride: diseño

> **Estado:** ⬜ DISEÑO (sin implementar). 2026-06-30.
> Diseño del ítem F4 de [`planning/MASS_FUTURE_WORK.md`](../planning/MASS_FUTURE_WORK.md):
> ajustar la masa de un perturbador que **no** está en `sb441-n16.bsp`.
> **Gate:** ≥ 1 perturbador fuera de los 16 ajustado con χ²_red ∈ [0.95, 1.05].

## Hallazgo central (corrige el scoping previo)

Un scoping inicial concluyó que F4 requería extender ASSIST con una API para
inyectar partículas masivas con órbita propia, lo cual no existe. **Eso es
incorrecto.** El motor ya integra perturbadores asteroidales arbitrarios desde
elementos; no usa la fuerza `ASTEROIDS` de la efeméride.

Evidencia en [`src/orbdet/dynamics_assist.py`](../src/orbdet/dynamics_assist.py)
`_build_sim` (líneas 121–160):

```python
for ap in asteroid_perturbers:
    r, v = _bary_eq_state_from_elements(ap.elements, sun_xyz, sun_vxyz)
    sim.add(m=ap.mass_msun * GM_SUN, x=r[0], ...)      # partícula masiva de rebound
...
extras.forces = ["SUN", "PLANETS", "GR_EIH"]           # asteroides NO en las fuerzas ASSIST
sim.gravity = "basic"; sim.N_active = len(asteroid_perturbers)
```

Es decir:
- Cada perturbador se agrega como **partícula masiva de rebound**, inicializada
  desde sus elementos en la época, e **integrada** bajo Sol+planetas+GR (vía
  ASSIST) + gravedad mutua entre perturbadores. **No** se lee de la efeméride
  durante la integración; sólo el estado inicial sale de ahí.
- `AsteroidPerturber(name, mass_msun, elements)` acepta **cualquier** órbita. Que
  hoy los 16 vengan de la efeméride es una decisión de
  `big_asteroid_perturbers` ([`dynamics_assist.py:209`](../src/orbdet/dynamics_assist.py)),
  no una restricción del integrador.
- La fuerza `ASTEROIDS` de ASSIST está **excluida** a propósito (líneas 12–16),
  así que **no hay doble conteo**: un 17º cuerpo que no está en `sb441-n16` no
  aparece en ninguna fuerza de la efeméride. Se agrega como una 17ª partícula
  masiva y listo.

**Conclusión:** F4 no toca el motor (`src/orbdet/`). Es wiring en la capa IO
([`scripts/mass/orbdet_fit_realdata.py`](../scripts/mass/orbdet_fit_realdata.py)).

## Qué bloquea hoy un 17º perturbador (todo en la capa IO)

1. `_ephem_name_for_perturber` (líneas 119–135) **lanza** si el número no está en
   `_EPHEM_NAME_BY_NUMBER` (los 16). Es la única barrera dura.
2. `studied = big_asteroid_perturbers(common_epoch, names=(pname,))[0]` (línea 497)
   obtiene órbita **y** masa-semilla del perturbador **de la efeméride**. Para un
   cuerpo fuera de los 16 no hay entrada → hay que tomar la órbita y la semilla de
   otra fuente.
3. El fondo `big_asteroid_perturbers(common_epoch, exclude=(pname,))` (línea 502)
   con `pname` fuera de los 16 ya devuelve los **16 completos** (el `exclude` no
   matchea) — que es exactamente lo correcto. **No requiere cambio.**

## Diseño

### Órbita del perturbador (fija; sólo la masa es libre)

El perturbador entra al ajuste con su órbita **fija** (igual que los 16; sus 6
elementos no son parámetros), sólo su masa es libre. Fuente de la órbita, en
orden de preferencia:

- **JPL Horizons** (recomendado): estado osculador en la época común vía
  `_horizons_elements(target, epoch)`
  ([`scripts/validate/validate_assist_horizons.py:43`](../scripts/validate/validate_assist_horizons.py)).
  Órbita de calidad JPL, integrada luego bajo nuestro modelo — el mismo
  tratamiento que los 16. Una sola query por corrida (barato).
- **MPCORB** (fallback offline): fila del snapshot vía `load_element_rows`
  ([`scripts/mass/fit_mass_gaia_loo.py:214`](../scripts/mass/fit_mass_gaia_loo.py))
  → `elements_from_mpcorb` → `propagate_elements` a la época común. Elementos
  osculadores de 2 cuerpos; menos exactos que Horizons sobre el arco.

### Masa-semilla

- `--seed-mass-kg` explícito (ya existe), **o**
- estimación desde la magnitud absoluta `H` de MPCORB con albedo+densidad
  asumidos, reusando `_mass_from_h`
  ([`scripts/mass/fit_mass_gaia_loo.py`](../scripts/mass/fit_mass_gaia_loo.py),
  usado en línea 736) y los helpers de
  [`src/characterize/physical.py`](../src/characterize/physical.py).

La semilla afecta sobre todo la convergencia, no el mínimo. Si el cuerpo es
débil, la identificabilidad la juzga F1/F2 (σ_jack + `mass_status`): un 17º con
deflexión bajo el ruido saldrá `not_identifiable`, igual que los 6 débiles de los
16. Eso es una salida válida del gate (medida **o** cota explícita).

### Fondo y no-doble-conteo

Fondo = los 16 de `sb441-n16` completos (no se excluye nada, porque el estudiado
no está entre ellos). El estudiado es la 17ª partícula masiva. `N_active = 17`.
ASSIST sigue aportando Sol/planetas/GR; ningún asteroide está en las fuerzas de
la efeméride → sin doble conteo.

### Objetivos

Sin cambios: `_read_targets_from_catalog` filtra encuentros < 0.05 AU donde
`number_1` o `number_2` == perturbador, para cualquier número.

## Cambios concretos (todos en `orbdet_fit_realdata.py`)

1. **`_custom_perturber(number, common_epoch, snapshot, args) -> (KeplerElements,
   seed_mass_msun, name)`**: resuelve órbita (Horizons o MPCORB según
   `--perturber-orbit-source`) y masa-semilla (`--seed-mass-kg` o `_mass_from_h`).
2. **`_run_perturber`**: envolver `_ephem_name_for_perturber` en try/except; si
   lanza (cuerpo fuera de los 16) → rama custom:
   - `perturber_elements, seed_mass_msun, pname = _custom_perturber(...)`
   - `background = big_asteroid_perturbers(common_epoch)` (los 16 completos)
   - el resto (fetch de objetivos, `_fit_with_rejection`, jackknife F1, salida)
     es idéntico.
3. **Flags nuevos**: `--perturber-orbit-source {horizons,mpcorb}` (default
   `horizons`), `--perturber-albedo`, `--perturber-density` para la semilla por H.
4. **Etiqueta de salida**: el JSON usa `perturber_name`; para customs sin nombre de
   efeméride, usar el nombre de MPCORB o el número.

Nada de esto toca `src/orbdet/`. El contrato de aislamiento del motor se mantiene.

## Selección del candidato para el gate

Hace falta **un** perturbador fuera de los 16 con suficientes objetivos y señal.
Criterio: grandes (D ≳ 150 km) del cinturón principal, no en `BIG_ASTEROIDS`, con
muchos encuentros < 0.05 AU contra objetivos pequeños en el catálogo congelado.

Candidatos clásicos de determinación de masa por encuentros (todos numerados,
fuera de los 16): **(24) Themis, (19) Fortuna, (29) Amphitrite, (354) Eleonora,
(532) Herculina, (45) Eugenia, (13) Egeria, (48) Doris**. Fortuna ya aparece como
referencia de masa de Goffin en
[`src/propagate/nbody.py`](../src/propagate/nbody.py).

Selección operativa: contar por número (no-16) los encuentros < 0.05 AU con
objetivos < 100k en `encounters_catalog_hybrid_stageb.parquet` y elegir el de
mayor `N` (mismo conteo que hace `_read_targets_from_catalog`). Puede ser un
helper chico `scripts/mass/find_extra_perturber_candidates.py` o una query ad-hoc.

## Validación de la órbita del perturbador (evidencia para el gate)

Antes de confiar en la masa, confirmar que la órbita fija del 17º integrada bajo
nuestro modelo sigue a Horizons a nivel mas sobre el arco (~900 d), reusando
`validate_assist_horizons` (propaga desde elementos y compara vs vectores
Horizons). Es el mismo gate de exactitud T2/T8 ya verde para los 16
(ASSIST vs Horizons 0.17 mas / 900 d).

## Limitaciones y trabajo posterior (F4b)

- **Órbita del perturbador fija.** Su incertidumbre **no** se propaga a σ(masa).
  Para los 16 es despreciable (órbitas de efeméride); para un 17º de MPCORB/Horizons
  es mayor. Mitigación inmediata: usar Horizons y validar la órbita; cota honesta
  vía F1/F2 si la masa no es identificable.
- **F4b (futuro): ajuste conjunto de la órbita del perturbador.** El 17º es un
  asteroide numerado **observado por Gaia**: se podría incluir su propia
  astrometría como un objetivo más cuya órbita se ajusta, a la vez que su masa
  perturba a los demás (tratamiento completo Fuentes-Muñoz). Agrega 6 parámetros y
  un bloque de observación para el perturbador; mayor rigor, fuera del MVP del gate.

## Resumen del gate

| Paso | Entregable | Verificación |
|------|-----------|--------------|
| 1 | Rama custom en `orbdet_fit_realdata.py` | corre un nº fuera de los 16 sin lanzar |
| 2 | Órbita del 17º vs Horizons | < 1 mas/900 d (reusa validate_assist_horizons) |
| 3 | Ajuste de masa del 17º | χ²_red ∈ [0.95, 1.05] |
| 4 | σ externa | σ_jack (F1) + `mass_status` (F2): medida o cota |

## Esfuerzo estimado

~1 sesión. Sin cambios de motor; reusa Horizons, `_mass_from_h`, `load_element_rows`,
`big_asteroid_perturbers`, la maquinaria de fit y el jackknife ya existentes. El
grueso es la rama custom + selección del candidato + correr/validar.
