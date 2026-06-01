# Gaia FPR SSO data model — reconocimiento (FPR_INGEST_PLAN Stage 0)

> Resultado del probe read-only contra el archivo Gaia en vivo
> ([scripts/dev/probe_gaia_fpr.py](../scripts/dev/probe_gaia_fpr.py),
> corrido 2026-06-01). Reporte crudo en `data/output/gaia_fpr_probe.json`
> (gitignored). **Este doc es el gate de salida del Stage 0**: fija lo que el
> resto del plan puede asumir del schema FPR.

## TL;DR

FPR es **viable para reabrir masas** y casi un drop-in:

- Tabla: **`gaiafpr.sso_observation`** (mismo patrón que `gaiadr3.sso_observation`).
- Época de referencia: **idéntica a DR3** — `epoch = JD_TCB − 2455197.5`
  (días desde J2010.0 TCB). Confirmado: la primera época de (1) Ceres coincide a
  9 decimales entre DR3 y FPR (`1705.73661579…`), o sea es la misma observación
  con la misma codificación.
- Baseline: **~65.8 meses** (2014-07-26 → 2020-01-20) vs ~34 meses en DR3 → casi
  **2×**. Es exactamente el leverage extra que motivó este plan.
- Las **14 columnas que usa la capa de masas están todas presentes.** La única
  ausente respecto a DR3 es `g_mag` (FPR no incluye fotometría en esta tabla),
  que la capa de masas **no usa**.

## Tablas del schema `gaiafpr`

12 tablas. La relevante:

| tabla | rol |
|-------|-----|
| `gaiafpr.sso_observation` | **per-transit astrometry de asteroides** (la que usamos) |
| `gaiafpr.sso_source` | resumen por objeto (no por tránsito) |
| (otras 10) | lentes, ISM, variables, crowded fields — no SSO |

## Ventana temporal y volumen (tabla completa)

| métrica | valor |
|---------|-------|
| `epoch` min / max (días desde J2010 TCB) | 1667.279 / 3671.500 |
| JD_TCB min / max | 2456864.779 / 2458869.000 |
| fecha UTC min / max | **2014-07-26 → 2020-01-20** |
| baseline | 2004.2 días (**65.8 meses**) |
| filas (tránsitos) | **46,264,083** |
| asteroides distintos (`number_mp`) | **156,793** |
| `max(number_mp)` | **399,961** |

⚠ **`mp_max` real = 399,961**, no 160,000. La ingesta FPR debe barrer hasta
~400,000 en `number_mp` (no todos existen; son 156,793 distintos, pero el rango
llega a ~400k). El `_DEFAULT_MP_MAX=160000` de DR3 **no sirve** para FPR.

## Diff de columnas vs DR3 (las que consumimos)

DR3 usa 18 columnas (14 en la capa de masas + denomination/solution_id/source_id/g_mag).
Contra `gaiafpr.sso_observation`:

| columna (uso actual) | en FPR? |
|----------------------|---------|
| `solution_id`, `source_id`, `denomination`, `number_mp`, `transit_id`, `observation_id` | ✅ |
| `epoch`, `epoch_utc` | ✅ (+ `epoch_err` nuevo) |
| `ra`, `dec` | ✅ |
| `ra_error_systematic`, `dec_error_systematic`, `ra_dec_correlation_systematic` | ✅ |
| `ra_error_random`, `dec_error_random`, `ra_dec_correlation_random` | ✅ |
| `position_angle_scan` | ✅ |
| `x_gaia`, `y_gaia`, `z_gaia` | ✅ |
| **`g_mag`** | ❌ **ausente** (sin fotometría en la tabla FPR) |

**Columnas nuevas en FPR (bonus, no requeridas):** `vx_gaia/vy_gaia/vz_gaia`
(+ variantes `*_geocentric`) — velocidad de Gaia; `epoch_err`; `is_rejected`
(flag de calidad); `astrometric_outcome_ccd/_transit`; `fov`.

### Implicaciones de diseño

1. **`columns_override` por release es necesario** (justifica el campo del
   Stage 1): FPR baja el set base **menos `g_mag`**. La columna no existe → un
   SELECT con `g_mag` falla.
2. **La capa de masas no se ve afectada por la ausencia de `g_mag`**:
   `fetch_gaia_full` ([fit_mass_gaia_loo.py:114-122](../scripts/mass/fit_mass_gaia_loo.py#L114))
   no lo selecciona. Las 14 columnas que sí pide están todas en FPR.
3. **Estimación de diámetros** (`characterize`, usa H de MPCORB, no `g_mag` de
   Gaia) no se ve afectada. Si en algún momento se quisiera fotometría por
   tránsito en FPR habría que ir a otra tabla — fuera de scope (masas).
4. **`is_rejected`**: FPR expone un flag de rechazo astrométrico. La ingesta/capa
   de masas FPR debería filtrar `is_rejected = false` para quedarse con
   astrometría limpia (DR3 no tenía este flag explícito). Anotar para Stage 2/3.
5. **`mp_max = 400,000`** para el barrido de batches FPR (ver arriba).

## Sanity check: (1) Ceres DR3 vs FPR

| | DR3 | FPR |
|--|-----|-----|
| nº de épocas | 186 | **382** (≈2×) |
| baseline | 916 días | **1840 días** (≈2×) |
| primera época (días) | 1705.7366158 | 1705.7366158 (idéntica) |

Confirma los dos supuestos centrales del plan: (a) FPR duplica densidad y
baseline, (b) la codificación de `epoch` es la misma → la conversión TCB→TDB
existente ([src/utils/time_utils.py](../src/utils/time_utils.py)) sirve sin
cambios.

## Config resultante para FPR (insumo del Stage 1)

```yaml
fpr:
  table: "gaiafpr.sso_observation"
  epoch_ref_jd_tcb: 2455197.5          # idéntica a DR3
  window_start: "2014-07-26T00:00:00"
  window_end:   "2020-01-21T00:00:00"
  mp_max: 400000
  columns_override:                     # base − g_mag (+ is_rejected opcional)
    drop: ["g_mag"]
```

**Gate Stage 0: PASADO.** Sin bloqueantes para Stage 1.
