# FPR_INGEST_PLAN — agregar Gaia FPR como fuente alternativa a DR3

> Plan para incorporar **Gaia FPR** (Focused Product Release, oct-2023) como
> fuente de astrometría de asteroides intercambiable con DR3, con el objetivo
> de **reabrir la capa de masas**. DR3 sigue siendo la fuente por default; FPR
> se selecciona con un flag. Ninguna etapa anterior del pipeline (detección,
> caracterización, catálogo) cambia su lógica: solo cambia de dónde salen las
> épocas y sobre qué ventana temporal opera.

Contexto: la capa de masas está cerrada **para DR3** (ver
[docs/mass_layer_track_a_closure.md](docs/mass_layer_track_a_closure.md) y
`project_mass_optimizer_bug` en memoria) — no por un bug sino por **leverage
físico insuficiente**: con ~34 meses de baseline la deflexión post-encuentro no
acumula señal sobre el ruido astrométrico (χ²(masa) multimodal, no
identificable, incluido Ceres). FPR ataca esa raíz: extiende el baseline a ~66
meses, que es exactamente la fuente que usa Fuentes-Muñoz (2024/2025) para sacar
231 masas — ya referenciada en la validación B1
(`project_literature_validation_b1`).

---

## Decisión de diseño: flag `release`, no fuente duplicada

**Opción elegida (A):** un selector `release: "dr3" | "fpr"` dentro de
`sources.gaia_sso`, con un sub-mapa `releases:` que guarda **solo** lo que
difiere entre versiones (nombre de tabla, época de referencia, ventana temporal,
mp_max, overrides de columnas). Los parámetros compartidos (archive_url,
batch_size, n_workers, columnas base) siguen siendo únicos.

**Por qué no la opción B (dos bloques `gaia_dr3` + `gaia_fpr` paralelos con un
`active_release` arriba):** duplica archive_url/batch/columnas, multiplica los
`_require(...)` en [config.py](src/utils/config.py) y obliga a tocar todos los
consumidores con `if release == ...`. El sub-mapa concentra la diferencia en un
solo lugar y deja el resto del código leyendo un único objeto "release activo".

**Por qué no un simple `table:` (lo que ya existe):** el campo `table` ya es
parametrizable, pero **no alcanza**. Cambiar de release arrastra cuatro cosas
acopladas que hoy están hardcodeadas en distintos lugares:
1. nombre de tabla TAP,
2. época de referencia del `epoch` (DR3 = días desde J2010.0 TCB = `2455197.5`),
3. ventana temporal de observación (DR3 = 2014-07-25 … 2017-05-28),
4. posibles diferencias de nombres de columnas.

Atarlas juntas bajo un nombre de release evita el estado inconsistente "tabla
FPR + ventana DR3".

### Aislamiento de artefactos (crítico)

Detección/caracterización/masas ingeridas bajo DR3 **no deben mezclarse** con
las de FPR. Todo path derivado se scopea por release:

- cache de ingesta: `data/cache/gaia_sso_chunks/{release}/`
- raw: `data/raw/gaia_sso_{release}.parquet`
- output de masas: `data/output/fits_{release}/fit_*.json`
- sidecar de provenance: registrar `gaia_release` + tabla + época de referencia.

Esto sigue la lección de provenance del FOLLOWUP_PLAN (sidecars que no capturan
el origen exacto generan ambigüedad irreproducible).

---

## Forma propuesta del config

```yaml
sources:
  gaia_sso:
    release: "dr3"                 # "dr3" | "fpr"  ← el flag
    archive_url: "https://gea.esac.esa.int/tap-server/tap"
    batch_size: 50000
    n_workers: 10
    max_retries: 10
    columns: [ ... ]               # columnas base compartidas
    releases:
      dr3:
        table: "gaiadr3.sso_observation"
        epoch_ref_jd_tcb: 2455197.5
        window_start: "2014-07-25T00:00:00"
        window_end:   "2017-05-28T00:00:00"
        mp_max: 160000
        columns_override: null     # usa las base
      fpr:
        table: "gaiafpr.sso_observation"   # ⚠ A VERIFICAR (Stage 0)
        epoch_ref_jd_tcb: 2455197.5        # ⚠ A VERIFICAR
        window_start: "2014-08-05T00:00:00"
        window_end:   "2020-01-21T00:00:00"  # ⚠ A VERIFICAR (~66 meses)
        mp_max: 160000                       # ⚠ A VERIFICAR
        columns_override: null               # ⚠ A VERIFICAR diffs de schema
```

Todos los valores de FPR marcados `⚠` son **hipótesis** hasta el Stage 0. No se
escribe ingesta sobre supuestos.

---

## Etapas

### Stage 0 — Reconocimiento del modelo de datos FPR (gate bloqueante)

**Por qué primero:** no asumir nada del schema FPR. Antes de escribir código hay
que confirmar contra el archivo Gaia en vivo.

Entregable: `scripts/dev/probe_gaia_fpr.py` (read-only, una corrida) que:
1. lista las tablas del schema `gaiafpr` vía TAP metadata,
2. dumpea las columnas de la tabla SSO de FPR y las diffea contra las 14 que hoy
   usamos en DR3 (`source_id, number_mp, epoch, ra, dec, x_gaia, y_gaia,
   z_gaia, *_error_systematic, *_error_random, *_correlation_*,
   position_angle_scan, g_mag`),
3. querea `MIN(epoch), MAX(epoch), COUNT(*), COUNT(DISTINCT number_mp)` para
   fijar empíricamente la ventana temporal real y `mp_max`,
4. trae las observaciones de **(1) Ceres** en ambos releases y compara nº de
   épocas y rango temporal (sanity: FPR debe tener más épocas y baseline más
   largo).

**Gate de salida:** documento `docs/gaia_fpr_data_model.md` con tabla exacta,
nombres de columnas, época de referencia confirmada y ventana temporal medida.
Si algún nombre de columna difiere, queda mapeado acá. **Sin este doc no se
avanza.** Decisiones de Stage 1–3 dependen de sus hallazgos.

### Stage 1 — Capa de config (release selector)

Archivos: [src/utils/config.py](src/utils/config.py),
[config.yaml](config.yaml).

- Nuevo dataclass `GaiaReleaseConfig` (table, epoch_ref_jd_tcb, window_start,
  window_end, mp_max, columns_override).
- `GaiaSSOSourceConfig` gana `release: str` y `releases: dict[str,
  GaiaReleaseConfig]`. Mantener `table` opcional/legacy para no romper configs
  viejos (deprecado, leído solo si no hay `releases`).
- Helper `GaiaSSOSourceConfig.active() -> GaiaReleaseConfig` que resuelve el
  release activo y valida que exista en el sub-mapa.
- `_build()` / `_require()` actualizados; tests de que un config sin `releases`
  (formato viejo) sigue validando.

**Gate:** `load_config()` resuelve `cfg.sources.gaia_sso.active()` para `dr3` y
`fpr`; test unitario nuevo en `tests/`.

### Stage 2 — Ingesta parametrizada

Archivos: [src/ingest/gaia_sso.py](src/ingest/gaia_sso.py),
[scripts/ingest/download_gaia_sso.py](scripts/ingest/download_gaia_sso.py).

- `download_gaia_sso(...)` recibe el `GaiaReleaseConfig` activo: tabla, época de
  referencia, ventana, mp_max, columnas efectivas (base + override).
- Cache y raw scopeados por release (ver "Aislamiento de artefactos"). Las
  chunks DR3 existentes no se invalidan (van a `…/dr3/`).
- El módulo deja de hardcodear `2455197.5` y `_DEFAULT_MP_MAX=160000`; los toma
  del release. Docstring actualizado (hoy afirma "DR3" en el header).
- `download_gaia_sso.py` agrega `--release` (default = el del config) para poder
  bajar ambos sin editar el yaml.

**Gate:** `download_gaia_sso --release fpr` baja Ceres+Vesta+Pallas+Hygiea y
escribe `data/raw/gaia_sso_fpr.parquet`; nº de épocas por cuerpo > que en DR3
(consistente con Stage 0).

### Stage 3 — Capa de masas sobre el release activo

Archivos: [scripts/mass/fit_mass_gaia_loo.py](scripts/mass/fit_mass_gaia_loo.py)
y los consumidores que querean Gaia directo
([fit_mass_gaia_joint.py](scripts/mass/fit_mass_gaia_joint.py),
[fit_mass_gaia_multitarget.py](scripts/mass/fit_mass_gaia_multitarget.py),
[analyze_mass_candidates.py](scripts/mass/analyze_mass_candidates.py)).

- `fetch_gaia_full(archive_url, target)` →
  `fetch_gaia_full(release_cfg, target)`: construye el ADQL con la tabla del
  release y la ventana del release; deja de hardcodear `gaiadr3.sso_observation`
  ([fit_mass_gaia_loo.py:118](scripts/mass/fit_mass_gaia_loo.py#L118)).
- Las constantes `_J2010_TCB_JD`, `_GAIA_START_JD_TCB`, `_GAIA_END_JD_TCB`
  ([:74-76](scripts/mass/fit_mass_gaia_loo.py#L74)) pasan a leerse del release
  activo (la conversión TCB→TDB de [time_utils.py](src/utils/time_utils.py) no
  cambia).
- Output `fit_*.json` scopeado a `data/output/fits_{release}/` y con
  `gaia_release` en el JSON.
- `--release` en los entry points de masas.

**Gate:** el fit corre end-to-end sobre FPR para un encuentro mid-mission y el
JSON registra el release correcto. (Este gate es de *plumbing*, no de física —
ver Stage 4.)

### Stage 4 — Re-validación de la física de masas (el gate que importa)

El plumbing funcionando **no** reabre nada por sí solo. Hay que re-correr el
diagnóstico que cerró DR3 y ver si FPR efectivamente da identificabilidad:

1. **Window-scan de χ²(masa)** sobre Ceres/Vesta/Pallas/Hygiea con datos FPR
   (mismo análisis que mostró el mínimo multimodal en DR3,
   `project_mass_optimizer_bug`). **Criterio de reapertura:** mínimo único y
   estable, ratio fit/lit dentro de ~2× para al menos los 4 cuerpos grandes.
2. **Specificity test** (`run_specificity_test.py`): la masa fiteada debe caer
   cuando se le asigna al perturber equivocado. En DR3 dio 0/27.
3. **Cruce con Fuentes-Muñoz 2025**: nuestras masas FPR vs sus 231. Es la
   comparación apples-to-apples (misma fuente).

**Gate de salida (honesto):** o bien las masas se vuelven defendibles bajo FPR
(y entonces se documenta la reapertura), o bien **no** y se documenta por qué
FPR tampoco alcanza con nuestra metodología — sin retractaciones tardías estilo
(111) Ate (`project_ate_mass_result`).

---

## Lo que NO toca este plan

- Detección, propagación, caracterización, dashboard: agnósticos a la fuente. Si
  se quisiera **detectar** sobre la ventana extendida de FPR, es solo cambiar
  `time_window` en config — fuera de scope acá (este plan es para masas).
- El catálogo congelado de 72.2M (DR3) queda intacto.

## Riesgos / gotchas específicos

- ⚠ **Propagación Kepler a 66 meses**: la época de MPCORB propagada con Kepler
  puro acumula más error en 5.5 años que en 3 (gotcha ya en CLAUDE.md). Para la
  capa de masas apoyarse en el refinamiento N-body (Track 1, ya publicable), no
  en Kepler puro.
- ⚠ **Frame/época de referencia FPR**: confirmar en Stage 0 que `epoch` usa la
  misma referencia J2010.0 TCB. Si difiere, es un error sistemático de segundos
  → sub-AU. No asumir.
- ⚠ **Schema drift**: FPR puede renombrar o agregar columnas de error/correlación
  respecto a DR3. El `columns_override` y el mapeo del Stage 0 lo absorben.
- ⚠ **No mezclar releases en un mismo análisis**: el aislamiento de paths lo
  fuerza, pero todo output debe llevar `gaia_release` para auditoría.

## Orden de ejecución

`Stage 0 (probe + doc)` → `Stage 1 (config)` → `Stage 2 (ingesta)` →
`Stage 3 (plumbing masas)` → `Stage 4 (re-validación física)`.

Stage 0 es bloqueante para todo lo demás. Cada stage es un PR separado hacia
`main` siguiendo el flujo habitual (branch → cambios → OK → commit → PR).
