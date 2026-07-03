# Dashboard status — Streamlit exploration app

> **Estado:** revisado y pulido (fase 7 del roadmap).
> **Fecha:** 2026-07-03.
> **Código:** `src/dashboard/app.py` (vistas Streamlit) · `src/dashboard/data.py`
> (capa de datos pura, sin Streamlit, testeable).
> **Arranque:** `docker compose up dashboard` → http://localhost:8501.

---

## Qué hace hoy el dashboard

Cuatro tabs sobre el catálogo Gaia DR3 de encuentros cercanos + la capa de masas:

1. **📦 Encounter Catalog** — catálogo caracterizado con filtros interactivos.
   - Carga memory-safe vía `src.dashboard.data`: prefiere
     `data/output/encounters_characterized_full.parquet` (72 M filas, ~5.8 GB);
     si no existe, cae a `encounters_characterized.parquet` (158 k).
   - Las métricas de cabecera (total, Gaia-observables, mínima aproximación)
     salen de una **agregación lazy en streaming** (`catalog_stats`) que nunca
     materializa el frame completo.
   - Las vistas interactivas cargan sólo los **300 000 encuentros más cercanos**
     (`load_catalog_display`, top-k lazy sobre `dist_au`); avisa cuando está
     capado. Columnas caracterizadas expuestas: `dist_au/dist_km`,
     `rel_vel_km_s`, `diameter_1/2_km`, `class_1/2`, `solar_elongation_deg`,
     `gaia_observable`.
   - Charts: histograma de distancias, velocidad vs separación (muestreo 10 k),
     encuentros por mes, torta por clase orbital. Tabla top-500 + descarga CSV.

2. **🔭 Novel Encounters** — subset candidato (`relevant_novel_encounters.csv`),
   ordenado por `deflection_score`; filtros por distancia/velocidad/observabilidad.

3. **📡 Gaia Coverage** — auditoría de cobertura de tránsitos por encuentro
   (`gaia_coverage_audit.csv`): pre/post tránsitos, viabilidad, gap al primer
   tránsito post-encuentro.

4. **⚖️ Asteroid Masses** — masas medidas + candidatos de seguimiento (ver abajo).

Todos los archivos referenciados **existen** en `data/output/` y el esquema del
catálogo caracterizado coincide con las columnas que la app consume. La app
**importa sin error** dentro de la imagen Docker (`streamlit`, `plotly`, `polars`
presentes; streamlit 1.58.0).

---

## Qué se arregló en esta revisión

1. **Tab de masas actualizada al estado real del proyecto.** Antes la app
   afirmaba "no publishable mass comes from this pipeline in DR3" y sólo mostraba
   dos fits LOO hardcodeados. Eso quedó obsoleto: el motor conjunto órbita+masa
   (`src/orbdet/`) produce masas medidas. Ahora la tab:
   - Carga `data/output/orbdet/mass_catalog_jack.csv` (16 perturbadores con σ por
     jackknife) vía las nuevas funciones puras `resolve_mass_catalog_path()` y
     `load_mass_catalog()` en `data.py` (fallback a `mass_catalog.csv`).
   - Muestra métricas (perturbadores ajustados, medidos con SNR_jack ≥ 3,
     calibradores), un scatter de `ratio_fit_over_ref` con barras de error
     jackknife coloreado por `mass_status`, tabla completa y descarga CSV.
   - Reencuadra los dos fits LOO viejos como **legacy / contexto**, no como
     resultado vigente (el método secuencial fue el límite, no la astrometría;
     ver `docs/mass_layer_track_a_closure.md`).

2. **Caption de cabecera corregida.** Ya no dice que la capa de masas "no es
   determinable en DR3"; ahora apunta a las masas medidas sobre Gaia FPR.

3. **Deprecación de Streamlit.** Migradas las 16 llamadas
   `use_container_width=True` → `width="stretch"` (el argumento viejo se elimina
   tras 2025-12-31 y ya emitía warning).

4. **Warning de serialización Arrow** en la tabla LOO legacy: la columna
   `ρ_lit` mezclaba `float` y `str` ("—"), disparando un fallback de pyarrow.
   Homogeneizada a string. Import ahora limpio, sin trazas.

---

## Mejoras propuestas (no implementadas, para no arriesgar)

- **Vista de encuentros notables.** `docs/notable_encounters.md` (jul 2026) tiene
  tablas grande-grande, NEA-MBA, y proxy de familia derivadas del catálogo
  caracterizado. Podría exponerse como una sub-sección o tab. No se implementó
  porque hoy son tablas Markdown curadas a mano; para hacerlo bien conviene
  primero volcarlas a un CSV/parquet estable (`data/output/notable_*.csv`) que la
  app pueda leer, en vez de parsear Markdown en runtime.
- **Cruce masas vs Fuentes-Muñoz (2025).** `docs/mass_crosscheck_jack.md` compara
  las 16 masas con la Tabla 5 de FM (z_formal vs z_jack). Si el cruce se
  materializa como CSV (`crosscheck_fuentes_munoz_jack.py`), añadir un panel de
  z-scores en la tab de masas sería directo y de alto valor científico.
- **Selector de catálogo.** Hoy `resolve_catalog_path()` elige automáticamente el
  mejor catálogo; un selector explícito (full 72 M vs 158 k vs hybrid stage-b)
  ayudaría a comparar runs, a costa de complejidad de memoria.
- **Column config de Streamlit** para formatear masas en notación científica y
  fijar unidades en los headers de tabla (km, km/s, AU) de forma nativa.

---

## Verificación

```bash
docker compose run --rm pipeline python -c "import streamlit, plotly; import src.dashboard.app"
# → IMPORT_OK, sin warnings de deprecación ni traza de Arrow
```

No se levantó el servidor (contexto de memoria ajustada); la verificación se
limitó a que el módulo importa y a que las funciones de carga apuntan a archivos
existentes con el esquema esperado.
