# FOLLOWUP_PLAN — correcciones de provenance y validación

> Plan acotado de tres ítems surgidos de una revisión crítica externa
> (2026-05-31). El `current_working_plan` (PRs #53–#59) está cerrado; estos son
> residuales de **reproducibilidad y documentación**, no de funcionalidad del
> pipeline. Ninguno reabre la capa de masas (cerrada definitivamente para DR3,
> ver [docs/mass_layer_track_a_closure.md](docs/mass_layer_track_a_closure.md)).

Contexto: la revisión reprodujo en su mayoría la autocrítica ya documentada en
[FROZEN_RUN.md](FROZEN_RUN.md) y los `docs/`. Estos tres ítems son lo
genuinamente accionable que quedó: dos inconsistencias de provenance y una
regresión de validación.

---

## Estado de implementación (2026-05-31) — los 3 ítems CERRADOS

| ítem | estado | entregable |
|------|--------|------------|
| 1 — prefilter declarado vs. efectivo | ✅ | `effective_prefilter_mode()` + writer registra `prefilter.effective`/`n_asteroids`/`max_n`; sidecar frozen parcheado a `skipped_large_n`; FROZEN_RUN.md §2/:92/:290 corregidos; 3 tests nuevos |
| 2 — techo del N-body | ✅ | `scripts/validate/measure_nbody_perturber_ceiling.py` + [docs/nbody_perturber_ceiling.md](docs/nbody_perturber_ceiling.md): medido \|Δdist\| mediana 1.3 μAU, máx **80 μAU** (≪ error Kepler 15.2 mAU) → perturber set **no** es el término dominante |
| 3 — regresión Fienga (804,733) | ✅ | `scripts/validate/diagnose_fienga_804_733.py`: hipótesis prefilter/cadencia **refutada**; es un **gap del artefacto** (código actual lo detecta a 0.013547 AU), recuperable re-corriendo; docs/literature_validation.md corregido |

**Hallazgo transversal (refuerza ítem 1):** el sidecar frozen es *backfilleado*
(`run_id="backfill_…"`, `git: {}`) y el parquet (2026-05-24) predata el commit
declarado `b1c4d9a` (2026-05-25) por ~24 h → el commit generador exacto no está
capturado. Anotado en FROZEN_RUN.md.

El detalle de cada ítem queda abajo como referencia del diseño original.

---

## Ítem 1 — Inconsistencia prefilter declarado vs. efectivo

**Severidad:** alta (afecta cómo se interpreta la completitud del catálogo).

**Problema.** Dos lugares afirman que el corte `|Δa| ≤ 0.5 AU` se aplicó al
freeze: el sidecar registra `prefilter.enabled = true` con
`semimajor_diff_max_au = 0.5`
([encounters_catalog_rebound_005au_provenance.json](data/output/encounters_catalog_rebound_005au_provenance.json)),
la fila de metadata de FROZEN_RUN.md lo declara `enabled`
([FROZEN_RUN.md:92](FROZEN_RUN.md#L92)), y un bullet afirma que pares con
`|Δa| > 0.5 AU` *"**are** missing"* del catálogo
([FROZEN_RUN.md:290](FROZEN_RUN.md#L290)). **Pero el código —pinneado a `main`
en commit `b1c4d9a`, [FROZEN_RUN.md:93](FROZEN_RUN.md#L93)— salta
`compatible_pairs` cuando `N > 5000`** y usa solo el filtro espacial KD-tree
([src/detect/pipeline.py:104-113](src/detect/pipeline.py#L104)):

```
_PREFILTER_MAX_N = 5_000
...
if n <= _PREFILTER_MAX_N:
    pairs = compatible_pairs(...)
else:
    # "skipping pair precomputation, KD-tree spatial filter only"
    pairs = None
```

Confirmado en log a `N=10000`:
[logs/kfn_full.log:8](logs/kfn_full.log) — *"N=10000 > 5000: skipping pair
precomputation, KD-tree spatial filter only"*. La corrida congelada es de
~150k cuerpos (72.2M encuentros, `tiered_mode = true`), muy por encima de 5000.

**Consecuencia.** El corte `|Δa| ≤ 0.5` **casi con certeza nunca se aplicó al
artefacto congelado**. El déficit de recall medido en
[docs/prefilter_recall.md](docs/prefilter_recall.md) (143,229 encuentros
perdidos, recall 76.38 %) cuantifica lo que el prefilter *perdería si se
aplicara*, no lo que el catálogo real perdió. A escala completa la completitud
se rige por el KD-tree espacial + la cadencia de grilla (12 h), ciego a `Δa`.

**Tareas.**
1. Confirmar el camino efectivo de la corrida congelada (revisar el log de
   `backfill_2026-05-25` o re-derivar el `N` efectivo del run).
2. En el writer de provenance ([src/catalog/writer.py:139](src/catalog/writer.py#L139)):
   registrar el **modo efectivo** del prefilter, no solo la config declarada.
   Añadir un campo `prefilter.effective` ∈ `{"applied", "skipped_large_n",
   "disabled"}` con el `N` y el umbral `_PREFILTER_MAX_N` usados.
3. Corregir FROZEN_RUN.md en sus tres puntos: la fila de metadata
   ([:92](FROZEN_RUN.md#L92)) y el bullet ([:290](FROZEN_RUN.md#L290)) deben
   decir que a `N > 5000` el corte `|Δa|` **no se aplica**; la afirmación
   *"are missing"* es contrafactual para el artefacto congelado. El caveat de
   recall de `prefilter_recall.md` mide daño potencial, no real. La completitud
   real está limitada por la cadencia del KD-tree (12 h), no por `Δa`.

**Criterio de aceptación.** El sidecar de cualquier corrida nueva distingue
prefilter declarado vs. efectivo; FROZEN_RUN.md §2 refleja que el catálogo de
72.2M no sufrió el corte `|Δa|`.

---

## Ítem 2 — Documentar el set de perturbers como techo de precisión

**Severidad:** media (precisión física; el budget de error ya está medido pero
contra el N-body propio, no contra una efeméride de referencia).

**Problema.** El scan N-body de la corrida congelada usa solo
**Sol + Júpiter + Saturno** como perturbers, sin Urano/Neptuno, planetas
terrestres ni asteroides mayores (`include_planets: [sun, jupiter, saturn]`,
`include_major_asteroids: false` en el sidecar). El budget de error de
[FROZEN_RUN.md](FROZEN_RUN.md) (Stage A/B) mide Kepler-refine **contra ese
N-body de 3 cuerpos**, no contra DE440/SPICE. Además
[src/propagate/nbody.py:12](src/propagate/nbody.py#L12) advierte que las
posiciones planetarias vienen del builtin de astropy como fallback interno.

**Consecuencia.** El "error N-body" reportado es relativo al modelo de 3
cuerpos, no al ground truth de una efeméride de precisión. El único chequeo
externo real contra DE44x es la validación JPL Horizons — que está
cadence-limited (caveat ya en [VALIDATION_SUMMARY.md](VALIDATION_SUMMARY.md)).

**Tareas.**
1. Documentar explícitamente en FROZEN_RUN.md y en el README del catálogo que
   el N-body es Sol+Júpiter+Saturno con efeméride builtin de astropy, y que el
   budget de error es *interno* (vs. ese modelo), no vs. DE440.
2. Cuantificar la cota superior del techo: re-refinar una muestra estratificada
   pequeña (~100 pares, cubriendo low-q / high-e / fast) con el set completo de
   planetas y comparar el desplazamiento de `dist_min`. Reportar el delta como
   incertidumbre sistemática del modelo de perturbers.

**Criterio de aceptación.** El catálogo declara su set de perturbers como
limitación de precisión nombrada, con una cota numérica del error sistemático
asociado.

---

## Ítem 3 — Regresión del par Fienga (804, 733)

**Severidad:** media (validación contra literatura).

**Problema.** El par Fienga `(804, 733)` está **ausente** del catálogo
frozen/hybrid aunque una corrida previa con efemérides JPL lo detectaba
([docs/literature_validation.md:54](docs/literature_validation.md#L54)).

**Tareas.**
1. Reproducir la detección del par con JPL Horizons para confirmar que el
   encuentro es real y obtener su `dist_min` y época de referencia.
2. Rastrear por qué cae del frozen/hybrid: ¿queda fuera por la cadencia de
   grilla (12 h) en el scan grueso, por el radio de query del KD-tree, o por el
   corte de umbral a 0.05 AU? Cruzar con el análisis de censura de
   [docs/kepler_threshold_bias_paper.md](docs/kepler_threshold_bias_paper.md).
3. Documentar la causa raíz. Si es censura de cadencia/umbral, registrarlo como
   falso negativo conocido; si es un bug del scan, abrir fix.

**Criterio de aceptación.** Causa raíz de la ausencia documentada en
`docs/literature_validation.md`, clasificada como censura conocida o como bug
con fix asociado.

---

## Fuera de alcance (explícito)

- **Capa de masas.** Cerrada para DR3; no se reabre. Cualquier crítica sobre
  no-identificabilidad, sesgos de calibradores o la ventana joint es correcta y
  ya está reflejada en el cierre de Track A. Sin acción.
- **Reframe candidato vs. catálogo final.** Ya hecho: el freeze se titula
  "candidate catalog", el híbrido lleva `refinement_method` por fila. Sin
  acción.
