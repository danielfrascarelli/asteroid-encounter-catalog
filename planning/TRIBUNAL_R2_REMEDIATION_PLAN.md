# Plan de remediación — Tribunal ronda 2 (informe Fable 2026-07-09)

> **Estado:** 🟡 CASI COMPLETO — 4/4 bloqueantes cerrados salvo depósito author-owned (B4); 6/6 mayores ✅; N1–N16 ✅. Falta solo Tarea 10 (ORCID + Zenodo + DOI/SHA-256) y la aprobación de commit.
> **Última actualización:** 2026-07-09
> Plan para cerrar los hallazgos de la segunda ronda del tribunal científico
> (`docs/tribunal/2026-07-09_tribunal_fable_informe.md`): **4 bloqueantes, 6 mayores,
> 16 menores**. Veredicto de esa ronda: *revisión mayor*. **Criterio de éxito:** los 4
> bloqueantes cerrados con su gate verificable, los 6 mayores resueltos, y un manuscrito
> donde todos los números provienen de un único freeze (b1fix + Stage B regenerado).
>
> **Fuente de verdad de los hallazgos:** `docs/tribunal/2026-07-09_tribunal_fable_informe.md`
> (generado con el prompt `docs/tribunal/2026-07-09_tribunal_fable_prompt.md`, auditoría
> solo-texto del manuscrito: el tribunal NO tuvo acceso al código — toda "hipótesis de
> ubicación" del informe debe confirmarse contra el código antes de actuar).
> Este tracker apunta a esa evidencia, no la duplica.
>
> La ronda 1 (`docs/tribunal/2026-07-04_tribunal_cientifico.md`) está **cerrada** —
> tracker: `planning/TRIBUNAL_REMEDIATION_PLAN.md` (✅ COMPLETO, PR #102). No repetir
> esa auditoría.

---

## Tabla de estado

| # | tarea | hallazgos | estado | depende de | entregable / evidencia |
|---|-------|-----------|--------|-----------|------------------------|
| 1 | Regenerar Stage B híbrido sobre el freeze b1fix + refrescar TODOS los números derivados | B3 | ✅ | — | **Stage B COMPLETO + híbrido ensamblado** (`..._hybrid_stageb_b1fix.parquet`, 80 072 774, gate n_missing=0, drift 5.6e-14). Números refrescados en el tex: §2.2 frágil 10 339 493 (12.91 %); §3.1 Stage B max **14.83 mAU** (era 15.2) + subset **10.34M** (era 8.73M), p99 1.99 confirmado; §3.2 upward **30 109** (era 25 283; 0.29 % s/cambio); §3.4 FM **13 770/40 176 = 34.3 %** (era 29.5 %), 26 406 no-recup (91.3 % domain / 8.7 % fuera / 1 >umbral) — TODO §3.4 eliminado. **Solo 2 TODOs vivos (ORCID, DOI)**. abstract/Concl. citan **Stage A** (12 µAU/2.5 mAU) → sin cambio (validación propagador, freeze-robusta; Stage A NO re-corrido). Fig.4 re-renderizada del freeze (07-11) y sincronizada a docs/figures; Tabla 6/Fig.3 cerradas en Tarea 2. **CERRADO** |
| 2 | Reparar join de datos físicos/elementos: Tabla 6, columna class, Fig. 3 | B1, M5 | ✅ | — | código listo (`_supplement_elements` MPCORB 100 %; Fig.3→MPCORB; tex caption/texto); gate ρ escrito (`scripts/validate/check_table6_density.py`, `--demo` caza B1: Nysa ρ=0.40, Fortuna ρ=6.90 FAIL). **CERRADO 2026-07-11**: recls promovido (80 072 774, class 100 % en numerados<1000), Tabla 6 regenerada (D medidos, clases MBA), gate ρ 15/16 PASS, Fig.3 100 % (449213/449213) — ver Tanda 3 |
| 3 | Auditar clamps de borde en Tablas 1–2 (fechas 2017-05-29) + guard de borde uniforme | B2 | 🟡 | 1 | Diagnóstico: las 2 filas están a JD_fin **+1.25 d** (2017-05-29, `near_boundary`=False, mínimos genuinos ~30 h fuera de la ventana; la grilla padea ~1.75 d). Tabla 1 → **10 filas** (quitadas Philomela×Nephele e Io×Adeona); "twelve"→"ten" en texto+caption+nota; Tabla 2 (dd50) verificada in-window; Tablas 3–5 con guard >130 d. Convención + flag `edge_censored` documentados en §2.1. **Falta:** materializar la columna `edge_censored` en el parquet depositado (con Tarea 10) |
| 4 | Reescribir §2.2: ventana de refinamiento realmente ejecutada + guard N-body ±12 h | M1 | ✅ | — | §2.2 reescrito (±6h=±Δt/2 + re-centrado iterativo + dedup entre bloques + modo fallo pre-fix); guard N-body `near_boundary` declarado; **injection-recovery 300/300 = 100 % PASS**; check "cero mínimos en el borde": los 12 eventos dd100 tienen `near_boundary`=False y Tablas 3–5 con guard >130 d → ningún evento tabulado en el borde de refinamiento ✅ |
| 5 | MC de σ de distancias con covarianza de época consistente | M2 | ✅ | 1 | vía (a): MC sobre `extreme_pairs_b1fix.json`. σ_d formal chica (1.9–7.5 % cercanas, <5e-4 % lentas, **ninguna >30 %**) PERO el corrimiento nominal MPCORB-2016→SBDB-2026 domina en la cola cercana: **25–382 %** (238813×391704: 1320→6358 km). Tabla 3 con **d₂₀₁₆/d₂₀₂₆/σ_d**; §3/§sec:extreme reescritos (época domina, no la σ formal; quitado el claim "el temor 10³–10⁴ km no aplica"); pasaje de familia alineado. Cola cercana = ilustrativa, no precisa; eventos anchos/lentos <3 % |
| 6 | Aclarar diseño del censo de umbral (estratificado vs aleatorio); rehacer extrapolación si corresponde | M4 | ✅ | — | §3.2 corregido (muestreo **aleatorio** uniforme, no estratificado); `scripts/validate/check_threshold_census_scaling.py` **PASS 0.46 %** (39 456×2020 = 7.97e7 vs 8.01e7) |
| 7 | Europa: discusión explícita, claim de determinabilidad acotado, cita SG22, ratios FM tabulados | M3, N9, N14a | ✅ | — | §5 reescrito: Europa z=−4.3 es vs **seed** DE441 (4.0e19, outlier); vs FM25 z=−0.8 (0.88), vs SG22 z=−1.7 (0.78, DR2+terrestre = independiente del FPR); 0.85 vs FM25 / 0.75 vs seed reconciliado; abstract N2 (Pallas bound); `siltala2022` en .bib; **`scripts/mass/verify_mass_status_rule.py` PASS 16/16**; compila 13 pp. Incidental: N7 y N3 cerrados |
| 8 | Cuarto término del presupuesto de incompletitud (membresía Kepler vs candidacy N-body) | M6 | ✅ | 1 | Sobre shards completos (10.34M): fracción \|Δd\|>7.2mAU en banda [0.045,0.05) = 4.6e-5; miembros near-threshold en el híbrido = 15.3M → **cota ≲704 pares** (orden 10², 2–3 órdenes bajo el ~10⁵ del censado, subdominante). §3.2 con **cuarto término** + script `measure_stageb_dd_distribution.py` |
| 9 | Pasada editorial única: N1–N16 (texto, bib, captions) | N1–N16 | ✅ | 1,2,3,5,7 (números finales) | adelantados ✅: **N3, N4, N6, N7, N8, N9, N10, N11, N12, N13, N14a, N14b, N14c, N15, N16a, N16b, N16c** (references.bib 19→25: +siltala2022, tedesco2002, usui2011, mainzer2011, carry2012, bernstein2025). **N1** cerrado (perturber-truncation vs Kepler: "~100×" → "~10× en mediana, ~40× en p99"). Solo queda **N5** (nombres .parquet del depósito, con Tarea 10). Suite verde (pytest 586) al 2026-07-10 |
| 10 | Submission: DOI Zenodo/VizieR, afiliación, ORCID, funding | B4 | 🟡 | 1–9 | Bloque de autor **completo en el tex** (2026-07-10): D.~Frascarelli, "Independent researcher, Montevideo, Uruguay", `dsanfra@gmail.com`; **sin financiamiento** (cómputo propio, en agradecimientos). Falta: **ORCID** (opcional, A&A-recomendado; el IEEE Author ID 37085466644 NO sirve para A&A/ADS — usuario evaluará crear ORCID en orcid.org), cuenta de depósito (Zenodo vs VizieR), y DOI + SHA-256 sobre el freeze final. TODOs vivos: `:30` (ORCID), `:1155` (DOI) |

**Regla de orden (del propio informe):** la Tarea 1 refresca los números de §2.2, §3.1,
§3.2 y §3.4 de los que dependen las tareas 3, 8 y el abstract — **ninguna edición de
números del paper antes de que la Tarea 1 cierre**. Las tareas 2, 4, 6 y 7 son
independientes y pueden arrancar en paralelo con la 1 (la 2 y la 6 tienen parte de
investigación de código que no toca números del catálogo).

---

## Contexto mínimo para el agente que tome esto

- **Manuscrito:** `docs/paper/aa_encounters.tex` (12 pp., compila limpio) +
  `docs/paper/references.bib` (19 entradas). Kit A&A (`aa.cls`/`aa.bst`) ya en
  `docs/paper/`. Se compila **en el host** (`docs/` no está montado en el contenedor):
  `cd docs/paper && pdflatex aa_encounters && bibtex aa_encounters && pdflatex aa_encounters && pdflatex aa_encounters`.
- **Branch de trabajo:** `fix/tribunal-remediation-core` (PR #102, abierta y MERGEABLE).
  Esta ronda 2 corrige el mismo paper de ese PR → **seguir commiteando ahí** mientras
  #102 siga abierto. Si el usuario lo mergea antes, crear branch nueva desde `main`
  (p. ej. `fix/tribunal-r2`) — **nunca** basar una branch en otra branch de PR abierta
  distinta de la propia.
- **Estado git pendiente al 2026-07-09:** `docs/tribunal/` está sin commitear (los dos
  prompts + el informe + el doc de ronda 1 movido desde `docs/tribunal_cientifico_2026-07-04.md`)
  y `config.local.yaml.bak_tribunal` está staged. Proponer al usuario commitear ese
  housekeeping primero, separado del trabajo científico.
- **Datos (verificado 2026-07-09):**
  - `data/output/encounters_catalog_rebound_005au_b1fix.parquet` — freeze B1-fixed,
    80.072.774 filas (2026-07-07). Es el catálogo base bueno.
  - `data/output/encounters_characterized_b1fix.parquet` — caracterizado b1fix (2026-07-08).
  - `data/output/encounters_catalog_hybrid_stageb.parquet` — **VIEJO** (2026-05-28,
    construido sobre `encounters_catalog_rebound_005au.parquet` pre-B1; su sidecar lo
    confirma). Es la raíz de B3. Los 874 shards de `data/output/stageb_nbody_shards/`
    son de esa corrida vieja — **no reutilizables** para la Tarea 1.
  - `config.local.yaml` (gitignored) apunta el output a `..._b1fix`; backup del overlay
    original en `config.local.yaml.bak_tribunal`.
  - Los parquet generados por el contenedor quedan **owned by root** — renombrar/mover
    dentro del contenedor, no con `mv` del host.

---

## Forma de presentar el progreso

Seguir `planning/TASK_TRACKING_GUIDE.md`. Concretamente:

1. **Este archivo es la única fuente de verdad** del avance de la ronda 2. Al avanzar:
   actualizar la tabla de estado (⬜→🟡→✅), con enlace a la evidencia (commit, script,
   número de gate) en la columna "entregable / evidencia". No duplicar la evidencia acá:
   apuntar a ella.
2. **Cada cierre de tarea cita el resultado de su gate** (el número exacto, no "pasa").
   Los gates de este plan son copia de los "criterios de validación de cierre" del
   informe — si un gate se relaja, marcar la decisión como tal (regla 6 de la guía:
   marcar, no borrar).
3. **Registro de sesión ("tandas"):** añadir al final de este archivo un apéndice
   `## Tanda N (YYYY-MM-DD)` por sesión de trabajo, estilo
   `planning/TRIBUNAL_REMEDIATION_PLAN.md` §A1–A6: qué se hizo, qué corrida quedó
   lanzada (log, output, container), qué quedó verificado con qué número.
4. **Fechas absolutas** siempre (`2026-07-09`, nunca "ayer").
5. **Respuestas al usuario:** prefijar cada respuesta con la fecha `YYYYMMDD-`
   (preferencia del usuario). Respuestas cortas y directas, sin resumen final de
   "lo que hice fue…".
6. **Git:** autor `Daniel Frascarelli <dsanfra@gmail.com>`; **sin `Co-Authored-By` ni
   footer de agente** en commits ni PRs. **No commitear ni pushear sin confirmación
   explícita del usuario** — presentar el diff y esperar el OK. No `--no-verify`.
   Al terminar, marcar este plan ✅ COMPLETO en la cabecera (no dejarlo leyéndose
   como activo).

---

## Docker: cómo ejecutar (engine NATIVO, memoria SIEMPRE topeada)

⚠️ **REGLA DURA: el host nunca puede morir.** El entorno usa el **Docker engine
nativo** (contexto `default`, verificado; ya no la VM de Docker Desktop). Sin límite
de memoria, un contenedor pesado usa toda la RAM + swap del host y **ya mató la PC
una vez** (2026-07-05). Host: 31 GB → cap 18 GB deja ≥10 GB al sistema.

- **Corridas normales:** `docker compose run --rm pipeline python -m <módulo> ...`
  — el `docker-compose.yml` ya trae `mem_limit: 18g` + `memswap_limit: 18g`
  (memswap = mem ⇒ sin swap). Es la vía preferida: el cap viene puesto.
- **Corridas largas (horas/días), detached** — patrón usado en la ronda 1:

  ```bash
  # Lanzar detached (compose hereda el cap de 18g). Si el shell no está en el
  # grupo docker, envolver cada comando en:  sg docker -c '<comando>'
  CID=$(docker compose run --rm -d pipeline \
        python -m scripts.validate.refine_stageb_nbody <args>)
  echo "$CID" > logs/stageb_regen_container_id.txt
  docker logs -f "$CID" >> logs/stageb_regen_$(date +%Y%m%d).log 2>&1 &

  # Waiter en background para enterarse del fin sin poll manual:
  docker wait "$CID" && echo "terminó con exit $(docker inspect -f '{{.State.ExitCode}}' "$CID")"
  ```

  Los runs detached sobreviven al cierre del shell/sesión.
- **Si se usa `docker run` crudo** (fuera de compose), es **obligatorio** replicar el
  cap a mano — nunca lanzarlo sin esto:

  ```bash
  docker run -d --name <nombre> --memory=18g --memory-swap=18g \
    -v "$PWD/data:/app/data" -v "$PWD/logs:/app/logs" \
    -v "$PWD/config.yaml:/app/config.yaml:ro" \
    -v "$PWD/config.local.yaml:/app/config.local.yaml:ro" \
    -v "$PWD/src:/app/src" -v "$PWD/scripts:/app/scripts" \
    gaia-asteroid-encounters python -m <módulo> <args>
  ```
- **Monitoreo:** `docker stats --no-stream` — MEM debe quedar bajo el cap.
- **Exit 137 = OOM del contenedor.** Es el resultado *aceptable* (host vivo). La
  respuesta correcta es **bajar workers** (probado: 2–4 workers; con el scan, 4 fue
  más rápido que 8) y/o usar el modo out-of-core, **no subir el cap**.
- **Tras editar código y antes de confiar en una corrida:** verificar
  `docker compose run --rm pipeline md5sum <archivo>` == `md5sum <archivo>` en el
  host. (Gotcha heredado de Docker Desktop; con engine nativo no debería pasar, pero
  el check cuesta segundos y una corrida con código viejo cuesta días.)

### Otros gotchas del entorno (evitan errores ya cometidos)

- **Docker only** — no hay Python local. `./src` y `./scripts` montados (sin rebuild
  al editar); rebuild solo si cambia `pyproject.toml`.
- **`docs/` NO está montado** en el contenedor: scripts que generan material para
  `docs/` escriben a `data/output/` y se mueve en el host (ojo con ownership root).
- **matplotlib NO está en la imagen** — instalar ad-hoc por corrida
  (`pip install --quiet matplotlib && python -m scripts.bench.make_paper_figures`) y
  todo import de matplotlib en `scripts/` debe ser **lazy** (dentro de funciones) o
  rompe `tests/test_script_entrypoints.py`.
- **Protección de `main` estricta**: checks requeridos `Tests` y `Docker Compose Tests`,
  `strict:true`, `enforce_admins:true` → serializar merges (update-branch → CI → merge).
- **TAP de Gaia** tira HTTP 500 transitorios; el motor tiene retry. Reintentar.
- **Tests:** `docker compose run --rm test` (suite completa) antes de proponer commit.
  Formato: `docker compose run --rm pipeline ruff check . --fix` y `black .`.

---

## Tareas

### Tarea 1 — Regenerar Stage B híbrido sobre b1fix y refrescar todos los números (B3)

**Problema (informe B3):** el manuscrito mezcla números de dos freezes. Prueba
aritmética: §2.2 dice "8 728 509 pairs (12.08 %)" pero 8 728 509 / 80 072 774 = 10.90 %
— el numerador es del catálogo viejo (72.2 M). Dos TODOs en el tex lo admiten
(`aa_encounters.tex:259-261` y `:539-540`).

**Entregable:** Stage B (subset frágil re-refinado N-body) + ensamble híbrido corridos
sobre `encounters_catalog_rebound_005au_b1fix.parquet`, con output a nombre nuevo
(p. ej. `encounters_catalog_hybrid_stageb_b1fix.parquet`) para **no pisar** el híbrido
viejo; sidecar completo. Cadena (verificar `--help` de cada script antes — los args
exactos no están documentados acá):

```bash
# 1. Selección del subset frágil sobre b1fix
docker compose run --rm pipeline python -m scripts.validate.select_stageb_nbody_subset <args>
# 2. Re-refinado N-body por shards — LA corrida pesada (8-9 M pares; días).
#    Detached + cap 18g + pocos workers (ver sección Docker). Shards a un dir NUEVO
#    (no mezclar con stageb_nbody_shards/ viejos).
docker compose run --rm -d pipeline python -m scripts.validate.refine_stageb_nbody <args>
# 3. Ensamble del híbrido
docker compose run --rm pipeline python -m scripts.validate.assemble_stageb_hybrid_catalog <args>
```

**Números a refrescar en el tex** (lista del informe, verificar uno por uno):
- §2.2: conteo y % del subset frágil (el % con denominador 80 072 774).
- §3.1: estadísticos Stage B (hoy p99 1.99 mAU, max 15.2 mAU).
- §3.2: upward crossings (hoy 25 283 / 0.29 %).
- §3.4: cross-match FM completo (hoy 40 176 / 11 842 / 91.6 % / 8.4 %; el TODO habla de
  25 962 no recuperados vs 28 334 implicados por el texto — incompatibles). Re-correr
  `scripts/validate/crosscheck_literature_encounters.py` (o
  `validate_fuentes_munoz_2025.py`) contra el híbrido nuevo.
- Fig. 4 si cambia (`scripts/bench/make_paper_figures.py`).
- Abstract y Conclusiones: 0.70 %, 12 µAU, 2.5 mAU deben ser los del freeze final.
- Eliminar los dos TODOs.

**Gate (informe):** `grep -c "TODO" docs/paper/aa_encounters.tex` = solo los
administrativos de B4 (idealmente 0); y todo porcentaje del paper con
numerador/denominador tabulados reproduce al redondeo declarado — en particular
frágil/total = el % impreso.

**Depende de:** nada. **Bloquea:** 3, 8, 9, 10.

### Tarea 2 — Join de datos físicos/elementos: Tabla 6, class, Fig. 3 (B1 + M5)

**Problema (informe B1):** la Tabla 6 (`tab:candperturbers`, tex:771-795) tiene
diámetros incompatibles con valores medidos y densidades no físicas con sus propias
masas FM25: Nysa D=140 (medido ~71→ρ 0.40 vs 3.1), Fortuna D=133 (medido ~210→ρ 6.9
vs 1.7), Lutetia D=120 (Rosetta 98), Angelina D=104 (medido ~50-60). Firma = D
derivado de H con albedo default 0.14 pese a existir diámetros medidos en SBDB.
7/10 filas con clase "—" (¡Massalia, Lutetia, Amphitrite!); Metis D=197 en Tabla 6 vs
190 en Tabla 5. **(M5):** Fig. 3 solo plotea 92 971/449 213 cuerpos (20.7 %) por usar
una "orbital-elements source" que no es MPCORB (candidato: `gaia_orbits.parquet`),
con caption y texto en contradicción directa.

**Raíz probable común** (hipótesis del informe, confirmar contra código): la Tabla 6 y
la columna class se generaron por un camino de datos distinto al de Tablas 1–2 (que sí
tienen diámetros correctos). Candidatos: el join en `src/characterize/physical.py`
contra `data/raw/sbdb_physical.parquet`, y el generador de la Tabla 6
(candidato: `scripts/mass/find_extra_perturber_candidates.py`); Fig. 3 en
`scripts/bench/make_paper_figures.py` (usa `gaia_orbits` — cambiarlo al snapshot
MPCORB del freeze, cobertura 100 %).

**Entregable:** (1) join reparado; (2) Tabla 6 regenerada con diámetros medidos
(SBDB) y clase poblada, ranking re-derivado (Fortuna con D≈210 probablemente sube) y
verificación de si el corte D≳100 sigue dando 10 candidatos; (3) D de Metis
reconciliado entre Tablas 5 y 6 (una sola fuente); (4) auditoría en el catálogo
caracterizado: fracción de cuerpos con diámetro medido en SBDB que quedó con
`diameter_source = default_albedo` — reportarla; (5) Fig. 3 regenerada desde MPCORB,
caption sin la disculpa del subsample, texto y caption ya no contradictorios;
(6) documentar qué era la fuente fallida.

**Gate (informe):** script que para cada fila de la Tabla 6 calcule ρ implicada por
(D tabulado, M_FM25) y verifique 0.8 < ρ < 4.5 g/cm³; assert de que todo cuerpo de las
tablas del paper con diámetro en SBDB tiene `diameter_source ∈ {measured,
albedo_measured}`; cero clases "—" entre numerados < 1000; caption de Fig. 3 reporta
cobertura completa (`N/N`).

**Depende de:** nada (los diámetros no dependen del Stage B). **Bloquea:** 9.

### Tarea 3 — Clamps de borde en Tablas 1–2 (B2)

**Problema:** ventana declarada "2014 July 25 – 2017 May 28" (tex:174), pero
(196) Philomela × (431) Nephele y (85) Io × (145) Adeona (Tabla 1, tex:602 y :607)
están fechados **2017-05-29**. Dos filas con la misma fecha en el borde = firma de
mínimo clampeado en el último nodo de la grilla. Las Tablas 3–5 tienen guard de
">130 d from either edge"; las Tablas 1–2 no tienen ninguno.

**Entregable:** (1) verificar en el catálogo (híbrido nuevo de la Tarea 1) si
`jd_tdb` de esas dos filas coincide con el último nodo de la grilla (código candidato:
`src/detect/refine.py`, `src/detect/pipeline.py`); (2) si son clamps: excluirlas o
marcarlas `edge_censored`, regenerar Tablas 1–2 (la Tabla 1 pasa de 12 a 10 filas →
corregir texto y abstract "twelve encounters"), y aplicar a todas las tablas de
eventos el mismo criterio de distancia al borde (o declarar cuáles no lo aplican y por
qué); (3) si son genuinas (p. ej. formateo TDB→UTC de un mínimo interior del 28):
documentar la convención de cierre de ventana (hasta qué instante, en qué escala).

**Gate (informe):** query sobre el catálogo publicado: cero filas de cualquier tabla
del paper con `jd_tdb` fuera de [JD_ini, JD_fin] de la provenance; assert de que
ninguna fila tabulada tiene su mínimo en el primer/último nodo de la grilla; flag
`edge_censored` por fila si se retienen.

**Depende de:** 1 (auditar sobre el catálogo final).

### Tarea 4 — Describir la ventana de refinamiento realmente ejecutada (M1)

**Problema:** §2.2 (tex:223-225) describe "±2 h window sampled at 120 s" con cadencia
gruesa de 12 h — matemáticamente el método roto pre-B1 (pierde el mínimo en 2/3 de los
casos). §4.1/§4.2 citan "the refinement-window fix (Sect. 2)" pero la Sect. 2 nunca lo
describe. **El código real ya está arreglado desde la ronda 1**: `src/detect/refine.py`
usa `window_hours: 6.0` (≥ Δt/2) + re-centrado iterativo (`_MAX_RECENTER = 4`) —
verificar contra el código actual y describir *eso*.

**Entregable:** párrafo de §2.2 reescrito con el algoritmo ejecutado en el freeze
b1fix (ventana ±6 h = ±Δt/2 + re-centrado iterativo si el argmin cae en borde), una
frase sobre el modo de fallo pre-fix (§4.1 lo usa como evidencia), y declaración del
guard de borde de la ventana N-body de ±12 h del híbrido (verificar qué hace
`refine_stageb_nbody.py` cuando el mínimo N-body cae en el borde: ¿re-centra? ¿flag?).

**Gate (informe):** el injection-recovery ya citado (tex:516-519) ejecutado con
mínimos uniformes en fase dentro del paso de 12 h muestra ≥99 % de recuperación con la
ventana **tal como quedó descrita en el texto** (mismos parámetros —
`scripts/validate/injection_recovery_detection.py` ya existe); cero mínimos del
catálogo final reportados en el borde de la ventana de refinamiento.

**Depende de:** nada (es texto + verificación de código existente).

### Tarea 5 — MC de σ con covarianza de época consistente (M2)

**Problema:** las distancias del paper se calcularon con MPCORB 2016-02-17, pero el MC
(`scripts/validate/mc_distance_uncertainty.py`, ronda 1 B4) muestrea covarianzas SBDB
de **2026** — soluciones con ~10 años más de astrometría. Peor donde más importa: los
cuerpos de la Tabla 3 son de número alto (arcos cortos en 2016); su σ de época 2016
puede ser 1–2 órdenes mayor y demoler el claim "significant to ~2–8 %". Además mezcla
nominal de una solución con covarianza de otra sin acotar el offset.

**Entregable:** elegir y declarar una de las dos vías del informe:
(a) evaluar significancia con la solución actual — recomputando también la distancia
nominal de los eventos de Tabla 3 con esa solución, y reportando el corrimiento
nominal MPCORB-2016→SBDB-2026 por evento como término empírico de error; o
(b) mantener el nominal 2016 y estimar covarianza de época 2016 (refit con astrometría
truncada, o escala por longitud de arco). Reescribir el pasaje "1–83 km… not the
10³–10⁴ km one would fear" según resulte.

**Gate (informe):** Tabla 3 regenerada con dos columnas — Δ nominal y σ_d consistente
en época — y el claim de significancia recalculado; si σ_d/d > 30 % para alguna fila,
el texto lo dice explícitamente en vez de "2–8 %".

**Depende de:** 1 (las filas de la Tabla 3 deben ser las del freeze final).

### Tarea 6 — Diseño del censo de umbral: ¿estratificado o aleatorio? (M4)

**Problema:** §3.2 (tex:429-434, 449-457) dice que la muestra de 10⁴ cuerpos fue
"stratified… rather than a density-weighted random draw" pero extrapola con
(449 454/10 000)² — válido solo para draw aleatorio. Dato duro del informe: el conteo
observado (17 469 pares en [0.05, 0.06)) coincide al 0.4 % con la predicción para draw
aleatorio. Una de las dos afirmaciones está mal escrita.

**Entregable:** determinar del registro de la corrida qué diseño se usó realmente
(candidatos: `scripts/validate/measure_threshold_false_negatives.py` y su sidecar/log;
buscar la semilla y el criterio de muestreo). Si fue aleatorio: quitar "stratified…"
y mantener la extrapolación. Si fue estratificado: rehacer la extrapolación con pesos
o re-correr con draw aleatorio (10⁴ cuerpos = barato) y recalcular el "~10⁵ censored"
de abstract y Conclusiones.

**Gate (informe):** el párrafo declara el diseño con una línea de procedencia (semilla,
criterio), y hay un check numérico publicado en el repo: pares en banda del sample vs
predicción N² del catálogo — si difieren >2σ, el texto usa el método corregido.

**Depende de:** nada (la re-corrida, si hace falta, es barata).

### Tarea 7 — Europa, determinabilidad, SG22, ratios FM (M3 + N9 + N14a)

**Problema (tres capas, tex:934, :841-858, :955-960):** (1) Europa es "measured" con
ratio 0.58 y z=−4.3 vs DE441 — incompatible con el "broadly consistent" del abstract y
jamás discutido. (2) El claim "true mass ≳2×10¹⁹ kg are recovered near unity ratio" lo
refutan Europa e Interamnia en la propia Tabla 8. (3) "geometric-mean ratio of 0.85"
no reproducible desde la tabla (da 0.75 con 0.77×0.95×0.58). Falta la determinación
independiente Gaia-based de Europa: **Siltala & Granvik 2022, A&A 658, A65** (existencia
verificada por el tribunal vía web).

**Entregable:** (1) discutir Europa explícitamente — o el fit subestima por regresión
masa↔órbita no capturada por el flag (→ reportar el falso positivo del flag y evaluar
endurecerlo, p. ej. exigir |z|<3 contra seed para "measured", o correr bootstrap en
Europa aunque lev=0.43), o la referencia DE441 es cuestionable (citar evidencia).
(2) Acotar la frase de determinabilidad a-priori a los casos que la cumplen, nombrando
excepciones. (3) Matizar el abstract. (4) Añadir cita SG22 al `.bib` (verificar contra
ADS) y comparar Europa/Eunomia/Amphitrite contra ella. (5) Tabular ratios vs FM25 o
corregir el 0.85 (N9).

**Gate (informe):** Tabla 8 con columna o párrafo de z vs FM25 para los seis
"measured"; regla de status reproducible por script (criterio publicado → mismos
flags); 0.85/0.75 reconciliado; SG22 citada y discutida.

**Depende de:** nada (usa los fits existentes; solo el bootstrap opcional computa).

### Tarea 8 — Cuarto término del presupuesto de incompletitud (M6)

**Problema:** la membresía es Kepler (<0.05 AU) pero la candidacy es N-body dentro de
r_q=0.0572 AU; margen 7.2 mAU < |Δd| máximo medido Kepler↔N-body (15.2 mAU en Stage B).
Existen miembros por definición que jamás fueron candidatos; el término no está en el
presupuesto.

**Entregable:** con la distribución empírica de Δd del **nuevo** Stage B (Tarea 1),
calcular P(Δd > r_q − d_samp) integrada sobre pares con d_Kep cerca del umbral, y
publicar la cota (esperable "≲10³ pares, subdominante frente a ~10⁵ del censoring").
Una frase en §3.2 o subsección corta.

**Gate (informe):** el presupuesto de §3 enumera **cuatro** términos pipeline-internos,
cada uno con número y método; el término nuevo tiene script reproducible como los otros
tres.

**Depende de:** 1 (usa la Δd del Stage B regenerado).

### Tarea 9 — Pasada editorial N1–N16

Una sola pasada sobre el tex/bib con los números ya finales. Checklist (ubicaciones y
soluciones exactas en el informe; acá el resumen accionable):

| N | fix | gate |
|---|-----|------|
| N1 | "~100×" → "~10× at the median (~40× at p99)" (tex:410-414) | ratio reproducible desde los percentiles de la misma frase |
| N2 | abstract: "recovers the four calibrators at \|z\|<3 (Pallas, with only eight encounters, as a bound)" | mismo vocabulario de status que Tabla 8 |
| N3 | fórmula de determinabilidad: añadir factor t (Δθ ~ 2GM·t/(b·v_rel·Δ)) y qué t se usó (tex:842-844) | análisis dimensional cierra en radianes |
| N4 | "~150,000 bodies" → 449 454 y ~10¹¹ pares (tex:193-194) | un solo N en todo el paper |
| N5 | unificar nombres .parquet entre §6 y Apéndice A (conexo B4) | grep de nombres consistente con el depósito |
| N6 | convención "subscript 1 = larger body" añadirla en §2.3 o quitar la cross-ref (tex:1144-1145) | la ref apunta a texto que existe |
| N7 | eliminar "(Sect.~below)" ×2 (tex:854, :858) | cero ocurrencias |
| N8 | abstract: "0.70 % downward-crossing rate in the discarded band, implying ~10⁵ censored (~0.25 % of the catalogue)" | ambas cantidades con denominador explícito |
| N9 | (con Tarea 7) tabular ratios vs FM25 o corregir 0.85 | número reproducible desde lo impreso |
| N10 | CIs binomiales: cluster bootstrap por cuerpo o frase reconociendo clustering (tex:433-434, :486) | CI recalculado o clustering declarado |
| N11 | declarar aproximaciones del flag gaia_observable: elongación geocéntrica (offset L2 ≤0.6°) y V como proxy de G (tex:270-275, :1178-1179) | caption/schema lo dicen |
| N12 | Conclusiones: quitar "selects individual events by hand"; alinear con la Intro; unificar candidate/real (tex:1085-1086 vs :116-118) | frases consistentes entre sí y con el abstract |
| N13 | `li2023`: añadir "and {Chen}, J." (references.bib:119-127) | entrada coincide con ADS |
| N14 | citas: SG22 (con Tarea 7); paper PSJ 2025 LSST como future work (verificar contenido antes de citar); fuentes primarias de diámetros/densidades (Tedesco, Usui, Mainzer; Carry 2012) | cada fuente física del pipeline con cita primaria verificada contra ADS |
| N15 | especificar: dataset de calibración de s_c; definición de "near-boundary" del piso 0.42 %; anti-double-counting del perturber en el force model para los 16 | las tres especificaciones presentes |
| N16 | cota de una línea para Yarkovsky/no-grav en el arco FPR; frase sobre manejo de mínimos en bordes de bloques paralelos (verificar `src/detect/parallel.py`); nota de caption para los z de Fortuna/Hygiea (¿incluyen σ de literatura?) | cotas impresas; z reproducibles desde el caption |

**Depende de:** 1, 2, 3, 5, 7 (números finales primero; si no, se edita dos veces).

### Tarea 10 — Submission (B4)

**Entregable:** depositar los tres parquet + sidecars en Zenodo/VizieR (scaffolding
listo: `docs/paper/DATA_AVAILABILITY.md`, `docs/paper/zenodo_data_deposit.json` —
faltan `creators` y DOI), insertar DOI en §6, completar afiliación/ORCID/funding
(tex:30-38, :1076-1077, :1111).

**⚠️ Parcialmente author-owned:** afiliación, ORCID, funding y la cuenta Zenodo/CDS
son del usuario — el agente prepara todo (SHA-256, metadata, texto) y **pide al
usuario** los datos personales y la ejecución del depósito. No inventar afiliación.

**Gate (informe):** el DOI resuelve públicamente y los SHA-256 listados en el registro
coinciden con los archivos citados en el texto.

**Depende de:** 1–9 (se deposita el freeze final, no uno intermedio).

---

## Fuera de scope

- **Re-auditar la ronda 1** — cerrada (`planning/TRIBUNAL_REMEDIATION_PLAN.md`).
- **N-body de todo el catálogo (~400 CPU-días)** — la Tarea 1 regenera el *subset
  frágil* (~8-9 M pares), no el catálogo entero. Sigue diferido.
- **F8 (Pallas con DR4)** — espera datos; `planning/MASS_FUTURE_WORK.md`.
- **Rediseño bayesiano del motor de masas** — M3 se resuelve con discusión +
  endurecimiento del flag, no reescribiendo el motor.

---

## Ruta crítica

```
Tarea 1 (Stage B b1fix, corrida pesada) ──→ 3, 5, 8 ──┐
Tarea 2 (join físico: Tabla 6, class, Fig. 3) ────────┤
Tarea 4 (texto refinamiento) ─────────────────────────┼──→ 9 (editorial) ──→ 10 (submission) → A&A
Tarea 6 (diseño del censo) ───────────────────────────┤
Tarea 7 (Europa/SG22) ────────────────────────────────┘
```

Lanzar la corrida de la Tarea 1 **primero** (días de wall-clock) y trabajar 2, 4, 6 y 7
en paralelo mientras corre. Nada de editar números del paper hasta que 1 cierre.

---

## Tanda 1 (2026-07-09)

**Pre-vuelo verificado:** engine Docker nativo (`context default`), cap 18g; md5 de los 4
scripts Stage B coincide host↔contenedor; host 28 CPU / 31 GB.

**Tarea 1 (B3) — LANZADA:**
- `select_stageb_nbody_subset --catalog ..._005au_b1fix.parquet` → subset frágil
  `data/cache/nbody_validation/stageb_selective_subset_b1fix.parquet` = **10 339 493 filas
  (12.91 % de 80 072 774)**; el conteo viejo (8 728 509 / 12.08 %) era del freeze pre-B1.
- Smoke (100 y 5 000 pares, arranque y medio del archivo): **0.005 s/row** con 24 workers,
  0 fallos, RAM 1.9 GiB/18 → estable. Total estimado ≈ **15 h**.
- Corrida real detached: `refine_stageb_nbody --input <subset b1fix> --output-dir
  data/output/stageb_nbody_shards_b1fix --shard-index 0 --shard-size 10000 --num-shards
  1034 --workers 24 --no-progress`. Container `1374d1d5feab…` (id en
  `logs/stageb_regen_container_id.txt`); log en `logs/stageb_regen_b1fix.log`; modo shards
  = reanudable saltando shards escritos. Waiter en background avisa al terminar.
- Pendiente al terminar: `assemble_stageb_hybrid_catalog` → `encounters_catalog_hybrid_stageb_b1fix.parquet`,
  y refrescar §3.1/§3.2/§3.4/Fig.4/abstract + borrar TODOs.

**Tarea 4 (M1) — texto cerrado, gate PASS:** §2.2 reescrito con el algoritmo real
(ventana ±6 h = ±Δt/2, coarse 12 h; re-centrado iterativo `_MAX_RECENTER=4`; vértice
cuadrático; dedup del mínimo por-par entre bloques → cierra N16b; modo de fallo pre-fix
±2 h ≈ 2/3 clipeados = el "refinement-window fix" que citan §4.1/§4.2). Guard N-body ±12 h
declarado (flag `near_boundary`, no re-centra). Gate: `injection_recovery_detection
--n-pairs 300` → **300/300 = 100 %**, |Δd|max 1.6e-10 AU, |Δt|max 0.1 s, PASS. Residual:
verificar "cero mínimos en el borde" sobre el catálogo final (junto a Tarea 3).

**Tarea 6 (M4) — CERRADA:** la corrida real del censo fue muestreo **aleatorio** uniforme
sembrado (código + `data/output/kepler_false_negatives/summary.json`: n_bodies=10000,
17 469 pares en banda, 0.6984 %), no estratificado como decía el texto. §3.2 corregido +
validación N² añadida. Script del gate: `scripts/validate/check_threshold_census_scaling.py`
→ **PASS 0.46 %** (39 456 × (449454/10000)² = 7.97e7 vs 8.01e7).

**Tarea 7 (M3+N9+N14a) — CERRADA:** hallazgo central: la z=−4.3 de Europa en la Tabla 8
es contra el **seed DE441** (4.0e19 kg), no contra masas medidas. Datos exactos
(`data/output/orbdet/mass_catalog_b1fix_targets.csv` + `.../fuentes_munoz_jack_mass_comparison.csv`):
Europa fit 2.33e19; FM25 2.66e19 (ratio 0.88, z=−0.8); SG22 2022 3.0e19 (ratio 0.78,
z=−1.7); seed 4.0e19 (ratio 0.58, z=−4.3). SG22 usa Gaia DR2 + terrestre → independiente
del FPR (que FM25 y este trabajo comparten). N9: "0.85" = geomean vs FM25 (0.75·0.95·0.88);
"0.75" = geomean vs seed (Tabla 8). Editado: §5 σ-externa (Europa + FM25/SG22 + z vs FM25
de los seis measured |z|<1.2), párrafo determinabilidad (acotado a determinaciones
independientes; N7 "(Sect.~below)"×2 eliminado; N3 fórmula con factor t), abstract (N2
Pallas bound + SG22). `siltala2022` en references.bib. Gate: `scripts/mass/verify_mass_status_rule.py`
→ **PASS 16/16** (criterios min_snr_jack=3, min_n_jack=10, max_leverage=0.5 → los flags de
Tabla 8). Bootstrap de Europa NO corrido: innecesario (lev=0.43; el flag mide precisión, y
Europa está bien constreñido — el problema es el seed, no el fit). Paper compila 13 pp.
**Adelanto de N-menores:** N3, N7 cerrados junto a M3; N9, N14a cerrados. Pendientes en
Tarea 9: N1, N4, N5, N6, N8, N10, N11, N12, N13, N15, N16.

**Tarea 2 (B1+M5) — código listo, corrida pesada pendiente:** `characterize_catalog.py::_supplement_elements`
ampliado (rellena a/e/i desde MPCORB para TODOS los cuerpos ausentes de gaia_orbits, no solo
[1,2,4,10]) → clase 100 % al re-caracterizar. Fig. 3 (`make_paper_figures.py`) apuntada al
snapshot MPCORB (`MPCORB_20160217.DAT`, 100 % cobertura); texto+caption del tex ya coherentes
(sin disculpa de subsample). Falta: re-caracterizar (tras Tarea 1) + regenerar Tabla 6 y
Fig. 3 + script de gate ρ + auditoría diameter_source.

**Nota de recursos:** un solo contenedor pesado a la vez (cada `docker compose run` recibe
su propio cap de 18 g; dos en simultáneo pueden sumar 36 g > 31 g del host). La
re-caracterización de la Tarea 2 (pesada) se corre cuando termine la Tarea 1, no en
paralelo.

## Tanda 2 (2026-07-10)

**Estado de la Tarea 1:** el host se **suspendió** durante el día (un gap de ~13.6 h en
`logs/stageb_regen_b1fix.log`), así que ~21 h de reloj = solo ~2 h de cómputo → 147/1034
shards al momento de esta tanda. Corrida sana (50 s/shard, 0 fallos, reanudable). Quedan
~12 h de **cómputo**; el reloj depende de cuánto siga despierta la máquina. Sin riesgo de
datos. (Si se quiere acelerar: mantener la máquina despierta.)

**Adelanto de N-menores independientes del freeze** (todos verificados contra código/datos,
paper compila 13 pp):
- **N8** (abstract + Conclusiones): 0.70 % aclarado como tasa dentro de la banda
  [0.05,0.06); censado ~0.25 % del catálogo con denominador explícito.
- **N10**: caveat de clustering en los CIs binomiales (censo + prefiltro) — declarados
  indicativos, no exactos.
- **N11**: aproximaciones del flag `gaia_observable` (elongación geocéntrica, offset L2
  ≤0.6°; V como proxy de G, ≲0.2 mag).
- **N12** (Conclusiones): quitado "selects individual events by hand" → alineado con la Intro;
  "candidate pairs" → "encounter pairs" (unificado con "real 3D close encounters").
- **N15** (verificado en `src/orbdet/`): (a) s_c se calibra **por-fit** por bisección sobre
  los residuos convergidos de cada perturbador (sin piso si χ²_red≤1 → explica los <1);
  (b) near-boundary = mínimo N-body en el borde de la ventana ±12 h; (c) fuerza `ASTEROIDS`
  de ASSIST **siempre** excluida, perturbadores como partículas masivas explícitas → sin
  double-counting para los 16.
- **N16a**: cota de Yarkovsky (~10⁻⁴ AU/Myr, ≲1 km ≲1 mas sobre el arco de ~6 yr, absorbido
  por el ajuste de órbita).
- **N16c**: caption de Tabla 8 — z = (M_fit−M_ref)/√(σ_tot²+σ_ref²) explícito.
- (N3, N4, N6, N7, N9, N13, N14a ya en Tanda 1.)

**Aceleración (2026-07-10, con caffeine activo):** 2º contenedor `refine_stageb_nbody`
lanzado sobre el rango alto (shards 900–1033, 4 workers, mismo `--output-dir`, reanudable),
usando los cores libres. Container `91551e53` (id en `logs/stageb_regen2_container_id.txt`,
log `logs/stageb_regen2_b1fix.log`). CPU total ~26/28 cores; RAM combinada ~1 GiB (host 31 GB,
sin riesgo). ETA de cómputo ~9.5–10.5 h. El job es CPU-bound: la RAM baja NO es headroom;
28 workers sobre 28 cores es el techo (más = oversubscription). Ambos waiters activos.

**Delegación a subagentes (working style nuevo, ver memoria):** Haiku validó la suite
(586 passed, ruff limpio, 3 gates PASS). Sonnet escribió `check_table6_density.py` (gate ρ,
Tarea 2) y `measure_stageb_dd_distribution.py` (Δd para M6, Tarea 8) — ambos ruff-limpios y
reusables. N14b resuelto por Sonnet (bernstein2025 citado).

**Pendiente al reanudar** (todo espera el híbrido de la Tarea 1): refresco §3.1/§3.2/§3.4/
Fig.4 + percentiles refine-error del abstract (N1); Tarea 2 re-caracterización + Tabla 6/Fig.3
+ correr gate ρ sobre la tabla corregida; Tareas 3/5; Tarea 8 sobre shards completos; N5 en la
pasada editorial final; Tarea 10 depósito.

## Tanda 3 (2026-07-11) — cierre de B1/M5 + N5, paper commit-ready

**Tarea 2 (B1+M5) CERRADA.** La re-caracterización (recls, escrita por corrida previa) se
**verificó** y **promovió** al canónico `encounters_characterized_b1fix.parquet` (9.18 GB;
pre-fix respaldado `_PREFIXBUG.parquet.bak`, reversible). Verificación (`data/cache/_verify_recls.py`):
80 072 774 filas; **class mala (null/vacío/"—"/"Other") entre numerados<1000 = 0** (de 291 557
en cuerpo 1 y 494 en cuerpo 2); `diameter_source_1` measured 27.9M / zone_albedo 52.2M. Fix B1
confirmado.

**Tabla 6 (`tab:candperturbers`) regenerada** contra el catálogo promovido (subagente Sonnet,
foreground-síncrono; el agente previo `a8d2bf5c` había fallado colgado en un Monitor sin correr
los pasos 2–6). Cambia la composición porque los E-type con D inflado por H·albedo caen del corte
D≥100 km al usar D **medidos**: Nysa 140→71, Angelina 104→58, Lutetia 120→98, Euterpe 141→96,
Urania 109→93, Amphitrite 240→190 (cae al 13.º). Entran Hestia(43), Pompeja(38), Hersilia(37),
Polyxo(29), Aegina(28), Beatrix(28). Nuevo líder **(46) Hestia** (Metis pasa a 38 útiles). Los 10
nuevos **todos tienen masa FM25** → "all ten already carry a published mass" ahora verdadero y
consistente (Angelina, que tenía "--", salió). Prosa §candperturbers reescrita (Hestia líder,
Metis el más lento, Pompeja el más cercano, Fortuna el mayor D≈200); §5 "top ranked"→"highly
ranked", Metis útiles 37→38. **Metis D reconciliado**: 190 en Tabla 5 y 6 (la vieja tenía 197 por
leer el catálogo `_full` equivocado).

**Gate ρ** (`check_table6_density.py`, subagente): **15/16 PASS, 0 FAIL, 1 SKIP** (Angelina sin
FM25). Nysa 0.40→**3.12**, Fortuna 6.90→**2.04** g/cm³, todos en 0.8–4.5. **Fig.3**: cobertura
**449 213/449 213 = 100 %** (fuente MPCORB snapshot). Figuras 1–4 re-renderizadas del freeze y
copiadas a `docs/figures/`. `mine_notable_encounters.py` `DEFAULT_CATALOG`→`_b1fix` (reproducibilidad).

**N5 CERRADO** (Tarea 9 completa): Apéndice A usaba nombres viejos (`encounters_characterized_full`,
`..._rebound_005au` sin `_b1fix`, `..._hybrid_stageb` sin `_b1fix`) — unificados al freeze
`_b1fix` en §6 y Apéndice A; el híbrido también nombrado en §6. Los tres productos existen en disco.

**Compila 13 pp, exit 0, sin refs/citas sin resolver.** 2 TODOs vivos, ambos author-owned (`:30`
ORCID, `:1184` DOI Zenodo). Cero cuerpos eliminados colgando en prosa (grep limpio).

**Housekeeping del árbol:** `.gitignore` +`*.fdb_latexmk`/`*.fls`; `config.local.yaml.bak_tribunal`
des-stageado (backup manual, untracked); `scripts/dev/_tribunal_b1fix_report.py` es scratch
("Safe to delete") → NO commitear. Docs del tribunal movidos a `docs/tribunal/`.

**Pendiente:** aprobación de commit del usuario (nada commiteado aún); Tarea 10 (ORCID + cuenta
Zenodo + DOI/SHA-256 sobre freeze final); opcional B2 columna `edge_censored` en el parquet
depositado (con Tarea 10).
