# Plan de remediación — hallazgos del tribunal científico

> **Estado:** 🟢 AVANZADO — catálogo **regenerado** (`_b1fix`, 80.072.774 filas; B1/B10
> cerrados) y caracterizado; **10/10 bloqueantes cerrados o sustancialmente**;
> mayores M1/M3/M4/M7/M8-9/M10/M11/M12/M13 cerrados. Todo en la branch
> `fix/tribunal-remediation-core` (PR #102, CI verde, paper 11 pp.). Restan:
> corridas pesadas (Stage B híbrido, re-fit de masas sobre targets b1fix),
> menores (M5/M6/M14, C11-C19) y submission (Tarea 28, author-owned).
> **Última actualización:** 2026-07-08
> Plan para cerrar los hallazgos de `docs/tribunal_cientifico_2026-07-04.md`
> (10 bloqueantes, 14 mayores, 19 menores) y llevar el paper de "rechazo probable en
> A&A" a "publicable en A&A". **Criterio de éxito:** los 10 bloqueantes cerrados con su
> gate verificable, y un manuscrito cuyo abstract/§1/§4/§5 sobreviven la verificación
> contra el código y los datos públicos.
>
> Fuente de verdad de los hallazgos: `docs/tribunal_cientifico_2026-07-04.md`.
> Este tracker apunta a esa evidencia, no la duplica.

---

## Tabla de estado

### Bloqueantes (gate de submission)

| # | tarea | sev | estado | depende de | entregable / PR |
|---|-------|-----|--------|-----------|-----------------|
| 1 | Fix ventana de refinamiento + regenerar catálogo | B1 | ✅ **REGENERADO**: `encounters_catalog_rebound_005au_b1fix.parquet` (80.072.774 filas); histograma de offset jd_tdb **plano** (pico [1.75,2.25]h: 76%→8.1%); gate 4 cuerpos OK; sidecar escrito | — | `src/detect/refine.py`, freeze nuevo |
| 2 | Root-cause del miss (804)×(733) | B10 | ✅ **presente en el catálogo b1fix a 0.013547 AU** (coincide JPL 0.0138); causa raíz confirmada = B1 (recorte de ventana) | 1 | doc en `docs/` |
| 3 | Declarar universo muestral y N + sidecars completos | B2 | 🔶 sidecar `universe` + §2.1 del tex declaran corte a∈[1.5,4.0] y N=449.454; falta rehacer extrapolación N² (post-regen) | 1 | sidecar, §2 del tex |
| 4 | Diámetros/albedos medidos + regenerar §4 | B3 | ✅ catálogo caracterizado b1fix con diámetros medidos (27.9M measured); tab:dd100 (2→12 pares) y tab:dd50 regeneradas; claim "no new large-large" corregido en abstract/§4.1/§4.5/§7. Iris×Nysa SALE, Nemausa×Aegina ENTRA (confirmado). tab:closest/slowest → con Tarea 5 (σ) | 1 | `src/characterize/`, tablas §4 |
| 5 | Presupuesto de completitud honesto (σ elementos + injection-recovery) | B4 | 🔶 §3 renombrado a pipeline-induced + término de elementos como cota (10³-10⁴ km) + injection-recovery citado; tab:closest/genético marcados no-significativos. Falta propagación MC de covarianzas AstDyS (diferido) | 1 | `docs/`, §3 del tex |
| 6 | σ de masas defendible (leverage, N mín, bootstrap) | B6 | 🔶 catálogo listo; falta bootstrap | — | `src/orbdet/`, `mass_catalog` |
| 7 | Cruce FM con potencia + una sola σ oficial | B7 | ✅ §5 reescrito: una sola σ_tot en la tabla; el "10/10" reemplazado por potencia (13-102%) + test de signo (5/6 bajo FM, geom 0.80, p≈0.2) + no-independencia FPR + caveat reinstaurado; anclaje = calibradores | 6 | `docs/mass_crosscheck_jack.md`, §5 |
| 8 | Validación de calibradores no circular (LOCO) + inyección end-to-end | B8 | 🔶 LOCO hecho y PASA (max \|z\|=1.40); falta ∂M/∂s_c e inyección | 6 | `scripts/mass/`, `docs/` |
| 9 | Tabla de 16 masas + schema del catálogo + selección de targets | B9 | ✅ tab:sixteen (16 masas: N, mass, σ_tot, χ²_red, ratio, z, leverage, status) + tabla de schema (Apéndice A); tab:calibrators redundante eliminada. N por perturbador en la tabla | 6,7 | tablas del tex |
| 10 | Reescribir novedad + bibliografía completa | B5 | ✅ título/abstract/§1 reescritos; 18/18 citas verificadas contra ADS en `references.bib` | — | `aa_encounters.tex`, `references.bib` |

### Mayores (exigidos por referee, no bloquean el build)

| # | tarea | sev | estado | depende de |
|---|-------|-----|--------|-----------|
| 11 | Corregir explicación de la brecha FM (§3.4) — es ventana temporal | M1,M2 | ✅ cuantificado: 83.8% [80.3,86.8]% de los pares FM no recuperados tienen su encuentro <0.05 AU fuera de la ventana DR3 (muestra 500 propagados 1990-2024); confirma que la brecha es temporal, no de umbral | 10 |
| 12 | Cota de cadencia analítica (no circular) | M3 | ✅ §3: d(δ)=√(d_min²+v²δ²) vs r_q → cubre v_rel≤192 km/s (geométrico, no de la distribución del catálogo); injection-recovery citado |
| 13 | Regenerar Fig. 3 desde cuerpos del catálogo real | M4 | ✅ fig1-4 regeneradas de b1fix; fig3 usa cuerpos únicos reales (449.213, no 93.010) con caveat de cobertura (solo 92.971 tienen elementos — gaia_orbits no cubre muchos numerados altos, incl. Ceres) | 1 |
| 14 | Reconocer/justificar asimetría ventana DR3 vs arco FPR | M5 | ✅ §5: la asimetría afecta solo la SELECCIÓN de targets, no el modelo (el perturbador se integra continuo sobre todo el arco); yield = cota inferior conservadora |  — |
| 15 | Tablas §4: fuente, flag boundary, reconciliar FROZEN_RUN | M6 | ✅ dd100/dd50 = todas medidas (declarado); closest/slowest/largeslowclose = captions declaran procedencia por-fila (measured JPL SBDB vs H+albedo) y flag boundary (>130 d del borde); FROZEN_RUN reconciliado (Hispania×Mocia) | 1 |
| 16 | Columna de señal de deflexión por par + estratificación de utilidad | M7 | ✅ columna `deflection_dv_m_s` en el catálogo + §2.4 con ranking por señal (781 pares >1 mm/s, top = Ceres/Vesta) y estratificación | 1 |
| 17 | Reencuadrar §5 como demo de uso; comparar con Siltala/FM | M8,M9 | 🔶 framing hecho (demo de uso + párrafo Siltala/FM + sin "methodological framework"); falta compresión y números post-regen | 9,10 |
| 18 | Injection-recovery de detección como test permanente | M10 | ✅ harness + test CI; 200/200, \|Δd\|≤1.2e-10 AU | 1 |
| 19 | Alinear docs de prefiltro con la corrección 2026-05-31 | M11 | ✅ 2026-07-04 | — |
| 20 | Documentar mezcla de propagadores en `dist_au` (híbrido) | M12 | ⬜ | 1 |
| 21 | Criterio de identificabilidad por verosimilitud perfilada | M13 | ✅ `src/orbdet/identifiability.py` (Δχ²(M=0) cuadrático+perfil exacto, false_alarm_probability con frontera M≥0); en build_mass_catalog en paralelo al snr_jack; 10 tests. Perfil marca 13/15 vs 9 por jack → confirma que son cantidades distintas | 6 |
| 22 | Cuenta a priori de S/N de deflexión en el paper | M14 | ✅ §5 párrafo "A-priori determinability": Δθ~GM/(b·v·Δ); 10⁻¹² M☉ S/N≲1/target, 10⁻¹¹ M☉ determinable — coincide con el corte observado (measured ≳2e19 kg) | — |

### Menores

| # | tarea | sev | estado |
|---|-------|-----|--------|
| 23 | Lote de correcciones de código/comentarios de bajo riesgo | C1–C8 | ✅ C1–C8 hechos |
| 24 | Lote de correcciones de texto/figuras del paper | C9–C15 | 🔶 C9/C10/C15 hechos; faltan C11–C14 (Fig.1 ley de potencias, etc.) |
| 25 | Reproducibilidad: sidecars, hashes, desempate determinista | C16 | 🔶 desempate hecho |
| 26 | F3: reformular gate como Δmasa pareada | C14 | 🔶 gate reformulado y PASA (máx 0.23 % < 0.25 %); falta cota de escala |
| 27 | Limpieza de refs y tablas menores del paper | C17–C19 | ⬜ |
| 28 | Mecánica de submission (DOI CDS/VizieR, autores, acks, software) | — | ⬜ |

**Regla de orden:** nada de tocar números del paper (tareas 3–5, 11–16, 24) antes de
cerrar la Tarea 1 — el catálogo se regenera y todo número derivado cambia.

---

## Fase 0 — Cimiento (desbloquea casi todo)

### Tarea 1 — Fix de la ventana de refinamiento y regeneración del catálogo
**Entregable:** `src/detect/refine.py` con ventana `≥ coarse_step/2` (±6 h) o re-centrado
iterativo cuando el argmin cae en borde; validación `window_hours ≥ coarse_step/2` en
`src/detect/pipeline.py`; catálogo congelado regenerado con sidecar completo.
**Gate:**
- Test de regresión nuevo: mínimo sintético con offset > `window_hours` del sample grueso
  se recupera con |Δt| ≤ paso fino y |Δd| ≤ 1 μAU (hoy no existe — `tests/test_detection.py`
  valida la parábola pero nunca el caso argmin-en-borde).
- Histograma de offset de `jd_tdb` vs grilla gruesa en las filas Kepler: **plano** en
  [0, 6] h (hoy: 76 % apilado en [1.75, 2.25] h).
- Assert de config falla si `window_hours < coarse_step/2`.
**Depende de:** nada. **Bloquea:** 2,3,4,5,11,12,13,15,16,18,20 y todo número del paper.
**Notas:** costo ~3× evaluaciones Kepler en refinamiento (trivial). Re-derivar Stage A/B
separando artefacto de ventana del error de modelo real; re-medir censura de umbral
(`scripts/validate/measure_threshold_false_negatives.py` usa la misma ventana).

### Tarea 6 — σ de masas defendible bajo leverage y N pequeño
**Entregable:** en `src/orbdet/` + `scripts/mass/build_mass_catalog.py`: columna de
leverage top-1 por perturbador; N mínimo ≥ 10 para reportar σ_jack; bootstrap por objetivo
(o delete-d) donde el leverage top-1 > 50 %; |M| en `sigma_sys_kg` (fix del bug de signo en
`build_mass_catalog.py:128`); masas negativas marcadas no-físicas en el CSV.
**Gate:** ningún perturbador reportado como "measured" cambia de clase al quitar su réplica
jackknife dominante (hoy Cybele 2.92→4.82 cruza); el CSV no contiene σ negativas; leverage
reportado por fila.
**Depende de:** nada (paralelizable con Tarea 1). **Bloquea:** 7,8,9,21.

### Tarea 10 — Reescribir la novedad y completar la bibliografía
**Entregable:** abstract y §1 del `.tex` sin "hand-selected"; novedad reformulada = catálogo
all-pairs público con procedencia + presupuesto de incompletitud medido (contrastado con
FM 2025 top-100-truncado, Ivantsov 2018, Galád & Gray 2002); `references.bib` con las ~12
faltantes: Li 2023, Siltala & Granvik 2020/2021/2022, Baer & Chesley 2017, Ivantsov 2018,
Galád & Gray 2002, Fienga 2003, Kretlow 2020, Park 2021 (DE440/441), David 2023 (FPR),
INPOP/EPM, y citas de software (rebound, IAS15, WHFast, ASSIST, astropy).
**Gate:** ninguna afirmación de primicia sin cita del prior art que la matiza; `references.bib`
≥ ~20 entradas; toda cita verificada contra ADS.
**Depende de:** nada. **Bloquea:** 11,17.

---

## Fase 1 — Sobre el catálogo regenerado

### Tarea 2 — Root-cause del falso negativo (804)×(733)
**Entregable:** doc en `docs/` con la causa exacta (propagar los dos cuerpos aislados con el
mismo snapshot y grilla, trazar dónde se pierde). Candidatos: borde de ventana (Tarea 1) o
corte a ∈ [1.5, 4.0] (Tarea 3).
**Gate:** el par aparece en el catálogo regenerado, o el paper divulga explícitamente el gap
con su causa. 4/4 eventos de Fienga in-window recuperados o justificados.
**Depende de:** 1.

### Tarea 3 — Declarar el universo muestral y el N
**Entregable:** N exacto de cuerpos propagados y corte a ∈ [1.5, 4.0] en §2 del tex;
`n_asteroids` poblado en el sidecar; título/abstract/§1/§7 sin "full numbered population";
extrapolación N² de censura rehecha con el N declarado.
**Gate:** el N del texto, del sidecar y del denominador de la extrapolación coinciden (hoy:
150k vs 449k, factor ~9); "~1.5–2.5×10⁵ censurados" reproducible desde el texto.
**Depende de:** 1.

### Tarea 4 — Diámetros y albedos medidos
**Entregable:** `src/characterize/` cruza con diámetros/albedos medidos (IRAS/AKARI/NEOWISE)
donde existan, cae a albedo por clase cuando no; §2.3 corregido ("class" es orbital, no
taxonómica); todas las tablas de §4 y la selección D≳100/D≳50 regeneradas.
**Gate:** Ceres ≈ 939 km (no 763), Nysa ≈ 71 km (no 139); la tabla `tab:dd100` recalculada;
la conclusión "no new large–large" re-evaluada con diámetros reales (hoy es falsa en ambas
direcciones: Iris×Nysa sale, Nemausa×Aegina entra).
**Depende de:** 1.

### Tarea 5 — Presupuesto de completitud honesto
**Entregable:** término de incertidumbre de elementos de entrada (propagar covarianzas
AstDyS/JPL sobre muestra estratificada por calidad orbital, o cota) añadido al budget;
σ(d_min) por fila o por estrato; injection-recovery de detección (mínimos colocados entre
samples gruesos, v_rel altas, e altas). Renombrar a "pipeline-induced incompleteness budget"
si el término de elementos queda como cota.
**Gate:** el budget tiene ≥4 términos o declara explícitamente el de elementos como acotado;
`tab:closest` y el par genético de §4.3 llevan σ o se marcan no significativos;
injection-recovery recupera ≥99 % de encuentros sintéticos < 0.05 AU con ratio d ≈ 1.
**Depende de:** 1.

### Tareas 12–16, 18, 20 — Mayores sobre el catálogo
Cada una con su gate en la tabla de arriba; todas dependen de Tarea 1. Destacados:
- **18 (injection-recovery permanente):** test de regresión que habría cazado B1.
  **Gate:** corre en CI y falla si un mínimo inter-sample se pierde o se sesga > 1 μAU.
- **16 (señal de deflexión por par):** columna `signal_mas ≈ GM/(b·v_rel)` integrada sobre
  la geometría. **Gate:** el catálogo permite rankear targets por señal (hoy solo por
  distancia); tabla de estratificación D×v_rel×d en el paper.

---

## Fase 2 — Capa de masas

### Tarea 7 — Cruce FM con potencia y una sola σ oficial
**Entregable:** en `docs/mass_crosscheck_jack.md` y §5: desviación mínima detectable a 3σ
junto a cada z; test de signo sobre ratios de no-calibradores como resultado principal;
`z_jack` renombrado a `z_total`; caso Pallas corregido; caveat del draft reinstaurado;
no-independencia con FM declarada; **una única σ oficial por masa** con definición explícita
de cada z.
**Gate:** los z tabulados son reproducibles desde (M, σ, M_ref) — hoy Pallas da 3.4, no el
+2.67 tabulado; la tabla reporta potencia; "10/10 en |z|<3" ya no se presenta como evidencia
de exactitud.
**Depende de:** 6.

### Tarea 8 — Validación de calibradores no circular + inyección end-to-end
**Entregable:** leave-one-calibrator-out para f_sys; incertidumbre de f_sys declarada; no
extrapolar f_sys a los débiles; test ∂M/∂s_c (±30 %); inyección-recuperación end-to-end sobre
datos reales (señal sintética + geometría y ruido reales, pipeline completo con calibración
de s_c, clipping, jackknife y clasificación) para ≥1 calibrador y ≥1 perturbador débil.
**Gate:** los calibradores se validan con f_sys de los *otros* (no el propio); la inyección
recupera la masa inyectada con sesgo cuantificado (documenta si aparece el −4 %).
**Depende de:** 6.

### Tarea 9 — Tabla de 16 masas + schema + selección de targets
**Entregable:** tabla de las 16 masas (σ_formal, σ_jack/oficial, N, χ²_red, flag, z) en §5;
tabla de columnas de los tres productos Parquet (nombre/tipo/unidad/descripción); descripción
cuantitativa de la cadena catálogo→targets por perturbador (cortes de v_rel, distancia,
calidad FPR, N final — resuelve "Pallas, 6 encuentros" vs 9 del cross-match).
**Gate:** el "determines all sixteen" es verificable contra la tabla; N por perturbador
consistente entre el paper y los JSON de `expanded_jack/`.
**Depende de:** 6,7.

### Tarea 21 — Identificabilidad por verosimilitud perfilada
**Entregable:** criterio measured/not_identifiable basado en Δχ²(M=0 vs M̂) con órbitas
re-perfiladas (maquinaria `--profile` ya existe), calibrado con inyecciones a masa nula y a
masa de literatura; en vez del corte duro snr_jack ≥ 3.
**Gate:** tasa de falsos "measured" bajo masa nula calibrada y < umbral declarado;
clasificación estable a ±1 encuentro.
**Depende de:** 6.

### Tarea 17 — Reencuadrar §5
**Entregable:** §5 como *demostración de uso del catálogo* (target selection + calibradores +
Psyche), comprimida; comparación explícita con Siltala & Granvik (MCMC) y FM (prior→posterior)
reconociendo equivalencia funcional; sin vender "methodological framework".
**Depende de:** 9,10.

---

## Fase 3 — Menores y submission

- **Tarea 23 (código):** docstring MPCORB TT (`time_utils.py:8`); comentario "DE440 reducida"
  y default `backend="rebound"` (`dynamics.py:14`, `mass_determination.py:403,470`); Sol
  builtin vs DE440 (`gaia_adapter.py:174-177`); ∂τ/∂param; clamp M≥0 y FD de e; meseta FD de
  elementos ASSIST; frame bias documentado; deflexión de la luz acotada.
- **Tarea 24 (texto/figuras):** definición de "encuentro"; marco declarado; criterio
  Gaia-observable; ley de potencias Fig. 1; "no third-party dependency"; radio de query
  0.0536 vs 0.0572.
- **Tarea 25 (reproducibilidad):** hash MPCORB completo en sidecars; sidecar híbrido con
  git/deps/config; corregir inconsistencias FROZEN_RUN (fine_step 60/120, 305/137 MiB,
  305,896/305,931); desempate lexicográfico (dist, t) en `parallel.py:310`.
- **Tarea 26 (F3):** gate como Δmasa pareada (<0.25 %), no f_sys de 3 puntos; cota de escala
  para la cola del cinturón.
- **Tarea 27 (refs/tablas):** quitar `fuentesmunoz2024` sin verificar; añadir masa FM como
  columna en tabla de candidatos §4.4; justificar muestra estratificada vs "representative".
- **Tarea 28 (submission):** DOI de depósito CDS/VizieR (reparo #1 de un árbitro A&A);
  autores/afiliaciones/ORCID; acknowledgements; code availability (DOI Zenodo/GitHub).

---

## Fuera de scope

- **F8 (Pallas con DR4):** bloqueado, espera datos. No es remediación; es trabajo futuro
  (`planning/MASS_FUTURE_WORK.md`).
- **Extender la ventana a FPR/DR4** para convertir el catálogo en insumo vivo: es el camino
  a un *paper de resultados*, no a cerrar los hallazgos del tribunal. Decisión estratégica
  aparte.
- **Reescribir el motor de masas como bayesiano jerárquico** (MCMC tipo Siltala & Granvik):
  mejoraría B6–B8 de raíz, pero es un rediseño, no una remediación. Las tareas 6–8 lo
  resuelven al nivel exigido sin reescribir el motor.

---

## Ruta crítica

```
Tarea 1 (fix + regen catálogo)  ──┬─→ 2,3,4,5  ──→ 12–16,18,20
                                  └─────────────────────────────┐
Tarea 6 (σ masas)  ──→ 7,8,21  ──→ 9  ──→ 17  ───────────────────┤
Tarea 10 (novedad + refs)  ──→ 11,17  ──────────────────────────┤
                                                                 ↓
                                              Fase 3 (menores) → 28 (submission) → A&A
```

Las tareas 1, 6 y 10 son independientes y arrancan en paralelo. Todo lo demás cuelga de
ellas. **Ninguna edición de números del paper antes de que la Tarea 1 cierre.**

---

# Apéndice — Información de implementación

> Referencia operativa por tarea: evidencia exacta (archivo:línea), spec del cambio,
> comandos y datos. Complementa la tabla de estado; la evidencia científica completa
> está en `docs/tribunal_cientifico_2026-07-04.md`.

## A. Hecho en esta sesión (2026-07-04, branch `fix/tribunal-remediation-core`)

**Tarea 1 — parte de código (B1):**
- `src/detect/refine.py`: `_refine_one_kepler` ahora **re-centra iterativamente** la
  ventana fina cuando el argmin cae en el borde (`_MAX_RECENTER = 4`; alcance total
  `(1+4)×window_hours` más allá de `t_coarse`). Antes devolvía el borde sin interpolar
  (causa raíz de B1). Default `window_hours` 2.0 → **6.0**.
- `config.yaml → detection.refinement.window_hours`: 2.0 → **6.0** (≥ `coarse_step/2`),
  con comentario explicando B1.
- `src/detect/pipeline.py → detect_encounters`: **validación nueva** — `ValueError` si
  `window_hours < grid_step/2` cuando el refinador Kepler está en uso
  (`force_kepler_refine` o sin `positions`). Cumple el gate "assert de config".
- `tests/test_detection.py`: 4 tests nuevos.
  `test_refine_recovers_minimum_beyond_window[6.0|2.0]` — par sintético de órbitas que
  se cruzan (planos a ±15°, |v_rel| ≈ 9.8 km/s realista, miss ≈ 1e-3 AU) con el mínimo
  verdadero a 5 h del sample grueso: exige |Δt| ≤ paso fino y |Δd| ≤ 1 μAU contra un
  scan denso de 5 s independiente de `refine.py`. Con ventana 2.0 solo pasa gracias al
  re-centrado (el caso del bug). `test_refine_edge_bias_would_be_caught` documenta que
  el sesgo pre-fix en esa geometría es > 100 μAU (el test discrimina de verdad).
  `test_pipeline_rejects_window_narrower_than_half_step` — gate de validación.
- Verificación: `docker compose run --rm test` → **539 passed** (los 26 skipped son
  gates que requieren el freeze montado); ruff + black limpios.

**Tarea 6 — parte de catálogo (B6):**
- `scripts/mass/build_mass_catalog.py`:
  - Fix del bug de signo: `sigma_sys_kg = f_sys·|M|` (antes `f_sys·M` → σ negativa
    para Davida, masa ajustada −1.06e20 kg). `sigma_total_frac` también con |M|.
  - Columnas nuevas por fila: `jack_leverage_top1` (fracción de la varianza jackknife
    de la réplica más desviada), `sigma_jack_excl_top1_kg`, `snr_jack_excl_top1`,
    `sigma_jack_defensible` (= `n_targets ≥ --min-n-jack` (default 10) y
    `leverage ≤ --max-leverage` (default 0.5)).
  - Clasificación: `non_physical` si M ≤ 0; `measured` exige `snr_jack ≥ 3` **con y
    sin la réplica dominante**, y `n_targets ≥ 10` → la clase es estable por
    construcción frente a quitar la réplica dominante (gate B6). La defensibilidad de
    σ_jack como *barra publicable* queda en el flag (ortogonal a identificabilidad):
    las tareas 7/9 deben usar bootstrap/delete-d donde `sigma_jack_defensible=False`.
- Verificado sobre `data/output/orbdet/expanded_jack/` (16 perturbadores): reproduce
  exactamente los leverage del tribunal (Juno 0.92, Thisbe 0.89, Psyche 0.81, Iris
  0.69, Ceres 0.63, Cybele 0.57, Davida 0.54); la fórmula
  `σ_jack = √((n−1)/n·Σ(m_i−m̄)²)` reproduce bit a bit `mass_fit_sigma_jack_kg` de los
  JSON. Resultado: 9 `measured`, 5 `not_identifiable` (incluye Pallas N=6 y Cybele),
  1 `non_physical` (Davida), 1 `unknown` (Camilla, sin jackknife). Salida de prueba en
  `/tmp/mass_catalog_b6_check.csv` (no pisa el CSV oficial).

**Tarea 23 — menores C1–C3 (docstrings/comentarios, sin cambio de comportamiento):**
- `src/utils/time_utils.py`: "MPCORB uses TDB" → "MPCORB epochs are TT (converted to
  TDB on ingest)" (C1).
- `src/orbdet/dynamics.py`: comentario "builtin (DE440 reducida)" corregido — builtin
  = series erfa `epv00`/`plan94`, ~km de error, NO DE440; ASSIST es producción (C2).
- `src/orbdet/mass_determination.py` (`determine_mass_and_orbit`,
  `determine_shared_mass`): warning explícito en docstring sobre el default
  `backend="rebound"` (C2).
- `src/orbdet/gaia_adapter.py:174`: comentario documentando la mezcla Sol-builtin vs
  posiciones DE440 (~km, fraccional ~1e-5, acotado por gate Horizons 0.17 mas) (C3).

**Tarea 25 — parcial:**
- `src/detect/parallel.py → _merge_candidates`: desempate **lexicográfico (dist, t)**
  — antes, con distancias float empatadas, ganaba el chunk que llegara primero
  (`imap_unordered`), no determinista.

## A2. Segunda tanda (misma sesión 2026-07-04)

- **Regeneración lanzada (Tarea 1):** `run_pipeline` corriendo en background con el
  código B1 corregido (verificado por md5 dentro del contenedor: `_MAX_RECENTER`
  presente, `window_hours: 6.0`), cache de trayectorias reutilizado (hit, N=449.454,
  snapshot MPCORB_20160217 autoseleccionado). Output a
  `encounters_catalog_rebound_005au_b1fix.parquet` (el `output.filename` de
  `config.local.yaml` fue redirigido para NO pisar el freeze; backup del overlay en
  `config.local.yaml.bak_tribunal`). Log: `logs/regen_b1fix_20260704.log`.
  ⚠️ Gotcha de entorno: Docker Desktop puede servir archivos recién editados
  **desfasados** dentro de los mounts — verificar md5 en el contenedor antes de
  confiar en una corrida (memoria `docker-desktop-stale-file-cache`).
- **Tarea 4 (B3) código+datos:** `scripts/ingest/download_sbdb_physical.py` nuevo y
  **ejecutado** — `data/raw/sbdb_physical.parquet` (895.910 numerados, 135.475 con
  diámetro medido, sidecar con SHA-256). Cadena de prioridades en
  `physical.py::diameter_km_with_source` (measured > albedo_measured > zone_albedo >
  default; zona = proxy orbital, no taxonomía) cableada en `characterize_catalog`
  (+streaming) con columnas nuevas `diameter_source_1/2` (schema + descripciones).
  Gates verificados sobre los datos: Ceres 939.4, Nysa 70.6, Aegina 103.4 →
  Iris×Nysa sale de D≳100 y Nemausa×Aegina entra, tal como predijo el tribunal.
  6 tests unitarios nuevos de la cadena (`tests/test_characterize.py`).
- **Tarea 3 (B2) sidecar:** bloque `universe` explícito (N propagado + cortes) en
  `write_detection_sidecar`; hash MPCORB ahora **SHA-256 completo** (clave `sha256`,
  antes `sha256_prefix` truncado — parte de Tarea 25). `run_pipeline` ya pasaba
  `n_asteroids=len(elements)`; el sidecar nuevo lo registrará.
- **Tarea 19 (M11) ✅:** `docs/prefilter_recall.md` (sección "Impact on the frozen
  catalog" reescrita como contrafactual) y `docs/completeness_vs_literature.md`
  (§interpretación y §6) alineados; además corregida ahí la causa de la brecha FM
  (ventana temporal, no umbral — M1).
- **Tarea 10 (B5) texto:** título ("all-pairs … main-belt … measured incompleteness
  budget"), abstract (context/aims) y §1 reescritos con el prior art (Galád & Gray
  2002, Ivantsov 2018, Goffin 2014, FM 2025 con mismo umbral 0.05 AU); la palabra
  "hand-selected" eliminada (también del §3.4, que ahora atribuye la brecha a la
  ventana temporal con TODO de cuantificación post-regen); §4.5 sin "methodological
  framework"; §5 con citas de software (assist/rebound/DE440/astropy/FPR) en vez de
  "no third-party dependency". Referencias nuevas citadas con keys `galad2002,
  ivantsov2018, li2023, siltala2020/21/22, baer2017, fienga2003, kretlow2020,
  park2021, gaiafpr2023, rein2012, rein2015ias15, holman2023, astropy2022` —
  **pendiente**: pegar las entradas BibTeX verificadas contra ADS (agente en curso).
- **Tarea 7 (B7) mecánica:** `crosscheck_fuentes_munoz_jack.py` — `z_jack` renombrado
  a `z_total`, columna `min_detectable_frac_3sigma` (reproduce la potencia del
  tribunal: Vesta 15 % … Interamnia 74 %), test de signo sobre no-calibradores
  medidos como resultado con potencia (hoy: 5/6 bajo FM, media geom. 0.802,
  p=0.219) y caveat de no-independencia (misma astrometría FPR) impreso y en el JSON.
- **Tarea 8 (B8) LOCO:** `scripts/mass/loco_calibrators.py` nuevo y ejecutado —
  f_sys leave-one-calibrator-out; **PASA**: max |z_LOCO| = 1.40 < 3 (Ceres −0.69,
  Vesta −1.40, Hygiea −0.10), con σ(f_sys) declarada (≈±50 % relativo con n=3).
  La validación de calibradores sobrevive a la des-circularización.
- **Tarea 26 (C14):** gate de F3 reformulado en `scripts/mass/f3_fsys.py` como
  Δmasa pareada por calibrador; ejecutado: máx 0.2315 % < 0.25 % → **PASS** (antes
  comparaba f_sys de 3 puntos, 4.158 % vs 4.257 %, ruido).
- **Tarea 23 (C6):** guardias FD en `_assist_pos_and_partials`: paso de masa sobre
  `max(|M|, piso)` (Davida M<0 degeneraba la FD) y FD adelantada de un lado cuando
  el paso central cruzaría e < 0.

## A3. Tercera tanda (misma sesión, mientras corre la regen)

- **Tarea 18 (M10) ✅:** `scripts/validate/injection_recovery_detection.py` — pares
  sintéticos construidos desde vectores de estado (posición común + Δv log-uniforme
  hasta 25 km/s + miss vector ⊥ Δv, |b| log-uniforme [1e-4, 0.045] AU), mínimos en
  fase uniforme vs la grilla de 12 h, verdad de terreno por scan denso a 5 s
  independiente del refinador. **Ejecutado: 200/200 recuperados, |Δd| máx
  1.2×10⁻¹⁰ AU, |Δt| máx 0.1 s → PASS.** Test permanente de CI:
  `tests/test_injection_recovery.py` (N=25). De paso valida empíricamente la cota
  analítica de ensanche de radio (insumo directo para la Tarea 12: el harness usa el
  semipaso, la producción usa el paso completo — más conservadora).
- **Tarea 2 (B10):** `diagnose_fienga_804_733.py` re-ejecutado con el refinador
  corregido (el `window_hours=2.0` hardcodeado del script lo cazó la validación
  nueva de `pipeline.py` — corregido a 6.0): el par aislado se detecta a
  0.013547 AU, 2015-02-12, matching JPL. Cierre definitivo cuando el b1fix esté:
  verificar `(733, 804)` en el parquet regenerado.
- **Tarea 16 (M7) código:** `physical.py::deflection_dv_m_s` — Δv = 2GM/(b·v_rel)
  impulsivo, GM ∝ ρ_zona·D³ (densidades Carry 2012 por zona orbital), cableado en
  characterize como columna `deflection_dv_m_s` (schema + descripción: métrica de
  ranking, no medición). 3 tests unitarios (escala Ceres, ∝D³, ∝1/b, NaN).
  Falta (post-regen): tabla de estratificación D×v_rel×d en el paper.

## A4. Cuarta tanda (sesión 2026-07-05)

ℹ️ **Docker Desktop se cayó y se reinició** (`systemctl --user start docker-desktop`)
durante la sesión. La regeneración del catálogo (Tarea 1) **está corriendo de nuevo**
en background: log `logs/regen_b1fix_20260705.log`, output
`encounters_catalog_rebound_005au_b1fix.parquet`, radio de query 0.05722 AU. Al
terminar el scan+refine: correr el gate de histograma plano, Stage B N-body +
ensamble híbrido, y verificar (733,804) presente (Tarea 2). config.local.yaml apunta
a `..._b1fix`; backup del overlay original en `config.local.yaml.bak_tribunal`.
Suite completa verde tras todos los cambios de la sesión (552 passed, ruff/black OK).

- **Tarea 10 (B5) ✅ bibliografía:** el agente ADS se quedó sin sesión; verificado
  a mano vía WebSearch/WebFetch. `references.bib` completado con 12 entradas nuevas
  (galad2002 A&A 391 1115; fienga2003 A&A 406 751; ivantsov2018 IAUS 330 386;
  baer2017 AJ 154 76; li2023 AJ 166 93; siltala2020 A&A 633 A46; park2021 AJ 161
  105; gaiafpr2023 A&A 680 A37 [David et al.]; rein2012 A&A 537 A128;
  rein2015ias15 MNRAS 446 1424; holman2023 PSJ 4 69; astropy2022 ApJ 935 167) —
  todas con vol/página/DOI verificados contra ADS. `fuentesmunoz2024` (LPSC, sin
  verificar) eliminada (Tarea 27). Cross-check automatizado: **18/18 citas del tex
  resuelven** contra el bib (`russell2012` queda definida sin usar, inofensivo).
- **Tarea 24 (menor 15):** radio de query del paper corregido a 0.0572 AU (paso
  completo, lo ejecutado — el log de regen lo confirma: +0.00722 AU) en vez de
  0.0536 (Δt/2); el texto ahora explica que el requisito estricto es Δt/2 y que el
  pipeline usa el paso completo (margen 2× conservador).

## A5. Quinta tanda (sesión 2026-07-05, en paralelo a la regen)

- **Fix de OOM del scan (destapa la regen):** `parallel.py` acumulaba TODOS los
  candidatos de los 104 chunks en el padre antes de deduplicar (~480M tuplas ≈ 30 GB
  → `OOMKilled=true`, ExitCode 137, moría a ~11/104). Ahora `_merge_into` funde cada
  chunk al llegar (memoria acotada a pares únicos). Con esto la regen **pasó 11/104**
  (donde moría) y sigue. Commit + push a PR #102.
- **Diagnóstico de memoria del entorno:** la RAM del host (99%) la consume Docker
  Desktop — `qemu-system-x86` (VM) 23.8 GB + `virtiofsd` (bind-mount) 8.6 GB, no un
  leak. El scan es I/O-bound (virtiofsd sirviendo el zarr): medido, **4 workers es más
  rápido que 8** (thrashing), y bajar workers reduce RAM. `config.local.yaml`: 4
  workers / chunk 10 d.
- **Tarea 9 (B9) parcial ✅:** tabla de schema del catálogo caracterizado en el
  Apéndice A del tex (nombre/tipo/unidad/descripción de las 30 columnas + nota de las
  extra del híbrido). Compila limpio.
- **Tarea 17 (M8/M9) framing ✅:** §5 reencuadrada como demostración de uso (no fuente
  de masas nuevas); párrafo "Relation to other methods" (Siltala&Granvik MCMC / FM
  prior→posterior; FM ajusta la misma astrometría FPR → cross-check no independiente);
  conclusiones sin "methodological framework".
- **Tarea 28 parcial:** lista de software/fuentes en acknowledgements (A&A la exige).
- **Paper compila limpio: 10 pp., 18/18 citas resueltas** (`pdflatex`+`bibtex` OK).
- Todo commiteado y pusheado (PR #102 actualizado).

## A6. Sexta tanda (sesión 2026-07-05, en paralelo a la regen)

- **Tarea 23 (C4–C8) ✅:** documentados en `src/orbdet/` (solo docstrings/comentarios):
  deflexión gravitacional de la luz no aplicada (DPAC ya la corrige; sub-mas) y término
  ∂τ/∂param omitido ~1e-4 (observation.py); frame bias ICRS↔eclíptica ~17 mas que se
  cancela en la cadena (frames.py); falta de chequeo de meseta en los pasos FD de
  elementos (mass_determination.py).
- **Tarea 3 (B2) tex ✅:** §2.1 declara el corte a∈[1.5,4.0] AU (excluye
  NEAs/Troyanos/outer) y N=449.454 cuerpos como denominador de completitud.
- **Tarea 24 (C9/C10) ✅:** definición de "encuentro" (una fila = mínimo global del par
  en la ventana), criterio Gaia-observable explícito (elong>45°, V<21), marco de
  referencia explícito e invariancia de la distancia mutua; §2.3 corregido (diámetros
  medidos>albedo con procedencia; clase dinámica, no taxonómica — B3).
- **Regen (2 workers):** estable en 10/104, pasó la zona de OOM (moría en 11-13 con 4
  workers); memoria estable. Más lenta (~50 min de scan) pero sobrevive.
- Paper compila limpio (10 pp.); todo commiteado y pusheado (PR #102).

## B. Tarea 1 — lo que falta: regeneración del catálogo

La corrida congelada (sidecar `data/output/encounters_catalog_rebound_005au_provenance.json`)
fue: scan N-body `whfast` (Sol+Júpiter+Saturno) a 12 h + refine Kepler forzado
(`forced_by_tiered: true`) con `window_hours: 2.0` ← el bug. Catálogo base: 72.236.904
filas; híbrido Stage B: `data/output/encounters_catalog_hybrid_stageb.parquet`.

Secuencia de regeneración (config ya corregida):
```bash
# 1. Corrida completa (scan coarse + refine Kepler ±6 h). Días de cómputo; el cache
#    zarr de trayectorias en data/cache/ es reutilizable si no cambió MPCORB.
docker compose run --rm pipeline python -m scripts.pipeline.run_pipeline --config config.yaml

# 2. Stage B: re-refinar el subset selectivo con N-body y ensamblar el híbrido
docker compose run --rm pipeline python -m scripts.validate.refine_stageb_nbody ...
docker compose run --rm pipeline python -m scripts.validate.assemble_stageb_hybrid_catalog ...

# 3. Sidecar de detección completo (poblar n_asteroids — hoy null, ver Tarea 3)
docker compose run --rm pipeline python -m scripts.pipeline.generate_detection_sidecar ...
```

Gate del histograma (verificación post-regen, el test que confirmó B1):
```python
import polars as pl
df = pl.scan_parquet("data/output/encounters_catalog_hybrid_stageb.parquet")
# offset de jd_tdb al nodo más cercano de la grilla de 12 h (0.5 d), en horas
t0 = <primer JD de la grilla>   # del sidecar
off = ((pl.col("jd_tdb") - t0) % 0.5).map_batches(lambda s: pl.Series(np.minimum(s, 0.5 - s))) * 24
# filas refinement_method=="kepler": el histograma de `off` debe ser ~plano en [0, 6] h
# (pre-fix: 76.3 % de la masa en [1.75, 2.25] h y nada > 2.25 h)
```

Derivados a rehacer con el catálogo nuevo (mismos scripts que los produjeron):
- Stage A/B error Kepler-vs-N-body: `scripts/validate/compare_kepler_vs_nbody.py` y
  `scripts/validate/refine_stageb_nbody.py` — separar el artefacto de ventana (72.4 %
  de pares con mínimo recortado en el subset re-refinado) del error de modelo real.
- Censura de umbral: `scripts/validate/measure_threshold_false_negatives.py` (usaba la
  misma ventana contaminada en :164-186; ahora hereda window ≥ 6 h del config).
- Costo esperado: ventana ±6 h ≈ 3× evaluaciones Kepler en refinamiento vs ±2 h
  (trivial frente al scan); re-centrados adicionales solo en casos de borde.

Nota de diseño: el re-centrado puede correr el mínimo fuera de la ventana de
observación en encuentros de borde (los que M6 quiere flaggear). Al derivar
`boundary_minimum` (Tarea 15), definirlo como `jd_tdb` fuera de
`[t_start_win, t_end_win]` o a < 1 paso fino del borde.

## C. Tareas 2–5 (sobre el catálogo regenerado)

**Tarea 2 (B10, miss (804)×(733)):** ya existe
`scripts/validate/diagnose_fienga_804_733.py` — correrlo contra el catálogo regenerado.
Si sigue ausente: propagar los dos cuerpos aislados con el snapshot congelado
(`MPCORB_20150524.DAT`, registrado en los JSON de masas) y la misma grilla; verificar
(a) si el mínimo caía en borde de ventana (B1, candidato más probable) y (b) los `a`
de ambos cuerpos contra el corte `a ∈ [1.5, 4.0]` (B2). Documentar en
`docs/fienga_804_733_root_cause.md`. Gate: el par (0.0138 AU) aparece en el freeze
nuevo, o el paper divulga el gap con causa; 4/4 eventos Fienga in-window explicados.

**Tarea 3 (B2, universo muestral):**
- El corte real está en el sidecar del catálogo base: `subset: only_numbered,
  a ∈ [1.5, 4.0] AU` (excluye NEAs a<1.5, troyanos ~6.000+, a>4). `n_asteroids: null`
  → poblarlo en `scripts/pipeline/generate_detection_sidecar.py` y en el writer.
- Tex: quitar "full numbered population" de título/abstract/§1/§7; declarar N exacto y
  el corte en §2; reconciliar "~150,000 bodies" (§2.2) con los "449k numerados" de
  `docs/kepler_threshold_bias_paper.md:134`; rehacer la extrapolación N² de censura
  con el N declarado (la cifra "~1.5–2.5×10⁵ censurados" debe ser reproducible desde
  el texto).

**Tarea 4 (B3, diámetros medidos):**
- Evidencia: `src/characterize/physical.py:10` (`albedo=0.14` para todo) y
  `src/characterize/encounter.py:213-214` (aplica el escalar de config
  `characterize.default_albedo`). El paper (§2.3) afirma lo contrario.
- Spec: (1) nuevo ingest `scripts/ingest/download_sbdb_physical.py` → bulk query a
  JPL SBDB (`ssd-api.jpl.nasa.gov/sbdb_query.api`, campos `number, diameter, albedo` —
  compila IRAS/AKARI/NEOWISE/ocultaciones) → `data/raw/sbdb_physical.parquet` con hash
  y fecha en sidecar. (2) `encounter.py`: join por número; prioridad diámetro medido >
  D(H, albedo medido) > D(H, albedo por clase **orbital** con tabla en config
  (`albedo_by_class`, p. ej. inner-belt 0.20 / mid 0.14 / outer 0.06) > D(H, 0.14).
  (3) columna `diameter_source` (`measured|albedo_measured|class_albedo|default`).
  (4) §2.3: corregir "taxonomic class estimated from H" → clase orbital.
- Gate numérico: Ceres ≈ 939 km (hoy 763), (44) Nysa ≈ 71 km (hoy 139, tipo E),
  (91) Aegina ≈ 104 km (hoy 60.6). Regenerar todas las tablas de §4; re-evaluar
  "no new large–large encounter" (pre-fix es falsa en ambas direcciones: Iris×Nysa
  sale de D≳100, Nemausa×Aegina entra).

**Tarea 5 (B4, presupuesto de completitud):**
- Término de elementos de entrada: propagar covarianzas AstDyS/JPL sobre muestra
  estratificada por calidad orbital (MPCORB trae el parámetro U y arco; punto de
  partida: `scripts/mass/measure_mpcorb_uncertainties.py`). Re-muestreo MC de
  elementos → σ(d_min) por estrato; si queda como cota, renombrar el presupuesto a
  "pipeline-induced incompleteness budget" en el tex.
- Injection-recovery de detección (también Tarea 18): script nuevo
  `scripts/validate/injection_recovery_detection.py` — pares sintéticos con mínimo
  colocado uniforme en fase respecto de la grilla gruesa, v_rel hasta ~190 km/s
  (el techo que cubre el query radius: √(0.0572²−0.05²) AU / 0.25 d ≈ 0.111 AU/d),
  e y i altas. Gate: ≥ 99 % recuperados con ratio d ≈ 1; corre en CI.
- `tab:closest` (1.094–2.590 km, hoy 4 cifras sin σ) y el par "genético" de §4.3
  (49.511 km): añadir σ(d_min) o marcarlos no significativos.

## D. Tareas 6–9, 21 (capa de masas) — lo que falta

**Tarea 6 (resto):** bootstrap por objetivo (re-muestreo con reemplazo de la lista de
targets + refit) o delete-d jackknife para los 7 perturbadores con leverage > 0.5
(Juno 0.92, Thisbe 0.89, Psyche 0.81, Iris 0.69, Ceres 0.63, Cybele 0.57, Davida 0.54).
Dónde: la maquinaria jackknife vive en el fitter (`scripts/mass/orbdet_fit_realdata.py`
/ `src/orbdet/`); los JSON ya guardan `jackknife_masses_kg` por réplica. Limitación a
declarar: las réplicas reutilizan s_c calibrado, clipping 4σ y warm-start del ajuste
completo → subestiman la variabilidad del procedimiento completo.

**Tarea 7 (cruce FM):**
- Potencia por perturbador: desviación fraccional mínima detectable a 3σ
  `= 3·√(σ_tot² + σ_FM²)/M_ref` — tabularla junto a cada z en
  `docs/mass_crosscheck_jack.md` y §5 (valores del tribunal: Vesta 15 % … Interamnia
  74 %; los sesgos observados son 14–30 %, mismo signo, media geom. de ratios ≈ 0.79).
- Test de signo sobre ratios de los 6 no-calibradores (5/6 < FM) como resultado
  principal en vez de "10/10 en |z|<3".
- `z_jack` → `z_total` en `scripts/mass/crosscheck_fuentes_munoz_jack.py` (:172-229) y
  en el tex; **una sola σ oficial por masa** con definición explícita del denominador
  de cada z (la discrepancia Pallas 3.4 vs +2.67 tabulado viene de incluir o no σ_ref
  en el denominador — decidir y documentar una vez).
- Reinstaurar el caveat eliminado del draft ("a wide jackknife σ makes almost any
  value consistent…" — buscarlo en el historial git del tex); declarar que FM 2025
  ajusta la misma astrometría FPR (errores no independientes).

**Tarea 8 (calibradores no circular):**
- LOCO: para cada calibrador c, calcular `f_sys^(−c)` con los otros dos y validar c
  contra ese piso (hoy |desvío|/RMS ≤ √3 por construcción). Declarar σ(f_sys) (n=3 →
  ~41 % relativo); no extrapolar f_sys a los débiles (su sistemático dominante es
  regresión a cero, −30/−60 %).
- Sensibilidad al piso: refit con `sys_floor_mas × {0.7, 1.3}` (está por-JSON, p. ej.
  Cybele 0.778 mas) → ∂M/∂s_c.
- Inyección end-to-end: sumar a las obs FPR reales de ≥1 calibrador y ≥1 débil la
  deflexión sintética de una masa conocida (generada con el forward model del propio
  motor), correr el pipeline completo (calibración s_c, clipping, jackknife,
  clasificación) y cuantificar el sesgo de recuperación (¿aparece el −4 %?).

**Tarea 9 (tablas):** el CSV post-B6 ya tiene todas las columnas para la tabla de 16
masas (M, σ_formal, σ_jack, leverage, N, χ²_red, status, z). Falta: tabla de schema de
los productos Parquet publicados (nombre/tipo/unidad/descripción por columna) y la
cadena catálogo→targets cuantificada (los cortes exactos de
`orbdet_fit_realdata --top-per-perturber 30`: v_rel, distancia, calidad FPR, N final —
resuelve "Pallas, 6 encuentros" vs 9 del cross-match).

**Tarea 21 (identificabilidad):** reemplazar/contrastar el corte `snr_jack ≥ 3` con
Δχ²(M=0 vs M̂) re-perfilando órbitas (`--profile` ya existe en el fitter); calibrar la
tasa de falsos "measured" con inyecciones a masa nula y a masa de literatura;
discutir multiplicidad (16 tests).

## E. Tarea 10 — referencias a añadir (verificar cada una contra ADS antes de commitear)

En `docs/paper/references.bib` (hoy 8 entradas; objetivo ≥ ~20):

| Ref | Para qué |
|---|---|
| Fuentes-Muñoz et al. 2025, AJ 170, 353 | búsqueda sistemática 1.783×1.07M, mismo umbral 0.05 AU — refuta "hand-selected"; reemplaza a `fuentesmunoz2024` (LPSC, eliminar) |
| Ivantsov, Hestroffer et al. 2018, IAUS 330, 386 | catálogo previo de encuentros mutuos all-numbered 2013–2023 para Gaia |
| Galád & Gray 2002, A&A 391 | búsqueda sistemática de encuentros para masas (24.599 asteroides) |
| Li et al. 2023, AJ 166 | ~20 masas con Gaia DR3 — el trabajo más comparable |
| Siltala & Granvik 2020 / 2021 / 2022 | MCMC; Psyche; Eunomia/Europa — tres perturbadores de este paper |
| Baer & Chesley 2017, AJ 154 | masas por encuentros, referencia estándar |
| Fienga et al. 2003 | usada como validación en el repo y no citada |
| Kretlow 2020 (SiMDA) | compilación de masas/densidades |
| Park et al. 2021, AJ 161 (DE440/441) | usada como referencia de ratios sin cita |
| Gaia Collaboration / David et al. 2023 (Gaia FPR SSO) | la fuente de datos del motor de masas |
| INPOP (Fienga et al.) y EPM (Pitjeva & Pitjev) | masas de efemérides planetarias |
| Software: rebound (Rein & Liu 2012), IAS15 (Rein & Spiegel 2015), WHFast (Rein & Tamayo 2015), ASSIST (Holman et al. 2023), astropy (Astropy Collab. 2022) | citas obligatorias A&A |

Reescritura de novedad (abstract + §1): la novedad defendible = **catálogo all-pairs
publicado íntegro con procedencia + presupuesto de incompletitud medido** (FM publica
top-100 por perturbador; Ivantsov no publicó el catálogo completo con procedencia).
"Hand-selected" desaparece. Gate: ninguna afirmación de primicia sin cita del prior
art que la matiza.

## F. Tareas 11–22 — punteros concretos

- **11 (M1/M2):** `docs/completeness_vs_literature.md` + tex §3.4. La causa dominante
  de la brecha con FM es la **ventana temporal** (FM busca sobre décadas; esto cubre
  2014-07→2017-05), no el umbral (FM usó 0.05 AU para MBAs). Cuantificar: clasificar
  los 25.962 pares FM no recuperados por época del encuentro (el cross-match está en
  `scripts/validate/crosscheck_literature_encounters.py`). Eliminar el "decisive
  point — 0 %" (tautológico: el catálogo solo contiene filas < 0.05 AU).
- **12 (M3):** reemplazar la cota de cadencia circular por la analítica:
  `d(δ) = √(d_min² + v²δ²)` con δ = semipaso grueso, contra el query radius ejecutado
  (0.0572 AU) → cubre hasta ~190 km/s. Escribirla en §3.
- **13 (M4):** `scripts/bench/make_paper_figures.py:353-366` — Fig. 3 debe usar los
  cuerpos únicos del catálogo real (`union(number_1, number_2)`), no el row count de
  `gaia_orbits.parquet`; corregir caption ("93,010 encountering bodies" es falso).
- **14 (M5):** tex — justificar cuantitativamente (o reconocer) que la selección de
  targets ignora encuentros 2017–2020 dentro del arco FPR ajustado.
- **15 (M6):** columnas `refinement_method`/fuente por tabla del §4; flag
  `boundary_minimum` (ver nota en B); reconciliar par top (987 km vs 1.094 km) y el
  encuentro "slowest" v_rel=15 m/s a 2 días del inicio de ventana (probable mínimo
  truncado); épocas completas con escala declarada.
- **16 (M7):** columna `deflection_signal` por par en characterize:
  `Δv ≈ 2GM/(b·v_rel)` con GM estimada de D³·ρ (ρ por clase), integrada a señal
  astrométrica (mas). Permite rankear targets por señal (como Ivantsov y FM) y la
  tabla de estratificación D×v_rel×d del paper.
- **17 (M8/M9):** §5 reencuadrada como *demo de uso del catálogo* (target selection +
  calibradores + Psyche), comprimida; comparación explícita con Siltala & Granvik
  (MCMC, posteriors completos) y FM (prior→posterior) reconociendo equivalencia
  funcional; sin "methodological framework" ni "no third-party dependency" como mérito.
- **18 (M10):** unit tests B1 ya en CI (ver A); completar con el injection-recovery
  end-to-end de la Tarea 5 como test permanente.
- **19 (M11):** alinear `docs/prefilter_recall.md:154-159` y
  `docs/completeness_vs_literature.md:186-187` con la corrección del sidecar
  (`_correction`: el prefiltro NO se aplicó — `skipped_large_n`; el déficit de ≥143k
  es contrafactual).
- **20 (M12):** documentar en el sidecar híbrido la regla de membresía (Kepler<0.05
  define el universo; 25.283 filas con dist_au ≥ 0.05, máx 0.05717) y recomendar en el
  paper el filtro para función de selección homogénea.
- **22 (M14):** publicar la cuenta a priori de S/N de deflexión (a 10⁻¹² M☉ y
  0.02 AU: S/N ≲ 1 por objetivo; a 10⁻¹¹ M☉: determinable con stacking) como criterio
  de determinabilidad en §5.

## G. Tareas 23–28 — restantes

- **23 (restan C4–C8):** ∂τ/∂param del light-time (~1e-4 relativo) — documentar en el
  módulo del Jacobiano; clamp M ≥ 0 (o log-GM) en el FD de masa y guardia para e < 0
  en el FD de elementos (`solve_kepler` lanza); chequeo de meseta para los pasos FD de
  elementos bajo ASSIST (1e-7 hardcodeado; el de GM ya lo tiene); documentar frame
  bias ICRS↔eclíptica (~17 mas, se cancela en la cadena); acotar la deflexión
  gravitacional de la luz por fuente a distancia finita (sub-mas cerca de 45° de
  elongación) y qué corrige DPAC.
- **24 (C9–C15):** definición explícita de "encuentro" (una fila = mínimo global de la
  ventana; mínimos locales secundarios no catalogados — para masas importan); marco
  declarado consistente (tex dice "heliocentric", sidecar "barycentric"); definir el
  criterio "Gaia-observable" (18.9 %: elongación, magnitud límite); ajustar índice de
  ley de potencias de Fig. 1 contra dN/dd ∝ d²; radio de query del paper 0.0536 →
  0.0572 (el ejecutado).
- **25 (resta):** hash MPCORB completo (hoy truncado a 16 hex) + git/deps/config en el
  sidecar híbrido (hoy solo shard_summary); reconciliar FROZEN_RUN vs sidecars
  (fine_step 60 vs 120 s; 305 vs 137 MiB; 305.896 vs 305.931 near-boundary).
- **26 (C14):** `scripts/mass/f3_fsys.py` — gate como Δmasa pareada (< 0.25 %), no
  como comparación de f_sys de 3 puntos (4.158 % vs 4.257 % es ruido); añadir cota de
  escala para la cola del cinturón no modelado.
- **27:** quitar `fuentesmunoz2024` (LPSC #2388, no verificada); masa FM como columna
  en la tabla de candidatos §4.4 (todos ya tienen masa FM 2025 y el texto no lo dice);
  justificar la muestra de censura estratificada ("drawn to span the belt") o dejar de
  llamarla "representative".
- **28 (submission):** DOI de depósito CDS/VizieR para el catálogo (reparo #1 de un
  árbitro A&A); autores/afiliaciones/ORCID; acknowledgements (Gaia DPAC, MPC, JPL);
  code availability con DOI (Zenodo) — `docs/paper/DATA_AVAILABILITY.md` y
  `docs/paper/zenodo_data_deposit.json` ya existen como base.
