# Continuación — dónde retomar el proyecto

> **Estado:** 🟢 PUNTO DE RETOMA — 2026-07-03
> Handoff para retomar en otra sesión. Resume qué está hecho, qué falta (con
> comandos concretos), y los gotchas del entorno para no repetir errores.
> Contexto científico en [`docs/dataset_paper_draft.md`](../docs/dataset_paper_draft.md),
> resultados de masas en [`docs/mass_determination_results.md`](../docs/mass_determination_results.md),
> plan de publicación en [`PUBLISH_PUSH_PLAN.md`](PUBLISH_PUSH_PLAN.md).

## 1. Estado actual (todo en `main`)

El grueso científico y de ingeniería está **cerrado y mergeado**. Resumen:

| Área | Estado |
|------|--------|
| Detección/caracterización de encuentros | ✅ Catálogo congelado de 72.2M encuentros 3D < 0.05 AU con **completitud medida** (ver `FROZEN_RUN.md`). |
| Motor de masas `src/orbdet/` | ✅ Ajuste conjunto órbita+masa, cov en bloques por FOV (ICC=0.32), 16 perturbadores. |
| F1 σ jackknife | ✅ `mass_catalog_jack.csv`; σ externa reemplaza la formal. |
| F2 identificabilidad | ✅ measured / not_identifiable / cota por curvatura de χ². |
| F3 sesgo −4 % | ✅ **hipótesis refutada** (fondo 16→35 no lo mueve). `docs/mass_f3_background_extension.md`. |
| F4 perturbadores fuera de los 16 | ✅ rama custom; Fortuna χ²_red=0.977, Metis 0.981. `docs/mass_layer_f4_design.md`. |
| F5 cruce Fuentes-Muñoz + σ jack | ✅ 10/10 medidas en \|z\|<3. `docs/mass_crosscheck_jack.md`. |
| F7 cruce Goffin 2014 | ✅ mergeado (PR #84). |
| Minería de encuentros notables | ✅ `docs/notable_encounters.md`; sin grande-grande nuevo. |
| Dataset paper | 🟡 **draft submittable** (`docs/dataset_paper_draft.md`): §1–§7, 4 figuras, completitud vs literatura, refs A&A. Falta lo de §2 abajo. |
| Dashboard (fase 7) | ✅ `src/dashboard/` con tab de masas; `docs/dashboard_status.md`. |

**Conclusión de novedad:** no hay masa nueva ni encuentro espectacular inédito.
Lo publicable = **dataset paper** (catálogo + completitud medida, nadie lo publica
así) + **methods** (FOV-block cov, σ jackknife, identificabilidad, F4 generaliza,
F3 negativo). Es dataset paper + methods paper, no paper de descubrimiento.

## 2. Qué falta (prioridad alta → baja)

### 2.1 Submission del dataset paper (mayor retorno; es trabajo de escritorio)
El draft está en `docs/dataset_paper_draft.md`. **Journal elegido: A&A.** El
paquete LaTeX A&A vive en `docs/paper/` (`aa_encounters.tex` + `references.bib`
+ `README.md`); compila limpio (probado con stub de `aa.cls`; el build real
necesita `aa.cls`/`aa.bst` del sitio A&A — ver `docs/paper/README.md`).

Hecho (2026-07-03, branch `docs/paper-submission-prep`):
- [x] **Referencias verificadas contra ADS**: Goffin 2014 (A&A 565 A56),
  Michalak 2000 (A&A 360 363), Park 2016 (Nature 537 515), Russell 2012
  (Science 336 684), Tanga 2023 (A&A 674 A12), Vernazza 2020 (Nat.Astron. 4 136),
  FM 2025 (AJ 170 353) — **todas OK salvo autores de FM 2025**: eran
  Fuentes-Muñoz, Farnocchia, Giorgini, Park (el draft había copiado mal
  Scheeres/Tanga → corregido en draft y `.bib`).
- [x] **Journal elegido: A&A**; paquete LaTeX generado en `docs/paper/`.
- [x] **Consistencia de números**: auditoría completa (abstract/cuerpo/tablas/
  captions vs docs de respaldo) → 0 contradicciones duras, aritmética de
  porcentajes OK. Corregida imprecisión del abstract ("sixteen … consistent" →
  "determines sixteen … of which the ten identifiable masses are consistent").
- [x] **Notas internas quitadas** en la versión LaTeX (el markdown conserva el
  bloque `<!-- Notas -->` y "Figures" a propósito, como copia de referencia).
- [x] Cosmético fig3: quitado el label duplicado y **fig3 regenerada** (2026-07-04,
  PDF/PNG en `docs/figures/` actualizados).
- [x] **Build A&A real verificado**: `aa.cls`/`aa.bst`/`linenoaa.sty` (kit oficial v9.4)
  incluidos en `docs/paper/`; PDF compilado `docs/paper/aa_encounters.pdf` (7 pp.).
- [x] **Números solo-de-figura** (93.010 cuerpos; 4.25 km/s) documentados en
  `docs/paper/figures_provenance.md`.
- [x] **Revisión tipo tribunal** completa en `docs/paper_referee_report.md` (panel de 4
  lentes + comparación con literatura verificada + veredicto). Cerró en el proceso el
  hueco de cadencia (radio ensanchado bracketea v_rel≤25 km/s; residual ≲10⁻³ %) y afinó
  ventana DR3, prefiltro, deflexión mutua, "231→232 masas". **Veredicto: minor–moderate
  revision** (ciencia sólida; resta lo mecánico/author-owned).

Falta (author-owned, marcado `TODO(author)` en el `.tex`):
- [ ] **Autores, afiliaciones, ORCID, e-mail de contacto** (`\author`,`\institute`).
- [ ] **Acknowledgements** (financiación, software, boilerplate Gaia/DPAC).
- [ ] **Data availability DOI** — **reparo #1 más probable de un árbitro A&A** (depósito
  CDS/VizieR mandatorio y contractualmente forzado). Scaffolding listo:
  `docs/paper/DATA_AVAILABILITY.md` + `zenodo_data_deposit.json`.
- [ ] **Injection–recovery de la completitud** — reparo científico más probable: un
  árbitro pedirá validación por inyección sintética (FN y FP), más allá de la medición
  por re-refinamiento. Ver `docs/paper_referee_report.md` §F.4.
- [ ] **Nº de abstract LPSC de FM 2024** (`55, 2388`): confirmar contra ADS.

**Regenerar figuras** (matplotlib NO está en la imagen — instalar ad-hoc):
```bash
docker compose run --rm pipeline bash -c \
  "pip install --quiet matplotlib && python -m scripts.bench.make_paper_figures"
mv data/output/figures/*.png data/output/figures/*.pdf docs/figures/   # docs/ no está montado
```

### 2.2 F6 — parcial ∂x/∂GM analítica ✅ (backend `rebound`; `assist` bloqueado)
**Hecho (2026-07-04).** `partial_wrt_gm_variational` (partícula variacional de masa
de REBOUND) integra ∂x/∂GM en una sola propagación por sentido; cableada en la rama
`rebound` de `mass_determination` (`gm_variational=True`). Gate cumplido: vs FD
< 1e-6 (`test_dgm_variational_matches_fd`), ahorra 2 props/Jacobiano.
- **Bloqueado en `assist`** (producción): las fuerzas de la efeméride no propagan
  partículas variacionales de REBOUND → allí ∂x/∂GM sigue por FD. Camino futuro:
  JVP `(∂a/∂r)·s` por diferencia finita direccional de la aceleración ASSIST.
- Detalle: `docs/mass_layer_f6_analytic_gm.md`, `planning/MASS_FUTURE_WORK.md` §F6.

### 2.3 F8 — Pallas con Gaia DR4 (bloqueado: espera datos)
Pallas tiene solo 6–7 objetivos < 0.05 AU en DR3 (target-limited). Cuando salga DR4:
```bash
docker compose run --rm pipeline python -m scripts.mass.orbdet_fit_realdata \
    --perturber 2 --release dr4 --from-catalog <catálogo DR4> \
    --top-per-perturber 40 --workers 24 --jackknife --out data/output/orbdet/pallas_dr4.json
```
Gate: N(Pallas) ≥ 20. No accionable hasta que DR4 esté disponible.

### 2.4 Cruce Goffin pair-level — ✅ verificado NO aplicable (2026-07-04)
El cruce de completitud (§3.4) usó FM 2025. Goffin 2014 a nivel de pares **no es
aplicable, no bloqueado**: se extrajo el PDF (`papers/goffin_2014_aa565_A56.pdf`,
pdfplumber ad-hoc) y se confirmó que Goffin hace un ajuste **simultáneo global** de
todas las masas contra los residuos astrométricos de toda la población — no hay listas
de objetivos por perturbador, así que no existe conjunto de pares (perturbador,
objetivo) que cruzar (table6 de VizieR es una compilación de masas de literatura, no
pares). El §3.4 del paper se afinó para decir esto explícitamente. Cerrado.

### 2.5 Submission del paper — scaffolding de datos listo (2026-07-04)
En `docs/paper/`: `DATA_AVAILABILITY.md` (checklist Zenodo+CDS/VizieR con tabla de
archivos/SHA), `zenodo_data_deposit.json` (metadata de depósito, faltan `creators` y el
DOI del paper), `figures_provenance.md` (los dos números solo-de-figura 93.010 y
4.25 km/s re-derivados). Sigue faltando lo author-owned: autores/afiliaciones/ORCID,
acknowledgements, y ejecutar el depósito para obtener el DOI.

## 3. Gotchas del entorno (LEER antes de operar — evitan errores ya cometidos)

- **Docker only**: todo corre en `docker compose run --rm pipeline python -m ...`.
  No hay Python local. `./src` y `./scripts` montados como volumen (sin rebuild al editar).
- **`docs/` NO está montado** en el contenedor (solo `data`, `logs`, `src`, `scripts`).
  Scripts que generan archivos para `docs/` deben escribir a `data/output/` y moverse
  con `mv` en el host.
- **matplotlib NO está en la imagen** (se instala ad-hoc por corrida). Cualquier módulo
  en `scripts/` debe importar matplotlib **lazy** (dentro de funciones, no a nivel de
  módulo) o rompe `tests/test_script_entrypoints.py::test_script_module_imports`.
- **Protección de `main` es ESTRICTA**: `strict:true` (branch debe estar up-to-date) +
  `enforce_admins:true` (ni `--admin` saltea) + checks requeridos `Tests` y
  `Docker Compose Tests`. Consecuencia: **serializar merges** (update-branch → esperar
  CI → merge, uno por uno). `Format (black)` y `Lint (ruff)` corren pero NO son
  requeridos (igual conviene dejarlos verdes: `docker compose run --rm pipeline black .`
  y `ruff check . --fix`).
- **Basar PRs en `main`, no en otra branch de PR abierta.** Mergear con `--delete-branch`
  la base de PRs stackeadas las **cierra** (no las reapunta). Error cometido esta sesión
  con #85–#88 (recreadas como #89–#92).
- **Git**: autor SIEMPRE `Daniel Frascarelli <dsanfra@gmail.com>`; **sin `Co-Authored-By`
  ni footer de agente** en commits ni PRs (ver `settings.local.json` attribution="").
  Nunca trabajar en `main` directo: branch descriptiva → PR → merge. No `--no-verify`.
- **Comandos Bash auto-aprobados** en esta sesión (`settings.local.json` allow: `"Bash"`).
- **TAP de Gaia** tira HTTP 500 transitorios; el motor ya tiene retry/backoff. Reintentar
  resuelve. No hay cache de tránsitos: cada corrida los baja de nuevo.
- **Memoria**: máquina con ~4 GB libres → OOM sobre 72M filas. Usar polars lazy +
  `engine="streaming"` y subsamplear scatters. Fits de masa: repartir `--workers`.

## 4. Archivos clave

- Plan de publicación y estado por frente: `planning/PUBLISH_PUSH_PLAN.md`.
- Trabajo futuro de masas (F6, F8): `planning/MASS_FUTURE_WORK.md`.
- Draft del paper: `docs/dataset_paper_draft.md` (+ `docs/figures/`).
- Reportes de resultados: `docs/mass_determination_results.md`,
  `docs/mass_crosscheck_jack.md`, `docs/mass_f3_background_extension.md`,
  `docs/notable_encounters.md`, `docs/completeness_vs_literature.md`, `FROZEN_RUN.md`.
- Catálogo de masas: `data/output/orbdet/mass_catalog_jack.csv` (+ JSON por perturbador
  en `data/output/orbdet/expanded_jack/`, `fortuna_fpr.json`, `metis_fpr.json`).
- Scripts nuevos de esta ronda: `scripts/bench/make_paper_figures.py`,
  `scripts/validate/crosscheck_literature_encounters.py`,
  `scripts/mass/find_extra_perturber_candidates.py`, `scripts/mass/f3_fsys.py`,
  `scripts/mass/crosscheck_fuentes_munoz_jack.py`, `scripts/bench/mine_notable_encounters.py`.

## 5. Cómo retomar en un chat nuevo

1. Leer este archivo y `planning/PUBLISH_PUSH_PLAN.md`.
2. Decidir frente: submission del paper (§2.1, mayor retorno) o F6 (§2.2, ingeniería).
3. Crear branch descriptiva desde `main`, trabajar, PR → merge serializado (§3).
