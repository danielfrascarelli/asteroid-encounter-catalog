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
El draft está en `docs/dataset_paper_draft.md`. Para someterlo:
- [ ] **Verificar referencias contra ADS**: volúmenes/páginas de Fuentes-Muñoz 2025
  (AJ 170, 353 — confirmar), Vernazza 2020 (Nat. Astron. 4, 136), Park 2016,
  Russell 2012, Tanga 2023, Goffin 2014, Michalak 2000. Sección "References" del draft.
- [ ] **Elegir journal**: A&A (encaja con Tanga/Goffin) o PSJ. Define template y formato.
- [ ] **Revisión de consistencia de números** entre secciones (72,236,904 encuentros;
  0.70 % censura; 76 % recall; 29.5 % recuperación FM; masas §5). Que ningún número
  se contradiga entre abstract, cuerpo y conclusiones.
- [ ] **Quitar notas internas** (el bloque `<!-- Notas de redacción -->` y la sección
  "Figures" interna) al pasar a template de journal.
- [ ] **Autores, afiliaciones, acknowledgements, data availability final** (DOIs/Zenodo
  para el catálogo si se publica el dato).
- [ ] Cosmético: fig3 tiene la etiqueta "Inclination i (deg)" duplicada (colorbar +
  marginal). Editar `scripts/bench/make_paper_figures.py` (`figure3_aei_map`) y regenerar.

**Regenerar figuras** (matplotlib NO está en la imagen — instalar ad-hoc):
```bash
docker compose run --rm pipeline bash -c \
  "pip install --quiet matplotlib && python -m scripts.bench.make_paper_figures"
mv data/output/figures/*.png data/output/figures/*.pdf docs/figures/   # docs/ no está montado
```

### 2.2 F6 — parcial ∂x/∂GM analítica (baja prioridad, única ingeniería sustancial que queda)
Hoy las parciales se calculan por diferencias finitas centrales sobre `propagate_assist`.
- **Acción**: implementar partícula variacional respecto a la masa (o esquema adjunto)
  sobre ASSIST, reduciendo nº de propagaciones por Jacobiano.
- **Gate**: parcial analítica vs FD < 1e-6 relativo; reducción medible de tiempo.
- Detalle en `planning/MASS_FUTURE_WORK.md` §F6. Toca `src/orbdet/`.

### 2.3 F8 — Pallas con Gaia DR4 (bloqueado: espera datos)
Pallas tiene solo 6–7 objetivos < 0.05 AU en DR3 (target-limited). Cuando salga DR4:
```bash
docker compose run --rm pipeline python -m scripts.mass.orbdet_fit_realdata \
    --perturber 2 --release dr4 --from-catalog <catálogo DR4> \
    --top-per-perturber 40 --workers 24 --jackknife --out data/output/orbdet/pallas_dr4.json
```
Gate: N(Pallas) ≥ 20. No accionable hasta que DR4 esté disponible.

### 2.4 Cruce Goffin pair-level (bloqueado: extracción de PDF)
El cruce de completitud (§3.4) usó FM 2025. Goffin 2014 a nivel de pares no se pudo:
los parquets VizieR solo traen masas por perturbador, no los pares, y el contenedor
no tiene extractor de PDF. Si se quiere: agregar `pdfplumber`/`PyMuPDF` y extraer las
tablas de pares de `papers/goffin_2014_aa565_A56.pdf`, luego extender
`scripts/validate/crosscheck_literature_encounters.py`.

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
