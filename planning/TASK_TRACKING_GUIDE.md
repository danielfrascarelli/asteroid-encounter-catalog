# Guía: cómo escribir un MD de tracking de tareas

> Convención del proyecto para planes y trackers de progreso. Todo MD de
> planificación o seguimiento de tareas **vive en `planning/`** y sigue esta
> estructura. Los `.md` de la raíz son charter/referencia/resultados, no
> trackers (ver abajo).

Esta guía existe porque el proyecto acumuló ~7 planes sueltos en la raíz, varios
ya terminados y nunca marcados como tales (p. ej. `FPR_INGEST_PLAN.md` quedó
leyéndose como "pendiente" cuando sus 5 etapas estaban hechas). Un tracker que no
declara su estado de forma inequívoca es deuda: el próximo lector —humano o
agente— no sabe si la tarea está viva, hecha o muerta.

---

## Reglas

### 1. Cabecera con estado inequívoco
Las primeras líneas deben responder, sin leer el resto: **qué es, para qué, en
qué estado, y desde cuándo**.

```markdown
# <título>

> **Estado:** 🟡 ACTIVO | ⬜ PENDIENTE | ✅ COMPLETO | ⚠️ SUPERADO
> **Última actualización:** 2026-06-01
> Plan de <una frase: qué resuelve y cuál es el criterio de éxito>.
```

Estados:
- **PENDIENTE** — planeado, sin arrancar.
- **ACTIVO** — trabajo en curso.
- **COMPLETO** — todos los ítems cerrados. **Marcarlo apenas cierra** (no dejarlo
  leyéndose como activo).
- **SUPERADO** — ya no es válido; decir por qué y enlazar a lo que lo reemplaza.

### 2. Tabla de estado arriba de todo
Un vistazo = el estado de cada ítem. Una fila por tarea:

```markdown
| # | tarea | estado | entregable / PR |
|---|-------|--------|-----------------|
| 1 | <qué> | ✅ | `path/al/archivo` · #42 |
| 2 | <qué> | 🟡 | — |
| 3 | <qué> | ⬜ | — |
```

### 3. Cada tarea define su "Done" (gate de aceptación)
No "implementar X" sino **cómo se sabe que X está bien**. Un criterio
verificable: un test que pasa, un número bajo un umbral, un artefacto que existe.

```markdown
### Tarea 3 — <título>
**Entregable:** <archivo/función/doc concreto>.
**Gate:** <condición verificable — "test Y pasa", "error ≤ Z", "recupera masa
inyectada ratio≈1.0">.
**Depende de:** Tarea 2.
```

### 4. Fechas absolutas, nunca relativas
"la semana que viene" envejece mal. Escribir `2026-06-08`.

### 5. Al cerrar un ítem, enlazar la evidencia
PR, commit, doc de resultados o test. El tracker apunta a la verdad, no la
duplica.

### 6. No reescribir la historia: marcar, no borrar
Si una decisión cambia, marcá el ítem viejo como SUPERADO con el motivo y el
enlace al nuevo. No edites en silencio (perdés el "por qué").

### 7. Una sola fuente de verdad
Si dos docs trackean lo mismo, uno manda y el otro enlaza. Evitar estados
divergentes.

### 8. Declarar lo que NO está en scope
Una sección "Fuera de scope" evita que el plan se infle y deja claro el borde.

### 9. Al terminar, cerrar el ciclo
Un plan COMPLETO se **borra** o se archiva — no se deja en la raíz pretendiendo
estar vivo. La info que sobrevive migra a `docs/` (resultados), `CHANGELOG.md` o
`ROADMAP.md` (referencia). Esta es la razón de ser de esta limpieza.

---

## Qué va en `planning/` y qué no

| va en `planning/` | NO (queda en raíz / `docs/`) |
|-------------------|------------------------------|
| Planes de acciones a implementar | Charter (`CLAUDE.md`), README, CONTRIBUTING |
| Trackers de progreso activos | Resultados/validación (`docs/*`, `VALIDATION*.md`) |
| Roadmaps de trabajo futuro | Referencia del catálogo (`FROZEN_RUN.md`) |
| Esta guía | Changelog, roadmap maestro histórico (`ROADMAP.md`) |

Regla mental: **si describe trabajo a hacer o en curso → `planning/`. Si
describe lo que el sistema es o lo que ya se midió → raíz o `docs/`.**
