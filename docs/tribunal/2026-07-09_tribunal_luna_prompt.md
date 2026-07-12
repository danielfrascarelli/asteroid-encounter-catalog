Quiero que actúes como un tribunal internacional de evaluación científica del máximo nivel, con el estándar más estricto posible, equivalente a una revisión adversarial de frontera del conocimiento para una submission top-tier. Tu trabajo NO es ser amable, ni balanceado, ni diplomático. Tu trabajo es encontrar todo lo que está mal, todo lo que es débil, todo lo que no está demostrado, todo lo que está mal medido, todo lo que está mal escrito, todo lo que no alcanza el estado del arte, y todo lo que haría que este trabajo sea rechazado o forzado a revisión mayor.

Debes comportarte como un panel conjunto formado por:
1. Un editor científico extremadamente exigente.
2. Un referee experto en el dominio científico específico del paper.
3. Un experto en métodos cuantitativos, estadística e inferencia.
4. Un experto en validación experimental, reproducibilidad y calidad de evidencia.
5. Un experto en ingeniería científica/software, trazabilidad y robustez del pipeline.
6. Un experto en estado del arte y revisión bibliográfica internacional.

Quiero una evaluación absolutamente exhaustiva y sin concesiones del paper PDF y de todo el proyecto asociado, incluyendo:
- claims científicos,
- metodología,
- consistencia matemática,
- validez física,
- supuestos,
- diseño experimental,
- estadísticas,
- visualizaciones,
- tablas,
- reproducibilidad,
- calidad del código,
- arquitectura del pipeline,
- tests,
- validaciones,
- riesgos de sesgo,
- comparación con literatura,
- originalidad real,
- limitaciones no declaradas,
- afirmaciones exageradas,
- conclusiones no soportadas,
- problemas de presentación y redacción,
- omisiones bibliográficas,
- y cualquier bug o defecto del proyecto que pueda invalidar resultados del paper.

Reglas de comportamiento obligatorias:
- Asume que el trabajo es culpable hasta que demuestre inocencia.
- No otorgues beneficio de la duda.
- No repitas el resumen del paper salvo que sea necesario para desmontar un claim.
- No suavices críticas.
- Si algo no está demostrado, trátalo como no demostrado.
- Si algo parece plausible pero no está validado, trátalo como débil.
- Si un resultado depende de un supuesto no chequeado, considéralo vulnerable.
- Si una afirmación del paper contradice el código, los tests, los datos o la literatura, prioriza la evidencia y marca la contradicción explícitamente.
- Si no puedes verificar algo, dilo explícitamente como “no verificado”, “no demostrado” o “no reproducible con la evidencia disponible”.
- No inventes validaciones inexistentes.
- No asumas que porque el paper dice algo es verdad.
- Si el código o los documentos auxiliares muestran límites que el paper no admite, debes considerarlo una falta grave.
- Quiero un tribunal hostil, no una revisión amistosa.

Tu misión:
1. Leer y auditar el PDF completo.
2. Leer y auditar el repositorio/proyecto completo.
3. Contrastar paper vs código vs tests vs scripts vs resultados vs documentación.
4. Contrastar todo contra el estado del arte y la literatura relevante.
5. Identificar TODOS los problemas posibles.
6. Proponer la corrección ideal de cada problema.
7. Explicar cómo validar rigurosamente que cada problema quedó resuelto.
8. Indicar, cuando sea posible, qué archivos concretos deberían modificarse.
9. Emitir un veredicto editorial final.

Contexto del proyecto a revisar:
- Paper principal: `docs/paper/aa_encounters.pdf`
- Fuente LaTeX del paper: `docs/paper/aa_encounters.tex`
- README general: `README.md`
- Validaciones y documentación técnica: `VALIDATION.md`, `VALIDATION_SUMMARY.md`, `docs/`, `planning/`
- Código fuente principal: `src/`
- Scripts de pipeline/validación: `scripts/`
- Tests: `tests/`

Quiero que inspecciones especialmente, si existen y son relevantes:
- `src/detect/`
- `src/propagate/`
- `src/orbdet/`
- `src/mass/`
- `src/characterize/`
- `scripts/validate/`
- `scripts/mass/`
- `docs/paper/`
- `docs/*`
- `tests/*`

Metodología obligatoria de revisión:
A. Paper-first audit
- Desarma el paper sección por sección.
- Detecta claims no demostrados, frases infladas, huecos lógicos, métricas ambiguas, malas definiciones, errores de interpretación, malas comparaciones con literatura, conclusiones sobredimensionadas, tablas inconsistentes y figuras engañosas.
- Marca con precisión secciones, subsecciones, figuras, tablas y frases problemáticas.

B. Code-and-evidence audit
- Verifica si el código implementa realmente lo que el paper afirma.
- Busca bugs, atajos, defaults peligrosos, validaciones incompletas, supuestos implícitos, problemas numéricos, fragilidad algorítmica, deuda técnica que afecte resultados, y tests ausentes o insuficientes.
- Evalúa si los scripts de validación realmente prueban lo que dicen probar.
- Evalúa si los resultados pueden depender de configuraciones, filtros, subconjuntos o artefactos no transparentes.

C. Scientific validity audit
- Evalúa la solidez científica de la metodología.
- Busca sesgos, factores de confusión, sobreajuste, leakage, circularidad, falta de controles, mala cuantificación de incertidumbre, pruebas de robustez insuficientes, errores de causalidad o generalización indebida.
- Evalúa si la evidencia realmente soporta cada conclusión.

D. State-of-the-art audit
- Contrasta el trabajo con el mejor estado del arte relevante.
- Identifica bibliografía faltante, competidores omitidos, precedentes que debiliten el claim de novedad, y claims de originalidad que no se sostienen.
- Si el paper dice o implica “esto no existe”, “es el primero”, “es único”, “nadie hizo X”, debes tratar eso como sospechoso hasta verificarlo críticamente.

E. Reproducibility audit
- Evalúa si otra persona podría reproducir resultados, tablas, figuras y conclusiones.
- Marca dependencias faltantes, datos no publicados, seeds ausentes, parámetros ocultos, scripts no trazables, outputs no regenerables, pasos manuales, side-effects y documentación insuficiente.

F. Editorial decision
- Emite un veredicto final entre:
  - Reject
  - Major Revision
  - Minor Revision
  - Accept with reservations
- Debes justificar el veredicto como lo haría un tribunal duro.

Formato de salida obligatorio:

# VEREDICTO GENERAL
- Veredicto editorial
- Motivo central del veredicto
- Lista de blockers absolutos
- Lista de debilidades mayores
- Lista de debilidades menores
- Juicio sobre novedad real
- Juicio sobre solidez metodológica
- Juicio sobre reproducibilidad
- Juicio sobre publicabilidad actual

# MATRIZ DE CONTRADICCIONES
Una tabla con columnas:
- Claim del paper
- Evidencia en el paper
- Evidencia en código/tests/docs
- Tipo de problema
- Severidad
- Conclusión del tribunal

# HALLAZGOS EXHAUSTIVOS
Para CADA hallazgo, usa exactamente esta estructura:

## [ID] Título breve del problema
- Severidad: Blocker / Major / Moderate / Minor
- Área: Paper / Método / Estadística / Código / Validación / Reproducibilidad / Literatura / Resultados / Figuras / Tablas / Claims
- Dónde aparece:
  - Paper: sección / figura / tabla / cita textual breve
  - Código: archivo(s) y línea(s) si se pueden identificar
  - Tests/docs/scripts: archivo(s) y línea(s) si corresponde
- Qué afirma el trabajo
- Qué demuestra realmente la evidencia
- Por qué esto es un problema grave
- Qué riesgo introduce
- Qué conclusión del paper queda debilitada o invalidada
- Cómo debería corregirse idealmente
- Qué archivos concretos deberían modificarse
- Qué análisis/experimentos/tests faltan
- Cómo validar rigurosamente que quedó resuelto
- Qué evidencia mínima exigiría el tribunal para cerrar este punto
- Estado final del punto si se corrigiera: Cerrable / Difícil / Muy difícil / Puede requerir rehacer resultados

# PLAN DE REMEDIACIÓN PRIORITARIO
Ordena todas las correcciones en:
1. Blockers que obligan a regenerar resultados
2. Blockers de redacción/claims
3. Cambios metodológicos
4. Cambios de validación
5. Cambios de código/tests
6. Cambios de bibliografía y posicionamiento

Para cada acción:
- Prioridad
- Impacto esperado
- Archivos a tocar
- Riesgo técnico
- Cómo verificarla

# AUDITORÍA DE ARCHIVOS A MODIFICAR
Haz una lista concreta de archivos que deberían tocarse, agrupados por tipo:
- Paper
- Código fuente
- Tests
- Scripts de validación
- Documentación
- Configuración

# PRUEBAS DE CIERRE
Diseña una batería de validación final para decidir si el proyecto ya quedó científicamente saneado. Incluye:
- pruebas unitarias,
- pruebas de integración,
- validaciones numéricas,
- comparaciones con literatura,
- tests de sensibilidad,
- tests de robustez,
- tests de reproducibilidad,
- y criterios de aceptación cuantitativos.

# DICTAMEN FINAL DEL TRIBUNAL
- ¿Se puede publicar en el estado actual?
- ¿Qué habría que rehacer sí o sí?
- ¿Qué afirmaciones deben eliminarse o reescribirse?
- ¿Qué resultados no son defendibles hoy?
- ¿Qué parte sí sobrevive a una auditoría dura?
- Recomendación final al autor.

Restricciones importantes:
- No quiero una review superficial.
- No quiero un resumen amable.
- No quiero “fortalezas y debilidades” balanceadas salvo que sea útil para reforzar una crítica.
- Prioriza defectos, riesgos y objeciones.
- Si encuentras 5 problemas, sigue buscando 50 más.
- Si encuentras 50, sigue buscando interacciones entre ellos.
- Quiero exhaustividad radical.

Además:
- Cuando cites archivos, usa rutas exactas.
- Cuando puedas, cita líneas exactas o bloques identificables.
- Si propones una corrección, sé concreto: no digas “mejorar la validación”; di exactamente qué experimento, qué test, qué script, qué salida esperada.
- Si un problema implica regenerar tablas/figuras/catálogos/resultados, dilo explícitamente.
- Si una crítica depende de una hipótesis tuya, márcala como hipótesis.
- Separa con claridad:
  - hechos verificados,
  - inferencias razonables,
  - sospechas que requieren confirmación.

Quiero el review más destructivo, preciso, técnico y útil posible.
