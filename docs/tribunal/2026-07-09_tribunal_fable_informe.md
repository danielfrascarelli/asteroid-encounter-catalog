Informe del tribunal — auditoría fresca del manuscrito aa_encounters (2026-07-09)
Estado de los materiales
aa_encounters.tex — leído completo (1189 líneas). Compila mentalmente; todas las \cite resuelven contra el .bib y todas las \ref internas resuelven, con una excepción semántica (ver M1: se cita un "fix" que la sección referenciada nunca describe). Contiene 4 bloques TODO en comentarios que afectan números publicados (ver B3, B4).
references.bib — leído completo (19 entradas, todas citadas en el texto, sin huérfanas).
aa_encounters.pdf — 12 páginas, compilado hoy 16:16. Sin referencias rotas (??) en el texto extraído; figuras 1–4 presentes, renderizadas y consistentes con sus captions.
Acceso a búsqueda web: SÍ. Verifiqué existencia/volumen/página de fuentesmunoz2025, li2023, ivantsov2018, siltala2020→2022, y el diámetro de (44) Nysa. El resto de la verificación bibliográfica es "según mi conocimiento, sin verificación externa" y lo marco así donde aplica.
HALLAZGOS BLOQUEANTES
[B1] BLOQUEANTE — Diámetros de la Tabla 6 (candidate perturbers) contradicen valores medidos y producen densidades no físicas con las propias masas FM25 tabuladas
Panel: Estado del arte + Física (verificación cruzada), Árbitro de manuscrito (inconsistencia interna).
Ubicación: Tabla 6 (tab:candperturbers), aa_encounters.tex:771-795; texto de Sect. 4.4.

Descripción del problema:
La tabla dice rankear "large bodies ($D\gtrsim100$ km)", pero sus diámetros son incompatibles con los valores medidos de asteroides extremadamente bien caracterizados, y la propia tabla trae la prueba: combinando el $D$ tabulado con la masa FM25 de la última columna, las densidades implicadas son no físicas.

(44) Nysa: $D=140$ km. El valor medido es ~70–79 km (radar Arecibo 79±10 km; ~71 km con albedo 0.48 — verificado por web). Con $M_{\rm FM25}=5.8\times10^{17}$ kg y $D=140$ km, $\rho = 0.40$ g/cm³ — absurdo para un E-type; con $D=71$ km da 3.1 g/cm³, correcto.
(19) Fortuna: $D=133$ km. Medido ~200–225 km (IRAS/ocultaciones). Con $M=8.5\times10^{18}$ kg y $D=133$: $\rho=6.9$ g/cm³ — absurdo para un C-type; con $D\approx210$ da 1.7 g/cm³, correcto.
(21) Lutetia: $D=120$ km. Visitada por Rosetta: $D_{\rm eff}\approx98$ km, $\rho=3.4$ g/cm³. Con el $D$ tabulado y $M=1.7\times10^{18}$ kg: $\rho=1.9$; con el $D$ de Rosetta: 3.4 — la masa FM25 reproduce exactamente la densidad Rosetta solo con el diámetro correcto.
(64) Angelina: $D=104$ km vs. medido ~50–60 km (E-type). Ni siquiera cumple el corte $D\gtrsim100$ km de la propia tabla. Lo mismo Nysa.
Patrón sistemático: los E-types (albedo ~0.5) aparecen ~2× demasiado grandes y los C-types oscuros demasiado chicos — exactamente la firma de $D$ derivado de $H$ con albedo default 0.14 ($D\propto p_V^{-1/2}$: $\sqrt{0.55/0.14}\approx2.0$), pese a que estos diez cuerpos tienen diámetros medidos en SBDB. Además, 7 de las 10 filas tienen clase dinámica "—" (¡para (20) Massalia, (21) Lutetia, (29) Amphitrite!), y (9) Metis figura con $D=197$ km aquí pero $D_1=190$ km en la Tabla 5 — dos valores para el mismo cuerpo en el mismo paper.
Esto invalida la selección de la Tabla 6 (el conjunto "$D\gtrsim100$ medido" está mal poblado), contradice la remediación previa de diámetros, y —dado que el pipeline de caracterización es el mismo— pone bajo sospecha la columna diameter_* y class_* del catálogo caracterizado de 80 M filas para algún subconjunto de cuerpos. Hipótesis de implementación asociada (requiere verificación contra el código): el join con el snapshot físico de SBDB falla para estas filas (candidato: src/characterize/physical.py), y la Tabla 6 se generó por un camino de datos distinto al de las Tablas 1–2 (que sí muestran diámetros correctos: Ceres 939, Thisbe 232, Hermione 209).

Solución exigida:
(1) Reparar el join de datos físicos y regenerar la Tabla 6 con diámetros medidos (IRAS/AKARI/NEOWISE/ocultaciones vía SBDB) y clase dinámica poblada; re-derivar el ranking (el orden puede cambiar: Fortuna con $D\approx210$ probablemente sube). (2) Reconciliar el $D$ de Metis entre Tablas 5 y 6 (una sola fuente). (3) Auditar en el catálogo caracterizado qué fracción de cuerpos con diámetro medido en SBDB quedó con diameter_source = default_albedo, y reportarla. (4) Verificar si el corte "$D\gtrsim100$ km" de la Sección 4.4 sigue dando 10 candidatos.

Criterio de validación de cierre:
Script que, para cada fila de la Tabla 6, calcule $\rho$ implicada por ($D$ tabulado, $M_{\rm FM25}$) y verifique $0.8 < \rho < 4.5$ g/cm³; y assert de que todo cuerpo de las tablas del paper con entrada de diámetro en el snapshot SBDB tiene diameter_source ∈ {measured, albedo_measured}. Cero filas con clase "—" entre los numerados < 1000.

[B2] BLOQUEANTE — Dos encuentros de la Tabla 1 están fechados FUERA de la ventana de observación declarada
Panel: Física/astrodinámica + Árbitro de manuscrito.
Ubicación: Tabla 1 (tab:dd100), filas 7 y 12: aa_encounters.tex:602 y aa_encounters.tex:607; ventana declarada en Sect. 2.1 (aa_encounters.tex:174).

Descripción del problema:
La ventana de propagación es "2014 July 25 – 2017 May 28". Sin embargo, (196) Philomela × (431) Nephele y (85) Io × (145) Adeona están fechados 2017-05-29 — un día después del fin de ventana. Dos filas independientes con exactamente la misma fecha, coincidente con el borde de la grilla, es la firma clásica de un mínimo clampeado en el último paso temporal (el mínimo verdadero cae después del fin de ventana y el pipeline reporta el valor de borde) — precisamente el "window-boundary truncation artefact" que el propio paper declara como modo de fallo eliminado (aa_encounters.tex:656-658). Nótese la asimetría editorial: las Tablas 3–5 garantizan ">130 d from either edge", pero las Tablas 1–2 —el escaparate del catálogo— no tienen ningún guard de borde. Si estas dos filas son clamps, sus distancias (0.0389 y 0.0470 AU) no son mínimos de encuentro reales y la Tabla 1 pasa de 12 a 10 entradas, con texto y abstract ("twelve encounters") a corregir.

Solución exigida:
(1) Verificar en el catálogo si jd_tdb de estas dos filas coincide con el último nodo de la grilla (candidatos: src/detect/refine.py, src/detect/pipeline.py — hipótesis de ubicación). (2) Si son clamps: excluirlas o marcarlas como censura de borde, regenerar Tablas 1–2, y aplicar a todas las tablas de eventos el mismo criterio de distancia al borde (o declarar explícitamente cuáles no lo aplican y por qué). (3) Si son genuinas (mínimo local interior el 28 de mayo que formatea a 29 por conversión TDB→UTC): documentar la convención de cierre de ventana (¿inclusive hasta qué instante, en qué escala?) — tal como está escrito, un lector no puede reconciliar la fecha con la ventana.

Criterio de validación de cierre:
Query sobre el catálogo publicado: cero filas de cualquier tabla del paper con jd_tdb fuera de $[\rm JD_{ini}, JD_{fin}]$ de la provenance, y assert de que ninguna fila tabulada tiene su mínimo en el primer/último nodo de la grilla; un flag edge_censored por fila si se decide retenerlas.

[B3] BLOQUEANTE — El manuscrito mezcla números de dos freezes distintos del catálogo (pre-B1 y B1-fixed), con TODOs internos que lo admiten
Panel: Algoritmos/completitud + Árbitro de manuscrito.
Ubicación: Sect. 2.2 (aa_encounters.tex:236-237 y TODO en aa_encounters.tex:259-261); Sect. 3.1–3.2; Sect. 3.4 (TODO en aa_encounters.tex:539-540).

Descripción del problema:
Hallazgo de texto verificable aritméticamente: la Sect. 2.2 dice que el subconjunto frágil son "8 728 509 pairs (12.08 % of the catalogue)", pero $8,728,509 / 80,072,774 = \mathbf{10.90,%}$. El 12.08 % solo cierra con un denominador de ~72.3 M filas — es decir, el conteo del subconjunto frágil proviene del catálogo pre-B1, mientras el denominador del paper es el catálogo B1-fixed. El comentario TODO(Stage B / b1fix) lo confirma: el híbrido "is being regenerated... will be refreshed... before submission". En cascada quedan bajo sospecha de freeze viejo: los 25 283 crossings upward (0.29 %) de Sect. 3.2, todo Stage B de Sect. 3.1 (p99 1.99 mAU, max 15.2 mAU), y el cross-match de Sect. 3.4 (cuyo TODO habla de "25 962 pares no recuperados" cuando el texto rendido implica $40,176-11,842 = 28,334$ — dos números incompatibles para la misma cantidad). El presupuesto de completitud —el claim distintivo del paper— está parcialmente calculado sobre un catálogo que ya no es el que se publica.

Solución exigida:
Completar la regeneración Stage B sobre el freeze B1-fixed y refrescar: conteo y porcentaje del subconjunto frágil, upward crossings, estadísticos Stage B, la Fig. 4 si cambia, y todo Sect. 3.4 (40 176 / 11 842 / 91.6 % / 8.4 %). Eliminar los dos TODOs. Verificar que Abstract y Conclusiones (que citan 0.70 %, 12 µAU, 2.5 mAU) sigan siendo los números del freeze final.

Criterio de validación de cierre:
grep -c "TODO" aa_encounters.tex = solo los administrativos de B4 (idealmente 0); y un check aritmético automatizable: todo porcentaje del paper con numerador/denominador tabulados reproduce al redondeo declarado (en particular frágil/total = el % impreso).

[B4] BLOQUEANTE (administrativo) — Metadatos de submission ausentes: afiliación, y DOI del archivo de datos
Panel: Árbitro de manuscrito.
Ubicación: aa_encounters.tex:30-38 (autor/afiliación); aa_encounters.tex:1076-1077 (TODO DOI); acknowledgements aa_encounters.tex:1111.

Descripción del problema:
"Affiliation --- TODO", DOI de Zenodo/VizieR pendiente, funding pendiente. Para un paper cuyo producto es el catálogo, la Sect. 6 sin identificador persistente es un rechazo de escritorio en A&A: el claim "published in full" es inverificable hasta que el depósito exista. Trivial pero bloqueante en sentido literal.

Solución exigida:
Depositar los tres parquet + sidecars en Zenodo/VizieR, insertar DOI en Sect. 6, completar afiliación/ORCID/funding.

Criterio de validación de cierre:
El DOI resuelve públicamente y los SHA-256 listados en el registro coinciden con los archivos citados en el texto.

HALLAZGOS MAYORES
[M1] MAYOR — La ventana de refinamiento descrita (±2 h) no puede bracketear un mínimo con cadencia gruesa de 12 h; el "refinement-window fix" se cita pero nunca se describe
Panel: Algoritmos/completitud.
Ubicación: Sect. 2.2 (aa_encounters.tex:223-225); citas al "fix" en aa_encounters.tex:583 y aa_encounters.tex:658.

Descripción del problema:
Hallazgo de texto: el método descrito dice "a $\pm2$ h window around the apparent minimum is sampled at 120 s" con cadencia gruesa $\Delta t = 12$ h. Bajo fase uniforme, el mínimo verdadero cae hasta $\Delta t/2 = 6$ h del sample más cercano; una ventana de ±2 h lo pierde en 2/3 de los casos. Es matemáticamente la descripción del método roto (el bug que sesgaba ~60 % del catálogo — 2/3 ≈ 67 %). Los §4.1/4.2 citan "the refinement-window fix (Sect. 2)", pero la Sect. 2 no contiene ninguna descripción de fix: la referencia interna cuelga y el lector queda con un método que, tal como está especificado, no puede haber producido el catálogo que el paper valida. Problema adicional del mismo tipo: el refinamiento N-body del híbrido usa "±12 h around the Kepler minimum" sin declarar qué pasa si el mínimo N-body cae en el borde de esa ventana (¿re-centrado? ¿flag?). Hipótesis de implementación: el código real usa ventana ≥ ±6 h o re-centrado iterativo (candidato: src/detect/refine.py) — no puedo confirmarlo con el texto.

Solución exigida:
Reescribir el párrafo de refinamiento describiendo el algoritmo realmente ejecutado en el freeze B1-fixed — p. ej.: "the refinement window spans $\pm\Delta t/2$ around the bracketing coarse sample" o "the window is iteratively re-centred until the minimum is interior". Declarar el guard de borde de la ventana N-body de ±12 h. Eliminar las menciones a "the refinement-window fix (Sect. 2)" o hacer que la Sect. 2 lo describa de verdad (una frase sobre el modo de fallo pre-fix bastaría, ya que el paper lo usa como evidencia en §4.1).

Criterio de validación de cierre:
El test de inyección-recuperación ya citado (aa_encounters.tex:516-519) ejecutado con mínimos colocados uniformemente en fase dentro del paso de 12 h, mostrando ≥99 % de recuperación con la ventana tal como quedó descrita en el texto (mismos parámetros), y cero mínimos reportados en el borde de la ventana de refinamiento en el catálogo final.

[M2] MAYOR — El Monte Carlo de incertidumbre de elementos usa covarianzas SBDB de HOY para distancias calculadas con elementos de 2016: subestima el error y mezcla soluciones
Panel: Estadística matemática + Física.
Ubicación: Sect. 3 (aa_encounters.tex:376-389); Tabla 3 y su caption; Sect. 4.2, 4.3.

Descripción del problema:
Las distancias tabuladas se calcularon con elementos del snapshot MPCORB 2016-02-17. El MC dibuja 3000 órbitas "from its full JPL SBDB covariance" — la covarianza de la solución orbital actual (2026), que incorpora ~10 años más de astrometría (incluida la propia Gaia) que los elementos usados. El σ relevante para "¿cuán lejos de la verdad puede estar el número impreso?" es el de la solución de 2016, no el de 2026. El sesgo es máximo exactamente donde más importa: los cuerpos de la Tabla 3 son de número alto — (238587), (391704), (393772) — recién numerados en la época del snapshot, con arcos cortos entonces; su incertidumbre 2016-era puede ser 1–2 órdenes de magnitud mayor que la actual, lo que puede llevar σ_d de 32–83 km a la escala de la separación misma (~10³ km) y demoler el claim "significant to ~2–8 %". Defecto metodológico adicional: se combina el nominal de MPCORB-2016 con la covarianza de SBDB-2026 — nominal y covarianza de soluciones distintas; el offset nominal MPCORB↔SBDB ni se menciona ni se acota. El descargo del texto ("formal covariances... optimistic by a factor of a few") no cubre este efecto, que no es de formalidad sino de época.

Solución exigida:
Elegir una de dos vías y declararla: (a) evaluar la significancia con la mejor solución actual — pero entonces recomputar también la distancia nominal de los eventos de Tabla 3 con esa solución (nominal y covarianza consistentes), reportando además el corrimiento nominal MPCORB-2016→SBDB-2026 por evento como término empírico de error; o (b) mantener el nominal 2016 y estimar la covarianza de época 2016 (refit con astrometría truncada a 2016, o escala aproximada por longitud de arco). Reescribir el pasaje "1–83 km... not the $10^3$–$10^4$ km one would fear" según el resultado — puede sobrevivir, pero hoy no está demostrado.

Criterio de validación de cierre:
Tabla 3 regenerada con dos columnas: Δ nominal (MPCORB-2016 vs solución actual) y σ_d consistente en época; el claim de significancia recalculado. Gate: para cada fila, |Δ nominal| reportado y σ_d re-derivado; si σ_d/d > 30 % para alguna fila, el texto lo dice explícitamente en vez de "2–8 %".

[M3] MAYOR — (52) Europa: "measured" con z = −4.3 contradice el "broadly consistent" del abstract; la narrativa de determinabilidad a-priori es refutada por la propia Tabla 8; y falta la determinación independiente publicada de Europa
Panel: Estadística matemática + Estado del arte.
Ubicación: Tabla 8 (aa_encounters.tex:934); abstract (results); parágrafo "A-priori determinability" (aa_encounters.tex:841-858); parágrafo "External σ" (aa_encounters.tex:955-960).

Descripción del problema:
Tres capas del mismo problema, todas verificables en el texto:

Europa tiene status measured con ratio 0.58 y $z=-4.3$ contra su referencia DE441. Un desvío de 4.3σ en una de las seis masas "identificables" es incompatible con el abstract ("these are broadly consistent with recent independent work") y nunca se discute — el texto analiza Thisbe y Juno en detalle pero calla el caso más discrepante que él mismo flaggea como medición.
El claim "perturbers whose true mass is $\gtrsim2\times10^{19}$ kg are recovered near unity ratio" es refutado por la propia tabla: Europa (masa verdadera $2.3/0.58\approx4\times10^{19}$ kg, ratio 0.58) e Interamnia ($\approx4.3\times10^{19}$ kg, ratio 0.47) están bien por encima del umbral y lejos de la unidad. La frontera de determinabilidad "coincides with what the engine delivers" no sobrevive su propia evidencia.
El test que el texto usa para calibrar la coherencia del pull (sign test, $p\approx0.25$) es el test de mínima potencia disponible; el z individual de Europa (−4.3) es significativo por sí solo. Elegir el test débil y reportar "suggestive but not significant" es exactamente el tipo de argumento que un árbitro de estadística marca. Además, "geometric-mean ratio of 0.85" no es reproducible desde la Tabla 8 (con los ratios tabulados 0.77 × 0.95 × 0.58 da 0.75); si los ratios contra FM25 difieren de los ratios contra seed tabulados, hay que tabularlos (ver N9).
Existe una determinación independiente Gaia-based de Europa (y Amphitrite, Eunomia): Siltala & Granvik 2022, A&A 658, A65 (verificada por web) — no citada. Es el cross-check externo natural para el caso más problemático del paper y para dos cuerpos más del estudio.
Solución exigida:
(1) Discutir Europa explícitamente: o el fit es una subestimación por regresión masa↔órbita no capturada por el flag (y entonces el flag de identificabilidad tiene un falso positivo conocido que hay que reportar y, quizá, endurecer — p. ej. exigir también $|z|<3$ contra el seed para el status "measured", o correr el bootstrap también en Europa aunque lev = 0.43), o la referencia DE441 de Europa es cuestionable (citar evidencia). (2) Reescribir la frase de determinabilidad a-priori acotándola a los casos que la cumplen y nombrando las excepciones. (3) Matizar el abstract. (4) Citar Siltala & Granvik (2022) y comparar Europa/Eunomia/Amphitrite contra ella.

Criterio de validación de cierre:
Tabla 8 acompañada de una columna o párrafo con $z$ vs. FM25 para los seis "measured"; regla de status reproducible por script (criterio publicado → mismos flags); el número 0.85/0.75 reconciliado; cita a SG22 presente y discutida.

[M4] MAYOR — La extrapolación catálogo-wide del censo de umbral usa escalado $N^2$ válido solo para muestra representativa, pero declara una muestra estratificada — y sus propios números delatan la contradicción
Panel: Estadística matemática.
Ubicación: Sect. 3.2 (aa_encounters.tex:429-434 y aa_encounters.tex:449-457).

Descripción del problema:
El texto dice que los 10⁴ cuerpos fueron "a stratified sample drawn to span the belt's $(a,e,i)$ range rather than a density-weighted random draw", y en el mismo aliento extrapola el conteo de pares en banda con $(449,454/10,000)^2$ — un escalado que solo es válido si la muestra es un draw aleatorio de la población (la densidad de pares depende del muestreo con pesos poblacionales; una muestra estirada a cubrir el rango tiene menos pares cercanos por construcción). Las dos afirmaciones no pueden ser ciertas a la vez. Y hay un dato duro: el conteo observado (17 469 pares en $[0.05,0.06)$) coincide al 0.4 % con la predicción para un draw aleatorio density-weighted ($80,072,774/2020 \times (1.2^{2.01}-1) \approx 17,546$, usando la pendiente 1.01 del propio paper). O sea: los números se comportan exactamente como un draw aleatorio, contradiciendo la descripción del diseño. Una de las dos cosas está mal escrita: o el sample fue efectivamente aleatorio (y entonces el argumento de "estimates the effect uniformly across the belt" es falso), o fue estratificado (y entonces la extrapolación $\sim3.5\times10^7$ pares en banda / $\sim10^5$ censados es inválida tal cual, y de paso la tasa 0.70 % no es la tasa poblacional que se multiplica).

Solución exigida:
Determinar del registro de la corrida qué diseño se usó realmente (candidato: script de scripts/validate/, hipótesis de ubicación) y reescribir el párrafo en consecuencia: si aleatorio, quitar "stratified... rather than density-weighted" y mantener la extrapolación; si estratificado, re-hacer la extrapolación con pesos (o re-correr con un draw aleatorio, que a 10⁴ cuerpos es barato) y recalcular el "~10⁵ censored" del abstract-adjacent y Conclusiones.

Criterio de validación de cierre:
El párrafo declara el diseño con una línea de procedencia (semilla, criterio de muestreo), y un check numérico publicado en el repo: conteo de pares en banda del sample vs. predicción $N^2$ del catálogo — si difieren >2σ, la extrapolación del texto usa el método corregido.

[M5] MAYOR — Fig. 3 usa una "orbital-elements source" que no es la del pipeline, dejando 79 % de los cuerpos afuera, con caption y texto en contradicción directa
Panel: Árbitro de manuscrito + Algoritmos.
Ubicación: Fig. 3 y caption (aa_encounters.tex:348-362); texto en aa_encounters.tex:309-313; conexo a Tabla 6 (clase "—").

Descripción del problema:
El catálogo entero se construyó propagando elementos MPCORB de los 449 454 cuerpos — los elementos osculantes existen por construcción para el 100 % de los cuerpos que encuentran. Sin embargo la Fig. 3 solo pudo plotear 92 971 de 449 213 (20.7 %) por elementos "absent from the orbital-elements source": esa fuente no puede ser MPCORB, y el texto jamás dice cuál es. Encima: (a) el texto afirma que la figura "confirm[s] the sample is the ordinary numbered population, not a biased subset" mientras su propio caption dice "this panel is a biased subsample of the parent population"; (b) el caption atribuye los faltantes a "mostly high-number, recently numbered objects", pero la Tabla 6 muestra el mismo síntoma (clase dinámica "—", derivada de elementos) en (20), (21), (27), (29), (30), (40), (44) — asteroides de dos dígitos. El diagnóstico del caption es, por evidencia interna del propio paper, incorrecto. Todo esto es innecesario: con los elementos MPCORB del propio freeze la figura se hace con el 100 % de los cuerpos.

Solución exigida:
Regenerar la Fig. 3 desde los elementos MPCORB del snapshot congelado (cobertura 100 %), eliminar del caption toda la disculpa por el subsample sesgado, y restaurar la afirmación del texto (que entonces sí será cierta). Reparar la misma fuente para class_1/class_2 (conexo con B1). Identificar y documentar qué era la "orbital-elements source" fallida (candidatos: src/ingest/gaia_orbits.py vs mpcorb.py — hipótesis de ubicación).

Criterio de validación de cierre:
Caption de Fig. 3 reporta "elements for 449,213/449,213 encountering bodies" (o la cifra completa real); texto y caption ya no se contradicen; cero clases "—" para numerados con elementos en el snapshot.

[M6] MAYOR — Término de incompletitud no presupuestado: la membresía es Kepler pero la candidacy es N-body, y el margen (7.2 mAU) es menor que la discrepancia máxima medida entre propagadores (15.2 mAU)
Panel: Algoritmos/completitud.
Ubicación: Sect. 2.2 (membership rule, aa_encounters.tex:249-258) contra Sect. 3.1 (aa_encounters.tex:401-404) y el presupuesto de Sect. 3.

Descripción del problema:
Hallazgo de texto (arquitectural): un par entra al catálogo si su distancia Kepler refinada es <0.05 AU, pero solo llega al refinador si el scan grueso N-body lo puso dentro de $r_q = 0.0572$ AU en algún sample. El margen entre umbral y radio de query es $0.0572-0.0500 = 7.2$ mAU (menos aún a mitad de paso, donde parte del margen la consume la velocidad). El propio paper mide $|\Delta d|$ Kepler↔N-body con máximos de 11.3 mAU (Stage A) y 15.2 mAU (Stage B) — mayores que el margen. Luego existen pares con distancia Kepler <0.05 AU cuya trayectoria N-body nunca entra en $r_q$: son miembros por definición que jamás fueron candidatos. La fracción es chica (cola > p99 de la distribución de Δd, condicionada a estar cerca del umbral; orden ≤10³ pares sobre 80 M), pero el paper vende un presupuesto de incompletitud medido con tres términos y este cuarto término pipeline-interno no aparece ni acotado — a diferencia del censoring de umbral, que sí está medido con esmero.

Solución exigida:
Añadir al presupuesto una cota de este término: con la distribución empírica de $\Delta d$ ya medida en Stage B, calcular $P(\Delta d > r_q - d_{\rm samp})$ integrada sobre los pares con $d_{\rm Kep}$ cerca del umbral, y publicar el número (probablemente "≲10³ pares, subdominante frente a los ~10⁵ del censoring"). Una frase en Sect. 3.2 o una nueva subsección corta basta.

Criterio de validación de cierre:
El presupuesto de Sect. 3 enumera cuatro términos pipeline-internos, cada uno con número y método; el nuevo término tiene un script reproducible como los otros tres.

HALLAZGOS MENORES
[N1] MENOR — "~100× below the Kepler two-body error" no se sigue de los propios números
Panel: Árbitro de manuscrito. Ubicación: aa_encounters.tex:410-414.
Descripción: Truncación de perturbadores: mediana 1.3 µAU vs. mediana Kepler 12 µAU → 9×; p99: 67 µAU vs 2.5 mAU → 37×. El "~100×" solo vale comparando máximos (80 µAU vs 11.3 mAU). Solución: escribir "~10× at the median (and ~40× at p99)". Cierre: el ratio impreso se reproduce desde los percentiles tabulados en la misma frase.

[N2] MENOR — El abstract dice "recovers all four calibrator masses" pero Pallas es status "bound"
Panel: Árbitro de manuscrito. Ubicación: abstract (results) vs Tabla 8 fila Pallas.
Descripción: Pallas: ratio 1.20, $z=+2.2$, status bound ($N=8$). "Recovers" sin calificación sobrevende: el engine la acota, no la mide. Solución: "recovers the four calibrators at $|z|<3$ (Pallas, with only eight encounters, as a bound)". Cierre: abstract y Tabla 8 usan el mismo vocabulario de status.

[N3] MENOR — Fórmula de determinabilidad dimensionalmente incorrecta
Panel: Física/astrodinámica. Ubicación: aa_encounters.tex:842-844.
Descripción: $\Delta\theta \sim GM/(b,v_{\rm rel},\Delta)$ tiene unidades de s⁻¹, no de ángulo: falta el tiempo de acumulación ($\Delta\theta \sim 2GM,t/(b,v_{\rm rel},\Delta)$, con $t$ el arco post-encuentro). Con $t\sim1$ yr el número final del párrafo (umbral $\sim10^{-11},M_\odot$) sí sale — la conclusión es correcta, la fórmula impresa no. Solución: añadir el factor $t$ y una frase de qué $t$ se usó. Cierre: análisis dimensional de la ecuación cierra en radianes.

[N4] MENOR — "$\sim$150,000 bodies" en Sect. 2.2 contradice el universo de 449,454 declarado en Sect. 2.1
Panel: Árbitro de manuscrito. Ubicación: aa_encounters.tex:193-194.
Descripción: El costo naive está estimado con un N obsoleto (150 k → "$10^{10}$ pairwise distances"); con el N real son ~$10^{11}$. Solución: actualizar a 449 454 y $10^{11}$. Cierre: un solo N en todo el paper.

[N5] MENOR — Nombres de archivo inconsistentes entre Sect. 6 y Apéndice A
Panel: Árbitro de manuscrito. Ubicación: aa_encounters.tex:1060-1065 vs aa_encounters.tex:1134-1140.
Descripción: Sect. 6: encounters_characterized_b1fix.parquet y ..._005au_b1fix.parquet; Apéndice: encounters_characterized_full.parquet y ..._005au.parquet (sin sufijo). Cuatro nombres para dos productos. Solución: unificar con los nombres del depósito final (conexo con B4). Cierre: grep de nombres .parquet en el tex devuelve un set consistente con el registro del archivo.

[N6] MENOR — El Apéndice atribuye a Sect. 2.3 una convención de ordenamiento que Sect. 2.3 no establece
Panel: Árbitro de manuscrito. Ubicación: aa_encounters.tex:1144-1145.
Descripción: "Bodies are ordered so that subscript 1 is the larger (perturber) body (Sect. 2.3)" — la Sect. 2.3 nunca lo dice. Solución: añadir la convención en 2.3 o quitar la cross-ref. Cierre: la referencia apunta a texto que existe.

[N7] MENOR — "(Sect.~below)" dos veces en Sect. 5
Panel: Árbitro de manuscrito. Ubicación: aa_encounters.tex:854 y aa_encounters.tex:858.
Descripción: Renderiza "Sect. below" — no es una referencia válida en estilo A&A. Solución: "see the External σ paragraph below" o \ref a párrafos numerables. Cierre: cero ocurrencias de "Sect.~below".

[N8] MENOR — El "0.70 %" del abstract/conclusiones se lee como fracción del catálogo, pero es la tasa dentro de la banda [0.05, 0.06)
Panel: Estadística matemática. Ubicación: abstract (results); aa_encounters.tex:1089-1090.
Descripción: La fracción de catálogo censada es ~0.2–0.3 % (~1.5–2.5×10⁵ de ~80 M + banda); el 0.70 % es la tasa de crossing en la banda descartada. Tal como está, un lector citará "0.70 % incompleteness" — número equivocado. Solución: en abstract: "a measured downward-crossing rate of 0.70 % in the discarded boundary band, implying ~10⁵ censored encounters (~0.25 % of the catalogue)". Cierre: ambas cantidades aparecen con sus denominadores explícitos.

[N9] MENOR — "geometric-mean ratio of 0.85" no reproducible desde la Tabla 8 (da 0.75) y base del ratio ambigua
Panel: Estadística matemática. Ubicación: aa_encounters.tex:956-958.
Descripción: Con los ratios tabulados de los tres measured non-calibrators (0.77, 0.95, 0.58) la media geométrica es 0.75; si el 0.85 usa ratios vs FM25 (no tabulados), el lector no puede verificarlo. Solución: tabular los ratios vs FM25 o corregir el número (conexo con M3). Cierre: el 0.85/0.75 se reproduce desde números impresos.

[N10] MENOR — CIs binomiales sobre pares tratan como independientes eventos que comparten cuerpos
Panel: Estadística matemática. Ubicación: aa_encounters.tex:433-434 (CI 0.59–0.83 %); aa_encounters.tex:486 (CI 76.27–76.49 %).
Descripción: Cada cuerpo participa en muchos pares (mismo cuerpo de alta-e genera crossings correlacionados); los CIs Wilson/binomiales reportados son anti-conservadores. Verifiqué que la aritmética de ambos CIs es correcta bajo independencia — el problema es el supuesto. Solución: CI por bootstrap por-cuerpo (cluster bootstrap), o una frase reconociendo el clustering y ensanchando. Cierre: CI recalculado con clustering declarado en el texto.

[N11] MENOR — Observabilidad "Gaia" calculada desde la Tierra y con límite en banda V para un límite instrumental definido en G
Panel: Física/astrodinámica. Ubicación: aa_encounters.tex:270-275; schema (aa_encounters.tex:1178-1179: "from Earth").
Descripción: Gaia está en L2 (~0.01 AU de la Tierra): cerca del corte de 45° la elongación geocéntrica misclasifica una capa de ~0.5°. Y el límite débil de Gaia es $G\approx20.7$–21; el flag usa $V<21$ sin declarar la conversión $V!-!G$ (~0.2 mag para asteroides). Ambas aproximaciones son razonables pero indeclaradas para un flag llamado gaia_observable. Solución: una frase declarando ambas aproximaciones y su efecto de borde. Cierre: caption/schema dicen "geocentric elongation (L2 offset ≤0.6°) and V as proxy for G".

[N12] MENOR — Conclusiones: "Unlike prior encounter-based work, which selects individual events by hand" contradice la Introducción
Panel: Árbitro de manuscrito + Estado del arte. Ubicación: aa_encounters.tex:1085-1086 vs aa_encounters.tex:116-118.
Descripción: La Intro describe (correctamente) a FM25 como "systematic, signal-ranked search" y a Goffin como fit global sin listas curadas; las Conclusiones descalifican todo el trabajo previo como selección a mano. Además "candidate pairs" (abstract) vs "real 3D close encounters" (conclusiones) para las mismas 80 M filas. Solución: alinear con la formulación de la Intro ("published truncated or without a characterised selection function") y unificar candidate/real. Cierre: ambas frases consistentes entre sí y con el abstract.

[N13] MENOR — Entrada li2023 del .bib omite al cuarto autor (Chen, J.)
Panel: Estado del arte (verificado por web). Ubicación: references.bib:119-127.
Descripción: El paper AJ 166, 93 es de Li, Yuan, Fu y Chen (Purple Mountain Observatory). Solución: añadir "and {Chen}, J.". Cierre: entrada coincide con ADS.

[N14] MENOR — Literatura relevante no citada
Panel: Estado del arte. Ubicación: Sect. 5 (relation to other methods); Sect. 2.3 (fuentes de diámetros/densidades).
Descripción: (a) Siltala & Granvik 2022, A&A 658, A65 (masas Gaia de Eunomia, Amphitrite, Europa, Edna) — verificada por web; imprescindible dado M3. (b) Un paper PSJ 2025 sobre masas por encuentros mutuos en LSST (existencia verificada solo por título/DOI PSJ/ade3de en resultados de búsqueda; contenido no leído) — candidato a "future work" junto a la mención de LSST implícita en el campo. (c) Las cadenas IRAS/AKARI/NEOWISE/ocultaciones y las "zone-average densities" se usan sin citar fuentes primarias (Tedesco et al., Usui et al., Mainzer et al.; densidades p. ej. Carry 2012) — estándar A&A las exige. Solución: añadir esas citas. Cierre: cada fuente de datos física del pipeline tiene cita primaria.

[N15] MENOR — Especificaciones estadísticas incompletas en la capa de masas
Panel: Estadística matemática + Algoritmos. Ubicación: aa_encounters.tex:868-880; aa_encounters.tex:434-435; aa_encounters.tex:994-998.
Descripción: (a) $s_c$ "calibrated by bisection to $\chi^2_{\rm red}=1$" — no se dice sobre qué dataset (¿global sobre calibradores? ¿por perturber?); con los $\chi^2_{\rm red}$ tabulados en 0.81–0.99 la respuesta importa: si es por-fit, la columna $\chi^2$ no es diagnóstica. (b) El "conservative floor of ~0.42 % after excluding near-boundary cases" no define "near-boundary". (c) Para los 16 perturbers del ephemeris no se declara cómo se evita el double-counting del propio perturber en el force model (sí se declara para los out-of-sixteen). Solución: una frase por ítem. Cierre: las tres especificaciones presentes; un lector puede re-derivar el flujo sin el código.

[N16] MENOR — Efectos no discutidos/no acotados y detalles de paralelización
Panel: Física + Algoritmos. Ubicación: Sect. 2.2–2.3, Sect. 5.
Descripción: (a) Yarkovsky y no-gravitacionales sobre los test bodies en el arco FPR (~6 yr): despreciables con seguridad, pero un árbitro esperará la cota de una línea — es sistemática conocida del campo de masas. (b) Manejo de mínimos que caen en el borde entre bloques temporales paralelos (¿solape entre bloques? ¿dedup del mínimo global?) — no especificado; hipótesis de implementación, candidato src/detect/parallel.py. (c) Aritmética menor no reproducible exactamente: z de Fortuna impreso 1.25 vs 1.27 desde sus propios números; z de Hygiea −0.2 vs −0.14 desde ratio y σ (¿incluyen σ de literatura? decláralo en el caption de Tabla 8). Solución: una línea de cota para (a); una frase para (b); nota de caption para (c). Cierre: cotas y convenciones impresas; los z reproducibles desde el caption.

Veredicto global de publicabilidad
Revisión mayor. La arquitectura del paper es sólida y honesta (el presupuesto de incompletitud medido es un aporte genuino y la capa estadística de masas es más autocrítica que la media del campo), pero hoy conviven en el texto números de dos freezes distintos (B3), una tabla insignia con fechas imposibles según su propia definición de ventana (B2), otra con diámetros refutables con datos públicos y con su propia columna de masas (B1), y un método de refinamiento descrito de forma que no puede haber producido el catálogo que se valida (M1) — nada de eso puede ir a un árbitro externo.

Plan de corrección priorizado
B3 (regenerar híbrido/Stage B sobre el freeze B1-fixed) — precondición: refresca los números de 2.2, 3.1, 3.2, 3.4 de los que dependen B2, M6 y el abstract.
B1 + M5 (misma raíz probable: join de datos físicos/elementos) — reparar la fuente, regenerar Tabla 6, columna class, Fig. 3.
B2 (auditar clamps de borde en Tablas 1–2 sobre el catálogo regenerado; aplicar guard de borde uniforme).
M1 (describir en Sect. 2.2 la ventana de refinamiento realmente ejecutada + guard de ±12 h del híbrido) — habilita que las citas al "fix" en 4.1/4.2 resuelvan.
M2 (MC con covarianza de época consistente) → regenera Tabla 3 y el claim de significancia de 4.2/4.3.
M4 (aclarar diseño del experimento de censura; recalcular extrapolación si corresponde) → toca abstract y conclusiones.
M3 + N9 + N14a (discusión de Europa, claim de determinabilidad, cita SG22, ratios FM tabulados).
M6 (cuarto término del presupuesto — barato: usa la distribución Δd ya medida en el nuevo Stage B).
N1–N16 (texto, bib, captions) en una sola pasada editorial.
B4 (depósito con DOI, afiliación) inmediatamente antes de someter.
Tabla de origen
Panel	Hallazgos
1. Árbitro de manuscrito	B2 (co), B3 (co), B4, M5 (co), N1, N2, N4, N5, N6, N7, N12 (co)
2. Física / astrodinámica	B2 (co), N3, N11, N16a
3. Estadística matemática	M2 (co), M3 (co), M4, N8, N9, N10, N15, N16c
4. Algoritmos / completitud	M1, M6, B3 (co), N15c, N16b
5. Estado del arte	B1 (verificación externa de diámetros), M3 (SG22), N12 (co), N13, N14
Verificaciones externas realizadas (todo lo no listado acá quedó "según mi conocimiento, sin verificación externa" — en particular, el claim de novedad del catálogo all-pairs sobrevivió mis búsquedas, pero una búsqueda negativa no es prueba de inexistencia):

Sources: Fuentes-Muñoz et al. 2025, AJ 170 (IOP) · Li et al. 2023, AJ 166, 93 (IOP) · Ivantsov, Hestroffer & Eggl 2018, IAUS 330, 386 (ADS) · Siltala & Granvik 2022, A&A 658, A65 · 44 Nysa (Wikipedia, diámetro radar/albedo) · Masas por encuentros mutuos en LSST, PSJ (IOP)