# Tribunal científico — evaluación adversarial del proyecto completo

> **Fecha:** 2026-07-04
> **Alcance:** paper (`docs/paper/aa_encounters.tex`), catálogo congelado (72.2M encuentros),
> capa de detección (`src/detect/`), motor de masas (`src/orbdet/`), presupuesto de
> completitud, y aporte científico frente al estado del arte 2002–2026.
> **Método:** cinco paneles independientes (árbitro de manuscrito, física/astrodinámica,
> estadística matemática, algoritmos/completitud, estado del arte con verificación web).
> Los hallazgos fueron verificados contra el **código y los datos reales**, no contra los
> docs; varios incluyen confirmación empírica sobre el parquet o recálculo desde los CSV.
> Este informe va **más allá** de `docs/paper_referee_report.md` (cuyo veredicto
> "minor–moderate revision" este tribunal considera demasiado indulgente).

---

## Veredicto global

**El paper NO es publicable en su estado actual. Revisión MAYOR, con regeneración del
catálogo incluida.** Hay un bug de implementación confirmado empíricamente que corrompe
la época y sesga la distancia de ~60 % de las filas del catálogo congelado (B1), el
universo muestral real no es el que el paper declara (B2), los diámetros tabulados están
mal y tumban una conclusión del abstract (B3), el claim central de novedad es falso
frente a literatura verificada que el paper no cita (B5), y la capa estadística de masas
tiene tres defectos estructurales (σ jackknife sin validez bajo leverage extremo, test de
consistencia sin potencia, validación circular de calibradores: B6–B8).

Lo que sí resiste auditoría (verificado independientemente): la cadena de tiempos
TCB→TDB y TT→TDB, light-time iterativo, paralaje de Gaia BCRS, proyección along-scan
exacta desde la covarianza publicada, el Jacobiano kepleriano (rederivado término a
término), la decisión FD-vs-variacional por backend, el techo de perturbadores
planetarios, la medición de censura de umbral como *diseño*, y la coherencia entre la
física de la señal de deflexión (estimación independiente de S/N) y lo que el proyecto
afirma poder medir. La honestidad interna del proyecto es real; el problema es que varias
de sus afirmaciones públicas no sobreviven la verificación contra su propio código y
datos, y su encuadre frente a la literatura es insostenible.

---

## A. Hallazgos BLOQUEANTES

### B1. Bug de ventana de refinamiento: ~60 % del catálogo tiene época recortada y distancia sesgada — **confirmado empíricamente**

El scan grueso corre a Δt = 12 h (`config.yaml → coarse_step_hours: 12.0`) pero el
refinador Kepler solo muestrea ±2 h alrededor del sample grueso
(`detection.refinement.window_hours: 2.0`; `src/detect/refine.py:163-186`). El mínimo
verdadero puede estar hasta 6 h del punto grueso; cuando el argmin cae en el borde de la
ventana, el código **devuelve el punto del borde sin interpolar** (`refine.py:184-186`).
Con fase uniforme, P(mínimo fuera de la ventana) = 2/3.

**Confirmación sobre los datos** (histograma del offset de `jd_tdb` respecto de la grilla
de 12 h, catálogo híbrido completo):

- Filas `refinement_method="kepler"` (63.5M, 87.9 % del catálogo): distribución
  **truncada en 2 h**, con 76.3 % de la masa apilada en [1.75, 2.25] h (esperado
  uniforme: ~8 %). Ninguna fila Kepler con offset > 2.25 h.
- Filas `nbody` (ventana ±12 h): distribución plana en [0, 6] h — descarta cualquier
  explicación dinámica.
- En el subset N-body re-refinado, el 72.4 % de los pares tenía el mínimo Kepler
  recortado: mediana |Δt_min| 2.33 h (vs 0.17 h en los no recortados), y 62.6 % de
  N-body-más-cerca (vs 51.6 % simétrico) — es decir, parte del "scatter Kepler vs
  N-body" de Stage A/B es este artefacto de implementación, no error de modelo.

**Consecuencias:** (1) `jd_tdb` incorrecto hasta ~4 h y `dist_au` sesgada hacia arriba
para ~43M filas; (2) el §refineerror del paper atribuye a "error de dos cuerpos" un
error que en parte es truncamiento algorítmico; (3) canal de falsos negativos no
presupuestado (pares descartados porque la distancia *en el borde* supera 0.05 AU);
(4) la descripción del método en el paper ("the minimum separation, encounter epoch,
and relative velocity are extracted") es falsa para ~60 % de las filas; (5) la medición
de censura del umbral (0.70 %) usó la misma ventana contaminada como referencia
(`scripts/validate/measure_threshold_false_negatives.py:164-186`).

**Exigido:** `window_hours ≥ coarse_step/2` (±6 h) o re-centrado iterativo; validación
`window_hours` vs paso grueso en `src/detect/pipeline.py` (hoy no existe); **regenerar el
catálogo congelado**; re-derivar Stage A/B separando artefacto de error de modelo;
re-medir la censura; corregir §Detection y §refineerror. Test de regresión: mínimo
sintético colocado entre samples gruesos con offset > window_hours (hoy no existe —
`tests/test_detection.py` valida la parábola pero nunca el caso argmin-en-borde).

### B2. El universo muestral del catálogo no es el que el paper afirma, y su tamaño no está declarado en ninguna parte

Título, abstract, §1 y §7 afirman exhaustividad "over the full numbered population". El
sidecar de provenance registra `subset: only_numbered, a ∈ [1.5, 4.0] AU`: el catálogo
**excluye NEAs con a < 1.5, troyanos de Júpiter (~6.000+ numerados) y todo a > 4 AU**.
Ese corte no aparece en el manuscrito como definición del universo (solo como "adverse
subset" del experimento de prefiltro, tex:370). El N de cuerpos propagados no se da en
ningún lado — ni en el sidecar (`n_asteroids: null`). §2.2 dice "~150,000 bodies"
mientras `docs/kepler_threshold_bias_paper.md:134` extrapola con "los 449k numerados":
un factor ~9 de ambigüedad en la extrapolación N² de censura (la cifra "~1.5–2.5×10⁵
censurados" del paper no es reproducible con el N que el propio paper declara).

**Exigido:** declarar N exacto y el corte en a; reescribir título/abstract/§1/§7;
registrar N en el sidecar; rehacer la extrapolación con el N declarado. Una
"completitud" sin denominador poblacional definido no es una completitud.

### B3. Diámetros falsos según el propio código; la tabla D≳100 km está mal en ambas direcciones y la conclusión "no new large–large encounter" no está soportada

§2.3 afirma "diameter … estimated from H with a class-dependent albedo **when no direct
measurement exists**". El código (`src/characterize/physical.py:10`,
`encounter.py:213-214`) aplica **albedo fijo 0.14 a todos los cuerpos**, sin albedo por
clase y sin usar jamás diámetros medidos. Verificable en las propias tablas: Ceres
aparece con 763 km (real: 939), (44) Nysa con 139 km (real: ~71, tipo E), (91) Aegina
con 60.6 km (real: ~104). Con diámetros medidos, **(7) Iris × (44) Nysa NO es un
encuentro both-D≳100 km y (51) Nemausa × (91) Aegina SÍ lo es**: la Tabla `tab:dd100`
está mal en ambas direcciones y la afirmación del abstract/§7 cae. Además "taxonomic
class estimated from H" es conceptualmente absurdo (lo que el código calcula es clase
*orbital*).

**Exigido:** cruzar con diámetros/albedos medidos (IRAS/AKARI/NEOWISE), regenerar todas
las tablas de §4, corregir la descripción del método y re-evaluar la conclusión.

### B4. El "completeness budget" omite el término dominante potencial (incertidumbre de los elementos de entrada) y ninguna distancia del catálogo lleva σ

Los tres términos medidos (Kepler-vs-N-body, censura, prefiltro) comparan el mismo
conjunto de elementos MPCORB bajo distintos modelos de fuerza: miden **consistencia
interna del pipeline**, no completitud física. La incertidumbre de los elementos
osculantes (along-track, creciente con |t − t_epoch| hasta ±1.4 años) no está propagada
ni acotada; para numerados de arco corto puede ser de miles a decenas de miles de km —
**del orden de las distancias de la Tabla `tab:closest` (1.094–2.590 km)**, publicadas
con cuatro cifras y sin barra de error. El par "genético" de §4.3 (49.511 km) tampoco
tiene σ.

**Exigido:** propagar covarianzas de elementos (o cota estadística estratificada por
calidad orbital vía re-muestreo) y añadir σ(d_min) por fila o por estrato; mientras
tanto, renombrar el presupuesto como "pipeline-induced incompleteness budget" y marcar
`tab:closest` y el par genético como no significativos.

### B5. El claim central de novedad es falso, y la bibliografía (8 entradas) omite justo los trabajos que lo refutan

El abstract y §1 fundan la novedad en que el trabajo previo usa eventos "hand-selected"
y que "no systematic catalogue of real 3D asteroid–asteroid encounters" existe.
Verificado contra las fuentes:

- **Fuentes-Muñoz et al. 2025 (AJ 170, 353)**: búsqueda sistemática 1.783 perturbadores ×
  1.07M asteroides, arcos completos hasta nov-2024, **mismo umbral 0.05 AU para MBAs**.
  Nada "a mano".
- **Ivantsov, Hestroffer et al. 2018 (IAUS 330, 386)**: catálogo de encuentros mutuos de
  todos los numerados 2013–2023, propagación relativista, construido específicamente
  para explotar la astrometría de Gaia — casi exactamente el objeto que el paper dice
  que no existe.
- **Galád & Gray 2002 (A&A 391)**: búsqueda sistemática de encuentros para masas sobre
  24.599 asteroides, hace 24 años.
- **Goffin 2014**: ajuste global contra toda la población numerada (el propio repo lo
  verificó en `CONTINUATION.md` §2.4).

Referencias ausentes que un árbitro experto conoce de memoria: **Li et al. 2023 (AJ 166
— 20 masas con Gaia DR3, el trabajo más directamente comparable)**, Siltala & Granvik
2020/2021/2022 (MCMC; Psyche; Eunomia/Europa — ¡tres perturbadores de este paper!),
Baer & Chesley 2017, Fienga et al. 2003 (¡usada como validación en el repo y no
citada!), Kretlow 2020, Park et al. 2021 (DE440/441 — usada como referencia de ratios
sin cita), David et al. 2023 (Gaia FPR SSO — la fuente de datos del motor de masas),
INPOP/EPM, y las citas de software (rebound, IAS15, ASSIST, astropy) que A&A exige.

**Exigido:** reescribir abstract y §1. La novedad defendible que queda: catálogo
**all-pairs** publicado íntegro con procedencia + presupuesto de incompletitud medido
(nadie publica eso; las listas de FM están truncadas al top-100 por perturbador). La
frase "hand-selected" debe desaparecer. Ronda completa de literatura: con 8 referencias
un editor de A&A lo devuelve sin enviarlo a referee.

### B6. σ jackknife sin validez bajo leverage extremo; la clasificación measured/not_identifiable es inestable al nivel de un encuentro

Recalculado desde los JSON de `expanded_jack/`: la réplica jackknife más influyente
aporta el **92 %** de la varianza en Juno, 89.5 % en Thisbe, 80.8 % en Psyche, 69.4 % en
Iris, 63 % en Ceres → σ_jack con ~1 grado de libertad efectivo. Quitando esa réplica:
snr_jack de Juno pasa de 3.42 → 14.95 (σ ×4.4), Thisbe 5.15 → 18.17, y **Cybele
2.92 → 4.82 (cruzaría el umbral y cambiaría de clase)**. Juno (3.42), Interamnia (3.07)
y Eunomia (3.85) están declarados "measured" con snr indistinguible del umbral 3 dada la
incertidumbre de la propia σ_jack. Además: el jackknife es **ciego a sistemáticos
común-modo** — Pallas lo demuestra dentro del propio dataset (6 réplicas en un rango del
5.6 %, σ_jack 4.1 %… y desvío real de +24 % vs Goffin); el código acepta N=3 (Davida,
masa ajustada **negativa** −1.06×10²⁰ kg, con `sigma_sys_kg` negativa por bug de signo
en `build_mass_catalog.py:128`); y las réplicas no repiten el pipeline completo
(reutilizan el s_c calibrado, el clipping a 4σ y warm-start del ajuste completo →
subestiman la variabilidad del procedimiento).

**Exigido:** reportar leverage top-1 por perturbador en JSON y catálogo; para leverage
> ~50 %, σ_jack no es defendible como σ (bootstrap o delete-d, y propagar SE(snr) al
criterio); N mínimo ≥ 10 para reportar σ_jack; |M| en σ_sys; masas negativas marcadas
no-físicas en el CSV; no presentar σ_jack como "σ externa" total (es dispersión
inter-encuentro; el común-modo va aparte y para N chico ninguna capa lo captura).

### B7. El cruce con Fuentes-Muñoz ("10/10 en |z|<3") es un test sin potencia frente al sesgo del 20–30 % que el propio proyecto documenta

Desviación fraccional mínima detectable a 3σ, calculada desde el CSV del cruce: Vesta
15 %, Ceres 17 %, Pallas 26 %, Hygiea 28 %, Psyche 29 %, Europa 33 %, Thisbe 48 %,
Eunomia 56 %, Juno 64 %, Interamnia 74 %. Para los 6 no-calibradores el test solo puede
rechazar desviaciones del 29–74 %; las desviaciones **observadas** son del 14–30 % y
sistemáticamente del mismo signo (5/6 por debajo de FM; media geométrica de ratios
≈ 0.79). "10/10 en |z|<3" no es evidencia de exactitud: es la afirmación de que barras
del 20–33 % no ven sesgos del 20–30 %. Agravantes: (i) 5 de los 10 "measured" son
cuerpos que `mass_determination_results.md` declara sesgados por regresión a masa cero
(ratio 0.39–0.72) — o son medidas (y hay que explicar z_total de −6σ/−8σ vs DE441) o son
sesgadas (y no cuentan como validación); (ii) la atribución del caso Pallas en el paper
es incorrecta (su z baja por el término f_sys·M, no por el jackknife — σ_jack de Pallas
es *menor* que la formal); (iii) la etiqueta `z_jack` en realidad usa
σ_total = √(max(σ_f,σ_j)² + (f_sys·M)²) — nombre engañoso; (iv) FM 2025 ajusta la misma
astrometría FPR: los errores no son independientes; (v) el caveat honesto del draft
("a wide jackknife σ makes almost any value consistent…") **fue eliminado del tex** —
regresión de honestidad exactamente donde el argumento es más frágil.

**Exigido:** reportar la potencia junto a cada z; test de signo sobre los ratios de
no-calibradores como resultado principal; renombrar z_jack → z_total; corregir Pallas;
reinstaurar el caveat; declarar la no-independencia con FM; una única σ oficial por masa
con definición explícita de cada z (hoy hay tres barras usadas según convenga y los z
tabulados no son reproducibles — p. ej. Pallas ratio 1.240 con σ_tot 7 % da z = 3.4, no
el +2.67 tabulado).

### B8. Validación circular de calibradores vía f_sys

f_sys = RMS(ratio−1) sobre 3 calibradores (Ceres/Vesta/Hygiea) = 4.16 %, y luego se
declara "los calibradores se recuperan con |z|<3" usando σ_total que incluye f_sys·M.
Con 3 puntos, |desvío|/RMS ≤ √3 ≈ 1.73 **por construcción**: los calibradores que
definen el piso no pueden salir mal. f_sys tiene además ~41 % de incertidumbre relativa
(n=3) y se aplica como constante universal a los 16 perturbadores, incluidos los débiles
cuyo sistemático dominante (regresión a cero, −30/−60 %) es de otra naturaleza. Conexo:
el piso s_c calibrado a χ²_red=1 con la masa libre puede absorber misfit de señal (la
deflexión es común-modo dentro de un cruce FOV, la misma forma que s_c²·11ᵀ) y nunca fue
estresado con inyecciones; la incertidumbre de s_c no se propaga a σ(M).

**Exigido:** leave-one-calibrator-out; declarar la incertidumbre de f_sys; no
extrapolarla a los débiles; test de sensibilidad ∂M/∂s_c (±30 %); e
**inyección-recuperación end-to-end sobre datos reales** (señal sintética + geometría y
ruido reales, pipeline completo incluyendo calibración de piso, clipping, jackknife y
clasificación) — hoy los mismos tránsitos FPR se usan para ajustar, calibrar s_c,
clippear, generar σ_jack, clasificar y validar, sin una sola partición ni inyección.

### B9. Falta la tabla de las 16 masas y el schema del catálogo

§5 afirma "The engine determines all sixteen" pero solo tabula 4 calibradores: las 16
masas con σ_formal, σ_jack, N, χ²_red, flag y z **no aparecen** — el resultado central
de media sección es inverificable. Y no existe tabla de columnas del catálogo
(nombre/tipo/unidad/descripción), obligatoria para un dataset paper A&A/CDS. La cadena
catálogo → targets tampoco está descrita (p. ej. "Pallas, limited to 6 encounters"
contradice los 9 pares del cross-match FM; los 6 son post-cortes de selección que el
paper nunca define).

**Exigido:** tabla completa de 16 masas + tabla de schema + descripción cuantitativa de
la selección de targets por perturbador.

### B10. Falso negativo externo conocido sin causa raíz, no divulgado en el paper

De 4 eventos de Fienga (2003) dentro de la ventana, el par (804)×(733) a 0.0138 AU —
bien adentro del umbral — **está ausente** del catálogo congelado. El diagnóstico
interno ("likely a prefilter issue") es inconsistente: el prefiltro se salta para
N > 5000 (`src/detect/pipeline.py:130-146`). Candidatos concretos: descarte por
distancia de borde de ventana (B1) o exclusión por el corte a ∈ [1.5, 4.0] (B2). Un
fallo de 1/4 en la única muestra de validación externa independiente, fuera de los tres
términos del presupuesto, contradice la tesis "medimos lo que falta" — y el paper no lo
menciona.

**Exigido:** root-cause reproduciendo el par aislado (propagar (804) y (733) con el
mismo snapshot y grilla); regenerar el freeze con el fix de B1 (opción correcta) o
divulgar el gap como cuarto término del presupuesto.

---

## B. Hallazgos MAYORES

**Paper / claims**
- **M1. Explicación factualmente errónea de la brecha con FM 2025 (§3.4):** atribuye el
  91.6 % de pares no recuperados a que FM "no se restringe a 0.05 AU" — falso: FM usó
  exactamente 0.05 AU para MBAs. La causa dominante es la **ventana temporal** (FM busca
  sobre décadas; este catálogo cubre 2014-07→2017-05). Cuantificar cuántos pares FM
  tienen su encuentro fuera de la ventana; corregir §3.4 y
  `docs/completeness_vs_literature.md`.
- **M2. El "decisive point — 0 %" de §3.4 es tautológico:** el catálogo solo contiene
  filas < ~0.05 AU, así que "0 % de pares que contenemos pero situamos más allá de
  0.05 AU" no puede ser otra cosa. La inferencia "every literature pair … is recovered"
  no se sigue: un par FM perdido caería en "never approached", indistinguible. Eliminar
  la inferencia o clasificar los 25.962 no recuperados por época.
- **M3. Cota de cadencia (≲10⁻³ %) circular:** medida sobre la distribución de v_rel del
  propio catálogo, que por construcción no contiene a los perdidos rápidos, y sobre un
  universo que excluye a < 1.5 AU. Derivar la cota geométrica analítica
  (d(t) = √(d_min² + v²δ²) contra r_q; el bracketing cubre hasta ~190 km/s — el
  argumento correcto es más fuerte que el publicado, pero hay que escribirlo bien).
- **M4. Caption de Fig. 3 falsa:** "93,010 encountering bodies" es el row count de
  `gaia_orbits.parquet` (fuente distinta del universo propagado) filtrado al marco del
  plot — no cuerpos con ≥1 encuentro (`make_paper_figures.py:353-366`). Regenerar desde
  los cuerpos distintos del catálogo real.
- **M5. Inconsistencia de diseño catálogo (ventana DR3) vs masas (arco FPR ~2020):**
  encuentros 2017–2020 que perturban la astrometría ajustada no entran a la selección
  de targets ni al modelo. Justificar cuantitativamente o reconocer la asimetría.
- **M6. Tablas de §4 sin fuente definida e inconsistentes con FROZEN_RUN:** el par top
  difiere (987 km vs 1.094 km), falta un par de FROZEN_RUN en `tab:closest`, y la tabla
  "slowest" contiene un encuentro a 2 días del inicio de ventana con v_rel = 15 m/s —
  con alta probabilidad un mínimo truncado por borde de ventana, no flaggeado. Añadir
  columna de método/catálogo, flag `boundary_minimum`, épocas completas con escala.
- **M7. Utilidad del 72.2M sobrevendida:** el paper solo extrae valor de subconjuntos
  minúsculos; la mayoría de los pares son sub-km × sub-km sin señal medible. Falta la
  columna que los usuarios reales necesitan: **señal de deflexión estimada por par**
  (∝ GM/(b·v_rel)) — Ivantsov y FM rankean por señal, este catálogo no. Añadirla (o al
  menos la estratificación D × v_rel × d) y reordenar el argumento de utilidad; sin eso,
  72.2M filas es volumen, no valor.
- **M8. La sección de masas (~⅓ del paper) no aporta ninguna masa competitiva y lo
  sabe:** los 16 ya tienen masas FM/Siltala/Li mejores; 6/12 no-calibradores salen con
  regresión a cero; Fortuna/Metis peores que lo publicado. Reencuadrar §5 explícitamente
  como *demostración de uso del catálogo* (target selection + calibradores + Psyche) y
  comprimirla, o separarla a un segundo paper cuando haya una masa nueva (Pallas DR4).
- **M9. "FOV-block cov + σ jackknife + identificabilidad" no es novedad metodológica:**
  frente al MCMC de Siltala & Granvik (posteriors completos, bimodalidades nativas) y al
  criterio prior→posterior de FM, es una reimplementación competente pero funcionalmente
  menos sofisticada. Venderlo como transparencia/reproducibilidad, con comparación
  explícita, no como "methodological framework".

**Detección / completitud**
- **M10. Sin injection-recovery end-to-end en detección:** todas las mediciones comparan
  propagadores del mismo pipeline sobre los mismos elementos. Un test de inyección
  sintética (mínimos entre samples gruesos, v_rel altas, e altas) habría detectado B1
  antes que esta auditoría. Añadirlo como regresión permanente. (Ya identificado
  internamente como reparo probable; este tribunal lo eleva a obligatorio.)
- **M11. `docs/prefilter_recall.md:154-159` y `completeness_vs_literature.md:186-187`
  contradicen la corrección de 2026-05-31:** afirman que al freeze le faltan ≥143k
  encuentros por el prefiltro, pero el prefiltro **no se aplicó** a esa corrida
  (`skipped_large_n`). El paper lo trata bien como contrafactual; los docs no fueron
  actualizados. Alinear.
- **M12. El catálogo híbrido mezcla propagadores en `dist_au`** (87.9 % Kepler /
  12.1 % N-body, membresía definida por Kepler<0.05, **25.283 filas ≥ 0.05 AU** dentro
  de un catálogo de umbral 0.05, máx 0.05717). Documentar la regla efectiva en el
  sidecar y recomendar en el paper cómo filtrar para análisis de función de selección
  homogénea.
- **M13. Clasificación measured/not_identifiable ad-hoc:** corte duro snr_jack ≥ 3 (no
  es "curvatura de χ²" como se describe), sin calibración de tasa de falsos "measured"
  bajo masa nula, con numerador sesgado a la baja y denominador de ~1 gdl. Reemplazar o
  contrastar con verosimilitud perfilada χ²(M) (la maquinaria `--profile` ya existe) y
  calibrar por simulación. Multiplicidad (16 tests) sin discusión.

**Física (segundo orden, pero exigibles)**
- **M14. Ausente del paper la cuenta a priori de S/N de deflexión** (señal cruda vs
  post-absorción por el ajuste orbital). La estimación independiente del panel de física
  (a 10⁻¹² M☉ y 0.02 AU: S/N ≲ 1 por objetivo; a 10⁻¹¹ M☉: determinable con stacking)
  coincide con lo que el pipeline reporta — publicarla como criterio a priori fortalece
  el paper y delimita honestamente la frontera de determinabilidad.

---

## C. Hallazgos MENORES (lista compacta)

1. Docstring falso "MPCORB uses TDB" (`src/utils/time_utils.py:8` — es TT; el código
   convierte bien).
2. Comentario falso "builtin = DE440 reducida" (`src/orbdet/dynamics.py:14`; es
   erfa epv00/plan94, 1400× peor) y default silencioso `backend="rebound"` en la API
   pública de masas (`mass_determination.py:403,470`).
3. Sol builtin mezclado con posiciones DE440 al extraer elementos heliocéntricos
   (`gaia_adapter.py:174-177`; ~km, entra directo en modo perturber-orbit mpcorb;
   fraccionalmente ~10⁻⁵ — documentar o unificar efeméride).
4. Deflexión gravitacional de la luz (residuo por fuente a distancia finita, hasta
   ~sub-mas cerca de 45° de elongación) no acotada ni discutida — verificar qué corrige
   DPAC para SSO y acotar.
5. Jacobiano sin término ∂τ/∂param del light-time (~10⁻⁴ relativo) — documentar.
6. FD de masa sin guardia para M ≤ 0 (LM sin positividad; Davida convergió a masa
   negativa) y FD de e puede cruzar e=0 (`solve_kepler` lanza). Clamp o log-GM.
7. Pasos FD de elementos bajo ASSIST (1e-7 hardcodeado) sin chequeo de meseta (la de
   GM sí lo tiene).
8. Frame bias ICRS↔eclíptica (~17 mas) no documentado — se cancela en la cadena de
   observación y está acotado por el gate Horizons 0.17 mas, pero los "elementos
   eclípticos" del motor son "ICRS rotado por ε".
9. Definición de "encuentro" implícita (una fila = mínimo global de la ventana; mínimos
   locales secundarios no catalogados) — definir; para masas importan.
10. Marco declarado inconsistente: tex dice "heliocentric", sidecar dice "barycentric"
    (irrelevante para distancias mutuas, pero el texto debe ser exacto).
11. Criterio "Gaia-observable" (18.9 %) sin definición (elongación, magnitud límite).
12. Fig. 1: ajustar el índice de la ley de potencias contra dN/dd ∝ d² — diagnóstico
    gratuito de censura.
13. "No third-party orbit-determination dependency" como mérito — engañoso
    (ASSIST/rebound/DE440 hacen el trabajo pesado); reescribir como declaración de
    software. Falta code availability (DOI de código).
14. F3 (refutación del −4 %): la refutación pareada es sólida, pero el gate formulado
    sobre f_sys (RMS de 3 puntos, ±41 %) compara 4.158 % vs 4.257 % — ruido de redondeo.
    Reescribir el gate como Δmasa pareada (<0.25 %) y añadir la cota de escala para la
    cola del cinturón no modelado.
15. Radio de query del paper (0.0536 AU) ≠ ejecutado (0.0572, paso completo) — corregir.
16. Reproducibilidad: sidecar Kepler backfilled con `git: {}` vacío y hash MPCORB
    truncado a 16 hex; sidecar híbrido sin commit/deps/config; inconsistencias numéricas
    FROZEN_RUN vs sidecars (fine_step 60 vs 120 s; 305 MiB vs 137 MiB; 305,896 vs
    305,931); empates float no deterministas en `_merge_candidates`
    (`src/detect/parallel.py:310` — desempate lexicográfico (dist, t)).
17. `fuentesmunoz2024` (LPSC #2388) sin verificar contra ADS — eliminarla y citar solo
    el paper 2025.
18. Tabla de "candidate perturbers" (§4.4): todos tienen ya masa FM 2025 y el texto no
    lo dice ahí — añadir la masa FM como columna.
19. Muestra de censura "drawn to span the belt in (a,e,i)" = estratificada, no
    aleatoria, presentada como "representative" — justificar o corregir.

---

## D. Veredicto de publicabilidad y plan de corrección

**Estado actual:** rechazo probable en primera ronda en A&A. Un referee experto en masas
conoce FM 2025 y Siltala & Granvik de memoria (B5), y cualquier verificación de las
tablas contra datos públicos expone B3. El bug B1 es peor: si se detecta después de
publicar, compromete el catálogo entero.

**Ruta mínima a "publicable en A&A" (orden de dependencia):**

1. **Fix B1 + regenerar el catálogo congelado** (ventana ±6 h o re-centrado; tests de
   regresión H8; re-derivar Stage A/B y censura; root-cause B10 en el camino). Todo lo
   demás depende de esto — no tiene sentido pulir números de un catálogo que se va a
   regenerar.
2. **Universo y provenance (B2 + menores 16):** declarar N y el corte en a; sidecars
   completos; rehacer extrapolaciones.
3. **Física de cuerpos (B3):** diámetros medidos; regenerar §4.
4. **Presupuesto honesto (B4, M3, M10):** término de elementos (o cota), σ(d_min) por
   estrato, cota de cadencia analítica, injection-recovery de detección.
5. **Capa de masas (B6–B9, M13):** leverage y potencia reportados, LOCO para f_sys,
   criterio de identificabilidad por perfil calibrado con inyecciones end-to-end, tabla
   de 16 masas, una sola σ oficial.
6. **Reescritura del framing (B5, M1, M2, M7–M9):** novedad = all-pairs + presupuesto
   medido; literatura completa; columna de señal de deflexión; §5 como demo de uso.
7. Mecánica de submission (ya identificada): DOI CDS/VizieR, autores/ORCID,
   acknowledgements, citas de software.

**Alternativas honestas** si no se invierte lo anterior: Planetary & Space Science /
Advances in Space Research (listón menor), o Zenodo+arXiv como salida de archivo digna.
**Para un aporte fuerte** (paper de resultados): extender la ventana a FPR/DR4, columna
de señal por par, y una masa genuinamente nueva o mejorada (Pallas con DR4 — ya
identificado como F8) — eso convertiría el conjunto en dos papers sólidos.

---

## E. Origen de los hallazgos

| Panel | Hallazgos principales |
|---|---|
| Árbitro de manuscrito (verificó tex vs código/datos) | B2, B3, B4, B9, M2, M4, M6, M7, menores 9–13, 15 |
| Física / astrodinámica (rederivó Jacobianos, S/N, escalas) | B10, M14, menores 1–8; verificaciones limpias de tiempos/marcos/observación |
| Estadística matemática (recalculó desde CSV/JSON) | B6, B7, B8, M13, menor 14 |
| Algoritmos / completitud (confirmación empírica en parquet) | **B1**, M3, M10–M12, menores 16, 19 |
| Estado del arte (verificación web 2002–2026) | B5, M1, M8, M9, menores 17–18 |
