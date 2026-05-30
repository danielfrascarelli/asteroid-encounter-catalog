# Track B Stage 1 — diagnóstico de los outliers de Stage 2 (Alkeste/57942, Industria/176865)

> Entregable de **Track B Stage 1** del follow-up post-deepwork (disuelto; ver [ROADMAP.md](../ROADMAP.md) § "Estado actual").
> Diagnostica los dos peores χ²_red del batch de Stage 2 con residuos
> por-observación descompuestos en along-scan (AL) y across-scan (AC).
> Herramienta: [scripts/mass/diagnose_stage2_outliers.py](../scripts/mass/diagnose_stage2_outliers.py).

## Contexto

El batch de Stage 2 (joint 7-param + likelihood Mahalanobis 2D) dejó dos pares
con χ²_red anómalo ([docs/mass_layer_stage2_diagnostic.md](mass_layer_stage2_diagnostic.md)):

- (124) Alkeste → 57942: χ²_red = **84.45**
- (389) Industria → 176865: χ²_red saltó de 4.86 (AL-only, Stage 1) a **31.87** (2D, Stage 2)

El diagnóstico original planteó tres hipótesis: (a) **perturbador secundario** no
modelado, (b) una **sistemática across-scan (AC)** del catálogo Gaia que el
likelihood 2D expone y el AL-only ocultaba, (c) un **transit individual malo**.
Este documento las discrimina.

## Método

Para cada par se reconstruye el set de observaciones joint **idéntico** a Stage 2
(ventana one-sided, `loo_window=180 d`, `blackout=7 d`, 20 perturbadores de fondo),
se re-corre el fit single-target Mahalanobis 2D, y en el óptimo se descompone cada
observación en:

- **AL**: residuo proyectado sobre along-scan `r_AL` con su `σ_AL` (la dirección
  de ~mas de precisión de Gaia).
- **AC**: residuo across-scan `r_AC` con su `σ_AC` (la dirección pobremente
  restringida, ~arcsec).
- **χ²_2D**: Mahalanobis 2D por observación (la métrica de Stage 2).

`χ²_red_AL` y `χ²_red_AC` son chequeos 1-D independientes por eje (no son
sumandos aditivos del χ²_2D, que incluye la correlación AL-AC). Se busca además
perturbador secundario en el catálogo híbrido (otros encuentros del target a
<0.3 AU, ±90 d) y se contrasta contra el registro de cuerpos masivos.

## Resultados

| Par | n_obs | χ²ᵣ_2D | χ²ᵣ_AL | χ²ᵣ_AC | \|AL pull\| med | \|AC pull\| med | σ_AL med | σ_AC med | n(χ²>25) | top-1 transit | sec. masivos |
|-----|------:|-------:|-------:|-------:|-----:|-----:|------:|------:|------:|-----:|----:|
| 124→57942 | 139 | 84.5 | **348.7** | **0.08** | 5.06 | 0.20 | 2.3 mas | 612 mas | 68 (49 %) | 4.7 % | 0/123 |
| 389→176865 | 272 | 31.9 | **65.8** | **0.11** | 2.89 | 0.21 | 4.0 mas | 612 mas | 48 (18 %) | 7.6 % | 0/173 |

(`χ²ᵣ_2D` reproduce Stage 2 exactamente: 84.45 y 31.87.)

### Hipótesis (b) — sistemática AC: **REFUTADA**

La AC **no** es el problema. `σ_AC` mediana ≈ 612 mas vs `σ_AL` ≈ 2–4 mas: el
across-scan está **150–270× menos restringido**. Los residuos AC son diminutos
frente a su error (|AC pull| ≈ 0.2; **χ²ᵣ_AC ≈ 0.08–0.11**, *mejor* que lo
esperado). Todo el exceso vive en el **along-scan**: |AL pull| mediana 3–5σ,
χ²ᵣ_AL = 66–349. El likelihood 2D no "expuso" una sistemática AC — el AC pesa
casi nada; el misfit AL siempre estuvo ahí.

### Hipótesis (c) — transit individual malo: **REFUTADA (es pervasivo)**

No hay un único transit que domine: el peor aporta sólo **4.7 % / 7.6 %** del
χ² total. El misfit está **distribuido**: en 57942, **68/139 (49 %)** de los
transits tienen χ²_2D > 25; en 176865, 48/272 (18 %), algo más concentrado
(top-5 % de transits → 36 % del χ²; 21 transits con χ² > 100 → un subconjunto
contaminado). RMS along-scan = **45 / 36 mas**, ~10–20× la precisión AL de Gaia.

### Hipótesis (a) — perturbador secundario: **REFUTADA**

123 (resp. 173) encuentros secundarios del target a <0.3 AU dentro de ±90 d;
**0** involucran un cuerpo del registro masivo. Los partners más cercanos
(95227 a 0.0028 AU; 256337 a 0.0064 AU) son cuerpos sub-km (números >95k,
masas ~10¹⁴ kg) → deflexión despreciable. No hay perturbador secundario.

## Veredicto

El χ²_red alto de ambos outliers es un **misfit along-scan pervasivo de pocas
decenas de mas** que el modelo joint (6 deltas + 1 masa) no logra absorber. No
es señal de masa, no es una sistemática AC (refuta la hipótesis de Stage 2), no
es un transit aislado, y no hay perturbador secundario. Es residuo de
mismodelado orbital/astrométrico along-scan — coherente con el hallazgo de A2.6
de que sobre datos reales el along-scan carga estructura más allá del modelo.
La precisión AL ~mas de Gaia convierte un residuo físico de decenas de mas en
un χ²_red enorme.

**Implicación práctica**: estos pares **no** deben promoverse como detecciones
de masa. Conviene un corte de calidad por residuo along-scan (p. ej. marcar
fits con χ²ᵣ_AL ≫ 1 o RMS_AL por encima de un umbral) al rankear candidatos.
No cambia la conclusión global de Stage 2 (2/27 ≈ 7 % de outliers); acota su
naturaleza.

## Reproducir

```bash
docker compose run --rm pipeline python -m scripts.mass.diagnose_stage2_outliers \
  --out-prefix data/output/stage2_outliers/diag
# escribe diag_perobs.csv (residuos AL/AC por obs) y diag_summary.json
```
