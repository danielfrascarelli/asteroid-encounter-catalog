"""Criterio de identificabilidad de masa por verosimilitud perfilada (M13 / T21).

Motivación (tribunal científico, M13)
--------------------------------------
La clasificación ``measured`` / ``not_identifiable`` del catálogo de masas usaba un
corte duro ``snr_jack = M̂/σ_jack ≥ 3``. El tribunal objetó que (a) eso NO es la
"curvatura de χ²" que se decía usar; (b) no hay calibración de la tasa de falsos
"measured" bajo masa nula; (c) el numerador está sesgado y el denominador tiene
~1 grado de libertad efectivo (una réplica jackknife domina).

Este módulo implementa el criterio pedido: la **verosimilitud perfilada** de la masa.
Se compara el ajuste completo (masa ``M̂`` + órbita) contra el ajuste con **masa
fijada a cero re-optimizando la órbita**, y se mide

    Δχ²(M=0)  =  χ²_prof(0) − χ²(M̂),

donde ``χ²_prof(M) = min_θ χ²(M, θ)`` es el χ² **perfilado** sobre los 6·N elementos
orbitales ``θ`` a masa fija. Bajo el teorema de Wilks Δχ²(0) es la estadística del
test de razón de verosimilitud para ``H0: M=0``; para un parámetro es ``χ²_1`` y el
umbral ``Δχ² > 9`` corresponde a ~3σ (con corrección de frontera ``M ≥ 0``, ver
:func:`false_alarm_probability`).

Dos caminos de cómputo
-----------------------
1. **Perfil exacto** (:func:`profile_chi2_curve`, :func:`delta_chi2_profiled`,
   :func:`refit_orbit_at_fixed_mass`): re-optimiza la órbita a cada masa fija con la
   maquinaria de ajuste N-cuerpos. Es lo que pide literalmente el tribunal, pero cada
   evaluación es un re-fit completo (caro). :func:`refit_orbit_at_fixed_mass` lo hace
   sobre :func:`orbdet.mass_determination.determine_shared_mass` fijando la masa.

2. **Aproximación cuadrática / Laplace** (:func:`delta_chi2_quadratic`,
   :func:`identifiability_from_covariance`): en el régimen lineal-Gaussiano el χ²
   perfilado es una parábola en ``M`` cuya curvatura es exactamente ``1/σ_formal²``,
   siendo ``σ_formal`` la σ **marginal** de la masa (índice 0 de la covarianza conjunta
   de Fisher, complemento de Schur del bloque de órbita). Por lo tanto

       Δχ²_prof(0)  ≈  (M̂ / σ_formal)²,

   que se calcula **directamente** de lo ya guardado en los JSON de ajuste
   (``mass_fit_kg`` y ``mass_fit_sigma_formal_kg``), sin re-correr fits. Es exacto en
   el límite lineal-Gaussiano y es la aproximación de Laplace en general; se aparta del
   perfil exacto sólo si el modelo es fuertemente no lineal en la órbita cerca de M=0.

Notas
-----
- ``σ_formal`` (no ``σ_jack``) es el denominador correcto para la curvatura del perfil:
  ``σ_jack`` mezcla el error de regresión masa↔órbita y adolece del problema de ~1 gdl
  (M13.c); la curvatura del perfil ya está en la covarianza conjunta vía Schur.
- El módulo respeta el contrato de aislamiento de ``orbdet`` (solo stdlib + numpy/scipy
  + imports relativos).
"""

from __future__ import annotations

import dataclasses
import math
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from scipy.stats import chi2 as _chi2

from .dynamics import DEFAULT_PERTURBERS, AsteroidPerturber
from .kepler import KeplerElements
from .least_squares import LeastSquaresResult, levenberg_marquardt
from .mass_determination import (
    TargetObservations,
    _make_config,
    _make_pool,
    _target_resid_and_blocks,
)

# Umbral por defecto: Δχ² > 9 ≈ 3σ para un parámetro (M̂/σ = 3 ⇒ Δχ² = 9).
DEFAULT_THRESHOLD: float = 9.0


@dataclass(frozen=True)
class IdentifiabilityResult:
    """Veredicto de identificabilidad de una masa por verosimilitud perfilada.

    Attributes
    ----------
    mass_hat:
        Masa ajustada ``M̂`` (mismas unidades con que se llamó; kg o M_sun).
    delta_chi2:
        Δχ²(M=0) = χ²_perfilado(0) − χ²(M̂). ≥ 0 salvo ruido numérico.
    threshold:
        Umbral de decisión sobre Δχ² (``9`` ≈ 3σ por defecto).
    identifiable:
        ``True`` si ``delta_chi2 >= threshold`` (masa distinguible de 0).
    method:
        ``"quadratic"`` (Laplace desde σ_formal) o ``"profiled"`` (re-fit exacto).
    sigma_used:
        σ de la masa usada en el camino cuadrático (``None`` en el exacto).
    n_sigma_equiv:
        Significancia equivalente ``√Δχ²`` en σ (para un parámetro).
    p_value:
        Probabilidad de falsa alarma bajo ``H0: M=0`` (ver
        :func:`false_alarm_probability`, con corrección de frontera ``M ≥ 0``).
    """

    mass_hat: float
    delta_chi2: float
    threshold: float
    identifiable: bool
    method: str
    sigma_used: float | None = None
    n_sigma_equiv: float = float("nan")
    p_value: float = float("nan")


def threshold_for_nsigma(n_sigma: float) -> float:
    """Umbral de Δχ² correspondiente a una significancia de ``n_sigma`` σ.

    Para un único parámetro (la masa), el perfil de χ² es una parábola y una
    desviación de ``k·σ`` respecto del óptimo eleva Δχ² en ``k²``. Así ``3σ → 9``.

    Parameters
    ----------
    n_sigma:
        Significancia deseada en desviaciones estándar (p. ej. ``3.0``).

    Returns
    -------
    float
        ``n_sigma**2``.
    """
    return float(n_sigma) ** 2


def false_alarm_probability(delta_chi2: float, *, boundary: bool = True) -> float:
    """Probabilidad de falsa alarma (p-valor) de ``Δχ²`` bajo ``H0: M=0``.

    Bajo Wilks, si ``M=0`` es cierto y la masa fuese un parámetro libre en toda la
    recta, ``Δχ² ~ χ²_1`` y ``p = P(χ²_1 > Δχ²)``. Pero la masa está acotada a
    ``M ≥ 0`` (una masa negativa no es física), y el estimador cae sobre la frontera
    la mitad de las veces bajo H0: la distribución nula es la **chi-barra-cuadrado**
    ``½ δ(0) + ½ χ²_1`` (Chernoff 1954), de modo que ``p = ½ P(χ²_1 > Δχ²)``. Esta es
    la tasa de falsos "measured" que el tribunal (M13.b) pedía calibrar analíticamente;
    la calibración empírica por inyecciones de masa nula debe reproducir esta curva.

    Parameters
    ----------
    delta_chi2:
        Estadística Δχ²(M=0) observada (≥ 0).
    boundary:
        Si ``True`` (default), aplica la corrección de frontera ``M ≥ 0`` (factor ½).
        Si ``False``, usa el ``χ²_1`` de dos colas (masa sin restricción de signo).

    Returns
    -------
    float
        p-valor en ``[0, 1]``.
    """
    d = max(float(delta_chi2), 0.0)
    p = float(_chi2.sf(d, df=1))
    return 0.5 * p if boundary else p


def _finalize(
    mass_hat: float,
    delta_chi2: float,
    threshold: float,
    method: str,
    *,
    sigma_used: float | None = None,
    boundary: bool = True,
) -> IdentifiabilityResult:
    """Ensambla el :class:`IdentifiabilityResult` con significancia y p-valor derivados."""
    d = float(delta_chi2)
    return IdentifiabilityResult(
        mass_hat=float(mass_hat),
        delta_chi2=d,
        threshold=float(threshold),
        identifiable=bool(np.isfinite(d) and d >= threshold),
        method=method,
        sigma_used=sigma_used,
        n_sigma_equiv=math.sqrt(d) if np.isfinite(d) and d >= 0.0 else float("nan"),
        p_value=false_alarm_probability(d, boundary=boundary) if np.isfinite(d) else float("nan"),
    )


def delta_chi2_quadratic(mass_hat: float, sigma_formal: float) -> float:
    """Δχ²(M=0) perfilado bajo la aproximación cuadrática / Laplace.

    En el régimen lineal-Gaussiano, perfilar los 6·N elementos orbitales a masa fija
    deja una parábola en ``M`` con curvatura ``1/σ_formal²``, siendo ``σ_formal`` la σ
    marginal de la masa (índice 0 de la covarianza conjunta ``(JᵀC⁻¹J)⁻¹``). Por tanto

        Δχ²(0) = (M̂ − 0)² / σ_formal² = (M̂ / σ_formal)².

    Este es el valor computable **directo de los JSON** (``mass_fit_kg`` /
    ``mass_fit_sigma_formal_kg``) sin re-correr ningún fit.

    Parameters
    ----------
    mass_hat:
        Masa ajustada ``M̂``.
    sigma_formal:
        σ marginal (Fisher) de la masa, en las mismas unidades que ``mass_hat``.

    Returns
    -------
    float
        ``(M̂/σ_formal)²``, o ``nan`` si ``sigma_formal`` no es positivo/finito.
    """
    s = float(sigma_formal)
    if not np.isfinite(s) or s <= 0.0 or not np.isfinite(mass_hat):
        return float("nan")
    return (float(mass_hat) / s) ** 2


def identifiability_from_covariance(
    mass_hat: float,
    sigma_formal: float,
    *,
    threshold: float = DEFAULT_THRESHOLD,
    boundary: bool = True,
) -> IdentifiabilityResult:
    """Criterio de identificabilidad por curvatura del perfil (camino cuadrático).

    Usa :func:`delta_chi2_quadratic`; es el criterio aplicable a partir de lo que ya
    guardan los JSON de ajuste, sin re-fits. Una masa ``M̂ ≤ 0`` (no física) nunca es
    identificable.

    Parameters
    ----------
    mass_hat:
        Masa ajustada ``M̂``.
    sigma_formal:
        σ marginal de Fisher de la masa (mismas unidades).
    threshold:
        Umbral sobre Δχ² (default ``9`` ≈ 3σ).
    boundary:
        Corrección de frontera ``M ≥ 0`` para el p-valor (default ``True``).

    Returns
    -------
    IdentifiabilityResult
        Con ``method="quadratic"``.
    """
    if not np.isfinite(mass_hat) or mass_hat <= 0.0:
        return IdentifiabilityResult(
            mass_hat=float(mass_hat),
            delta_chi2=float("nan"),
            threshold=float(threshold),
            identifiable=False,
            method="quadratic",
            sigma_used=float(sigma_formal) if np.isfinite(sigma_formal) else None,
        )
    d = delta_chi2_quadratic(mass_hat, sigma_formal)
    return _finalize(
        mass_hat, d, threshold, "quadratic", sigma_used=float(sigma_formal), boundary=boundary
    )


def profile_chi2_curve(
    chi2_at_fixed_mass: Callable[[float], float],
    mass_values: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Perfil ``χ²_prof(M)`` sobre una grilla de masas fijas.

    Evalúa el χ² **perfilado** (mínimo sobre la órbita a masa fija) en cada valor de la
    grilla, delegando el re-ajuste de la órbita al callable ``chi2_at_fixed_mass``. El
    callable desacopla este perfilado de la maquinaria N-cuerpos concreta, de modo que
    la lógica es testeable con modelos sintéticos baratos; en producción se cablea a
    :func:`refit_orbit_at_fixed_mass`.

    Parameters
    ----------
    chi2_at_fixed_mass:
        ``M → χ²`` que fija la masa a ``M``, re-optimiza los elementos orbitales y
        devuelve el χ² mínimo resultante.
    mass_values:
        Grilla de masas ``(K,)`` a evaluar.

    Returns
    -------
    (masses, chi2s)
        ``masses`` ``(K,)`` (copia de la grilla) y ``chi2s`` ``(K,)`` los χ² perfilados.
    """
    masses = np.asarray(mass_values, dtype=float)
    chi2s = np.array([float(chi2_at_fixed_mass(float(m))) for m in masses], dtype=float)
    return masses, chi2s


def delta_chi2_profiled(
    chi2_at_fixed_mass: Callable[[float], float],
    chi2_hat: float,
    *,
    mass_null: float = 0.0,
) -> float:
    """Δχ²(M=mass_null) exacto perfilando la órbita a masa fija.

    Calcula ``χ²_prof(mass_null) − χ²(M̂)`` re-optimizando la órbita a
    ``M = mass_null`` (típicamente 0). Requiere un re-fit por evaluación → caro.

    Parameters
    ----------
    chi2_at_fixed_mass:
        ``M → χ²`` perfilado (ver :func:`profile_chi2_curve`).
    chi2_hat:
        χ² del ajuste conjunto en el óptimo ``M̂`` (p. ej. ``result.chi2``).
    mass_null:
        Masa de la hipótesis nula (default ``0.0``).

    Returns
    -------
    float
        Δχ² = χ²_prof(mass_null) − χ²(M̂).
    """
    return float(chi2_at_fixed_mass(float(mass_null))) - float(chi2_hat)


def identifiability_from_profile(
    mass_hat: float,
    chi2_hat: float,
    chi2_at_null: float,
    *,
    threshold: float = DEFAULT_THRESHOLD,
    boundary: bool = True,
) -> IdentifiabilityResult:
    """Criterio de identificabilidad a partir del perfil exacto ya evaluado.

    Parameters
    ----------
    mass_hat:
        Masa ajustada ``M̂``.
    chi2_hat:
        χ² en el óptimo conjunto.
    chi2_at_null:
        χ² perfilado a ``M=0`` (órbita re-optimizada a masa nula).
    threshold:
        Umbral sobre Δχ² (default ``9`` ≈ 3σ).
    boundary:
        Corrección de frontera ``M ≥ 0`` para el p-valor.

    Returns
    -------
    IdentifiabilityResult
        Con ``method="profiled"``.
    """
    d = float(chi2_at_null) - float(chi2_hat)
    return _finalize(mass_hat, d, threshold, "profiled", boundary=boundary)


def refit_orbit_at_fixed_mass(
    targets: list[TargetObservations],
    mass_fixed_msun: float,
    initial_elements: list[KeplerElements],
    perturber_elements: KeplerElements,
    epoch_jd_tdb: float,
    *,
    perturber_name: str = "perturber",
    background_perturbers: tuple[AsteroidPerturber, ...] = (),
    perturbers: tuple[str, ...] = DEFAULT_PERTURBERS,
    integrator: str = "ias15",
    dt_days: float = 1.0,
    n_lighttime_iter: int = 3,
    gm_rel_delta: float = 1e-3,
    backend: str = "assist",
    gr: bool = True,
    sys_floor_mas: float = 0.0,
    gm_variational: bool = True,
    max_iter: int = 80,
    n_workers: int = 1,
    **lm_kwargs,
) -> tuple[list[KeplerElements], LeastSquaresResult]:
    """Re-optimiza los 6·N elementos orbitales con la **masa fijada** a un valor.

    Espejo de :func:`orbdet.mass_determination.determine_shared_mass` pero SIN la
    columna de masa en el vector de parámetros: el Jacobiano es puramente en bloques
    (6 elementos por objetivo, sin columna densa de masa). Es la evaluación
    ``χ²_prof(M)`` que consume :func:`profile_chi2_curve` para el perfil **exacto**.

    .. warning::
        Cada llamada es un ajuste N-cuerpos completo (caro). Para el perfil basta con
        una evaluación a ``M=0``; úsese con criterio.

    Parameters
    ----------
    targets, initial_elements:
        Los ``N`` objetivos y su semilla orbital (mismo orden y longitud). Conviene
        semillar con los elementos ya ajustados del fit conjunto (warm start).
    mass_fixed_msun:
        Masa fija del perturbador (M_sun); para la hipótesis nula, ``0.0``.
    perturber_elements, epoch_jd_tdb, y demás:
        Idénticos a :func:`determine_shared_mass`.

    Returns
    -------
    (fitted_elements, result)
        ``result.chi2`` es ``χ²_prof(mass_fixed_msun)``. ``result.covariance`` es
        ``(6N, 6N)`` (solo órbita).
    """
    if not targets:
        raise ValueError("refit_orbit_at_fixed_mass requiere al menos un objetivo")
    if len(targets) != len(initial_elements):
        raise ValueError("targets e initial_elements deben tener la misma longitud")
    cfg = _make_config(
        epoch_jd_tdb,
        perturber_elements,
        perturber_name,
        background_perturbers,
        perturbers,
        integrator,
        dt_days,
        n_lighttime_iter,
        gm_rel_delta,
        backend=backend,
        gr=gr,
        sys_floor_mas=sys_floor_mas,
        gm_variational=gm_variational,
    )
    n_t = len(targets)
    n_par = 6 * n_t
    mass = float(mass_fixed_msun)
    warm = [dataclasses.replace(t, initial_elements=el) for t, el in zip(targets, initial_elements)]

    def _assemble(results) -> tuple[np.ndarray, np.ndarray]:
        resids: list[np.ndarray] = []
        blocks: list[np.ndarray] = []
        for k, (resid_k, _jac_mass_k, jac_elem_k) in enumerate(results):
            block = np.zeros((resid_k.size, n_par), dtype=float)
            block[:, 6 * k : 6 * k + 6] = jac_elem_k
            resids.append(resid_k)
            blocks.append(block)
        return np.concatenate(resids), np.vstack(blocks)

    pool = _make_pool(n_workers, n_t)
    try:

        def residual_and_jac(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
            args = [(mass, KeplerElements(*x[6 * k : 6 * k + 6]), warm[k], cfg) for k in range(n_t)]
            if pool is not None:
                results = pool.starmap(_target_resid_and_blocks, args)
            else:
                results = [_target_resid_and_blocks(*a) for a in args]
            return _assemble(results)

        x0 = np.concatenate([np.asarray(el.as_array()) for el in initial_elements])
        result = levenberg_marquardt(residual_and_jac, x0, max_iter=max_iter, **lm_kwargs)
    finally:
        if pool is not None:
            pool.close()
            pool.join()
    fitted = [KeplerElements(*result.x[6 * k : 6 * k + 6]) for k in range(n_t)]
    return fitted, result


__all__ = [
    "DEFAULT_THRESHOLD",
    "IdentifiabilityResult",
    "threshold_for_nsigma",
    "false_alarm_probability",
    "delta_chi2_quadratic",
    "identifiability_from_covariance",
    "profile_chi2_curve",
    "delta_chi2_profiled",
    "identifiability_from_profile",
    "refit_orbit_at_fixed_mass",
]
