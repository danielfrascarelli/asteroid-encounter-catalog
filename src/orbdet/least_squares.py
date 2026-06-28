"""Mínimos cuadrados no lineales: corrector diferencial Levenberg-Marquardt.

Solver genérico y autocontenido para el problema ``min_x ½‖r(x)‖²`` con ``r`` el
vector de residuos **ya blanqueados** (whitened, adimensionales) y Jacobiano
``J = ∂r/∂x`` analítico. Es el motor del corrector diferencial de órbitas (T5) y,
extendido con el parámetro de masa, del ajuste conjunto órbita+masa (T6).

Levenberg-Marquardt interpola entre Gauss-Newton (rápido cerca del mínimo) y
descenso de gradiente (robusto lejos), con amortiguamiento ``λ`` sobre la
diagonal de ``JᵀJ`` (escalado de Marquardt). La covarianza de los parámetros es
``(JᵀJ)⁻¹`` evaluada en la solución (los residuos ya están blanqueados, así que
``JᵀJ`` es la matriz de información de Fisher).
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np

# Firma del modelo: x (n,) → (residuos (m,), jacobiano (m, n)).
ResidualAndJac = Callable[[np.ndarray], tuple[np.ndarray, np.ndarray]]


@dataclass(frozen=True)
class LeastSquaresResult:
    """Resultado del corrector diferencial.

    Attributes
    ----------
    x:
        Vector de parámetros en la solución ``(n,)``.
    residuals:
        Residuos blanqueados en la solución ``(m,)``.
    cost:
        ``½‖r‖²`` en la solución.
    chi2:
        ``‖r‖²`` (χ² total, residuos ya blanqueados).
    dof:
        Grados de libertad ``m − n`` (≥1).
    chi2_reduced:
        ``chi2 / dof``.
    covariance:
        Covarianza de los parámetros ``(JᵀJ)⁻¹`` ``(n, n)`` en la solución.
    n_iter:
        Iteraciones consumidas.
    converged:
        ``True`` si paró por tolerancia (no por tope de iteraciones).
    """

    x: np.ndarray
    residuals: np.ndarray
    cost: float
    chi2: float
    dof: int
    chi2_reduced: float
    covariance: np.ndarray
    n_iter: int
    converged: bool


def levenberg_marquardt(
    residual_and_jac: ResidualAndJac,
    x0: np.ndarray,
    *,
    max_iter: int = 100,
    ftol: float = 1e-12,
    xtol: float = 1e-12,
    gtol: float = 1e-12,
    lambda0: float = 1e-3,
    lambda_up: float = 5.0,
    lambda_down: float = 0.2,
    lambda_min: float = 1e-12,
    lambda_max: float = 1e12,
) -> LeastSquaresResult:
    """Minimiza ``½‖r(x)‖²`` por Levenberg-Marquardt con Jacobiano analítico.

    Parameters
    ----------
    residual_and_jac:
        Callable ``x → (r, J)`` con ``r`` ``(m,)`` y ``J`` ``(m, n)``.
    x0:
        Punto inicial ``(n,)``.
    max_iter:
        Tope de iteraciones externas (cada una acepta un paso).
    ftol, xtol, gtol:
        Tolerancias de convergencia: reducción relativa de costo, norma relativa
        del paso, y norma del gradiente ``Jᵀr``, respectivamente.
    lambda0, lambda_up, lambda_down, lambda_min, lambda_max:
        Estado y actualización del amortiguamiento de Marquardt.

    Returns
    -------
    LeastSquaresResult
    """
    x = np.array(x0, dtype=float)
    r, jac = residual_and_jac(x)
    r = np.asarray(r, dtype=float)
    jac = np.asarray(jac, dtype=float)
    cost = 0.5 * float(r @ r)
    n = x.size
    m = r.size

    lam = float(lambda0)
    converged = False
    used_iters = 0
    for it in range(1, max_iter + 1):
        used_iters = it
        jtj = jac.T @ jac
        grad = jac.T @ r
        if np.linalg.norm(grad, ord=np.inf) < gtol:
            converged = True
            break

        diag = np.diag(jtj).copy()
        diag[diag <= 0.0] = 1.0  # escalado de Marquardt seguro
        # Bucle interno: subir λ hasta que el paso reduzca el costo.
        accepted = False
        while lambda_min <= lam <= lambda_max:
            a_mat = jtj + lam * np.diag(diag)
            try:
                dx = np.linalg.solve(a_mat, -grad)
            except np.linalg.LinAlgError:
                lam *= lambda_up
                continue
            x_new = x + dx
            r_new, jac_new = residual_and_jac(x_new)
            r_new = np.asarray(r_new, dtype=float)
            cost_new = 0.5 * float(r_new @ r_new)
            if cost_new < cost:
                step_norm = float(np.linalg.norm(dx))
                x_scale = float(np.linalg.norm(x)) + xtol
                cost_drop = (cost - cost_new) / max(cost, 1e-300)
                x, r, jac = x_new, r_new, np.asarray(jac_new, dtype=float)
                cost = cost_new
                lam = max(lam * lambda_down, lambda_min)
                accepted = True
                if cost_drop < ftol or step_norm < xtol * x_scale:
                    converged = True
                break
            lam *= lambda_up
        if converged or not accepted:
            break

    jtj = jac.T @ jac
    try:
        covariance = np.linalg.inv(jtj)
    except np.linalg.LinAlgError:
        covariance = np.full((n, n), np.nan)
    dof = max(m - n, 1)
    chi2 = float(r @ r)
    return LeastSquaresResult(
        x=x,
        residuals=r,
        cost=cost,
        chi2=chi2,
        dof=dof,
        chi2_reduced=chi2 / dof,
        covariance=covariance,
        n_iter=used_iters,
        converged=converged,
    )


__all__ = ["LeastSquaresResult", "levenberg_marquardt", "ResidualAndJac"]
