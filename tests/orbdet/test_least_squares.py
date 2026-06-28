"""Tests del solver Levenberg-Marquardt genérico (src/orbdet/least_squares.py)."""

from __future__ import annotations

import numpy as np

from src.orbdet.least_squares import levenberg_marquardt


def test_linear_least_squares_exact() -> None:
    """Para un modelo lineal r = A x − b, LM recupera la solución normal en 1 paso."""
    rng = np.random.default_rng(0)
    a_mat = rng.normal(size=(20, 3))
    x_true = np.array([1.5, -2.0, 0.7])
    b = a_mat @ x_true

    def rj(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        return a_mat @ x - b, a_mat

    res = levenberg_marquardt(rj, np.zeros(3), lambda0=1e-6)
    assert res.converged
    assert np.allclose(res.x, x_true, atol=1e-8)
    assert res.chi2 < 1e-16


def test_rosenbrock_residual_form() -> None:
    """Mínimo de Rosenbrock en forma de residuos: r=[10(y−x²), 1−x] → (1,1)."""

    def rj(p: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        x, y = p
        r = np.array([10.0 * (y - x * x), 1.0 - x])
        jac = np.array([[-20.0 * x, 10.0], [-1.0, 0.0]])
        return r, jac

    res = levenberg_marquardt(rj, np.array([-1.2, 1.0]), max_iter=200)
    assert res.converged
    assert np.allclose(res.x, [1.0, 1.0], atol=1e-6)


def test_covariance_matches_linear_theory() -> None:
    """Con residuos blanqueados, la covarianza es (AᵀA)⁻¹ para el caso lineal."""
    rng = np.random.default_rng(1)
    a_mat = rng.normal(size=(50, 4))
    x_true = rng.normal(size=4)
    b = a_mat @ x_true

    def rj(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        return a_mat @ x - b, a_mat

    res = levenberg_marquardt(rj, np.zeros(4))
    assert np.allclose(res.covariance, np.linalg.inv(a_mat.T @ a_mat), rtol=1e-6)
