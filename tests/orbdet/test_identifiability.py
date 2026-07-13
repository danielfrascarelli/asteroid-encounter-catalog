"""Tests del criterio de identificabilidad por verosimilitud perfilada (M13 / T21).

Valida, con modelos lineal-Gaussianos sintéticos y con los JSON de ajuste reales,
que:

- El Δχ²(M=0) por curvatura (aproximación cuadrática) coincide con el perfil exacto
  ``χ²_prof(0) − χ²(M̂)`` en el régimen lineal-Gaussiano (donde ambos deben ser
  idénticos), y equivale a ``(M̂/σ_formal)²``.
- Masa fuerte (alta SNR) → identificable (Δχ² ≫ 9); masa nula/ruido → no
  identificable (Δχ² ≪ 9).
- El p-valor de falsa alarma bajo H0 tiene la corrección de frontera ``M ≥ 0`` y el
  umbral 9 corresponde a ~3σ.

Los tests son livianos: NO corren fits N-cuerpos. El perfil exacto se ejercita con un
problema de mínimos cuadrados lineal explícito (la órbita se "re-optimiza" resolviendo
las ecuaciones normales), para verificar la equivalencia con la curvatura.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

from src.orbdet.identifiability import (
    DEFAULT_THRESHOLD,
    delta_chi2_profiled,
    delta_chi2_quadratic,
    false_alarm_probability,
    identifiability_from_covariance,
    identifiability_from_profile,
    profile_chi2_curve,
    threshold_for_nsigma,
)

_JACK_DIR = Path(__file__).resolve().parents[2] / "data" / "output" / "orbdet" / "expanded_jack"


class _LinearGaussianModel:
    """Modelo lineal-Gaussiano ``y = A·[M, θ] + ruido`` con blanqueo unidad.

    ``A`` es el diseño ``(m, 1 + p)`` (columna 0 = masa, resto = órbita). Los residuos
    ``r = A·x − y`` ya están blanqueados (σ=1). Sirve para verificar que el perfil
    exacto (re-optimizar θ a M fija resolviendo mínimos cuadrados) reproduce la
    curvatura ``(M̂/σ_formal)²``, porque en este régimen ambos son exactos.
    """

    def __init__(self, design: np.ndarray, y: np.ndarray) -> None:
        self.A = np.asarray(design, dtype=float)
        self.y = np.asarray(y, dtype=float)
        # Ajuste conjunto (masa + órbita) por ecuaciones normales.
        self.x_hat, *_ = np.linalg.lstsq(self.A, self.y, rcond=None)
        r = self.A @ self.x_hat - self.y
        self.chi2_hat = float(r @ r)
        cov = np.linalg.inv(self.A.T @ self.A)
        self.mass_hat = float(self.x_hat[0])
        self.sigma_formal = float(math.sqrt(cov[0, 0]))

    def chi2_at_fixed_mass(self, mass: float) -> float:
        """χ² perfilado: fija la masa y re-optimiza θ (ecuaciones normales)."""
        A_theta = self.A[:, 1:]
        resid_target = self.y - self.A[:, 0] * float(mass)
        theta, *_ = np.linalg.lstsq(A_theta, resid_target, rcond=None)
        r = A_theta @ theta - resid_target
        return float(r @ r)


def _make_model(
    mass_true: float, noise: float, seed: int, m_obs: int = 200
) -> _LinearGaussianModel:
    rng = np.random.default_rng(seed)
    # Columna de masa + 3 columnas de "órbita", no degeneradas.
    mass_col = np.linspace(-1.0, 1.0, m_obs)
    theta_cols = np.column_stack(
        [np.ones(m_obs), np.sin(np.linspace(0, 6, m_obs)), np.cos(np.linspace(0, 4, m_obs))]
    )
    design = np.column_stack([mass_col, theta_cols])
    theta_true = np.array([0.3, -0.5, 0.2])
    y = mass_col * mass_true + theta_cols @ theta_true + rng.normal(0.0, noise, m_obs)
    return _LinearGaussianModel(design, y)


def test_threshold_and_nsigma() -> None:
    """3σ ↦ Δχ²=9 y es el default."""
    assert threshold_for_nsigma(3.0) == pytest.approx(9.0)
    assert threshold_for_nsigma(5.0) == pytest.approx(25.0)
    assert DEFAULT_THRESHOLD == pytest.approx(9.0)


def test_false_alarm_boundary_correction() -> None:
    """El p-valor con frontera M≥0 es la mitad del χ²_1 de dos colas y decrece con Δχ²."""
    from scipy.stats import chi2 as _chi2

    d = 9.0
    two_sided = float(_chi2.sf(d, df=1))
    assert false_alarm_probability(d, boundary=True) == pytest.approx(0.5 * two_sided)
    assert false_alarm_probability(d, boundary=False) == pytest.approx(two_sided)
    # Δχ²=9 (~3σ) es un p pequeño; Δχ²=1 (~1σ) mucho mayor.
    assert false_alarm_probability(9.0) < 1e-2
    assert false_alarm_probability(1.0) > 0.1
    assert false_alarm_probability(-5.0) == pytest.approx(0.5)  # clamp a 0


def test_quadratic_matches_exact_profile_linear_gaussian() -> None:
    """GATE: en el régimen lineal-Gaussiano, Δχ² cuadrático == perfil exacto."""
    model = _make_model(mass_true=5.0, noise=1.0, seed=1)
    d_quad = delta_chi2_quadratic(model.mass_hat, model.sigma_formal)
    d_prof = delta_chi2_profiled(model.chi2_at_fixed_mass, model.chi2_hat, mass_null=0.0)
    assert d_prof == pytest.approx(d_quad, rel=1e-8)
    # Y ambos == (M̂/σ_formal)².
    assert d_quad == pytest.approx((model.mass_hat / model.sigma_formal) ** 2, rel=1e-8)


def test_profile_curve_minimum_at_mhat() -> None:
    """El perfil χ²(M) es parabólico con mínimo en M̂ e igual a χ²(M̂) allí."""
    model = _make_model(mass_true=5.0, noise=1.0, seed=2)
    grid = np.linspace(model.mass_hat - 3.0, model.mass_hat + 3.0, 25)
    masses, chi2s = profile_chi2_curve(model.chi2_at_fixed_mass, grid)
    assert chi2s.shape == masses.shape
    k_min = int(np.argmin(chi2s))
    assert masses[k_min] == pytest.approx(model.mass_hat, abs=grid[1] - grid[0])
    assert chi2s.min() == pytest.approx(model.chi2_hat, rel=1e-6)


def test_strong_mass_is_identifiable() -> None:
    """Masa fuerte (alta SNR) → Δχ² ≫ 9 → identificable, por ambos caminos."""
    model = _make_model(mass_true=8.0, noise=1.0, seed=3)
    res_q = identifiability_from_covariance(model.mass_hat, model.sigma_formal)
    res_p = identifiability_from_profile(
        model.mass_hat, model.chi2_hat, model.chi2_at_fixed_mass(0.0)
    )
    assert res_q.identifiable and res_p.identifiable
    assert res_q.delta_chi2 > DEFAULT_THRESHOLD
    assert res_q.delta_chi2 == pytest.approx(res_p.delta_chi2, rel=1e-6)
    assert res_q.n_sigma_equiv == pytest.approx(math.sqrt(res_q.delta_chi2))
    assert res_q.p_value < 1e-3


def test_weak_mass_not_identifiable() -> None:
    """Masa débil positiva ahogada en ruido (SNR ~2) → Δχ² < 9 → NO identificable."""
    model = _make_model(mass_true=2.0, noise=6.0, seed=4)
    assert model.mass_hat > 0.0  # positiva pero de baja SNR
    res = identifiability_from_covariance(model.mass_hat, model.sigma_formal)
    assert not res.identifiable
    assert math.isfinite(res.delta_chi2) and res.delta_chi2 < DEFAULT_THRESHOLD


def test_nonphysical_mass_never_identifiable() -> None:
    """Masa ajustada ≤ 0 (no física) nunca es identificable."""
    res = identifiability_from_covariance(-1.0e20, 5.0e18)
    assert not res.identifiable
    assert math.isnan(res.delta_chi2)


def test_zero_sigma_gives_nan() -> None:
    """σ_formal no positiva → Δχ² nan (sin curvatura definida)."""
    assert math.isnan(delta_chi2_quadratic(1.0e20, 0.0))
    assert math.isnan(delta_chi2_quadratic(1.0e20, -1.0))


@pytest.mark.parametrize("name", ["ceres_fpr", "pallas_fpr"])
def test_real_json_calibrators_identifiable(name: str) -> None:
    """Los calibradores Big-4 (Ceres, Pallas) son identificables desde el JSON real.

    Usa ``mass_fit_kg`` y ``mass_fit_sigma_formal_kg`` guardados — demuestra que el
    criterio de perfil (camino cuadrático) se computa sin re-correr fits.
    """
    p = _JACK_DIR / f"{name}.json"
    if not p.exists():
        pytest.skip(f"sin JSON de ejemplo {p}")
    d = json.loads(p.read_text())
    res = identifiability_from_covariance(d["mass_fit_kg"], d["mass_fit_sigma_formal_kg"])
    assert res.method == "quadratic"
    assert res.identifiable, f"{name}: Δχ²={res.delta_chi2:.1f} debería superar 9"
    assert res.delta_chi2 > DEFAULT_THRESHOLD
