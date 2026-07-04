"""Tests de src/orbdet/variational.py — ecuaciones variacionales.

Gates de T3:
  (i)  ∂x/∂elementos (analítico = Φ·J_elem) coincide con la diferencia finita de
       la propagación completa a < 1e-6 relativo.
  (ii) ∂x/∂GM por diferencias finitas centrales es estable bajo refinamiento de δ
       (meseta de Richardson).

Todas las pruebas que integran con REBOUND van marcadas ``slow``.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from src.orbdet.constants import GM_SUN
from src.orbdet.dynamics import AsteroidPerturber, propagate
from src.orbdet.kepler import KeplerElements
from src.orbdet.variational import (
    partial_wrt_gm,
    partial_wrt_gm_variational,
    partials_wrt_elements,
    propagate_with_stm,
    richardson_convergence_dgm,
)

_EPOCH = 2_457_000.5  # JD TDB, era Gaia
_EL = KeplerElements(
    a=2.7,
    e=0.15,
    i=math.radians(10.0),
    Omega=math.radians(80.0),
    omega=math.radians(60.0),
    M=math.radians(45.0),
)
_PERTURBERS = ("sun", "jupiter", "saturn")


# --- Matriz de transición de estado -----------------------------------------


@pytest.mark.slow
def test_stm_identity_at_epoch() -> None:
    """En dt=0 la matriz de transición de estado es la identidad."""
    st = propagate_with_stm(_EL, _EPOCH, np.array([_EPOCH]), perturbers=_PERTURBERS)
    assert np.allclose(st.stm[0], np.eye(6), atol=1e-10)


@pytest.mark.slow
def test_stm_symplectic_determinant() -> None:
    """El flujo hamiltoniano preserva el volumen de fase: det(Φ) = 1."""
    out = _EPOCH + np.array([150.0, -200.0])
    st = propagate_with_stm(_EL, _EPOCH, out, perturbers=_PERTURBERS)
    for k in range(out.size):
        assert np.linalg.det(st.stm[k]) == pytest.approx(1.0, abs=1e-7)


@pytest.mark.slow
def test_partials_positions_match_propagate() -> None:
    """Las posiciones del cómputo variacional coinciden con dynamics.propagate."""
    out = _EPOCH + np.array([-100.0, 0.0, 250.0])
    pos, _vel, _ds = partials_wrt_elements(_EL, _EPOCH, out, perturbers=_PERTURBERS)
    pos_ref = propagate(_EL, _EPOCH, out, perturbers=_PERTURBERS)
    assert np.allclose(pos, pos_ref, atol=1e-12)


def _fd_dstate_delements(
    el: KeplerElements, epoch: float, out: np.ndarray, steps: np.ndarray
) -> np.ndarray:
    """Diferencia finita central de ∂[r,v](t)/∂elementos a través de la dinámica."""
    base = np.array(el.as_array(), dtype=float)
    fd = np.zeros((out.size, 6, 6))
    for j in range(6):
        plus, minus = base.copy(), base.copy()
        plus[j] += steps[j]
        minus[j] -= steps[j]
        pp, vp = propagate(
            KeplerElements(*plus), epoch, out, perturbers=_PERTURBERS, return_velocity=True
        )
        pm, vm = propagate(
            KeplerElements(*minus), epoch, out, perturbers=_PERTURBERS, return_velocity=True
        )
        fd[:, :, j] = (np.hstack([pp, vp]) - np.hstack([pm, vm])) / (2.0 * steps[j])
    return fd


@pytest.mark.slow
def test_partials_wrt_elements_match_full_propagation_fd() -> None:
    """GATE (i): ∂x/∂elementos analítico vs diferencia finita de la propagación
    completa (a través de la dinámica N-cuerpos), < 1e-6 relativo.

    La referencia FD se extrapola con Richardson (pasos h y h/2 → orden h⁴) para
    eliminar la truncación O(h²): sobre un arco de ±300 d el partial respecto a
    ``a`` enrolla con el período y una FD simple a un solo paso deja residuo
    ~3e-6 puramente de truncación, no del partial analítico."""
    out = _EPOCH + np.array([-300.0, 300.0])
    _pos, _vel, dstate = partials_wrt_elements(_EL, _EPOCH, out, perturbers=_PERTURBERS)

    h = np.array([1e-6 * _EL.a, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6])
    fd_h = _fd_dstate_delements(_EL, _EPOCH, out, h)
    fd_h2 = _fd_dstate_delements(_EL, _EPOCH, out, h / 2.0)
    # Richardson: D ≈ (4·D(h/2) − D(h)) / 3, truncación O(h⁴).
    dstate_fd = (4.0 * fd_h2 - fd_h) / 3.0

    # Error relativo por columna (vector de 6 componentes de estado por elemento):
    # normalizar por la norma de la columna evita que una entrada casi nula —donde
    # el ruido de redondeo del integrador domina— infle el relativo punto a punto.
    diff = np.linalg.norm(dstate - dstate_fd, axis=1)  # (N, 6) por columna/elemento
    base = np.linalg.norm(dstate_fd, axis=1)
    rel = diff / np.where(base > 0.0, base, 1.0)
    assert rel.max() < 1e-6, f"max rel error por columna {rel.max():.2e}"


# --- Parcial respecto a GM del perturbador -----------------------------------

# Perturbador masivo (~4× Ceres en M_sun) en órbita cercana al objetivo: deflexión
# claramente medible y, a esta masa, respuesta lineal (curvatura en masa < 1e-3).
_PERT_AST = AsteroidPerturber(
    name="test_pert",
    mass_msun=2e-9,
    elements=KeplerElements(
        a=2.71,
        e=0.16,
        i=math.radians(10.2),
        Omega=math.radians(80.0),
        omega=math.radians(60.0),
        M=math.radians(46.5),
    ),
)
_OUT_GM = _EPOCH + np.array([-400.0, 400.0])


@pytest.mark.slow
def test_dgm_nonzero() -> None:
    """La parcial respecto a GM es no nula (el perturbador deflecta al objetivo)."""
    dgm = partial_wrt_gm(
        _EL,
        _EPOCH,
        _OUT_GM,
        perturber_index=0,
        perturbers=_PERTURBERS,
        asteroid_perturbers=(_PERT_AST,),
    )
    assert np.isfinite(dgm).all()
    assert np.linalg.norm(dgm) > 0.0


@pytest.mark.slow
def test_dgm_requires_positive_mass() -> None:
    bad = (_PERT_AST.__class__(name="z", mass_msun=0.0, elements=_PERT_AST.elements),)
    with pytest.raises(ValueError, match="masa nominal"):
        partial_wrt_gm(
            _EL,
            _EPOCH,
            _OUT_GM,
            perturber_index=0,
            perturbers=_PERTURBERS,
            asteroid_perturbers=bad,
        )


@pytest.mark.slow
def test_dgm_richardson_plateau() -> None:
    """GATE (ii): la parcial ∂x/∂GM es estable bajo refinamiento del paso δ.

    La deflexión es lineal en la masa a primer orden, así que la diferencia
    central es prácticamente exacta y el barrido de δ describe una meseta: los
    cambios relativos sucesivos quedan muy por debajo de 1."""
    diag = richardson_convergence_dgm(
        _EL,
        _EPOCH,
        _OUT_GM,
        perturber_index=0,
        perturbers=_PERTURBERS,
        asteroid_perturbers=(_PERT_AST,),
    )
    max_rel_change = diag["max_rel_change"]
    # Meseta: todos los pasos coinciden entre sí muy por debajo de 1.
    assert max_rel_change.max() < 1e-4, f"sin meseta: {max_rel_change}"


@pytest.mark.slow
def test_dgm_central_difference_is_linear() -> None:
    """A masa pequeña la respuesta es lineal: ∂x/∂GM por diferencia central debe
    coincidir con la pendiente secante (x(m)-x(0))/(GM_SUN·m)."""
    dgm = partial_wrt_gm(
        _EL,
        _EPOCH,
        _OUT_GM,
        perturber_index=0,
        perturbers=_PERTURBERS,
        asteroid_perturbers=(_PERT_AST,),
        rel_delta=1e-3,
    )
    # Secante respecto a masa nula (sin perturbador) vs masa nominal.
    pos0, vel0 = propagate(_EL, _EPOCH, _OUT_GM, perturbers=_PERTURBERS, return_velocity=True)
    posm, velm = propagate(
        _EL,
        _EPOCH,
        _OUT_GM,
        perturbers=_PERTURBERS,
        asteroid_perturbers=(_PERT_AST,),
        return_velocity=True,
    )
    secant = (np.hstack([posm, velm]) - np.hstack([pos0, vel0])) / (GM_SUN * _PERT_AST.mass_msun)
    scale = np.maximum(np.abs(dgm), np.abs(secant))
    rel = np.abs(dgm - secant) / np.where(scale > 1e-30, scale, 1.0)
    # Linealidad a ~1e-8 de masa solar: acuerdo a mejor que 1e-3.
    assert rel.max() < 1e-3, f"max rel {rel.max():.2e}"


# --- Parcial ∂x/∂GM analítica (partícula variacional de masa, F6) -------------


@pytest.mark.slow
def test_dgm_variational_matches_fd() -> None:
    """GATE F6: la parcial ∂x/∂GM analítica (partícula variacional de masa)
    coincide con la diferencia finita central < 1e-6 relativo por época.

    La referencia FD se extrapola con Richardson (δ y δ/2 → orden δ⁴) para eliminar
    la truncación O(δ²) de la diferencia central; así el residuo mide la partícula
    variacional contra la derivada verdadera, no contra el sesgo de paso de la FD.
    """
    kw = dict(perturber_index=0, perturbers=_PERTURBERS, asteroid_perturbers=(_PERT_AST,))
    dgm_var = partial_wrt_gm_variational(_EL, _EPOCH, _OUT_GM, **kw)

    fd_h = partial_wrt_gm(_EL, _EPOCH, _OUT_GM, rel_delta=2e-3, **kw)
    fd_h2 = partial_wrt_gm(_EL, _EPOCH, _OUT_GM, rel_delta=1e-3, **kw)
    dgm_fd = (4.0 * fd_h2 - fd_h) / 3.0

    diff = np.linalg.norm(dgm_var - dgm_fd, axis=1)  # (N,) sobre las 6 componentes
    base = np.linalg.norm(dgm_fd, axis=1)
    rel = diff / np.where(base > 0.0, base, 1.0)
    assert rel.max() < 1e-6, f"max rel error por época {rel.max():.2e}"


@pytest.mark.slow
def test_dgm_variational_finite_and_nonzero() -> None:
    """La parcial variacional es finita y no nula (el perturbador deflecta)."""
    dgm = partial_wrt_gm_variational(
        _EL,
        _EPOCH,
        _OUT_GM,
        perturber_index=0,
        perturbers=_PERTURBERS,
        asteroid_perturbers=(_PERT_AST,),
    )
    assert dgm.shape == (_OUT_GM.size, 6)
    assert np.isfinite(dgm).all()
    assert np.linalg.norm(dgm) > 0.0
