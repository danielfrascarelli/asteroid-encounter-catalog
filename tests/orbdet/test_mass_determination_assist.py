"""Closing-loop del ajuste de masa con el backend ASSIST (T8b).

Gate: inyectar una masa de perturbador en observaciones sintéticas generadas con
el **mismo modelo de fuerzas state-of-the-art** (DE440 + GR + perturbadores) y
recuperarla con el ajuste conjunto cuyas parciales ∂x/∂elementos y ∂x/∂masa salen
por diferencias finitas sobre :func:`propagate_assist`. Es el equivalente del gate
T6 (rebound) pero sobre el backend que se usa en datos reales.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from src.orbdet.dynamics import AsteroidPerturber, planet_state_ecliptic
from src.orbdet.dynamics_assist import propagate_assist
from src.orbdet.frames import ecliptic_to_equatorial
from src.orbdet.kepler import KeplerElements, elements_to_state, state_to_elements
from src.orbdet.mass_determination import determine_mass_and_orbit
from src.orbdet.observation import light_time_correct, radec_from_positions

_EPOCH = 2_457_000.5
_TARGET = KeplerElements(
    a=2.70, e=0.15, i=math.radians(10.0),
    Omega=math.radians(80.0), omega=math.radians(60.0), M=math.radians(45.0),
)
_MASS_TRUE = 3e-9  # M_sun, ~6× Ceres → deflexión clara


def _perturber_elements() -> KeplerElements:
    r_t, v_t = elements_to_state(_TARGET)
    r_p = r_t + np.array([0.004, 0.0, 0.0])
    v_p = v_t + np.array([0.0, 5e-4, 0.0])
    return state_to_elements(r_p, v_p)


def _predict_assist(el, obs_jd, gaia_icrs, ast_perts):
    def bary(jd):
        return propagate_assist(
            el, _EPOCH, np.atleast_1d(jd), asteroid_perturbers=ast_perts, gr=True
        )

    _jd_ret, ast_icrs = light_time_correct(bary, obs_jd, gaia_icrs, n_iter=3)
    return radec_from_positions(ast_icrs, gaia_icrs)


@pytest.mark.slow
def test_closing_loop_assist_noiseless() -> None:
    """GATE: sin ruido, el backend ASSIST recupera la masa inyectada a ratio ≈ 1."""
    n = 20
    obs = _EPOCH + np.linspace(-200.0, 200.0, n)
    gaia_ecl = np.array([planet_state_ecliptic("earth", float(t))[0] for t in obs])
    gaia_icrs = ecliptic_to_equatorial(gaia_ecl)
    pert_el = _perturber_elements()
    pert_true = AsteroidPerturber("pert", _MASS_TRUE, pert_el)

    ra_t, dec_t = _predict_assist(_TARGET, obs, gaia_icrs, (pert_true,))
    pa = np.linspace(15.0, 165.0, n)
    sigma_al = np.full(n, 1.0)

    mass_fit, _el_fit, res = determine_mass_and_orbit(
        _TARGET, 0.4 * _MASS_TRUE, pert_el, _EPOCH,
        obs, ra_t, dec_t, pa, sigma_al, gaia_icrs,
        perturber_name="pert", backend="assist", gr=True, max_iter=40,
    )
    ratio = mass_fit / _MASS_TRUE
    assert res.converged
    assert abs(ratio - 1.0) < 0.05, f"ratio masa fit/true = {ratio:.4f}"
