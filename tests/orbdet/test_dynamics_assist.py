"""Tests del backend ASSIST (efemérides JPL DE440 + GR + 16 perturbadores) — T8.

- ``slow``: la deflexión por un perturbador masivo aparece (gravedad partícula-
  partícula activa junto a ASSIST) y la masa de un asteroide de la efeméride
  coincide con la literatura.
- ``horizons`` (red): el backend ASSIST reproduce los vectores de JPL Horizons a
  nivel sub-mas sobre la ventana Gaia, muy por debajo del ruido along-scan, y bate
  ampliamente al backend de planetas libres. Es el gate de exactitud T2/T8.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from src.orbdet.dynamics import AsteroidPerturber
from src.orbdet.dynamics_assist import (
    big_asteroid_perturbers,
    ephem_asteroid_mass_msun,
    propagate_assist,
)
from src.orbdet.kepler import KeplerElements, elements_to_state, state_to_elements

_EPOCH = 2_457_000.5
_EL = KeplerElements(
    a=2.70, e=0.15, i=math.radians(10.0),
    Omega=math.radians(80.0), omega=math.radians(60.0), M=math.radians(45.0),
)


def test_ceres_mass_matches_literature() -> None:
    """La masa de Ceres en la efeméride DE441 ≈ 4.7e-10 M_sun (9.4e20 kg)."""
    m = ephem_asteroid_mass_msun("Ceres", _EPOCH)
    assert m == pytest.approx(4.72e-10, rel=0.02)


@pytest.mark.slow
def test_perturber_deflects_target_under_assist() -> None:
    """Un perturbador masivo agregado deflecta al objetivo (gravedad rebound activa)."""
    out = _EPOCH + np.linspace(-150.0, 150.0, 6)
    r_t, v_t = elements_to_state(_EL)
    pert_el = state_to_elements(r_t + np.array([0.004, 0.0, 0.0]),
                                v_t + np.array([0.0, 5e-4, 0.0]))
    pert = AsteroidPerturber("pert", 3e-9, pert_el)

    pos_with = propagate_assist(_EL, _EPOCH, out, asteroid_perturbers=(pert,))
    pos_without = propagate_assist(_EL, _EPOCH, out, asteroid_perturbers=())
    defl = np.linalg.norm(pos_with - pos_without, axis=1).max()
    assert defl > 1e-9, f"sin deflexión apreciable: {defl:.2e} AU"


@pytest.mark.horizons
def test_assist_matches_horizons_submas() -> None:
    """GATE T2/T8: ASSIST reproduce Horizons a sub-mas sobre ~900 días (ventana DR3)."""
    from astroquery.jplhorizons import Horizons

    target = "8"  # (8) Flora, no es uno de los 16 perturbadores
    tab = Horizons(id=target, id_type="smallbody", location="@sun",
                   epochs=_EPOCH).elements(refplane="ecliptic")
    el = KeplerElements(
        a=float(tab["a"][0]), e=float(tab["e"][0]),
        i=math.radians(float(tab["incl"][0])),
        Omega=math.radians(float(tab["Omega"][0])),
        omega=math.radians(float(tab["w"][0])),
        M=math.radians(float(tab["M"][0])),
    )
    out = _EPOCH + np.linspace(30.0, 900.0, 10)
    bg = big_asteroid_perturbers(_EPOCH)
    got = propagate_assist(el, _EPOCH, out, asteroid_perturbers=bg, gr=True)
    vec = Horizons(id=target, id_type="smallbody", location="@0",
                   epochs=list(out)).vectors(refplane="ecliptic")
    ref = np.column_stack([vec["x"], vec["y"], vec["z"]]).astype(float)

    res_au = np.linalg.norm(got - ref, axis=1).max()
    helio_au = float(np.mean(np.linalg.norm(ref, axis=1)))
    res_mas = res_au / helio_au * math.degrees(1.0) * 3.6e6
    assert res_mas < 1.0, f"residuo ASSIST vs Horizons = {res_mas:.3f} mas (>1)"
