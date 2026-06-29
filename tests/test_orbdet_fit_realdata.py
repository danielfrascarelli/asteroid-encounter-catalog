"""Unit tests for the pure (network-free) helpers of ``orbdet_fit_realdata``.

The end-to-end fit needs the Gaia TAP and ASSIST ephemerides, but the IO glue —
residual splitting, transit masking, perturber-name resolution and reading the
literature metadata from the validation CSV — is pure Python and tested here.
"""

from __future__ import annotations

import numpy as np
import pytest

from scripts.mass.orbdet_fit_realdata import (
    _ephem_name_for_perturber,
    _mask_target,
    _read_perturber_meta,
    _split_by_targets,
)
from src.orbdet.kepler import KeplerElements
from src.orbdet.mass_determination import TargetObservations


def _toy_target(n: int) -> TargetObservations:
    el = KeplerElements(a=2.5, e=0.1, i=0.1, Omega=0.2, omega=0.3, M=0.4)
    return TargetObservations(
        initial_elements=el,
        obs_jd_tdb=np.arange(n, dtype=float),
        ra_obs_deg=np.full(n, 10.0),
        dec_obs_deg=np.full(n, 5.0),
        pa_scan_deg=np.full(n, 45.0),
        sigma_al_mas=np.full(n, 0.5),
        gaia_bary_icrs=np.tile(np.array([1.0, 0.0, 0.0]), (n, 1)),
    )


def test_split_by_targets_respects_lengths() -> None:
    targets = [_toy_target(3), _toy_target(5), _toy_target(2)]
    values = np.arange(10, dtype=float)
    parts = _split_by_targets(values, targets)
    assert [p.size for p in parts] == [3, 5, 2]
    np.testing.assert_array_equal(parts[0], [0, 1, 2])
    np.testing.assert_array_equal(parts[1], [3, 4, 5, 6, 7])
    np.testing.assert_array_equal(parts[2], [8, 9])


def test_mask_target_keeps_selected_transits() -> None:
    t = _toy_target(6)
    keep = np.array([True, False, True, True, False, True])
    masked = _mask_target(t, keep)
    assert masked.obs_jd_tdb.size == 4
    np.testing.assert_array_equal(masked.obs_jd_tdb, [0, 2, 3, 5])
    assert masked.gaia_bary_icrs.shape == (4, 3)
    # otras columnas se recortan en paralelo
    assert masked.sigma_al_mas.size == 4
    # los elementos iniciales se preservan
    assert masked.initial_elements.a == t.initial_elements.a


def test_ephem_name_for_big4() -> None:
    assert _ephem_name_for_perturber(1, None) == "Ceres"
    assert _ephem_name_for_perturber(2, "pallas") == "Pallas"
    assert _ephem_name_for_perturber(4, None) == "Vesta"
    assert _ephem_name_for_perturber(10, None) == "Hygiea"


def test_ephem_name_from_csv_name() -> None:
    # Eunomia está en BIG_ASTEROIDS pero no en el mapa Big-4 explícito.
    assert _ephem_name_for_perturber(15, "eunomia") == "Eunomia"


def test_ephem_name_rejects_unknown_perturber() -> None:
    with pytest.raises(ValueError, match="efeméride DE441"):
        _ephem_name_for_perturber(99999, "NotAnEphemAsteroid")


def test_read_perturber_meta(tmp_path) -> None:
    csv = tmp_path / "summary.csv"
    csv.write_text(
        "perturber,perturber_name,target,mass_lit_kg,mass_lit_sigma_kg,literature_source\n"
        "1,Ceres,18937,4.71e20,4e18,DAWN (Park+ 2016)\n"
        "1,Ceres,24836,4.71e20,4e18,DAWN (Park+ 2016)\n"
        "10,Hygiea,1234,8.3e19,1e18,Some Ref\n"
    )
    meta = _read_perturber_meta(csv, 1)
    assert meta["name"] == "Ceres"
    assert meta["mass_lit_kg"] == pytest.approx(4.71e20)
    assert meta["mass_lit_sigma_kg"] == pytest.approx(4e18)
    assert meta["source"] == "DAWN (Park+ 2016)"

    meta10 = _read_perturber_meta(csv, 10)
    assert meta10["name"] == "Hygiea"
    assert meta10["mass_lit_kg"] == pytest.approx(8.3e19)


def test_read_perturber_meta_missing_file(tmp_path) -> None:
    meta = _read_perturber_meta(tmp_path / "nope.csv", 1)
    assert meta == {"name": None, "mass_lit_kg": None, "mass_lit_sigma_kg": None, "source": None}
