"""Tests for release-aware Gaia fetch in scripts/mass/fit_mass_gaia_loo.py.

The TAP call is mocked; we assert the ADQL targets the right table/window and
applies the FPR rejection filter only when the release exposes one.
"""

from unittest.mock import MagicMock, patch

import pandas as pd

from scripts.mass.fit_mass_gaia_loo import fetch_gaia_full
from src.utils.config import GaiaReleaseConfig

_COLS = [
    "number_mp",
    "epoch",
    "ra",
    "dec",
    "ra_error_systematic",
    "dec_error_systematic",
    "ra_dec_correlation_systematic",
    "ra_error_random",
    "dec_error_random",
    "ra_dec_correlation_random",
    "position_angle_scan",
    "x_gaia",
    "y_gaia",
    "z_gaia",
]


def _tap_mock() -> MagicMock:
    job = MagicMock()
    job.get_results.return_value.to_pandas.return_value = pd.DataFrame(
        {c: [0.0, 1.0] for c in _COLS}
    )
    tap = MagicMock()
    tap.launch_job_async.return_value = job
    return tap


def _captured_adql(tap: MagicMock) -> str:
    return tap.launch_job_async.call_args[0][0]


@patch("scripts.mass.fit_mass_gaia_loo.TapPlus")
def test_fetch_defaults_to_dr3(mock_tap_cls: MagicMock) -> None:
    mock_tap_cls.return_value = tap = _tap_mock()
    fetch_gaia_full("http://fake", 1)
    adql = _captured_adql(tap)
    assert "gaiadr3.sso_observation" in adql
    assert "is_rejected" not in adql


@patch("scripts.mass.fit_mass_gaia_loo.TapPlus")
def test_fetch_fpr_table_window_and_reject_filter(mock_tap_cls: MagicMock) -> None:
    mock_tap_cls.return_value = tap = _tap_mock()
    fpr = GaiaReleaseConfig(
        table="gaiafpr.sso_observation",
        epoch_ref_jd_tcb=2455197.5,
        window_start="2014-07-26T00:00:00",
        window_end="2020-01-21T00:00:00",
        mp_max=400000,
        columns_drop=["g_mag"],
        reject_flag_column="is_rejected",
    )
    fetch_gaia_full("http://fake", 1, fpr)
    adql = _captured_adql(tap)
    assert "gaiafpr.sso_observation" in adql
    assert "gaiadr3" not in adql
    # FPR window upper bound ~ JD 2458869 - 2455197.5 ≈ 3671 days since J2010
    assert "AND is_rejected = 'false'" in adql
    assert "3671" in adql or "3672" in adql  # window_end maps to ~3671 days


@patch("scripts.mass.fit_mass_gaia_loo.TapPlus")
def test_fetch_dr3_release_cfg_has_no_reject_filter(mock_tap_cls: MagicMock) -> None:
    mock_tap_cls.return_value = tap = _tap_mock()
    dr3 = GaiaReleaseConfig(
        table="gaiadr3.sso_observation",
        epoch_ref_jd_tcb=2455197.5,
        window_start="2014-07-25T00:00:00",
        window_end="2017-05-28T00:00:00",
        mp_max=160000,
    )
    fetch_gaia_full("http://fake", 1, dr3)
    adql = _captured_adql(tap)
    assert "gaiadr3.sso_observation" in adql
    assert "is_rejected" not in adql
