"""Tests for the joint orbit-drift + mass residual model."""

from __future__ import annotations

import numpy as np

from src.mass import forward_model_joint as fmj


def _elements() -> dict:
    return {
        "a_au": 2.5,
        "e": 0.1,
        "i_deg": 5.0,
        "Omega_deg": 359.8,
        "omega_deg": 0.2,
        "M_deg": 359.0,
        "epoch_jd": 2457200.5,
    }


def _obs(n: int = 4) -> fmj.GaiaObservationBundle:
    return fmj.GaiaObservationBundle(
        jd_tdb=np.linspace(2457200.5, 2457210.5, n),
        gaia_xyz_bary=np.zeros((n, 3)),
        ra_deg=np.full(n, 10.0),
        dec_deg=np.zeros(n),
        position_angle_scan_deg=np.full(n, 90.0),
        ra_error_systematic_mas=np.full(n, 1.0),
        dec_error_systematic_mas=np.full(n, 1.0),
        ra_dec_correlation_systematic=np.zeros(n),
        ra_error_random_mas=np.zeros(n),
        dec_error_random_mas=np.zeros(n),
        ra_dec_correlation_random=np.zeros(n),
    )


def test_apply_target_deltas_uses_relative_a_and_wraps_angles() -> None:
    params = np.array([18.0, 1e-3, 5e-4, -0.1, 0.5, -0.5, 2.0])

    adjusted = fmj.apply_target_deltas(_elements(), params)

    np.testing.assert_allclose(adjusted["a_au"], 2.5025)
    np.testing.assert_allclose(adjusted["e"], 0.1005)
    np.testing.assert_allclose(adjusted["i_deg"], 4.9)
    np.testing.assert_allclose(adjusted["Omega_deg"], 0.3)
    np.testing.assert_allclose(adjusted["omega_deg"], 359.7)
    np.testing.assert_allclose(adjusted["M_deg"], 1.0)


def test_residuals_joint_zero_when_forward_model_matches_observations(monkeypatch) -> None:
    obs = _obs()

    def fake_forward_model(*args, **kwargs):
        return obs.ra_deg.copy(), obs.dec_deg.copy()

    monkeypatch.setattr(fmj, "forward_model", fake_forward_model)
    params = np.array([18.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    residuals = fmj.residuals_joint(params, _elements(), _elements(), obs, include_planets=("sun",))

    assert residuals.shape == (len(obs.jd_tdb) + 6,)
    np.testing.assert_allclose(residuals, 0.0)


def test_residuals_joint_appends_normalized_orbit_priors(monkeypatch) -> None:
    obs = _obs()

    def fake_forward_model(*args, **kwargs):
        return obs.ra_deg.copy(), obs.dec_deg.copy()

    monkeypatch.setattr(fmj, "forward_model", fake_forward_model)
    params = np.array([18.0, 2e-4, 5e-4, 0.05, 0.2, 0.2, 0.5])

    residuals = fmj.residuals_joint(params, _elements(), _elements(), obs, include_planets=("sun",))

    np.testing.assert_allclose(residuals[: len(obs.jd_tdb)], 0.0)
    np.testing.assert_allclose(residuals[len(obs.jd_tdb) :], 1.0)


def test_fit_joint_recovers_synthetic_log_mass(monkeypatch) -> None:
    obs = _obs(n=6)
    true_log_mass = 18.5
    scale_deg = 1e-6

    def fake_forward_model(
        target_elements,
        perturber_elements,
        perturber_mass_kg,
        obs_jd_tdb,
        gaia_xyz_bary,
        **kwargs,
    ):
        del target_elements, perturber_elements, gaia_xyz_bary, kwargs
        log_mass = np.log10(perturber_mass_kg)
        ra_pred = np.full(len(obs_jd_tdb), 10.0 + (log_mass - true_log_mass) * scale_deg)
        return ra_pred, np.zeros(len(obs_jd_tdb))

    monkeypatch.setattr(fmj, "forward_model", fake_forward_model)

    result = fmj.fit_joint(
        _elements(),
        _elements(),
        obs,
        initial_log10_mass=17.0,
        include_planets=("sun",),
        max_nfev=50,
    )

    assert result.success
    assert abs(result.x[0] - true_log_mass) < 1e-4
    np.testing.assert_allclose(result.x[1:], 0.0, atol=1e-8)
