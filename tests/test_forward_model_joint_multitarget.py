"""Tests for the multi-target joint orbit-drift + mass residual model."""

from __future__ import annotations

import numpy as np

from src.mass import forward_model_joint_multitarget as fmm
from src.mass.forward_model_joint import DEFAULT_PRIORS, GaiaObservationBundle


def _elements(a_au: float = 2.5) -> dict:
    return {
        "a_au": a_au,
        "e": 0.1,
        "i_deg": 5.0,
        "Omega_deg": 359.8,
        "omega_deg": 0.2,
        "M_deg": 359.0,
        "epoch_jd": 2457200.5,
    }


def _obs(n: int = 4, ra_deg: float = 10.0) -> GaiaObservationBundle:
    return GaiaObservationBundle(
        jd_tdb=np.linspace(2457200.5, 2457210.5, n),
        gaia_xyz_bary=np.zeros((n, 3)),
        ra_deg=np.full(n, ra_deg),
        dec_deg=np.zeros(n),
        position_angle_scan_deg=np.full(n, 90.0),
        ra_error_systematic_mas=np.full(n, 1.0),
        dec_error_systematic_mas=np.full(n, 1.0),
        ra_dec_correlation_systematic=np.zeros(n),
        ra_error_random_mas=np.zeros(n),
        dec_error_random_mas=np.zeros(n),
        ra_dec_correlation_random=np.zeros(n),
    )


def _bundles(n_targets: int = 3, n_obs: int = 4) -> list[fmm.TargetBundle]:
    return [
        fmm.TargetBundle(
            target_number=1000 + i,
            elements=_elements(a_au=2.5 + 0.05 * i),
            obs=_obs(n_obs, ra_deg=10.0 + i),
        )
        for i in range(n_targets)
    ]


def test_target_slice_extracts_mass_and_per_target_deltas() -> None:
    params = np.concatenate([[18.0], np.arange(6 * 3, dtype=float)])

    s0 = fmm._target_slice(params, 0)
    s1 = fmm._target_slice(params, 1)
    s2 = fmm._target_slice(params, 2)

    np.testing.assert_allclose(s0, [18.0, 0, 1, 2, 3, 4, 5])
    np.testing.assert_allclose(s1, [18.0, 6, 7, 8, 9, 10, 11])
    np.testing.assert_allclose(s2, [18.0, 12, 13, 14, 15, 16, 17])


def test_make_bounds_replicates_delta_bounds_per_target() -> None:
    lo, hi = fmm.make_bounds(3, DEFAULT_PRIORS)

    assert lo.shape == (1 + 18,)
    assert hi.shape == (1 + 18,)
    assert lo[0] == DEFAULT_PRIORS.log10_mass_bounds[0]
    assert hi[0] == DEFAULT_PRIORS.log10_mass_bounds[1]
    for i in range(3):
        s = 1 + 6 * i
        np.testing.assert_allclose(lo[s : s + 6], DEFAULT_PRIORS.lower_bounds[1:])
        np.testing.assert_allclose(hi[s : s + 6], DEFAULT_PRIORS.upper_bounds[1:])


def test_prior_residuals_multitarget_normalises_each_target_block() -> None:
    sigma = DEFAULT_PRIORS.sigma_vector
    deltas = np.concatenate([sigma, 2.0 * sigma, -sigma])
    params = np.concatenate([[18.0], deltas])

    out = fmm.prior_residuals_multitarget(params, 3, DEFAULT_PRIORS)

    np.testing.assert_allclose(out[0:6], 1.0)
    np.testing.assert_allclose(out[6:12], 2.0)
    np.testing.assert_allclose(out[12:18], -1.0)


def test_residuals_joint_multitarget_zero_when_forward_model_matches(monkeypatch) -> None:
    bundles = _bundles(n_targets=3, n_obs=4)

    def fake_forward_model(target_elements, *args, **kwargs):
        del args, kwargs
        n = len(bundles[0].obs.jd_tdb)
        ra = float(target_elements.get("_ra_marker", 10.0))
        return np.full(n, ra), np.zeros(n)

    for bundle in bundles:
        bundle.elements["_ra_marker"] = float(bundle.obs.ra_deg[0])

    monkeypatch.setattr(fmm, "forward_model", fake_forward_model)
    params = np.concatenate([[18.0], np.zeros(6 * len(bundles))])

    res = fmm.residuals_joint_multitarget(
        params,
        bundles,
        _elements(),
        include_planets=("sun",),
    )

    n_obs_total = sum(len(b.obs.jd_tdb) for b in bundles)
    n_priors = 6 * len(bundles)
    assert res.shape == (n_obs_total + n_priors,)
    np.testing.assert_allclose(res, 0.0)


def test_fit_joint_multitarget_recovers_synthetic_log_mass(monkeypatch) -> None:
    """Two synthetic targets with the same true mass: the joint fit must recover it."""
    bundles = _bundles(n_targets=2, n_obs=6)
    true_log_mass = 18.5
    scale_deg_per_log10 = 1e-6

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
        offset = (log_mass - true_log_mass) * scale_deg_per_log10
        ra = np.full(len(obs_jd_tdb), 10.0 + offset)
        return ra, np.zeros(len(obs_jd_tdb))

    for bundle in bundles:
        bundle.obs.ra_deg[:] = 10.0
    monkeypatch.setattr(fmm, "forward_model", fake_forward_model)

    result = fmm.fit_joint_multitarget(
        bundles,
        _elements(),
        initial_log10_mass=17.0,
        include_planets=("sun",),
        max_nfev=80,
    )

    assert result.success
    assert abs(result.x[0] - true_log_mass) < 1e-4
    np.testing.assert_allclose(result.x[1:], 0.0, atol=1e-8)


def test_fit_joint_multitarget_recovers_mass_with_per_target_da(monkeypatch) -> None:
    """Per-target time-constant orbit offset (mimics da_rel) + time-varying mass signal.

    The mass signal is shaped as a Heaviside step centred at each target's
    encounter epoch — emulating a gravitational deflection — so the optimiser
    can distinguish it from the time-constant per-target orbit offset.
    """
    bundles = _bundles(n_targets=3, n_obs=10)
    true_log_mass = 19.0
    base_a = [b.elements["a_au"] for b in bundles]
    true_da_rel = [0.0, 5e-6, -3e-6]
    encounter_jd = [b.obs.jd_tdb[len(b.obs.jd_tdb) // 2] for b in bundles]
    scale_M = 1e-6
    scale_a = 1.0
    base_ra = 10.0

    def fake_forward_model(
        target_elements,
        perturber_elements,
        perturber_mass_kg,
        obs_jd_tdb,
        gaia_xyz_bary,
        **kwargs,
    ):
        del perturber_elements, gaia_xyz_bary, kwargs
        a = float(target_elements["a_au"])
        idx = int(np.argmin([abs(a - b * (1.0 + 5e-6)) for b in base_a]))
        a0 = base_a[idx]
        log_mass = np.log10(perturber_mass_kg)
        kick = np.where(obs_jd_tdb > encounter_jd[idx], 1.0, 0.0)
        mass_signal = (log_mass - true_log_mass) * scale_M * kick
        orbit_signal = (a / a0 - 1.0) * scale_a
        ra = base_ra + mass_signal + orbit_signal
        return ra, np.zeros(len(obs_jd_tdb))

    for i, bundle in enumerate(bundles):
        bundle.obs.ra_deg[:] = base_ra + true_da_rel[i] * scale_a

    monkeypatch.setattr(fmm, "forward_model", fake_forward_model)

    result = fmm.fit_joint_multitarget(
        bundles,
        _elements(),
        initial_log10_mass=17.5,
        include_planets=("sun",),
        max_nfev=400,
    )

    assert result.success
    assert abs(result.x[0] - true_log_mass) < 1e-3
    for i, da_truth in enumerate(true_da_rel):
        da_fit = result.x[1 + 6 * i]
        assert abs(da_fit - da_truth) < 5e-7, f"target {i}: da_fit={da_fit} vs truth={da_truth}"


def test_fit_joint_multitarget_profiled_recovers_injected_mass(monkeypatch) -> None:
    """Profiled optimiser recovers a shared mass from a clean time-varying signal.

    Mirrors the closing-loop test on real data: the perturber mass enters as a
    Heaviside deflection after each encounter; the profiled outer/inner split
    must land on the injected log10_M. Returns a ``ProfiledFitResult`` whose
    ``x`` shares the ``(1 + 6N)`` layout of the joint fit.
    """
    bundles = _bundles(n_targets=3, n_obs=10)
    true_log_mass = 19.0
    encounter_jd = [b.obs.jd_tdb[len(b.obs.jd_tdb) // 2] for b in bundles]
    scale_M = 1.0
    base_ra = 10.0

    def fake_forward_model(
        target_elements,
        perturber_elements,
        perturber_mass_kg,
        obs_jd_tdb,
        gaia_xyz_bary,
        **kwargs,
    ):
        del perturber_elements, gaia_xyz_bary, kwargs
        a = float(target_elements["a_au"])
        idx = int(np.argmin([abs(a - b.elements["a_au"]) for b in bundles]))
        log_mass = np.log10(perturber_mass_kg)
        kick = np.where(obs_jd_tdb > encounter_jd[idx], 1.0, 0.0)
        ra = base_ra + (log_mass - true_log_mass) * scale_M * kick
        return ra, np.zeros(len(obs_jd_tdb))

    # Truth observations: generated at the true mass with zero deltas.
    for bundle in bundles:
        bundle.obs.ra_deg[:] = base_ra

    monkeypatch.setattr(fmm, "forward_model", fake_forward_model)

    result = fmm.fit_joint_multitarget_profiled(
        bundles,
        _elements(),
        include_planets=("sun",),
        inner_max_nfev=400,
    )

    assert result.x.shape == (1 + 6 * len(bundles),)
    assert abs(result.log10_mass - true_log_mass) < 1e-2
    assert abs(result.x[0] - true_log_mass) < 1e-2
