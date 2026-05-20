"""Tests for src/astrometry/forward_model.py.

Covers:
- residuals_mas: RA wrap at 0/360, sign convention, Dec residuals
- forward_model: shape validation (no N-body needed), RA/Dec range
- forward_model with zero perturber mass: smoke test with real N-body
"""

from __future__ import annotations

import numpy as np
import pytest

from src.astrometry.forward_model import forward_model, residuals_mas


# ---------------------------------------------------------------------------
# residuals_mas — pure geometry, no N-body
# ---------------------------------------------------------------------------

class TestResidualsMas:
    def test_zero_residual(self):
        ra = np.array([10.0, 180.0, 359.0])
        dec = np.array([-30.0, 0.0, 45.0])
        dra, ddec = residuals_mas(ra, dec, ra, dec)
        np.testing.assert_allclose(dra,  0.0, atol=1e-9)
        np.testing.assert_allclose(ddec, 0.0, atol=1e-9)

    def test_dec_residual_sign(self):
        """obs_dec > pred_dec → positive ddec."""
        ra  = np.array([45.0])
        obs_dec  = np.array([10.5])
        pred_dec = np.array([10.0])
        _, ddec = residuals_mas(ra, obs_dec, ra, pred_dec)
        assert ddec[0] > 0.0

    def test_dec_residual_magnitude(self):
        """1° Dec difference → 3600000 mas."""
        ra  = np.array([0.0])
        obs_dec  = np.array([1.0])
        pred_dec = np.array([0.0])
        _, ddec = residuals_mas(ra, obs_dec, ra, pred_dec)
        np.testing.assert_allclose(ddec, [3_600_000.0], rtol=1e-9)

    def test_ra_wrap_near_zero(self):
        """RA wrap: obs=0.001°, pred=359.999° → small positive residual ~2mas×cos(dec)."""
        dec = np.array([0.0])
        obs_ra   = np.array([0.001])
        pred_ra  = np.array([359.999])
        dra, _ = residuals_mas(obs_ra, dec, pred_ra, dec)
        # (0.001 - 359.999 + 360) = 0.002° * cos(0) * 3600000 = 7200 mas
        np.testing.assert_allclose(dra, [7200.0], rtol=1e-6)

    def test_ra_wrap_near_360(self):
        """Symmetric RA wrap: obs near 360, pred near 0."""
        dec = np.array([0.0])
        obs_ra  = np.array([359.999])
        pred_ra = np.array([0.001])
        dra, _ = residuals_mas(obs_ra, dec, pred_ra, dec)
        np.testing.assert_allclose(dra, [-7200.0], rtol=1e-6)

    def test_ra_residual_cos_dec_scaling(self):
        """RA residual is projected by cos(dec)."""
        obs_ra  = np.array([10.1, 10.1])
        pred_ra = np.array([10.0, 10.0])
        dec0    = np.array([0.0, 60.0])
        dra0, _ = residuals_mas(obs_ra, dec0[:1], pred_ra, dec0[:1])
        dra60, _ = residuals_mas(obs_ra, dec0[1:], pred_ra, dec0[1:])
        np.testing.assert_allclose(dra60[0] / dra0[0], np.cos(np.radians(60.0)), rtol=1e-6)

    def test_output_shape(self):
        n = 20
        ra  = np.random.default_rng(0).uniform(0, 360, n)
        dec = np.random.default_rng(1).uniform(-90, 90, n)
        dra, ddec = residuals_mas(ra, dec, ra + 0.01, dec + 0.01)
        assert dra.shape == (n,)
        assert ddec.shape == (n,)


# ---------------------------------------------------------------------------
# forward_model — shape/type validation (no full N-body call needed)
# ---------------------------------------------------------------------------

_DUMMY_ELEMENTS = {
    "a_au": 2.767, "e": 0.079, "i_deg": 10.6,
    "Omega_deg": 80.3, "omega_deg": 73.1, "M_deg": 77.4,
    "epoch_jd": 2457200.5,
}


class TestForwardModelValidation:
    def test_wrong_gaia_xyz_shape_raises(self):
        """gaia_xyz_bary with shape (N, 2) instead of (N, 3) raises ValueError."""
        obs_jd  = np.array([2457200.5, 2457201.5])
        gaia_bad = np.zeros((2, 2))   # wrong: should be (N, 3)
        with pytest.raises(ValueError, match="gaia_xyz_bary shape"):
            forward_model(
                _DUMMY_ELEMENTS, _DUMMY_ELEMENTS, 0.0,
                obs_jd, gaia_bad,
            )

    def test_mismatched_lengths_raise(self):
        """obs_jd length != gaia_xyz rows raises ValueError."""
        obs_jd   = np.array([2457200.5, 2457201.5, 2457202.5])
        gaia_ok  = np.zeros((2, 3))   # 2 rows, but 3 epochs
        with pytest.raises(ValueError, match="gaia_xyz_bary shape"):
            forward_model(
                _DUMMY_ELEMENTS, _DUMMY_ELEMENTS, 0.0,
                obs_jd, gaia_ok,
            )


# ---------------------------------------------------------------------------
# forward_model — integration (requires REBOUND)
# ---------------------------------------------------------------------------

@pytest.mark.slow
class TestForwardModelIntegration:
    """Integration-level tests requiring REBOUND.  Slow but hermetic (no network)."""

    # Ceres-like elements at epoch 2457200.5
    _CERES = {
        "a_au": 2.767, "e": 0.079, "i_deg": 10.6,
        "Omega_deg": 80.3, "omega_deg": 73.1, "M_deg": 77.4,
        "epoch_jd": 2457200.5,
    }
    # Vesta-like perturber
    _VESTA = {
        "a_au": 2.361, "e": 0.089, "i_deg": 7.1,
        "Omega_deg": 103.9, "omega_deg": 151.2, "M_deg": 122.6,
        "epoch_jd": 2457200.5,
    }
    # Gaia at 1 AU from Sun along +X
    _GAIA_XYZ = np.array([[1.0, 0.0, 0.0]])

    def test_output_ra_in_range(self):
        """Predicted RA is always in [0, 360)."""
        obs_jd = np.linspace(2457200.5, 2457210.5, 5)
        gaia   = np.tile(self._GAIA_XYZ, (5, 1))
        ra, dec = forward_model(self._CERES, self._VESTA, 0.0, obs_jd, gaia)
        assert np.all(ra >= 0.0)
        assert np.all(ra < 360.0)

    def test_output_dec_in_range(self):
        """Predicted Dec is always in [-90, 90]."""
        obs_jd = np.linspace(2457200.5, 2457210.5, 5)
        gaia   = np.tile(self._GAIA_XYZ, (5, 1))
        ra, dec = forward_model(self._CERES, self._VESTA, 0.0, obs_jd, gaia)
        assert np.all(dec >= -90.0)
        assert np.all(dec <= 90.0)

    def test_output_shape(self):
        n = 7
        obs_jd = np.linspace(2457200.5, 2457207.5, n)
        gaia   = np.tile(self._GAIA_XYZ, (n, 1))
        ra, dec = forward_model(self._CERES, self._VESTA, 0.0, obs_jd, gaia)
        assert ra.shape == (n,)
        assert dec.shape == (n,)

    def test_zero_vs_nonzero_mass_differ(self):
        """Non-zero perturber mass shifts predicted positions relative to M=0."""
        obs_jd = np.linspace(2457200.5, 2457230.5, 10)
        gaia   = np.tile(self._GAIA_XYZ, (10, 1))
        ra0,  dec0  = forward_model(self._CERES, self._VESTA, 0.0,    obs_jd, gaia)
        ra1,  dec1  = forward_model(self._CERES, self._VESTA, 1.0e20, obs_jd, gaia)
        # With a tiny mass most positions should be identical; with 1e20 kg there
        # should be at least some difference somewhere
        assert not (np.allclose(ra0, ra1) and np.allclose(dec0, dec1))

    def test_large_mass_changes_trajectory(self):
        """An unphysically large mass (10^24 kg) must produce large residuals."""
        obs_jd = np.linspace(2457200.5, 2457300.5, 20)
        gaia   = np.tile(self._GAIA_XYZ, (20, 1))
        ra0, dec0 = forward_model(self._CERES, self._VESTA, 0.0,    obs_jd, gaia)
        ra1, dec1 = forward_model(self._CERES, self._VESTA, 1.0e24, obs_jd, gaia)
        dra, ddec = residuals_mas(ra0, dec0, ra1, dec1)
        rms = float(np.sqrt(np.mean(dra**2 + ddec**2)))
        assert rms > 1.0   # must shift by > 1 mas somewhere
