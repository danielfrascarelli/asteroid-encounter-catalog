"""Tests for the 2D Mahalanobis residual implementation."""

from __future__ import annotations

import numpy as np
import pytest

from src.mass.likelihood_al import mahalanobis_residuals_2d


def _const(value: float, n: int) -> np.ndarray:
    return np.full(n, value, dtype=float)


def test_diagonal_no_correlation_chi2_matches_closed_form() -> None:
    """sigma_RA = sigma_Dec, rho = 0: chi2 = (dra/sigma)^2 + (ddec/sigma)^2."""
    n = 4
    sigma = 5.0
    rng = np.random.default_rng(seed=0)
    dra = rng.normal(0.0, sigma, size=n)
    ddec = rng.normal(0.0, sigma, size=n)

    whitened, chi2 = mahalanobis_residuals_2d(
        dra,
        ddec,
        ra_err_sys=_const(sigma, n),
        dec_err_sys=_const(sigma, n),
        corr_sys=_const(0.0, n),
        ra_err_rand=_const(0.0, n),
        dec_err_rand=_const(0.0, n),
        corr_rand=_const(0.0, n),
    )
    expected = (dra / sigma) ** 2 + (ddec / sigma) ** 2
    np.testing.assert_allclose(chi2, expected, rtol=1e-10)
    # Whitened vector squared norm should equal chi2.
    pairs = whitened.reshape(-1, 2)
    np.testing.assert_allclose(np.sum(pairs * pairs, axis=1), chi2, rtol=1e-10)


def test_explicit_chi2_value_known() -> None:
    """delta=(3,4), sigma_RA=1, sigma_Dec=2, rho=0 -> chi2 = 9 + 4 = 13."""
    dra = np.array([3.0])
    ddec = np.array([4.0])
    _, chi2 = mahalanobis_residuals_2d(
        dra,
        ddec,
        ra_err_sys=np.array([1.0]),
        dec_err_sys=np.array([2.0]),
        corr_sys=np.array([0.0]),
        ra_err_rand=np.array([0.0]),
        dec_err_rand=np.array([0.0]),
        corr_rand=np.array([0.0]),
    )
    np.testing.assert_allclose(chi2, [13.0], rtol=1e-12)


def test_correlation_matches_closed_form() -> None:
    """delta=(1,1), sigma_RA=sigma_Dec=1, rho=0.5 -> chi2 = (delta^T Sigma^-1 delta)."""
    dra = np.array([1.0])
    ddec = np.array([1.0])
    rho = 0.5
    sigma = 1.0
    # Sigma = [[1, 0.5],[0.5, 1]]; det = 0.75; inv = (1/0.75)[[1,-0.5],[-0.5,1]]
    # chi2 = (1/0.75) * (1 - 2*0.5 + 1) = (1/0.75) * 1 = 1.333...
    expected = (1.0 - 2.0 * rho + 1.0) / (1.0 - rho * rho)
    _, chi2 = mahalanobis_residuals_2d(
        dra,
        ddec,
        ra_err_sys=np.array([sigma]),
        dec_err_sys=np.array([sigma]),
        corr_sys=np.array([rho]),
        ra_err_rand=np.array([0.0]),
        dec_err_rand=np.array([0.0]),
        corr_rand=np.array([0.0]),
    )
    np.testing.assert_allclose(chi2, [expected], rtol=1e-12)


def test_systematic_plus_random_adds_in_quadrature_diagonal() -> None:
    """Two diagonal contributions sum: sigma_total^2 = sigma_sys^2 + sigma_rand^2."""
    dra = np.array([6.0])
    ddec = np.array([0.0])
    sigma_sys = 3.0
    sigma_rand = 4.0
    # total sigma^2 in RA = 9 + 16 = 25 -> chi2 = 36 / 25 = 1.44
    _, chi2 = mahalanobis_residuals_2d(
        dra,
        ddec,
        ra_err_sys=np.array([sigma_sys]),
        dec_err_sys=np.array([sigma_sys]),
        corr_sys=np.array([0.0]),
        ra_err_rand=np.array([sigma_rand]),
        dec_err_rand=np.array([sigma_rand]),
        corr_rand=np.array([0.0]),
    )
    np.testing.assert_allclose(chi2, [36.0 / 25.0], rtol=1e-12)


def test_degenerate_rho_falls_back_safely(caplog: pytest.LogCaptureFixture) -> None:
    """rho clipped to (-0.9999, 0.9999) but values near 1 must not produce NaNs."""
    dra = np.array([2.0])
    ddec = np.array([1.0])
    with caplog.at_level("WARNING", logger="src.mass.likelihood_al"):
        whitened, chi2 = mahalanobis_residuals_2d(
            dra,
            ddec,
            ra_err_sys=np.array([1.0]),
            dec_err_sys=np.array([1.0]),
            # The clamp at +/-0.9999 plus the determinant floor protects us:
            # det = 1 - 0.9999^2 ~ 2e-4 > 1e-8 floor, so no fallback triggers,
            # but the result must remain finite.
            corr_sys=np.array([0.99999]),
            ra_err_rand=np.array([0.0]),
            dec_err_rand=np.array([0.0]),
            corr_rand=np.array([0.0]),
        )
    assert np.all(np.isfinite(whitened))
    assert np.all(np.isfinite(chi2))


def test_whitened_norm_equals_chi2_with_correlation() -> None:
    """||whitened_i||^2 must equal chi2_per_obs[i] even with non-zero rho."""
    n = 5
    rng = np.random.default_rng(seed=42)
    dra = rng.normal(0.0, 3.0, size=n)
    ddec = rng.normal(0.0, 5.0, size=n)
    sigma_ra = rng.uniform(1.0, 4.0, size=n)
    sigma_dec = rng.uniform(1.0, 4.0, size=n)
    rho = rng.uniform(-0.7, 0.7, size=n)
    whitened, chi2 = mahalanobis_residuals_2d(
        dra,
        ddec,
        ra_err_sys=sigma_ra,
        dec_err_sys=sigma_dec,
        corr_sys=rho,
        ra_err_rand=_const(0.0, n),
        dec_err_rand=_const(0.0, n),
        corr_rand=_const(0.0, n),
    )
    pairs = whitened.reshape(-1, 2)
    norms = np.sum(pairs * pairs, axis=1)
    np.testing.assert_allclose(norms, chi2, rtol=1e-10)


def test_mahalanobis_reduces_to_al_when_ac_axis_dominates() -> None:
    """When the across-scan eigenvalue blows up, Mahalanobis 2D ~ AL-only."""
    # Build a covariance whose principal axis is aligned with RA and is the
    # noisy (across-scan) direction; the Dec axis is then the precise AL one.
    # That happens when the scan direction PA = 90 deg (e_AL = (1, 0)) is
    # rotated by 90 deg w.r.t. the eigenvector frame -- the easiest equivalent
    # is to compare against the explicit AL projection used in Stage 1.
    n = 3
    rng = np.random.default_rng(seed=7)
    dra = rng.normal(0.0, 1.0, size=n)
    ddec = rng.normal(0.0, 1.0, size=n)
    sigma_al = 0.3
    sigma_ac = 30.0  # 100x noisier
    # Diagonal covariance in the (AL, AC) frame -> set RA = AL axis, Dec = AC
    # axis: PA = 90 deg in the al_residuals_and_weights convention,
    # so r_AL = dra and Sigma = diag(sigma_AL^2, sigma_AC^2).
    _, chi2_2d = mahalanobis_residuals_2d(
        dra,
        ddec,
        ra_err_sys=_const(sigma_al, n),
        dec_err_sys=_const(sigma_ac, n),
        corr_sys=_const(0.0, n),
        ra_err_rand=_const(0.0, n),
        dec_err_rand=_const(0.0, n),
        corr_rand=_const(0.0, n),
    )
    # AL-only chi2 ignoring the across-scan axis
    chi2_al = (dra / sigma_al) ** 2
    extra = (ddec / sigma_ac) ** 2  # tiny because sigma_ac is large
    np.testing.assert_allclose(chi2_2d, chi2_al + extra, rtol=1e-12)
    # The 2D result should be within ~1% of the AL-only result because the
    # AC term contributes (1/100)^2 * O(1) ~ 1e-4 per obs.
    assert np.max(np.abs(chi2_2d - chi2_al)) < 1e-2
