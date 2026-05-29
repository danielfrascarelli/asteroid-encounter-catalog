"""Per-observation 2D Mahalanobis residuals for Gaia DR3 SSO astrometry.

The Gaia DR3 ``sso_observation`` table provides, for each transit, both
random and systematic uncertainties on (RA, Dec) with a correlation
coefficient. Stage 1 of the mass layer projected those errors onto the
along-scan axis and discarded the across-scan component (see
:mod:`src.mass.forward_model_joint`). Stage 2 keeps **both** tangential
components and weights them with the full 2x2 covariance matrix, giving
a proper Mahalanobis chi-squared per observation.

Notation
--------
* ``dra``, ``ddec`` -- tangential residuals in mas
  (``(RA_obs - RA_pred) cos Dec_pred`` and ``Dec_obs - Dec_pred``).
* ``sigma_ra``, ``sigma_dec`` -- 1-sigma uncertainties in RA*, Dec
  (mas) for one of the two error families (systematic or random).
* ``rho`` -- correlation coefficient in [-1, 1] for the same family.

The total covariance is the sum of the systematic and random matrices.
"""

from __future__ import annotations

import logging

import numpy as np

logger = logging.getLogger(__name__)

# Numerical floor for the determinant of the per-observation covariance
# before we fall back to a diagonal approximation. The unit is mas^4.
_DET_FLOOR = 1e-8
# Hard clamp for |rho|; the Gaia archive guarantees |rho| <= 1 but bad
# exports occasionally overshoot, which would yield a non-positive-definite
# covariance and break the Cholesky factorisation downstream.
_RHO_CLAMP = 0.9999


def _build_covariance(
    sigma_ra: np.ndarray,
    sigma_dec: np.ndarray,
    rho: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (Sigma_xx, Sigma_xy, Sigma_yy) component arrays for one family."""
    rho_clipped = np.clip(rho, -_RHO_CLAMP, _RHO_CLAMP)
    sxx = sigma_ra * sigma_ra
    syy = sigma_dec * sigma_dec
    sxy = rho_clipped * sigma_ra * sigma_dec
    return sxx, sxy, syy


def _invert_2x2(
    sxx: np.ndarray,
    sxy: np.ndarray,
    syy: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Invert a stack of symmetric 2x2 matrices.

    Returns
    -------
    inv_xx, inv_xy, inv_yy, fallback_mask
        ``inv_*`` are the entries of Sigma^-1 (symmetric). ``fallback_mask``
        marks observations where the closed-form inverse was numerically
        unsafe and a diagonal fallback was used instead.
    """
    det = sxx * syy - sxy * sxy
    fallback_mask = det < _DET_FLOOR
    safe_det = np.where(fallback_mask, 1.0, det)
    inv_xx = np.where(fallback_mask, 1.0 / np.maximum(sxx, _DET_FLOOR), syy / safe_det)
    inv_xy = np.where(fallback_mask, 0.0, -sxy / safe_det)
    inv_yy = np.where(fallback_mask, 1.0 / np.maximum(syy, _DET_FLOOR), sxx / safe_det)
    if np.any(fallback_mask):
        logger.warning(
            "Mahalanobis 2D: %d/%d observations fell back to a diagonal covariance " "(det < %g).",
            int(np.sum(fallback_mask)),
            int(fallback_mask.size),
            _DET_FLOOR,
        )
    return inv_xx, inv_xy, inv_yy, fallback_mask


def _cholesky_upper_2x2(
    a: np.ndarray,
    b: np.ndarray,
    c: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return the upper-triangular Cholesky factor L of [[a, b], [b, c]].

    Uses the convention
        L = [[L11, L12],
             [0,   L22]]
    so that ``L^T L = Sigma^-1`` when ``[[a, b], [b, c]] = Sigma^-1``.
    """
    a_safe = np.maximum(a, _DET_FLOOR)
    L11 = np.sqrt(a_safe)
    L12 = b / L11
    diag = c - L12 * L12
    L22 = np.sqrt(np.maximum(diag, _DET_FLOOR))
    return L11, L12, L22


def mahalanobis_residuals_2d(
    dra_mas: np.ndarray,
    ddec_mas: np.ndarray,
    ra_err_sys: np.ndarray,
    dec_err_sys: np.ndarray,
    corr_sys: np.ndarray,
    ra_err_rand: np.ndarray,
    dec_err_rand: np.ndarray,
    corr_rand: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Whiten per-observation 2D residuals using the Gaia (RA, Dec) covariance.

    Parameters
    ----------
    dra_mas, ddec_mas
        Tangential residuals (mas), shape ``(N,)``.
        ``dra_mas`` is expected to be ``(RA_obs - RA_pred) * cos(Dec_pred)``.
    ra_err_sys, dec_err_sys, corr_sys
        Systematic uncertainties and correlation, shape ``(N,)``.
    ra_err_rand, dec_err_rand, corr_rand
        Random uncertainties and correlation, shape ``(N,)``.

    Returns
    -------
    whitened : ndarray, shape ``(2N,)``
        Vector ``L @ delta`` per observation, stacked. The sum of squared
        entries equals the total Mahalanobis chi-squared.
    chi2_per_obs : ndarray, shape ``(N,)``
        Per-observation Mahalanobis chi-squared values.
    """
    dra = np.asarray(dra_mas, dtype=float)
    ddec = np.asarray(ddec_mas, dtype=float)
    n = dra.shape[0]
    if ddec.shape != (n,):
        raise ValueError(f"ddec_mas shape {ddec.shape} != ({n},)")
    arrays = {
        "ra_err_sys": ra_err_sys,
        "dec_err_sys": dec_err_sys,
        "corr_sys": corr_sys,
        "ra_err_rand": ra_err_rand,
        "dec_err_rand": dec_err_rand,
        "corr_rand": corr_rand,
    }
    casted: dict[str, np.ndarray] = {}
    for name, arr in arrays.items():
        a = np.asarray(arr, dtype=float)
        if a.shape != (n,):
            raise ValueError(f"{name} shape {a.shape} != ({n},)")
        casted[name] = a

    sxx_sys, sxy_sys, syy_sys = _build_covariance(
        casted["ra_err_sys"], casted["dec_err_sys"], casted["corr_sys"]
    )
    sxx_rand, sxy_rand, syy_rand = _build_covariance(
        casted["ra_err_rand"], casted["dec_err_rand"], casted["corr_rand"]
    )
    sxx = sxx_sys + sxx_rand
    sxy = sxy_sys + sxy_rand
    syy = syy_sys + syy_rand

    inv_xx, inv_xy, inv_yy, _ = _invert_2x2(sxx, sxy, syy)

    # chi^2 = delta^T Sigma^-1 delta with delta = (dra, ddec)
    chi2_per_obs = dra * dra * inv_xx + 2.0 * dra * ddec * inv_xy + ddec * ddec * inv_yy

    # Whitened residual via the upper-triangular Cholesky of Sigma^-1:
    #   r1 = L11 * dra + L12 * ddec
    #   r2 = L22 * ddec
    # so that r1^2 + r2^2 == chi2_per_obs.
    L11, L12, L22 = _cholesky_upper_2x2(inv_xx, inv_xy, inv_yy)
    r1 = L11 * dra + L12 * ddec
    r2 = L22 * ddec
    whitened = np.empty(2 * n, dtype=float)
    whitened[0::2] = r1
    whitened[1::2] = r2
    return whitened, chi2_per_obs


__all__ = [
    "mahalanobis_residuals_2d",
]
