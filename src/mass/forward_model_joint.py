"""Joint orbit-drift + perturber-mass residual model.

This module keeps the expensive astrometric forward model unchanged and
changes only the optimisation surface: the perturber mass and six target-orbit
deltas are evaluated in one residual vector.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from scipy.optimize import OptimizeResult, least_squares

from src.astrometry.forward_model import forward_model, residuals_mas
from src.mass.likelihood_al import mahalanobis_residuals_2d

LikelihoodKind = Literal["al", "mahalanobis2d"]


@dataclass(frozen=True)
class GaiaObservationBundle:
    """Gaia SSO observations and covariance columns needed by the fit."""

    jd_tdb: np.ndarray
    gaia_xyz_bary: np.ndarray
    ra_deg: np.ndarray
    dec_deg: np.ndarray
    position_angle_scan_deg: np.ndarray
    ra_error_systematic_mas: np.ndarray
    dec_error_systematic_mas: np.ndarray
    ra_dec_correlation_systematic: np.ndarray
    ra_error_random_mas: np.ndarray
    dec_error_random_mas: np.ndarray
    ra_dec_correlation_random: np.ndarray


@dataclass(frozen=True)
class JointFitPriors:
    """Gaussian priors and hard bounds for the seven-parameter joint fit."""

    sigma_da_rel: float = 2e-4
    sigma_de: float = 5e-4
    sigma_di_deg: float = 0.05
    sigma_dOmega_deg: float = 0.2  # noqa: N815
    sigma_domega_deg: float = 0.2
    sigma_dM_deg: float = 0.5  # noqa: N815
    log10_mass_bounds: tuple[float, float] = (14.0, 23.0)
    da_rel_bounds: tuple[float, float] = (-1e-3, 1e-3)
    de_bounds: tuple[float, float] = (-2e-3, 2e-3)
    di_deg_bounds: tuple[float, float] = (-0.25, 0.25)
    dOmega_deg_bounds: tuple[float, float] = (-1.0, 1.0)  # noqa: N815
    domega_deg_bounds: tuple[float, float] = (-1.0, 1.0)
    dM_deg_bounds: tuple[float, float] = (-2.5, 2.5)  # noqa: N815

    @property
    def sigma_vector(self) -> np.ndarray:
        return np.array(
            [
                self.sigma_da_rel,
                self.sigma_de,
                self.sigma_di_deg,
                self.sigma_dOmega_deg,
                self.sigma_domega_deg,
                self.sigma_dM_deg,
            ],
            dtype=float,
        )

    @property
    def lower_bounds(self) -> np.ndarray:
        return np.array(
            [
                self.log10_mass_bounds[0],
                self.da_rel_bounds[0],
                self.de_bounds[0],
                self.di_deg_bounds[0],
                self.dOmega_deg_bounds[0],
                self.domega_deg_bounds[0],
                self.dM_deg_bounds[0],
            ],
            dtype=float,
        )

    @property
    def upper_bounds(self) -> np.ndarray:
        return np.array(
            [
                self.log10_mass_bounds[1],
                self.da_rel_bounds[1],
                self.de_bounds[1],
                self.di_deg_bounds[1],
                self.dOmega_deg_bounds[1],
                self.domega_deg_bounds[1],
                self.dM_deg_bounds[1],
            ],
            dtype=float,
        )


DEFAULT_PRIORS = JointFitPriors()


def _as_float_array(values: np.ndarray, *, name: str) -> np.ndarray:
    out = np.asarray(values, dtype=float)
    if out.ndim != 1:
        raise ValueError(f"{name} must be 1-D, got shape {out.shape}")
    return out


def _wrap_degrees(angle: float) -> float:
    return float(angle % 360.0)


def validate_observations(obs: GaiaObservationBundle) -> int:
    """Validate observation array shapes and return the number of observations."""
    jd = _as_float_array(obs.jd_tdb, name="jd_tdb")
    n_obs = len(jd)
    xyz = np.asarray(obs.gaia_xyz_bary, dtype=float)
    if xyz.shape != (n_obs, 3):
        raise ValueError(f"gaia_xyz_bary shape {xyz.shape} != ({n_obs}, 3)")
    for field_name in (
        "ra_deg",
        "dec_deg",
        "position_angle_scan_deg",
        "ra_error_systematic_mas",
        "dec_error_systematic_mas",
        "ra_dec_correlation_systematic",
        "ra_error_random_mas",
        "dec_error_random_mas",
        "ra_dec_correlation_random",
    ):
        arr = _as_float_array(getattr(obs, field_name), name=field_name)
        if len(arr) != n_obs:
            raise ValueError(f"{field_name} length {len(arr)} != {n_obs}")
    return n_obs


def apply_target_deltas(target_elements: dict, params: np.ndarray) -> dict:
    """Apply joint-fit orbital deltas to a target MPCORB element dict."""
    p = np.asarray(params, dtype=float)
    if p.shape != (7,):
        raise ValueError(f"joint params must have shape (7,), got {p.shape}")

    _, da_rel, de, di_deg, dOmega_deg, domega_deg, dM_deg = p
    adjusted = dict(target_elements)
    adjusted["a_au"] = float(target_elements["a_au"]) * (1.0 + float(da_rel))
    adjusted["e"] = float(target_elements["e"]) + float(de)
    adjusted["i_deg"] = float(target_elements["i_deg"]) + float(di_deg)
    adjusted["Omega_deg"] = _wrap_degrees(float(target_elements["Omega_deg"]) + dOmega_deg)
    adjusted["omega_deg"] = _wrap_degrees(float(target_elements["omega_deg"]) + domega_deg)
    adjusted["M_deg"] = _wrap_degrees(float(target_elements["M_deg"]) + dM_deg)
    return adjusted


def al_residuals_and_weights(
    dra_mas: np.ndarray,
    ddec_mas: np.ndarray,
    pa_scan_deg: np.ndarray,
    ra_err_sys: np.ndarray,
    dec_err_sys: np.ndarray,
    corr_sys: np.ndarray,
    ra_err_rand: np.ndarray,
    dec_err_rand: np.ndarray,
    corr_rand: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Project tangent-plane residuals onto Gaia along-scan with sigma_AL."""
    pa = np.radians(pa_scan_deg)
    e_al_ra = np.sin(pa)
    e_al_dec = np.cos(pa)
    r_al = dra_mas * e_al_ra + ddec_mas * e_al_dec

    def _projected_var(s_ra: np.ndarray, s_dec: np.ndarray, rho: np.ndarray) -> np.ndarray:
        out = (
            e_al_ra**2 * s_ra**2
            + 2.0 * e_al_ra * e_al_dec * rho * s_ra * s_dec
            + e_al_dec**2 * s_dec**2
        )
        return np.asarray(out, dtype=float)

    var_al = _projected_var(ra_err_sys, dec_err_sys, corr_sys) + _projected_var(
        ra_err_rand, dec_err_rand, corr_rand
    )
    sigma_al = np.sqrt(np.maximum(var_al, 1e-6))
    return r_al, sigma_al


def prior_residuals(params: np.ndarray, priors: JointFitPriors = DEFAULT_PRIORS) -> np.ndarray:
    """Return normalised Gaussian prior residuals for the six orbit deltas."""
    p = np.asarray(params, dtype=float)
    if p.shape != (7,):
        raise ValueError(f"joint params must have shape (7,), got {p.shape}")
    return np.asarray(p[1:] / priors.sigma_vector, dtype=float)


def residuals_joint(
    params: np.ndarray,
    target_elements: dict,
    perturber_elements: dict,
    obs: GaiaObservationBundle,
    *,
    priors: JointFitPriors = DEFAULT_PRIORS,
    background_elements: dict[str, dict] | None = None,
    include_planets: tuple[str, ...] = (
        "sun",
        "mercury",
        "venus",
        "earth",
        "mars",
        "jupiter",
        "saturn",
        "uranus",
        "neptune",
    ),
    dt_days: float = 1.0,
    integrator: str = "whfast",
    likelihood: LikelihoodKind = "al",
) -> np.ndarray:
    """Evaluate astrometric residuals plus orbital priors for the joint fit.

    ``likelihood`` selects how the per-observation (RA*, Dec) residual is
    weighted: ``"al"`` projects onto along-scan only (Stage 1 baseline),
    ``"mahalanobis2d"`` whitens the full 2D tangential residual with the
    Gaia (RA, Dec) covariance, keeping the across-scan component too
    (Stage 2).
    """
    validate_observations(obs)
    p = np.asarray(params, dtype=float)
    if p.shape != (7,):
        raise ValueError(f"joint params must have shape (7,), got {p.shape}")

    adjusted_target = apply_target_deltas(target_elements, p)
    mass_kg = float(10.0 ** p[0])
    ra_pred, dec_pred = forward_model(
        adjusted_target,
        perturber_elements,
        mass_kg,
        obs.jd_tdb,
        obs.gaia_xyz_bary,
        include_planets=include_planets,
        include_background=bool(background_elements),
        background_elements=background_elements,
        dt_days=dt_days,
        integrator=integrator,
    )
    dra, ddec = residuals_mas(obs.ra_deg, obs.dec_deg, ra_pred, dec_pred)
    if likelihood == "al":
        r_al, sigma_al = al_residuals_and_weights(
            dra,
            ddec,
            obs.position_angle_scan_deg,
            obs.ra_error_systematic_mas,
            obs.dec_error_systematic_mas,
            obs.ra_dec_correlation_systematic,
            obs.ra_error_random_mas,
            obs.dec_error_random_mas,
            obs.ra_dec_correlation_random,
        )
        astrometric = r_al / sigma_al
    elif likelihood == "mahalanobis2d":
        astrometric, _ = mahalanobis_residuals_2d(
            dra,
            ddec,
            obs.ra_error_systematic_mas,
            obs.dec_error_systematic_mas,
            obs.ra_dec_correlation_systematic,
            obs.ra_error_random_mas,
            obs.dec_error_random_mas,
            obs.ra_dec_correlation_random,
        )
    else:
        raise ValueError(f"unknown likelihood {likelihood!r}")
    return np.concatenate([astrometric, prior_residuals(p, priors)])


def fit_joint(
    target_elements: dict,
    perturber_elements: dict,
    obs: GaiaObservationBundle,
    *,
    initial_log10_mass: float,
    initial_deltas: np.ndarray | None = None,
    priors: JointFitPriors = DEFAULT_PRIORS,
    background_elements: dict[str, dict] | None = None,
    include_planets: tuple[str, ...] = (
        "sun",
        "mercury",
        "venus",
        "earth",
        "mars",
        "jupiter",
        "saturn",
        "uranus",
        "neptune",
    ),
    dt_days: float = 1.0,
    integrator: str = "whfast",
    max_nfev: int = 800,
    likelihood: LikelihoodKind = "al",
) -> OptimizeResult:
    """Run the seven-parameter least-squares joint fit."""
    if initial_deltas is None:
        deltas = np.zeros(6, dtype=float)
    else:
        deltas = np.asarray(initial_deltas, dtype=float)
        if deltas.shape != (6,):
            raise ValueError(f"initial_deltas must have shape (6,), got {deltas.shape}")
    x0 = np.concatenate([[float(initial_log10_mass)], deltas])
    return least_squares(
        residuals_joint,
        x0,
        args=(target_elements, perturber_elements, obs),
        kwargs={
            "priors": priors,
            "background_elements": background_elements,
            "include_planets": include_planets,
            "dt_days": dt_days,
            "integrator": integrator,
            "likelihood": likelihood,
        },
        method="trf",
        bounds=(priors.lower_bounds, priors.upper_bounds),
        max_nfev=max_nfev,
        ftol=1e-10,
        xtol=1e-10,
        gtol=1e-10,
    )
