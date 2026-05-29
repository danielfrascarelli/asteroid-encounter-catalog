"""Multi-target joint orbit-drift + perturber-mass residual model.

Shares one perturber-mass parameter across multiple Gaia close encounters of
the same perturber and keeps six free orbital deltas per target. The aim is
to break the M ↔ deltas degeneracy that Stage 4 of the deepwork exposed
(structural fit/lit bias, see ``docs/mass_layer_validation.md``).

Parameter vector layout
-----------------------
``params`` has length ``1 + 6 * N`` where ``N == len(target_bundles)``::

    params[0]              = log10_M_perturber
    params[1 + 6*i : 1 + 6*(i+1)] = (da_rel, de, di, dOmega, domega, dM) for target i

Residual vector layout
----------------------
``residuals`` concatenates astrometric residuals of every target followed by
the 6*N normalised Gaussian priors. Astrometric whitening reuses the same
likelihood functions as the single-target joint fit (along-scan or 2D
Mahalanobis).
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from scipy.optimize import OptimizeResult, least_squares

from src.astrometry.forward_model import forward_model, residuals_mas
from src.mass.forward_model_joint import (
    DEFAULT_PRIORS,
    GaiaObservationBundle,
    JointFitPriors,
    LikelihoodKind,
    al_residuals_and_weights,
    apply_target_deltas,
    validate_observations,
)
from src.mass.likelihood_al import mahalanobis_residuals_2d


@dataclass(frozen=True)
class TargetBundle:
    """One target's MPCORB elements + Gaia observations for the multi-target fit."""

    target_number: int
    elements: dict
    obs: GaiaObservationBundle


def _validate_params(params: np.ndarray, n_targets: int) -> np.ndarray:
    p = np.asarray(params, dtype=float)
    expected = 1 + 6 * n_targets
    if p.shape != (expected,):
        raise ValueError(
            f"multitarget params must have shape ({expected},) for {n_targets} targets, "
            f"got {p.shape}"
        )
    return p


def _target_slice(params: np.ndarray, idx: int) -> np.ndarray:
    """Return the 7-vector (log10_M, 6 deltas) for the i-th target."""
    start = 1 + 6 * idx
    out = np.empty(7, dtype=float)
    out[0] = params[0]
    out[1:] = params[start : start + 6]
    return out


def make_bounds(
    n_targets: int, priors: JointFitPriors = DEFAULT_PRIORS
) -> tuple[np.ndarray, np.ndarray]:
    """Stack the 7-parameter bounds N times (mass + N×6 deltas)."""
    lo = np.empty(1 + 6 * n_targets, dtype=float)
    hi = np.empty_like(lo)
    lo[0] = priors.log10_mass_bounds[0]
    hi[0] = priors.log10_mass_bounds[1]
    delta_lo = priors.lower_bounds[1:]
    delta_hi = priors.upper_bounds[1:]
    for i in range(n_targets):
        s = 1 + 6 * i
        lo[s : s + 6] = delta_lo
        hi[s : s + 6] = delta_hi
    return lo, hi


def prior_residuals_multitarget(
    params: np.ndarray, n_targets: int, priors: JointFitPriors = DEFAULT_PRIORS
) -> np.ndarray:
    """Return the 6*N normalised Gaussian prior residuals."""
    p = _validate_params(params, n_targets)
    sigma = priors.sigma_vector
    out = np.empty(6 * n_targets, dtype=float)
    for i in range(n_targets):
        s = 1 + 6 * i
        out[6 * i : 6 * (i + 1)] = p[s : s + 6] / sigma
    return out


def _astrometric_residuals_one_target(
    params7: np.ndarray,
    target_elements: dict,
    perturber_elements: dict,
    obs: GaiaObservationBundle,
    background_elements: dict[str, dict] | None,
    include_planets: tuple[str, ...],
    dt_days: float,
    integrator: str,
    likelihood: LikelihoodKind,
) -> np.ndarray:
    validate_observations(obs)
    adjusted_target = apply_target_deltas(target_elements, params7)
    mass_kg = float(10.0 ** params7[0])
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
        return np.asarray(r_al / sigma_al, dtype=float)
    if likelihood == "mahalanobis2d":
        whitened, _ = mahalanobis_residuals_2d(
            dra,
            ddec,
            obs.ra_error_systematic_mas,
            obs.dec_error_systematic_mas,
            obs.ra_dec_correlation_systematic,
            obs.ra_error_random_mas,
            obs.dec_error_random_mas,
            obs.ra_dec_correlation_random,
        )
        return np.asarray(whitened, dtype=float)
    raise ValueError(f"unknown likelihood {likelihood!r}")


def residuals_joint_multitarget(
    params: np.ndarray,
    target_bundles: Sequence[TargetBundle],
    perturber_elements: dict,
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
    """Evaluate astrometric + prior residuals for the multi-target joint fit.

    Astrometric residuals are concatenated in the order of ``target_bundles``;
    the 6*N prior residuals are appended at the end. The single ``log10_M`` is
    shared across all targets.
    """
    n_targets = len(target_bundles)
    if n_targets == 0:
        raise ValueError("target_bundles must contain at least one target")
    p = _validate_params(params, n_targets)

    chunks: list[np.ndarray] = []
    for i, bundle in enumerate(target_bundles):
        p_i = _target_slice(p, i)
        chunks.append(
            _astrometric_residuals_one_target(
                p_i,
                bundle.elements,
                perturber_elements,
                bundle.obs,
                background_elements,
                include_planets,
                dt_days,
                integrator,
                likelihood,
            )
        )
    chunks.append(prior_residuals_multitarget(p, n_targets, priors))
    return np.concatenate(chunks)


def fit_joint_multitarget(
    target_bundles: Sequence[TargetBundle],
    perturber_elements: dict,
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
    max_nfev: int = 1200,
    likelihood: LikelihoodKind = "al",
) -> OptimizeResult:
    """Run the (1 + 6*N)-parameter least-squares fit with shared perturber mass."""
    n_targets = len(target_bundles)
    if n_targets == 0:
        raise ValueError("target_bundles must contain at least one target")
    if initial_deltas is None:
        deltas = np.zeros(6 * n_targets, dtype=float)
    else:
        deltas = np.asarray(initial_deltas, dtype=float)
        if deltas.shape != (6 * n_targets,):
            raise ValueError(
                f"initial_deltas must have shape ({6 * n_targets},), got {deltas.shape}"
            )
    x0 = np.concatenate([[float(initial_log10_mass)], deltas])
    lo, hi = make_bounds(n_targets, priors)
    return least_squares(
        residuals_joint_multitarget,
        x0,
        args=(target_bundles, perturber_elements),
        kwargs={
            "priors": priors,
            "background_elements": background_elements,
            "include_planets": include_planets,
            "dt_days": dt_days,
            "integrator": integrator,
            "likelihood": likelihood,
        },
        method="trf",
        bounds=(lo, hi),
        max_nfev=max_nfev,
        ftol=1e-10,
        xtol=1e-10,
        gtol=1e-10,
    )
