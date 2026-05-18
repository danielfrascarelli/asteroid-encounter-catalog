"""Unit tests for src/propagate/nbody.py — propagate_grid_nbody."""

from __future__ import annotations

import numpy as np
import polars as pl
import pytest

from src.propagate.kepler import kepler_to_cartesian
from src.propagate.nbody import propagate_grid_nbody

_EPOCH_JD = 2451545.0  # J2000


def _elements(n: int = 2) -> pl.DataFrame:
    """Synthetic circular belt asteroids, all sharing epoch J2000."""
    rng = np.random.default_rng(0)
    return pl.DataFrame(
        {
            "number": np.arange(1, n + 1, dtype=np.int32),
            "designation": [str(k) for k in range(1, n + 1)],
            "a_au": rng.uniform(2.2, 3.2, size=n),
            "e": np.zeros(n),
            "i_deg": np.zeros(n),
            "Omega_deg": np.zeros(n),
            "omega_deg": np.zeros(n),
            "M_deg": rng.uniform(0.0, 360.0, size=n),
            "epoch_jd": np.full(n, _EPOCH_JD),
        },
        schema_overrides={"number": pl.Int32, "designation": pl.Utf8},
    )


# ---------------------------------------------------------------------------
# Shape and contract tests
# ---------------------------------------------------------------------------


def test_output_shape_and_dtype() -> None:
    """Returns (T, N, 3) float32."""
    t_grid = np.linspace(_EPOCH_JD, _EPOCH_JD + 10.0, 7)
    out = propagate_grid_nbody(_elements(3), t_grid, include_planets=["sun"])
    assert out.shape == (7, 3, 3)
    assert out.dtype == np.float32


def test_out_parameter_writes_in_place() -> None:
    """out= is populated and returned; no copy is made."""
    t_grid = np.linspace(_EPOCH_JD, _EPOCH_JD + 5.0, 4)
    buf = np.zeros((4, 2, 3), dtype=np.float32)
    result = propagate_grid_nbody(_elements(2), t_grid, include_planets=["sun"], out=buf)
    assert result is buf
    # Positions must be non-zero (asteroids are in the belt, not at the origin)
    assert np.any(buf != 0.0)


def test_empty_elements_returns_empty() -> None:
    """Empty elements DF → shape (T, 0, 3)."""
    t_grid = np.linspace(_EPOCH_JD, _EPOCH_JD + 2.0, 3)
    empty = pl.DataFrame(
        {
            "number": np.array([], dtype=np.int32),
            "designation": pl.Series([], dtype=pl.Utf8),
            "a_au": np.array([], dtype=np.float64),
            "e": np.array([], dtype=np.float64),
            "i_deg": np.array([], dtype=np.float64),
            "Omega_deg": np.array([], dtype=np.float64),
            "omega_deg": np.array([], dtype=np.float64),
            "M_deg": np.array([], dtype=np.float64),
            "epoch_jd": np.array([], dtype=np.float64),
        },
        schema_overrides={"number": pl.Int32, "designation": pl.Utf8},
    )
    out = propagate_grid_nbody(empty, t_grid, epoch_jd=_EPOCH_JD)
    assert out.shape == (3, 0, 3)


def test_nonmonotonic_time_grid_raises() -> None:
    """Non-strictly-increasing time grid raises ValueError."""
    t_grid = np.array([_EPOCH_JD, _EPOCH_JD + 2.0, _EPOCH_JD + 1.0])
    with pytest.raises(ValueError, match="monotonically"):
        propagate_grid_nbody(_elements(1), t_grid)


def test_out_shape_mismatch_raises() -> None:
    """Passing out= with wrong shape raises ValueError."""
    t_grid = np.linspace(_EPOCH_JD, _EPOCH_JD + 5.0, 3)
    wrong = np.empty((5, 2, 3), dtype=np.float32)  # T=5 but grid has T=3
    with pytest.raises(ValueError, match="shape"):
        propagate_grid_nbody(_elements(2), t_grid, include_planets=["sun"], out=wrong)


# ---------------------------------------------------------------------------
# Physical accuracy test
# ---------------------------------------------------------------------------


def test_sun_only_matches_kepler() -> None:
    """With Sun as the sole massive body, N-body positions agree with Kepler
    2-body to within 1e-3 AU over a 30-day window.

    For circular orbits and a Sun-only force field the equations of motion are
    identical to Keplerian, so the only error source is the numerical integrator.
    IAS15 (adaptive) reduces this to near floating-point precision, well inside
    the float32 output precision (~1e-7 relative).
    """
    t_grid = np.linspace(_EPOCH_JD, _EPOCH_JD + 30.0, 30)
    elems = _elements(2)

    pos_nbody = propagate_grid_nbody(
        elems,
        t_grid,
        include_planets=["sun"],
        integrator="ias15",
    )  # (30, 2, 3) float32

    a_arr = elems["a_au"].to_numpy()
    m0_arr = elems["M_deg"].to_numpy() * (np.pi / 180.0)

    for k in range(len(elems)):
        pos_kepler = kepler_to_cartesian(
            a_au=a_arr[k],
            e=0.0,
            i_rad=0.0,
            Omega_rad=0.0,
            omega_rad=0.0,
            M0_rad=m0_arr[k],
            epoch_jd=_EPOCH_JD,
            t_jd=t_grid,
        )  # (30, 3)

        diff = np.linalg.norm(pos_nbody[:, k, :].astype(np.float64) - pos_kepler, axis=1)
        assert (
            diff.max() < 1e-3
        ), f"asteroid {k} (a={a_arr[k]:.2f} AU): max |Δr| = {diff.max():.2e} AU"
