"""Unit tests for AL projection, sigma_AL, and MPCORB auto-selection.

These tests cover the core along-scan projection logic that is the backbone
of the LOO mass-fit methodology.  All tests are pure Python — no REBOUND,
no network, no real MPCORB data required.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

from scripts.fit_mass_gaia_loo import (
    _best_mpcorb_snapshot,
    _mass_from_h,
    al_residuals_and_weights,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_al_args(
    dra: np.ndarray,
    ddec: np.ndarray,
    pa_deg: np.ndarray,
    s_ra_sys: float = 0.5,
    s_dec_sys: float = 0.5,
    rho_sys: float = 0.0,
    s_ra_rand: float = 0.3,
    s_dec_rand: float = 0.3,
    rho_rand: float = 0.0,
) -> tuple:
    n = len(dra)
    return (
        dra,
        ddec,
        pa_deg,
        np.full(n, s_ra_sys),
        np.full(n, s_dec_sys),
        np.full(n, rho_sys),
        np.full(n, s_ra_rand),
        np.full(n, s_dec_rand),
        np.full(n, rho_rand),
    )


# ---------------------------------------------------------------------------
# al_residuals_and_weights — analytical cases
# ---------------------------------------------------------------------------


class TestALProjection:
    def test_pa_zero_is_pure_dec(self):
        """PA=0° (scan toward North): AL = Dec component only."""
        dra = np.array([5.0, -3.0])
        ddec = np.array([2.0, 4.0])
        pa = np.zeros(2)
        r_al, _ = al_residuals_and_weights(*_make_al_args(dra, ddec, pa))
        # sin(0)=0, cos(0)=1 → r_al = ddec
        np.testing.assert_allclose(r_al, ddec, atol=1e-12)

    def test_pa_90_is_pure_ra(self):
        """PA=90° (scan toward East): AL = RA component only."""
        dra = np.array([7.0, -1.0])
        ddec = np.array([2.0, 9.0])
        pa = np.full(2, 90.0)
        r_al, _ = al_residuals_and_weights(*_make_al_args(dra, ddec, pa))
        # sin(90)=1, cos(90)=0 → r_al = dra
        np.testing.assert_allclose(r_al, dra, atol=1e-12)

    def test_pa_45_equal_weight(self):
        """PA=45°: AL = (dRA + dDec) / sqrt(2)."""
        dra = np.array([3.0])
        ddec = np.array([3.0])
        pa = np.full(1, 45.0)
        r_al, _ = al_residuals_and_weights(*_make_al_args(dra, ddec, pa))
        expected = (dra + ddec) / math.sqrt(2)
        np.testing.assert_allclose(r_al, expected, atol=1e-12)

    def test_zero_residual_gives_zero_al(self):
        n = 5
        pa = np.linspace(0, 90, n)
        r_al, _ = al_residuals_and_weights(*_make_al_args(np.zeros(n), np.zeros(n), pa))
        np.testing.assert_allclose(r_al, 0.0, atol=1e-12)

    def test_sign_preserved(self):
        """A negative Dec residual at PA=0 gives a negative AL residual."""
        dra = np.array([0.0])
        ddec = np.array([-5.0])
        pa = np.zeros(1)
        r_al, _ = al_residuals_and_weights(*_make_al_args(dra, ddec, pa))
        assert r_al[0] < 0.0

    def test_pa_180_is_minus_dec(self):
        """PA=180° reverses the Dec component (scan toward South)."""
        dra = np.array([0.0])
        ddec = np.array([3.0])
        pa0 = np.zeros(1)
        pa180 = np.full(1, 180.0)
        r0, _ = al_residuals_and_weights(*_make_al_args(dra, ddec, pa0))
        r180, _ = al_residuals_and_weights(*_make_al_args(dra, ddec, pa180))
        np.testing.assert_allclose(r0, -r180, atol=1e-12)


class TestSigmaAL:
    def test_floor_when_all_errors_zero(self):
        """sigma_AL has a floor of sqrt(1e-6) ≈ 0.001 mas."""
        n = 4
        pa = np.zeros(n)
        zeros = np.zeros(n)
        _, sigma_al = al_residuals_and_weights(
            zeros, zeros, pa, zeros, zeros, zeros, zeros, zeros, zeros
        )
        assert np.all(sigma_al >= math.sqrt(1e-6) - 1e-15)

    def test_systematic_dominates(self):
        """Large systematic >> random: sigma_AL ≈ systematic projection."""
        pa = np.array([0.0])  # scan = Dec direction
        s_dec_sys = 10.0  # 10 mas systematic in Dec
        n = 1
        _, sigma_al = al_residuals_and_weights(
            np.zeros(n),
            np.zeros(n),
            pa,
            np.zeros(n),
            np.full(n, s_dec_sys),
            np.zeros(n),
            np.zeros(n),
            np.zeros(n),
            np.zeros(n),
        )
        # σ²_AL = cos²(0)*s_dec² = 100 → σ_AL = 10
        np.testing.assert_allclose(sigma_al, [s_dec_sys], atol=1e-10)

    def test_combined_variances_add(self):
        """σ²_AL = σ²_AL_sys + σ²_AL_rand (independent components)."""
        pa = np.array([0.0])
        s_sys = 3.0  # mas
        s_rand = 4.0  # mas
        _, sigma_al = al_residuals_and_weights(
            np.zeros(1),
            np.zeros(1),
            pa,
            np.zeros(1),
            np.full(1, s_sys),
            np.zeros(1),
            np.zeros(1),
            np.full(1, s_rand),
            np.zeros(1),
        )
        expected = math.sqrt(s_sys**2 + s_rand**2)  # 5.0 mas
        np.testing.assert_allclose(sigma_al, [expected], atol=1e-10)

    def test_correlation_increases_sigma(self):
        """Positive correlation between RA/Dec increases sigma_AL at PA=45°."""
        pa = np.full(1, 45.0)
        s = 1.0
        n = 1
        _, s_uncorr = al_residuals_and_weights(
            np.zeros(n),
            np.zeros(n),
            pa,
            np.full(n, s),
            np.full(n, s),
            np.zeros(n),
            np.zeros(n),
            np.zeros(n),
            np.zeros(n),
        )
        _, s_corr = al_residuals_and_weights(
            np.zeros(n),
            np.zeros(n),
            pa,
            np.full(n, s),
            np.full(n, s),
            np.full(n, 1.0),  # rho=+1
            np.zeros(n),
            np.zeros(n),
            np.zeros(n),
        )
        # With rho=+1: σ²_AL = (sin+cos)² = 2 > uncorrelated = 1
        assert s_corr[0] > s_uncorr[0]

    def test_output_shapes_match_input(self):
        n = 17
        pa = np.linspace(0, 360, n)
        dra = np.random.default_rng(0).normal(0, 1, n)
        ddec = np.random.default_rng(1).normal(0, 1, n)
        r_al, sigma_al = al_residuals_and_weights(*_make_al_args(dra, ddec, pa))
        assert r_al.shape == (n,)
        assert sigma_al.shape == (n,)


# ---------------------------------------------------------------------------
# _best_mpcorb_snapshot — MPCORB auto-selection
# ---------------------------------------------------------------------------


def _write_snapshot(directory: Path, name: str, epoch_jd: float) -> None:
    dat = directory / f"{name}.DAT"
    dat.write_text("placeholder")
    sidecar = directory / f"{name}.json"
    sidecar.write_text(json.dumps({"epoch_jd_tdb": epoch_jd}))


class TestBestMpcorbSnapshot:
    def test_picks_closest_epoch(self, tmp_path):
        _write_snapshot(tmp_path, "MPCORB_A", 2456000.0)  # far
        _write_snapshot(tmp_path, "MPCORB_B", 2457000.0)  # close  ← expected
        _write_snapshot(tmp_path, "MPCORB_C", 2458000.0)  # far
        result = _best_mpcorb_snapshot(tmp_path, 2457100.0)
        assert "MPCORB_B" in result.name

    def test_exact_match(self, tmp_path):
        _write_snapshot(tmp_path, "MPCORB_X", 2457000.5)
        result = _best_mpcorb_snapshot(tmp_path, 2457000.5)
        assert "MPCORB_X" in result.name

    def test_tie_broken_alphabetically(self, tmp_path):
        """Two snapshots with identical epoch — alphabetically first wins."""
        _write_snapshot(tmp_path, "MPCORB_20160217", 2457400.5)
        _write_snapshot(tmp_path, "MPCORB_20160303", 2457400.5)
        result = _best_mpcorb_snapshot(tmp_path, 2457547.5)
        assert "20160217" in result.name

    def test_returns_dat_path(self, tmp_path):
        _write_snapshot(tmp_path, "MPCORB_Z", 2457200.5)
        result = _best_mpcorb_snapshot(tmp_path, 2457200.5)
        assert result.suffix == ".DAT"

    def test_skips_json_without_dat(self, tmp_path):
        """JSON sidecar without matching DAT file is silently skipped."""
        # Only JSON, no DAT
        (tmp_path / "MPCORB_GHOST.json").write_text(json.dumps({"epoch_jd_tdb": 2457000.0}))
        # Valid pair further away
        _write_snapshot(tmp_path, "MPCORB_REAL", 2458000.0)
        result = _best_mpcorb_snapshot(tmp_path, 2457000.0)
        assert "MPCORB_REAL" in result.name

    def test_skips_malformed_json(self, tmp_path):
        """Malformed JSON sidecar is silently skipped."""
        (tmp_path / "MPCORB_BAD.DAT").write_text("x")
        (tmp_path / "MPCORB_BAD.json").write_text("not valid json {{")
        _write_snapshot(tmp_path, "MPCORB_GOOD", 2457000.0)
        result = _best_mpcorb_snapshot(tmp_path, 2457000.0)
        assert "MPCORB_GOOD" in result.name

    def test_skips_missing_epoch_key(self, tmp_path):
        """JSON without epoch_jd_tdb key is silently skipped."""
        (tmp_path / "MPCORB_NOKEY.DAT").write_text("x")
        (tmp_path / "MPCORB_NOKEY.json").write_text(json.dumps({"other": 123}))
        _write_snapshot(tmp_path, "MPCORB_OK", 2457000.0)
        result = _best_mpcorb_snapshot(tmp_path, 2457000.0)
        assert "MPCORB_OK" in result.name

    def test_empty_archive_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="No valid MPCORB snapshots"):
            _best_mpcorb_snapshot(tmp_path, 2457000.0)

    def test_gaia_era_encounters_pick_contemporaneous(self, tmp_path):
        """For Gaia-era encounters the 2015-2016 snapshots win over 2025."""
        _write_snapshot(tmp_path, "MPCORB_20150524", 2457200.5)
        _write_snapshot(tmp_path, "MPCORB_20160217", 2457400.5)
        _write_snapshot(tmp_path, "MPCORB_20251201", 2461000.5)
        # Ate encounter 2016-06-08 ≈ JD 2457547.5 → 20160217 wins (147 days)
        result = _best_mpcorb_snapshot(tmp_path, 2457547.5)
        assert "20160217" in result.name


# ---------------------------------------------------------------------------
# _mass_from_h — H-magnitude to mass estimator
# ---------------------------------------------------------------------------


class TestMassFromH:
    def test_none_returns_default(self):
        """H=None → default 10^18 kg (placeholder for unknown mass)."""
        assert _mass_from_h(None) == 1.0e18

    def test_brighter_gives_larger_mass(self):
        """Lower H (brighter) implies larger diameter and hence more mass."""
        m_bright = _mass_from_h(5.0)
        m_faint = _mass_from_h(15.0)
        assert m_bright > m_faint

    def test_ceres_like_mass(self):
        """H=3.4 (Ceres) → mass in plausible range for ~950 km body."""
        m = _mass_from_h(3.4)
        # Expected ~10^21 kg order of magnitude; our formula uses ρ=1500 kg/m³
        assert 1e20 < m < 1e23

    def test_positive_output(self):
        for h in [1.0, 5.0, 10.0, 20.0]:
            assert _mass_from_h(h) > 0
