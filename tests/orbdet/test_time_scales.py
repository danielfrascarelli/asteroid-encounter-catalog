"""Tests de src/orbdet/time_scales.py."""

from __future__ import annotations

import numpy as np
import pytest

from src.orbdet import time_scales as ts


def test_tcb_tdb_roundtrip() -> None:
    jd = 2_457_000.0
    back = ts.tdb_to_tcb(ts.tcb_to_tdb(jd))
    assert float(back[0]) == pytest.approx(jd, abs=1e-9)


def test_tcb_minus_tdb_is_seconds_scale() -> None:
    """TCB y TDB difieren del orden de segundos en la era Gaia (no horas)."""
    jd = 2_457_000.0  # ~2014-12
    diff_days = float(ts.tcb_to_tdb(jd)[0]) - jd
    diff_s = abs(diff_days) * 86_400.0
    assert 1.0 < diff_s < 60.0  # ~14 s en esta época


def test_utc_tdb_roundtrip() -> None:
    jd = 2_457_000.5
    back = ts.tdb_to_utc(ts.utc_to_tdb(jd))
    assert float(back[0]) == pytest.approx(jd, abs=1e-8)


def test_utc_tdb_offset_about_minute() -> None:
    """TDB ≈ UTC + (32.184 s + leap seconds) ≈ ~67 s en 2015."""
    jd = 2_457_000.5
    diff_s = (float(ts.utc_to_tdb(jd)[0]) - jd) * 86_400.0
    assert 60.0 < diff_s < 75.0


def test_vectorized() -> None:
    jds = np.array([2_456_900.0, 2_457_000.0, 2_457_100.0])
    out = ts.tcb_to_tdb(jds)
    assert out.shape == (3,)


def test_gaia_epoch_to_jd_tdb_reference() -> None:
    """epoch=0 → JD_TCB = J2010 ref → TDB cercano a esa fecha."""
    jd_tdb = float(ts.gaia_epoch_to_jd_tdb(0.0)[0])
    # Debe quedar a < 1 minuto del JD de referencia (solo difiere por TCB→TDB).
    assert abs(jd_tdb - ts.J2010_TCB_JD) * 86_400.0 < 60.0


def test_iso_utc_to_jd_tdb() -> None:
    jd = ts.iso_utc_to_jd_tdb("2015-06-04T00:00:00")
    assert 2_457_177.0 < jd < 2_457_178.0
