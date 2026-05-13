"""Tests for src/utils/time_utils.py.

Reference values computed with astropy independently to guard against
regressions if the implementation changes.
"""

import pytest
from astropy.time import Time

from src.utils.time_utils import (
    jd_tdb_to_time,
    tcb_to_tdb,
    tdb_to_tcb,
    tdb_to_utc,
    utc_to_tdb,
)

# J2000.0 in JD
J2000_JD = 2451545.0

# Gaia DR3 observation window endpoints (UTC strings)
GAIA_START_UTC = "2014-07-25T00:00:00"
GAIA_END_UTC = "2017-05-28T00:00:00"


# ---------------------------------------------------------------------------
# utc_to_tdb
# ---------------------------------------------------------------------------


def test_utc_to_tdb_returns_tdb_scale() -> None:
    t = utc_to_tdb(GAIA_START_UTC)
    assert t.scale == "tdb"


def test_utc_to_tdb_from_string() -> None:
    t = utc_to_tdb(GAIA_START_UTC)
    # TDB - UTC ≈ leap_seconds(35 in 2014) + 32.184 s (TAI→TT) ≈ 67 s
    t_utc = Time(GAIA_START_UTC, format="isot", scale="utc")
    delta_seconds = abs((t.jd - t_utc.jd) * 86400)
    assert 60.0 < delta_seconds < 75.0


def test_utc_to_tdb_from_time_object() -> None:
    t_utc = Time(GAIA_START_UTC, format="isot", scale="utc")
    t = utc_to_tdb(t_utc)
    assert t.scale == "tdb"


def test_utc_to_tdb_preserves_instant() -> None:
    """Converting to TDB and back to UTC must recover the original instant."""
    t_orig = Time(GAIA_START_UTC, format="isot", scale="utc")
    t_tdb = utc_to_tdb(t_orig)
    t_back = tdb_to_utc(t_tdb)
    # Tolerance: 1 microsecond in JD units
    assert abs(t_back.jd - t_orig.jd) < 1e-11


# ---------------------------------------------------------------------------
# tdb_to_utc
# ---------------------------------------------------------------------------


def test_tdb_to_utc_returns_utc_scale() -> None:
    t_tdb = Time(J2000_JD, format="jd", scale="tdb")
    assert tdb_to_utc(t_tdb).scale == "utc"


# ---------------------------------------------------------------------------
# tcb_to_tdb
# ---------------------------------------------------------------------------


def test_tcb_to_tdb_returns_float() -> None:
    result = tcb_to_tdb(J2000_JD)
    assert isinstance(result, float)


def test_tcb_to_tdb_direction() -> None:
    """TDB < TCB for epochs after T0 (TCB ticks faster)."""
    result = tcb_to_tdb(J2000_JD)
    assert result < J2000_JD


def test_tcb_to_tdb_magnitude() -> None:
    """At J2000, accumulated TDB-TCB offset ≈ L_B * (J2000 - T0) days.

    T0 = 2443144.5003725, J2000 = 2451545.0
    Δ ≈ 1.55e-8 * (2451545 - 2443144.5) * 86400 s ≈ 11.25 s
    """
    result = tcb_to_tdb(J2000_JD)
    delta_seconds = (J2000_JD - result) * 86400.0
    assert 10.0 < delta_seconds < 13.0


def test_tcb_to_tdb_at_epoch_is_identity() -> None:
    """At the defining epoch T0, TCB == TDB by definition."""
    tcb_epoch = 2443144.5003725
    result = tcb_to_tdb(tcb_epoch)
    assert abs(result - tcb_epoch) < 1e-10


def test_tcb_to_tdb_gaia_window() -> None:
    """Within the Gaia window, offset grows ~1.6 ms/day — check sign and order."""
    t_start = utc_to_tdb(GAIA_START_UTC).jd
    t_end = utc_to_tdb(GAIA_END_UTC).jd
    # Gaia epoch col is in TCB; converting to TDB should give slightly smaller JD
    t_start_from_tcb = tcb_to_tdb(t_start)
    t_end_from_tcb = tcb_to_tdb(t_end)
    assert t_start_from_tcb < t_start
    assert t_end_from_tcb < t_end


# ---------------------------------------------------------------------------
# tdb_to_tcb (inverse)
# ---------------------------------------------------------------------------


def test_tdb_to_tcb_is_inverse_of_tcb_to_tdb() -> None:
    """Round-trip TCB → TDB → TCB must recover the original value."""
    jd_tcb = J2000_JD + 1234.5
    recovered = tdb_to_tcb(tcb_to_tdb(jd_tcb))
    assert abs(recovered - jd_tcb) < 1e-10


def test_tdb_to_tcb_direction() -> None:
    """TCB > TDB for epochs after T0."""
    result = tdb_to_tcb(J2000_JD)
    assert result > J2000_JD


# ---------------------------------------------------------------------------
# jd_tdb_to_time
# ---------------------------------------------------------------------------


def test_jd_tdb_to_time_scale() -> None:
    t = jd_tdb_to_time(J2000_JD)
    assert t.scale == "tdb"


def test_jd_tdb_to_time_value() -> None:
    t = jd_tdb_to_time(J2000_JD)
    assert t.jd == pytest.approx(J2000_JD)
