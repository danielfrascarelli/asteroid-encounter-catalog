"""Time scale conversions for the pipeline.

All internal computations use Julian Date in TDB (Barycentric Dynamical Time).
Conversions from/to other scales happen only at I/O boundaries.

Scale relationships (approximate, for context):
- TDB vs TCB: TCB runs faster by L_B = 1.550519768e-8. Over one year the
  accumulated difference is ~1.8 s (≈ 1.6 ms/day). MPCORB epochs are TT
  (converted to TDB on ingest); Gaia DR3 epoch column uses TCB.
- TDB vs TT: differ by periodic terms only, max amplitude ~1.7 ms. For the
  purposes of this pipeline (detection threshold 0.01 AU) TDB ≈ TT is fine,
  but we use the proper astropy conversion regardless.
- TDB vs UTC: UTC has leap seconds; the offset grows over time.

References
----------
- IAU 2006 resolutions on time scales
- Lindegren & Dravins (2003), A&A 401, 1185 (TCB definition)
- astropy.time documentation
"""

from __future__ import annotations

from astropy.time import Time

# IAU 2006 defining constant for TCB–TDB conversion
_LB = 1.550519768e-8  # dimensionless rate difference TCB/TDB - 1
_TCB_JD0 = 2443144.5003725  # JD of TCB/TDB epoch (1977-01-01 00:00:32.184 TAI)


def utc_to_tdb(t: str | Time) -> Time:
    """Convert a UTC time to TDB scale.

    Parameters
    ----------
    t:
        Input time. If a string, must be ISO 8601 format (e.g.
        ``"2014-07-25T00:00:00"``). If a :class:`~astropy.time.Time`, its
        current scale is respected and converted to TDB.

    Returns
    -------
    astropy.time.Time
        Same instant expressed in TDB scale.

    Examples
    --------
    >>> t = utc_to_tdb("2014-07-25T00:00:00")
    >>> t.scale
    'tdb'
    """
    if isinstance(t, str):
        t = Time(t, format="isot", scale="utc")
    return t.tdb


def tdb_to_utc(t: Time) -> Time:
    """Convert a TDB time to UTC scale.

    Parameters
    ----------
    t:
        Input time in TDB (or any scale recognised by astropy).

    Returns
    -------
    astropy.time.Time
        Same instant expressed in UTC scale.
    """
    return t.utc


def tcb_to_tdb(jd_tcb: float) -> float:
    """Convert a Julian Date in TCB to TDB.

    Uses the linear approximation defined by IAU 2006:
        JD_TDB ≈ JD_TCB - L_B * (JD_TCB - T_0)

    where L_B = 1.550519768e-8 and T_0 = 2443144.5003725 JD (the TCB/TDB
    epoch of 1977-01-01T00:00:32.184 TAI).

    The residual nonlinear periodic terms are at the sub-microsecond level,
    negligible for this pipeline.

    Parameters
    ----------
    jd_tcb:
        Julian Date in TCB scale.

    Returns
    -------
    float
        Julian Date in TDB scale.

    Examples
    --------
    >>> abs(tcb_to_tdb(2451545.0) - 2451545.0) < 0.001
    True
    """
    return jd_tcb - _LB * (jd_tcb - _TCB_JD0)


def tdb_to_tcb(jd_tdb: float) -> float:
    """Convert a Julian Date in TDB to TCB (inverse of :func:`tcb_to_tdb`).

    Parameters
    ----------
    jd_tdb:
        Julian Date in TDB scale.

    Returns
    -------
    float
        Julian Date in TCB scale.
    """
    # Invert: jd_tdb = jd_tcb - L_B*(jd_tcb - T0)
    #          jd_tdb = jd_tcb*(1 - L_B) + L_B*T0
    #          jd_tcb = (jd_tdb - L_B*T0) / (1 - L_B)
    return (jd_tdb - _LB * _TCB_JD0) / (1.0 - _LB)


def jd_tdb_to_time(jd: float) -> Time:
    """Wrap a bare JD float (TDB) into an astropy Time object.

    Convenience function to make the scale explicit when constructing
    Time objects from raw JD values obtained from MPCORB or internal
    propagation grids.

    Parameters
    ----------
    jd:
        Julian Date in TDB.

    Returns
    -------
    astropy.time.Time
    """
    return Time(jd, format="jd", scale="tdb")
