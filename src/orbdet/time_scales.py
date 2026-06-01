"""Conversiones entre escalas de tiempo (JD), apoyadas en astropy.

El motor trabaja internamente en **JD TDB**. Gaia reporta épocas en **TCB**
(como días desde J2010.0 TCB, ver más abajo) y los inputs de usuario suelen ser
**UTC ISO**. Estas funciones hacen las conversiones explícitas; nunca mezclar
floats de JD sin saber su escala.

A diferencia de la aproximación lineal de ``src/utils/time_utils.py``, acá se usa
``astropy.time.Time`` para máxima exactitud (objetivo: resultados publicables).
"""

from __future__ import annotations

import numpy as np
from astropy.time import Time

# Época de referencia que usa la tabla SSO de Gaia (DR3 y FPR): el campo ``epoch``
# es días desde J2010.0 TCB, i.e. epoch = JD_TCB - J2010_TCB_JD.
J2010_TCB_JD: float = 2_455_197.5


def _to_array(x: float | np.ndarray) -> np.ndarray:
    return np.atleast_1d(np.asarray(x, dtype=float))


def tcb_to_tdb(jd_tcb: float | np.ndarray) -> np.ndarray:
    """JD en TCB → JD en TDB."""
    return np.asarray(Time(_to_array(jd_tcb), format="jd", scale="tcb").tdb.jd, dtype=float)


def tdb_to_tcb(jd_tdb: float | np.ndarray) -> np.ndarray:
    """JD en TDB → JD en TCB."""
    return np.asarray(Time(_to_array(jd_tdb), format="jd", scale="tdb").tcb.jd, dtype=float)


def utc_to_tdb(jd_utc: float | np.ndarray) -> np.ndarray:
    """JD en UTC → JD en TDB."""
    return np.asarray(Time(_to_array(jd_utc), format="jd", scale="utc").tdb.jd, dtype=float)


def tdb_to_utc(jd_tdb: float | np.ndarray) -> np.ndarray:
    """JD en TDB → JD en UTC."""
    return np.asarray(Time(_to_array(jd_tdb), format="jd", scale="tdb").utc.jd, dtype=float)


def iso_utc_to_jd_tdb(iso: str) -> float:
    """String ISO 8601 en UTC → JD TDB (escalar)."""
    return float(Time(iso, scale="utc").tdb.jd)


def gaia_epoch_to_jd_tdb(
    epoch_days: float | np.ndarray,
    epoch_ref_jd_tcb: float = J2010_TCB_JD,
) -> np.ndarray:
    """Época de Gaia (días desde J2010.0 TCB) → JD TDB.

    ``jd_tcb = epoch_days + epoch_ref_jd_tcb``; luego TCB → TDB.
    ``epoch_ref_jd_tcb`` es idéntico en DR3 y FPR (ver docs/gaia_fpr_data_model.md).
    """
    jd_tcb = _to_array(epoch_days) + epoch_ref_jd_tcb
    return tcb_to_tdb(jd_tcb)
