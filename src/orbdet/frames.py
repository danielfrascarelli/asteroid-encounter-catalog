"""Rotaciones y transformaciones de marco de referencia.

Convenciones:
- Vectores como arrays ``(3,)`` o ``(N, 3)`` (fila = vector).
- Rotaciones activas: ``R @ v`` rota el vector dentro del mismo marco.
- ``rotation_x/y/z(theta)`` rotan **alrededor** del eje indicado, en radianes.

El marco interno del propagador kepleriano es **eclíptico** (la eclíptica J2000
es natural para elementos orbitales); Gaia reporta en **ecuatorial ICRS**. Las
conversiones eclíptica↔ecuatorial son una rotación alrededor del eje X por la
oblicuidad.
"""

from __future__ import annotations

import numpy as np

from .constants import OBLIQUITY_J2000_RAD


def rotation_x(theta: float) -> np.ndarray:
    """Matriz de rotación 3×3 alrededor del eje X por *theta* (rad)."""
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])


def rotation_y(theta: float) -> np.ndarray:
    """Matriz de rotación 3×3 alrededor del eje Y por *theta* (rad)."""
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]])


def rotation_z(theta: float) -> np.ndarray:
    """Matriz de rotación 3×3 alrededor del eje Z por *theta* (rad)."""
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def _apply(matrix: np.ndarray, vec: np.ndarray) -> np.ndarray:
    """Aplica *matrix* (3×3) a un vector ``(3,)`` o a un lote ``(N, 3)``."""
    v = np.asarray(vec, dtype=float)
    result = matrix @ v if v.ndim == 1 else v @ matrix.T
    return np.asarray(result, dtype=float)


def ecliptic_to_equatorial(
    vec: np.ndarray, obliquity_rad: float = OBLIQUITY_J2000_RAD
) -> np.ndarray:
    """Eclíptica J2000 → ecuatorial ICRS (rotación +obliquity alrededor de X)."""
    return _apply(rotation_x(obliquity_rad), vec)


def equatorial_to_ecliptic(
    vec: np.ndarray, obliquity_rad: float = OBLIQUITY_J2000_RAD
) -> np.ndarray:
    """Ecuatorial ICRS → eclíptica J2000 (rotación -obliquity alrededor de X)."""
    return _apply(rotation_x(-obliquity_rad), vec)
