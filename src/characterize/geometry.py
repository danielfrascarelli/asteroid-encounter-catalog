"""Unit conversion helpers for encounter geometry."""

from __future__ import annotations

import numpy as np

# 1 AU = 149 597 870.7 km  (IAU 2012)
AU_KM: float = 149_597_870.7
# 1 AU/day in km/s
_AU_DAY_TO_KM_S: float = AU_KM / 86_400.0


def dist_au_to_km(d: float | np.ndarray) -> np.ndarray:
    """Convert distance from AU to km."""
    return np.asarray(d, dtype=float) * AU_KM


def vel_au_per_day_to_km_s(v: float | np.ndarray) -> np.ndarray:
    """Convert velocity from AU/day to km/s."""
    return np.asarray(v, dtype=float) * _AU_DAY_TO_KM_S


def vel_au_per_day_to_m_s(v: float | np.ndarray) -> np.ndarray:
    """Convert velocity from AU/day to m/s."""
    return np.asarray(v, dtype=float) * _AU_DAY_TO_KM_S * 1_000.0
