"""Tests de src/orbdet/frames.py."""

from __future__ import annotations

import math

import numpy as np
import pytest

from src.orbdet import frames as fr
from src.orbdet.constants import OBLIQUITY_J2000_RAD


@pytest.mark.parametrize("rot", [fr.rotation_x, fr.rotation_y, fr.rotation_z])
def test_rotation_matrices_are_orthonormal(rot) -> None:
    R = rot(0.7)
    assert np.allclose(R @ R.T, np.eye(3), atol=1e-14)
    assert np.linalg.det(R) == pytest.approx(1.0, abs=1e-14)


def test_rotation_z_known_angle() -> None:
    # Rz(90°) lleva x̂ → ŷ
    out = fr.rotation_z(math.pi / 2.0) @ np.array([1.0, 0.0, 0.0])
    assert np.allclose(out, [0.0, 1.0, 0.0], atol=1e-12)


def test_rotation_x_known_angle() -> None:
    # Rx(90°) lleva ŷ → ẑ
    out = fr.rotation_x(math.pi / 2.0) @ np.array([0.0, 1.0, 0.0])
    assert np.allclose(out, [0.0, 0.0, 1.0], atol=1e-12)


def test_ecliptic_equatorial_roundtrip() -> None:
    v = np.array([0.3, -0.7, 0.5])
    back = fr.equatorial_to_ecliptic(fr.ecliptic_to_equatorial(v))
    assert np.allclose(back, v, atol=1e-14)


def test_x_axis_invariant_under_obliquity() -> None:
    # La rotación es alrededor de X → x̂ no cambia.
    out = fr.ecliptic_to_equatorial(np.array([1.0, 0.0, 0.0]))
    assert np.allclose(out, [1.0, 0.0, 0.0], atol=1e-14)


def test_ecliptic_pole_maps_by_obliquity() -> None:
    # El polo eclíptico ẑ_ecl pasa a (0, -sinε, cosε) en ecuatorial.
    out = fr.ecliptic_to_equatorial(np.array([0.0, 0.0, 1.0]))
    eps = OBLIQUITY_J2000_RAD
    assert np.allclose(out, [0.0, -math.sin(eps), math.cos(eps)], atol=1e-12)


def test_batch_application() -> None:
    vs = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    out = fr.ecliptic_to_equatorial(vs)
    assert out.shape == (2, 3)
    assert np.allclose(out[0], [1.0, 0.0, 0.0], atol=1e-14)
