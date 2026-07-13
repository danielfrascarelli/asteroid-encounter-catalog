"""Test de inyección-recuperación del estimador de masas (B8).

Inyecta una masa de perturbador conocida, genera la astrometría sintética que esa
masa produce sobre varios objetivos (con ruido along-scan) y verifica que el ajuste
conjunto masa+órbita la recupera sin sesgo (|z| acotado) con χ²_red sano. Es el gate
que detecta un sesgo del estimador (distinto del anclaje de exactitud vía
calibradores, que valida el modelo de fuerzas).
"""

from __future__ import annotations

import argparse
import math

import pytest

from scripts.validate.injection_recovery_mass import run


@pytest.mark.slow
def test_injection_recovery_unbiased() -> None:
    """GATE B8: masa inyectada recuperada dentro de 4σ con χ²_red sano."""
    args = argparse.Namespace(
        n_targets=10,
        m_inj_kg=2e19,
        noise_mas=2.0,
        n_sigma=4.0,
        seed=42,
        jackknife=False,
    )
    assert run(args) == 0


@pytest.mark.slow
def test_injection_recovery_second_seed() -> None:
    """Insesgado también en otra semilla/masa (no un artefacto de una realización)."""
    args = argparse.Namespace(
        n_targets=10,
        m_inj_kg=5e19,
        noise_mas=1.5,
        n_sigma=4.0,
        seed=7,
        jackknife=False,
    )
    assert run(args) == 0


def test_injection_recovery_signature() -> None:
    """Sanity barato: la firma esperada existe (no ejecuta el ajuste)."""
    assert callable(run)
    assert math.isfinite(2e19)
