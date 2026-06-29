"""Helper de tests: skip cuando las efemérides DE440 de ASSIST no están presentes.

Los archivos ``linux_p1550p2650.440`` (~100 MB) y ``sb441-n16.bsp`` (~650 MB) no se
versionan ni viven en CI. Los tests que ejercitan el backend ASSIST se saltan
automáticamente si faltan, y corren donde sí están (local con ``ORBDET_EPHEM_DIR``).
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

_PLANETS_FILE = "linux_p1550p2650.440"
_ASTEROIDS_FILE = "sb441-n16.bsp"


def ephem_available() -> bool:
    base = Path(os.environ.get("ORBDET_EPHEM_DIR", os.path.join("data", "raw", "ephem")))
    return (base / _PLANETS_FILE).exists() and (base / _ASTEROIDS_FILE).exists()


requires_ephem = pytest.mark.skipif(
    not ephem_available(),
    reason="efemérides DE440 ausentes (definir ORBDET_EPHEM_DIR con .440 + sb441-n16.bsp)",
)
