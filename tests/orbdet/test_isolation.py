"""Verifica el contrato de aislamiento de orbdet.

orbdet no debe importar ningún otro módulo ``src.*`` del proyecto (solo
stdlib + numpy/scipy/astropy/rebound + imports relativos dentro de orbdet).
Esto mantiene el motor testeable de forma independiente del pipeline.
"""

from __future__ import annotations

import ast
from pathlib import Path

_ORBDET_DIR = Path(__file__).resolve().parents[2] / "src" / "orbdet"


def _orbdet_py_files() -> list[Path]:
    return sorted(_ORBDET_DIR.glob("*.py"))


def test_orbdet_dir_exists() -> None:
    assert _ORBDET_DIR.is_dir(), f"no existe {_ORBDET_DIR}"
    assert _orbdet_py_files(), "no hay módulos .py en orbdet"


def test_no_cross_src_imports() -> None:
    """Ningún archivo de orbdet importa src.<algo distinto de orbdet>."""
    offenders: list[str] = []
    for path in _orbdet_py_files():
        tree = ast.parse(path.read_text(), filename=str(path))
        for node in ast.walk(tree):
            mods: list[str] = []
            if isinstance(node, ast.Import):
                mods = [alias.name for alias in node.names]
            elif isinstance(node, ast.ImportFrom):
                # Imports relativos (level>0) son internos de orbdet: OK.
                if node.level and node.level > 0:
                    continue
                if node.module:
                    mods = [node.module]
            for mod in mods:
                if mod == "src" or (mod.startswith("src.") and not mod.startswith("src.orbdet")):
                    offenders.append(f"{path.name}: import {mod}")

    assert (
        not offenders
    ), "orbdet rompe el aislamiento importando módulos del pipeline:\n" + "\n".join(offenders)


def test_only_whitelisted_third_party() -> None:
    """Las dependencias de terceros se limitan a la whitelist documentada."""
    allowed_top = {"numpy", "scipy", "astropy", "rebound"}
    stdlib_ok = {  # subconjunto de stdlib que el motor puede usar
        "__future__",
        "math",
        "dataclasses",
        "typing",
        "ast",
        "pathlib",
        "functools",
        "itertools",
        "collections",
        "logging",
        "warnings",
        "enum",
    }
    offenders: list[str] = []
    for path in _orbdet_py_files():
        tree = ast.parse(path.read_text(), filename=str(path))
        for node in ast.walk(tree):
            names: list[str] = []
            if isinstance(node, ast.Import):
                names = [alias.name for alias in node.names]
            elif isinstance(node, ast.ImportFrom) and not (node.level and node.level > 0):
                if node.module:
                    names = [node.module]
            for name in names:
                top = name.split(".")[0]
                if top in allowed_top or top in stdlib_ok or top == "src":
                    continue
                offenders.append(f"{path.name}: import {name}")
    assert not offenders, "orbdet importa dependencias fuera de whitelist:\n" + "\n".join(offenders)
