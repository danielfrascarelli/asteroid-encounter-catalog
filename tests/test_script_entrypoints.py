"""Smoke tests that every CLI script under ``scripts/<subpkg>/`` *imports*.

Goal: catch the class of regression we just hit, where two scripts still
referenced the pre-reorganisation flat path ``scripts.step_model_test``
after the reorg moved it to ``scripts.dev.step_model_test``.  Those bugs
only fired when the user actually ran the entrypoint — the unit-test suite
never imported them.

We test ``importlib.import_module`` instead of ``python -m … --help`` so the
check works uniformly for scripts that don't expose argparse (running them
with ``--help`` falls through to ``main()`` and may hang on network calls
or run heavy bench code).  Importing only resolves module-level statements,
which is exactly what catches broken imports without doing any work.
"""

from __future__ import annotations

import importlib
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent
_SCRIPTS_ROOT = _REPO_ROOT / "scripts"

_SUBPKGS = ("ingest", "pipeline", "mass", "validate", "bench", "dev")


def _discover_modules() -> list[str]:
    mods: list[str] = []
    for sub in _SUBPKGS:
        d = _SCRIPTS_ROOT / sub
        if not d.is_dir():
            continue
        for p in sorted(d.glob("*.py")):
            if p.name == "__init__.py":
                continue
            mods.append(f"scripts.{sub}.{p.stem}")
    return mods


@pytest.mark.parametrize("module", _discover_modules())
def test_script_module_imports(module: str) -> None:
    """``import <module>`` must succeed.

    Catches stale flat imports (``scripts.foo`` instead of
    ``scripts.<subpkg>.foo``), missing dependencies, syntax errors at module
    scope, and circular imports.  Does not exercise the script's main logic.
    """
    importlib.import_module(module)
