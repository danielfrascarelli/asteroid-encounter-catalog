"""Tests for src/utils/config.py."""

import textwrap
from pathlib import Path

import pytest

from src.utils.config import PipelineConfig, load_config

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def config_path(tmp_path: Path) -> Path:
    """Copy the real config.yaml to a temp directory."""
    import shutil

    src = Path(__file__).parent.parent / "config.yaml"
    dst = tmp_path / "config.yaml"
    shutil.copy(src, dst)
    return dst


# ---------------------------------------------------------------------------
# Happy-path tests
# ---------------------------------------------------------------------------


def test_load_returns_pipeline_config(config_path: Path) -> None:
    cfg = load_config(base_path=config_path, local_path=None)
    assert isinstance(cfg, PipelineConfig)


def test_default_values(config_path: Path) -> None:
    cfg = load_config(base_path=config_path, local_path=None)

    assert cfg.run.seed == 42
    assert cfg.detection.threshold_au == pytest.approx(0.01)
    assert cfg.time_window.start == "2014-07-25T00:00:00"
    assert cfg.time_window.scale == "utc"
    assert cfg.propagation.method == "kepler"
    assert cfg.subset.only_numbered is True
    assert cfg.subset.max_asteroids is None
    assert cfg.subset.semimajor_axis_au.min == pytest.approx(1.5)
    assert cfg.subset.semimajor_axis_au.max == pytest.approx(4.0)


def test_known_pairs(config_path: Path) -> None:
    cfg = load_config(base_path=config_path, local_path=None)
    perturbers = {kp.perturber for kp in cfg.validation.known_pairs}
    assert 1 in perturbers  # Ceres
    assert 4 in perturbers  # Vesta
    assert 2 in perturbers  # Pallas


def test_sources_config(config_path: Path) -> None:
    cfg = load_config(base_path=config_path, local_path=None)
    assert "minorplanetcenter" in cfg.sources.mpcorb.url
    assert cfg.sources.gaia_sso.table == "gaiadr3.sso_observation"
    assert "epoch" in cfg.sources.gaia_sso.columns


# ---------------------------------------------------------------------------
# Local override tests
# ---------------------------------------------------------------------------


def test_local_override_merges_deeply(config_path: Path, tmp_path: Path) -> None:
    local = tmp_path / "config.local.yaml"
    local.write_text(textwrap.dedent("""\
        detection:
          threshold_au: 0.05
        run:
          name: "test-run"
    """))
    cfg = load_config(base_path=config_path, local_path=local)
    # Overridden values
    assert cfg.detection.threshold_au == pytest.approx(0.05)
    assert cfg.run.name == "test-run"
    # Non-overridden values still present
    assert cfg.run.seed == 42
    assert cfg.propagation.method == "kepler"


def test_missing_local_override_is_ignored(config_path: Path, tmp_path: Path) -> None:
    cfg = load_config(
        base_path=config_path,
        local_path=tmp_path / "nonexistent.yaml",
    )
    assert isinstance(cfg, PipelineConfig)


# ---------------------------------------------------------------------------
# Error-path tests
# ---------------------------------------------------------------------------


def test_missing_base_config_raises(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        load_config(base_path=tmp_path / "nope.yaml", local_path=None)


def test_missing_required_key_raises(tmp_path: Path) -> None:
    broken = tmp_path / "config.yaml"
    broken.write_text("run:\n  name: x\n  description: y\n  seed: 1\n")
    with pytest.raises(ValueError, match="Missing required config key"):
        load_config(base_path=broken, local_path=None)
