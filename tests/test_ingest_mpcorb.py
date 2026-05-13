"""Tests for src/ingest/mpcorb.py."""

from pathlib import Path

import pytest

from src.ingest.mpcorb import parse_mpcorb, unpack_epoch
from tests.fixtures.generate_mpcorb import (
    CERES_A,
    HEADER,
    LINES,
    PALLAS_A,
    VESTA_A,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def mpcorb_file(tmp_path: Path) -> Path:
    """Write a minimal MPCORB.DAT with Ceres, Pallas, Vesta + extras."""
    content = (
        HEADER
        + LINES["ceres"]
        + "\n"
        + LINES["pallas"]
        + "\n"
        + LINES["vesta"]
        + "\n"
        + LINES["provisional"]
        + "\n"
        + LINES["trojan"]
        + "\n"
    )
    p = tmp_path / "MPCORB.DAT"
    p.write_text(content, encoding="ascii")
    return p


# ---------------------------------------------------------------------------
# unpack_epoch
# ---------------------------------------------------------------------------


def test_unpack_epoch_returns_float() -> None:
    assert isinstance(unpack_epoch("K205A"), float)


def test_unpack_epoch_k205a() -> None:
    # K205A = 2020-May-10 TT ≈ JD 2458979.5
    jd = unpack_epoch("K205A")
    assert 2458970.0 < jd < 2458990.0


def test_unpack_epoch_j991a() -> None:
    # J991A = 1999-Jan-10 TT ≈ JD 2451188.5
    jd = unpack_epoch("J991A")
    assert 2451188.0 < jd < 2451200.0


def test_unpack_epoch_month_oct() -> None:
    # Month 'A' = October
    jd_oct = unpack_epoch("K20A1")
    jd_sep = unpack_epoch("K2091")
    assert jd_oct > jd_sep  # October is after September


# ---------------------------------------------------------------------------
# parse_mpcorb — core
# ---------------------------------------------------------------------------


def test_parse_returns_dataframe(mpcorb_file: Path) -> None:
    import polars as pl

    df = parse_mpcorb(mpcorb_file)
    assert isinstance(df, pl.DataFrame)


def test_parse_expected_columns(mpcorb_file: Path) -> None:
    df = parse_mpcorb(mpcorb_file)
    expected = {
        "number", "designation", "H", "G", "epoch_jd",
        "M_deg", "omega_deg", "Omega_deg", "i_deg", "e", "a_au",
    }
    assert expected.issubset(set(df.columns))


# ---------------------------------------------------------------------------
# Ceres / Vesta / Pallas presence and accuracy
# ---------------------------------------------------------------------------


def test_ceres_present(mpcorb_file: Path) -> None:
    df = parse_mpcorb(mpcorb_file)
    assert 1 in df["number"].to_list()


def test_vesta_present(mpcorb_file: Path) -> None:
    df = parse_mpcorb(mpcorb_file)
    assert 4 in df["number"].to_list()


def test_pallas_present(mpcorb_file: Path) -> None:
    df = parse_mpcorb(mpcorb_file)
    assert 2 in df["number"].to_list()


def test_ceres_semimajor_axis(mpcorb_file: Path) -> None:
    df = parse_mpcorb(mpcorb_file)
    row = df.filter(df["number"] == 1)
    assert row["a_au"][0] == pytest.approx(CERES_A, abs=0.001)


def test_vesta_semimajor_axis(mpcorb_file: Path) -> None:
    df = parse_mpcorb(mpcorb_file)
    row = df.filter(df["number"] == 4)
    assert row["a_au"][0] == pytest.approx(VESTA_A, abs=0.001)


def test_pallas_semimajor_axis(mpcorb_file: Path) -> None:
    df = parse_mpcorb(mpcorb_file)
    row = df.filter(df["number"] == 2)
    assert row["a_au"][0] == pytest.approx(PALLAS_A, abs=0.001)


# ---------------------------------------------------------------------------
# Filtering
# ---------------------------------------------------------------------------


def test_only_numbered_excludes_provisional(mpcorb_file: Path) -> None:
    df = parse_mpcorb(mpcorb_file, only_numbered=True)
    # "J98S00A" is a provisional designation — should not appear
    designations = df["designation"].to_list()
    assert not any("2014 AA1" in d for d in designations)


def test_only_numbered_false_includes_provisional(mpcorb_file: Path) -> None:
    df = parse_mpcorb(mpcorb_file, only_numbered=False)
    designations = df["designation"].to_list()
    assert any("2014 AA1" in d for d in designations)


def test_semimajor_filter_excludes_trojan(mpcorb_file: Path) -> None:
    # Trojan at a≈5.2 AU should be outside the main belt filter
    df = parse_mpcorb(mpcorb_file, semimajor_min_au=1.5, semimajor_max_au=4.0)
    assert 588 not in df["number"].to_list()


def test_semimajor_filter_retains_main_belt(mpcorb_file: Path) -> None:
    df = parse_mpcorb(mpcorb_file, semimajor_min_au=1.5, semimajor_max_au=4.0)
    numbers = df["number"].to_list()
    assert 1 in numbers  # Ceres
    assert 2 in numbers  # Pallas
    assert 4 in numbers  # Vesta


# ---------------------------------------------------------------------------
# Error paths
# ---------------------------------------------------------------------------


def test_missing_file_raises(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        parse_mpcorb(tmp_path / "nonexistent.DAT")


def test_empty_file_returns_empty_df(tmp_path: Path) -> None:
    p = tmp_path / "empty.DAT"
    p.write_text("")
    df = parse_mpcorb(p)
    assert len(df) == 0
