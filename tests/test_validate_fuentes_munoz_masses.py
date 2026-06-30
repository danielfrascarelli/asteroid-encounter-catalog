"""Unit tests for the Fuentes-Muñoz 2025 mass cross-check parser.

Self-contained: builds a tiny synthetic MRT snippet so the test never depends on
the gitignored ``data/raw/fuentes_munoz_2025/`` download and runs in CI.
"""

from __future__ import annotations

from pathlib import Path

from scripts.validate.validate_fuentes_munoz_masses import (
    _G_KM3_KG_S2,
    parse_table5_masses,
)

# A minimal MRT: the parser keys off the LAST '---' divider, then reads
# fixed-width byte columns. Columns (1-indexed): Asteroid 1-22, GMfin 66-76,
# e_GMfin 78-88. We pad each row to put GMfin/e_GMfin at the right bytes.
_HEADER = """Title: synthetic
--------------------------------------------------------------------------------
   Bytes Format Units     Label
--------------------------------------------------------------------------------
"""


def _row(asteroid: str, gm_fin: str, e_gm_fin: str) -> str:
    line = list(" " * 88)
    line[0 : len(asteroid)] = asteroid
    line[65 : 65 + len(gm_fin)] = gm_fin
    line[77 : 77 + len(e_gm_fin)] = e_gm_fin
    return "".join(line)


def test_parse_and_unit_conversion(tmp_path: Path) -> None:
    mrt = tmp_path / "t5.txt"
    mrt.write_text(
        _HEADER
        + _row("1 Ceres", "62.63917624", "0.01854053")
        + "\n"
        + _row("16 Psyche", "1.5980", "0.05")
        + "\n"
        + _row("2013 KY18", "9.9", "0.1")
        + "\n"  # provisional → dropped
        + _row("999 Blank", "", "")
        + "\n"  # no GMfin → dropped
    )

    df = parse_table5_masses(mrt)

    # Provisional designation and the blank-GM row are excluded.
    assert sorted(df["perturber"].to_list()) == [1, 16]

    # GM → mass conversion: M = GM / G.
    ceres = df.filter(df["perturber"] == 1).to_dicts()[0]
    expected_ceres_kg = 62.63917624 / _G_KM3_KG_S2
    assert abs(ceres["fm_mass_kg"] - expected_ceres_kg) / expected_ceres_kg < 1e-9
    # Ceres mass parameter lands near the known ~9.38e20 kg.
    assert 9.3e20 < ceres["fm_mass_kg"] < 9.5e20

    psyche = df.filter(df["perturber"] == 16).to_dicts()[0]
    assert abs(psyche["fm_mass_kg"] - 1.5980 / _G_KM3_KG_S2) / psyche["fm_mass_kg"] < 1e-9
    assert abs(psyche["fm_sigma_kg"] - 0.05 / _G_KM3_KG_S2) / psyche["fm_sigma_kg"] < 1e-9
