"""Download Goffin (2014) asteroid **mass** determinations from VizieR.

Reference:
    Goffin E. (2014)
    "New determination of asteroid masses from close encounters"
    A&A 565, A56 — DOI: 10.1051/0004-6361/201322766
    VizieR catalog: J/A+A/565/A56

This is the companion to ``download_goffin_2014.py`` (which fetches the
*encounter* catalog). Here we fetch the **derived masses** that Goffin obtained
from those ground-based close encounters, for use as an independent ground-based
cross-check of our Gaia-FPR orbdet mass catalog.

The mass table is VizieR ``table5`` of the catalog. Its byte-by-byte description
gives, per asteroid:

    Seq    sequence number = MPC number (1=Ceres, 2=Pallas, ...)
    Name   asteroid name
    Nd     number of close encounters used
    M      mass             [units: 1e-10 solar masses]
    e_M    1-sigma error    [units: 1e-10 solar masses]
    Diam   diameter (km), Dens (g/cm^3), Cl (taxonomic class), ...

UNITS — CRITICAL
----------------
Goffin reports masses in units of **1e-10 solar masses**. We convert to kg via

    mass_kg = M[1e-10 Msun] * 1e-10 * M_SUN_KG

using ``M_SUN_KG`` from ``src.orbdet.constants`` (= GM_sun_SI / G ≈ 1.98841e30 kg)
so the conversion is consistent with the rest of the mass layer. Sanity check:
Ceres M=4.748 → 9.44e20 kg, Vesta M=1.30 → 2.58e20 kg.

Usage:
    docker compose run --rm pipeline python -m scripts.ingest.download_goffin_2014_masses
"""

from __future__ import annotations

import argparse
import json
import logging
from datetime import UTC, datetime
from pathlib import Path

import polars as pl

from src.orbdet.constants import M_SUN_KG

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)

_CATALOG_ID = "J/A+A/565/A56"
_MASS_TABLE = "J/A+A/565/A56/table5"
_OUT_PARQUET = Path("data/raw/goffin_2014_masses.parquet")
_OUT_META = Path("data/raw/goffin_2014_masses_metadata.json")

# Goffin table5 reports M and e_M in units of 1e-10 solar masses.
_GOFFIN_MASS_UNIT_MSUN = 1e-10


def download_goffin_2014_masses(
    out_parquet: Path = _OUT_PARQUET,
    out_meta: Path = _OUT_META,
) -> Path:
    """Fetch Goffin (2014) table5 masses from VizieR and write a parquet.

    Parameters
    ----------
    out_parquet : Path
        Destination parquet for the parsed mass table.
    out_meta : Path
        Destination JSON sidecar with provenance and unit convention.

    Returns
    -------
    Path
        Path to the written parquet file.
    """
    from astroquery.vizier import Vizier

    out_parquet.parent.mkdir(parents=True, exist_ok=True)

    logger.info("Querying VizieR for mass table %s …", _MASS_TABLE)
    v = Vizier(row_limit=-1, columns=["Seq", "Name", "Nd", "M", "e_M", "Diam", "Dens"])
    table_list = v.get_catalogs(_MASS_TABLE)
    if len(table_list) == 0:
        raise RuntimeError(f"VizieR returned no tables for {_MASS_TABLE}")
    tab = table_list[0]
    logger.info("Retrieved table with %d rows, cols=%s", len(tab), tab.colnames)

    rows: list[dict[str, object]] = []
    for r in tab:
        seq = r["Seq"]
        m_unit = r["M"]
        e_m_unit = r["e_M"]
        if seq is None or m_unit is None:
            continue
        mass_kg = float(m_unit) * _GOFFIN_MASS_UNIT_MSUN * M_SUN_KG
        sigma_kg = (
            float(e_m_unit) * _GOFFIN_MASS_UNIT_MSUN * M_SUN_KG if e_m_unit is not None else None
        )
        name = r["Name"]
        if isinstance(name, bytes):
            name = name.decode("ascii", errors="replace")
        rows.append(
            {
                "perturber": int(seq),
                "name": str(name).strip(),
                "n_encounters": int(r["Nd"]) if r["Nd"] is not None else None,
                "goffin_mass_1e10_msun": float(m_unit),
                "goffin_e_mass_1e10_msun": float(e_m_unit) if e_m_unit is not None else None,
                "goffin_mass_kg": mass_kg,
                "goffin_sigma_kg": sigma_kg,
                "goffin_diam_km": float(r["Diam"]) if r["Diam"] is not None else None,
                "goffin_dens_g_cm3": float(r["Dens"]) if r["Dens"] is not None else None,
            }
        )

    df = pl.DataFrame(rows).unique(subset=["perturber"], keep="first").sort("perturber")
    df.write_parquet(out_parquet, compression="zstd")
    logger.info("Wrote %d Goffin masses to %s", df.height, out_parquet)

    # Sanity check on the two canonical calibrators.
    for num, expect in ((1, 9.4e20), (4, 2.6e20)):
        sub = df.filter(pl.col("perturber") == num)
        if sub.height:
            got = sub["goffin_mass_kg"][0]
            ok = 0.5 * expect < got < 2.0 * expect
            logger.info(
                "Sanity: perturber %d mass = %.3e kg (expect ~%.1e) %s",
                num,
                got,
                expect,
                "OK" if ok else "FAIL — check unit scale!",
            )

    metadata = {
        "vizier_catalog": _CATALOG_ID,
        "vizier_table": _MASS_TABLE,
        "downloaded_at": datetime.now(UTC).isoformat(),
        "filename": out_parquet.name,
        "n_rows": df.height,
        "unit_convention": (
            "Goffin table5 reports M and e_M in units of 1e-10 solar masses; "
            f"converted to kg via mass_kg = M * 1e-10 * M_SUN_KG with "
            f"M_SUN_KG = {M_SUN_KG:.6e} kg (src.orbdet.constants)."
        ),
        "m_sun_kg": M_SUN_KG,
        "reference": "Goffin E. (2014) A&A 565, A56 — DOI: 10.1051/0004-6361/201322766",
    }
    out_meta.write_text(json.dumps(metadata, indent=2))
    logger.info("Metadata → %s", out_meta)
    return out_parquet


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Download Goffin (2014) asteroid mass determinations from VizieR"
    )
    parser.add_argument("--out-parquet", type=Path, default=_OUT_PARQUET)
    parser.add_argument("--out-meta", type=Path, default=_OUT_META)
    args = parser.parse_args()
    download_goffin_2014_masses(args.out_parquet, args.out_meta)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
