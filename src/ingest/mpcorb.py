"""Parser for the MPC Orbit Database (MPCORB.DAT) fixed-width format.

MPCORB.DAT column layout (1-based, per MPC documentation):
  1-  7  A7    Number or packed provisional designation
  9- 13  F5.2  H (absolute magnitude)
 15- 19  F5.3  G (slope parameter)
 21- 25  A5    Epoch (MPC packed date, .0 TT)
 27- 35  F9.5  M  — mean anomaly at epoch (degrees)
 37- 46  F10.5 ω  — argument of perihelion, J2000 (degrees)
 48- 57  F10.5 Ω  — longitude of ascending node, J2000 (degrees)
 59- 68  F10.5 i  — inclination to ecliptic, J2000 (degrees)
 70- 79  F10.7 e  — orbital eccentricity
 81- 91  F11.8 n  — mean daily motion (deg/day)
 93-103  F11.7 a  — semimajor axis (AU)
105      A1    U  — uncertainty parameter
167-194  A28   Readable designation
"""

from __future__ import annotations

import logging
from pathlib import Path

import polars as pl
from astropy.time import Time

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Packed epoch decoding
# ---------------------------------------------------------------------------

_MONTH: dict[str, int] = {str(i): i for i in range(1, 10)}
_MONTH.update({"A": 10, "B": 11, "C": 12})

_DAY: dict[str, int] = {str(i): i for i in range(1, 10)}
_DAY.update({c: 10 + j for j, c in enumerate("ABCDEFGHIJKLMNOPQRSTUV")})

_CENTURY: dict[str, int] = {"I": 1800, "J": 1900, "K": 2000}


def unpack_epoch(packed: str) -> float:
    """Convert a 5-char MPC packed epoch string to JD (TDB).

    Format: <century><YY><month><day>
    where century ∈ {I,J,K}, month and day use alphanumeric encoding.

    Parameters
    ----------
    packed:
        5-character packed epoch, e.g. ``"K205A"`` → 2020-May-10 TT.

    Returns
    -------
    float
        Julian Date in TDB scale.
    """
    packed = packed.strip()
    year = _CENTURY[packed[0]] + int(packed[1:3])
    month = _MONTH[packed[3]]
    day = _DAY[packed[4]]
    t = Time(f"{year:04d}-{month:02d}-{day:02d}", format="iso", scale="tt")
    return float(t.tdb.jd)


# ---------------------------------------------------------------------------
# Line classification helpers
# ---------------------------------------------------------------------------


def _is_numbered(no_field: str) -> bool:
    """Return True if the number field contains a plain integer."""
    try:
        int(no_field.strip())
        return True
    except ValueError:
        return False


def _is_data_line(line: str) -> bool:
    """Return True if the line has enough columns and a parseable semimajor axis."""
    if len(line) < 103:
        return False
    try:
        float(line[92:103])
        return True
    except ValueError:
        return False


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def parse_mpcorb(
    path: str | Path,
    *,
    only_numbered: bool = True,
    semimajor_min_au: float = 0.0,
    semimajor_max_au: float = float("inf"),
) -> pl.DataFrame:
    """Parse MPCORB.DAT and return a typed Polars DataFrame.

    Parameters
    ----------
    path:
        Path to the MPCORB.DAT file (plain text, not compressed).
    only_numbered:
        If True (default), skip provisional/unnumbered objects.
    semimajor_min_au:
        Minimum semimajor axis to include (AU).
    semimajor_max_au:
        Maximum semimajor axis to include (AU).

    Returns
    -------
    polars.DataFrame
        Columns: ``number`` (Int32), ``designation`` (Utf8), ``H`` (Float64),
        ``G`` (Float64), ``epoch_jd`` (Float64), ``M_deg``, ``omega_deg``,
        ``Omega_deg``, ``i_deg``, ``e``, ``a_au`` (all Float64).
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"MPCORB file not found: {path}")

    logger.info("Parsing MPCORB from %s", path)

    rows: list[dict] = []
    skipped_unnumbered = 0

    with path.open(encoding="ascii", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not _is_data_line(line):
                continue

            no_field = line[0:7]
            if only_numbered and not _is_numbered(no_field):
                skipped_unnumbered += 1
                continue

            a = float(line[92:103])
            if not (semimajor_min_au <= a <= semimajor_max_au):
                continue

            try:
                epoch_jd = unpack_epoch(line[20:25])
            except (KeyError, ValueError, IndexError):
                logger.debug("Skipping line with unparseable epoch: %.40s…", line)
                continue

            number = int(no_field.strip()) if _is_numbered(no_field) else None
            designation = line[166:194].strip() if len(line) >= 194 else no_field.strip()

            rows.append(
                {
                    "number": number,
                    "designation": designation,
                    "H": float(line[8:13]) if line[8:13].strip() else None,
                    "G": float(line[14:19]) if line[14:19].strip() else None,
                    "epoch_jd": epoch_jd,
                    "M_deg": float(line[26:35]),
                    "omega_deg": float(line[36:46]),
                    "Omega_deg": float(line[47:57]),
                    "i_deg": float(line[58:68]),
                    "e": float(line[69:79]),
                    "a_au": a,
                }
            )

    if skipped_unnumbered:
        logger.debug("Skipped %d unnumbered objects", skipped_unnumbered)

    df = pl.DataFrame(
        rows,
        schema={
            "number": pl.Int32,
            "designation": pl.Utf8,
            "H": pl.Float64,
            "G": pl.Float64,
            "epoch_jd": pl.Float64,
            "M_deg": pl.Float64,
            "omega_deg": pl.Float64,
            "Omega_deg": pl.Float64,
            "i_deg": pl.Float64,
            "e": pl.Float64,
            "a_au": pl.Float64,
        },
    )

    logger.info(
        "Loaded %d objects from MPCORB (filters: numbered=%s, a=[%.2f, %.2f])",
        len(df),
        only_numbered,
        semimajor_min_au,
        semimajor_max_au,
    )
    return df
