"""Check bulk-density plausibility of Table 6 candidate perturbers (tribunal B1 gate).

For each candidate perturber row ``(designation, D_km, M_kg)`` this script computes
the implied bulk density assuming a sphere of diameter ``D_km``:

    rho = M / ( (4/3) * pi * (D/2)^3 )      [g/cm^3]

(``D`` is converted km -> cm by x1e5; ``M`` is converted kg -> g by x1e3.) A row
passes the gate iff ``0.8 < rho < 4.5`` g/cm^3, the physically plausible range for
main-belt asteroid bulk densities (icy/rubble-pile to metallic). This is the
reproducible check for tribunal finding B1: mass/diameter pairs published for a
handful of Table 6 candidate perturbers implied densities well outside this
range, indicating a mass- or diameter-source bug.

Rows are read from a CSV via ``--in`` (columns: ``designation``, ``D_km``,
``M_kg``; a blank/missing ``M_kg`` is treated as "no mass available" and the
row is skipped rather than failed). A ``--demo`` mode is also provided with the
*current* (pre-fix) Table 6 values hardcoded, so the gate can be exercised
immediately without a CSV; it is expected to show several FAILs until Table 6
is corrected upstream, at which point this script should be re-run against the
corrected values and pass.

Usage
-----
    docker compose run --rm pipeline python -m scripts.validate.check_table6_density --demo
    docker compose run --rm pipeline python -m scripts.validate.check_table6_density \\
        --in data/output/table6_candidates.csv
"""

from __future__ import annotations

import argparse
import logging
import math
import sys
from pathlib import Path

import polars as pl

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
logger = logging.getLogger(__name__)

_KM_TO_CM = 1.0e5
_KG_TO_G = 1.0e3

_RHO_MIN_G_CM3 = 0.8
_RHO_MAX_G_CM3 = 4.5

# Current (pre-fix) Table 6 values. This is the exact input the B1 finding
# flagged as buggy: some implied densities fall outside the plausible range
# (e.g. Nysa too low, Fortuna too high). Angelina has no FM25 mass and is
# skipped rather than failed.
_DEMO_ROWS: list[tuple[str, float, float | None]] = [
    ("(9) Metis", 197.0, 6.5e18),
    ("(30) Urania", 109.0, 8.2e17),
    ("(40) Harmonia", 141.0, 2.0e18),
    ("(19) Fortuna", 133.0, 8.5e18),
    ("(21) Lutetia", 120.0, 1.7e18),
    ("(64) Angelina", 104.0, None),
    ("(44) Nysa", 140.0, 5.8e17),
    ("(20) Massalia", 178.0, 4.3e18),
    ("(29) Amphitrite", 240.0, 1.4e19),
    ("(27) Euterpe", 141.0, 1.5e18),
]


def bulk_density_g_cm3(d_km: float, m_kg: float) -> float:
    """Compute bulk density in g/cm^3 from diameter (km) and mass (kg).

    Parameters
    ----------
    d_km : float
        Diameter in kilometers, assuming a spherical body.
    m_kg : float
        Mass in kilograms.

    Returns
    -------
    float
        Bulk density in g/cm^3.
    """
    r_cm = (d_km * _KM_TO_CM) / 2.0
    volume_cm3 = (4.0 / 3.0) * math.pi * r_cm**3
    m_g = m_kg * _KG_TO_G
    return m_g / volume_cm3


def check_rows(rows: list[tuple[str, float, float | None]]) -> tuple[list[dict], bool]:
    """Compute density and PASS/FAIL/SKIP status for each candidate perturber row.

    Parameters
    ----------
    rows : list of (designation, D_km, M_kg)
        ``M_kg`` may be ``None`` (or NaN) when no mass estimate is available,
        in which case the row is marked ``SKIP`` and does not affect the
        overall pass/fail outcome.

    Returns
    -------
    (results, all_pass)
        ``results`` is a list of per-row dicts with keys ``designation``,
        ``D_km``, ``M_kg``, ``rho_g_cm3``, ``status``. ``all_pass`` is
        ``True`` iff no row is ``FAIL`` (``SKIP`` rows are ignored).
    """
    results: list[dict] = []
    all_pass = True
    for designation, d_km, m_kg in rows:
        if m_kg is None or (isinstance(m_kg, float) and math.isnan(m_kg)):
            logger.warning("%s: no mass available; skipping density check", designation)
            results.append(
                {
                    "designation": designation,
                    "D_km": d_km,
                    "M_kg": None,
                    "rho_g_cm3": None,
                    "status": "SKIP",
                }
            )
            continue

        rho = bulk_density_g_cm3(d_km, m_kg)
        passed = _RHO_MIN_G_CM3 < rho < _RHO_MAX_G_CM3
        all_pass = all_pass and passed
        results.append(
            {
                "designation": designation,
                "D_km": d_km,
                "M_kg": m_kg,
                "rho_g_cm3": rho,
                "status": "PASS" if passed else "FAIL",
            }
        )
    return results, all_pass


def _load_csv_rows(path: Path) -> list[tuple[str, float, float | None]]:
    """Read candidate perturber rows from a CSV with columns designation, D_km, M_kg."""
    df = pl.read_csv(path, null_values=["", "NA", "NaN", "None", "none"])
    required = {"designation", "D_km", "M_kg"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Input CSV {path} is missing required column(s): {sorted(missing)}")

    rows: list[tuple[str, float, float | None]] = []
    for row in df.iter_rows(named=True):
        m_kg = row["M_kg"]
        m_kg = None if m_kg is None else float(m_kg)
        rows.append((str(row["designation"]), float(row["D_km"]), m_kg))
    return rows


def print_table(results: list[dict]) -> None:
    """Print an aligned designation/D/M/rho/status table to stdout."""
    header = f"{'designation':<18} {'D_km':>8} {'M_kg':>12} {'rho_g_cm3':>12} {'status':>6}"
    print(header)
    print("-" * len(header))
    for r in results:
        m_str = f"{r['M_kg']:.3e}" if r["M_kg"] is not None else "n/a"
        rho_str = f"{r['rho_g_cm3']:.3f}" if r["rho_g_cm3"] is not None else "n/a"
        print(
            f"{r['designation']:<18} {r['D_km']:>8.1f} {m_str:>12} {rho_str:>12} "
            f"{r['status']:>6}"
        )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--in",
        dest="in_csv",
        type=Path,
        help="Input CSV with columns designation, D_km, M_kg.",
    )
    group.add_argument(
        "--demo",
        action="store_true",
        help="Use the hardcoded current (pre-fix) Table 6 values instead of a CSV.",
    )
    args = parser.parse_args()

    if args.demo:
        logger.info("Running in --demo mode with current (pre-fix) Table 6 values")
        rows = _DEMO_ROWS
    else:
        logger.info("Loading candidate perturber rows from %s", args.in_csv)
        rows = _load_csv_rows(args.in_csv)

    logger.info("Checking bulk density plausibility (%.1f < rho < %.1f g/cm^3) for %d rows",
                _RHO_MIN_G_CM3, _RHO_MAX_G_CM3, len(rows))

    results, all_pass = check_rows(rows)
    print_table(results)

    n_fail = sum(1 for r in results if r["status"] == "FAIL")
    n_skip = sum(1 for r in results if r["status"] == "SKIP")
    logger.info(
        "Result: %d/%d PASS, %d FAIL, %d SKIP",
        len(results) - n_fail - n_skip,
        len(results),
        n_fail,
        n_skip,
    )

    if not all_pass:
        logger.error("Bulk-density gate FAILED: %d row(s) outside plausible range", n_fail)
        return 1

    logger.info("Bulk-density gate PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
