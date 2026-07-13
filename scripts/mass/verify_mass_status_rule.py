"""Reproduce the Table 8 mass-status flags from the published identifiability rule.

The paper assigns each perturber a status --- ``measured`` / ``bound`` /
``non-physical`` / ``unconstrained`` --- from three published thresholds applied
to the jackknife (or, under high leverage, bootstrap) diagnostics:

    min_snr_jack = 3.0   (signal-to-noise of the mass over its external sigma)
    min_n_jack   = 10    (test bodies)
    max_leverage = 0.5   (top-1 jackknife leverage above which sigma_jack is not
                          defensible and the bootstrap sigma is used instead)

Rule (mirrors ``scripts/mass/build_mass_catalog.py`` main loop):

    mass <= 0                                  -> non_physical  (paper: non-physical)
    sigma_jack missing                         -> unknown       (paper: unconstrained)
    leverage > max_leverage and sigma_boot set -> measured if snr_boot >= 3 and N >= 10
                                                  else not_identifiable
    otherwise                                  -> measured if snr_jack >= 3 and
                                                  snr_jack_excl >= 3 and N >= 10
                                                  else not_identifiable

This script recomputes that status from the catalogue's diagnostic columns and
asserts it matches the stored ``mass_status``. It is the reproducible gate for
the identifiability finding (tribunal R2, M3): the published criteria must
regenerate the tabulated flags. The z-score against the seed/literature is
*deliberately not* part of the rule (it tests accuracy, not identifiability),
which is why (52) Europa is ``measured`` despite z = -4.3 against its DE441 seed.

Usage
-----
docker compose run --rm pipeline python -m scripts.mass.verify_mass_status_rule \
    --catalog data/output/orbdet/mass_catalog_b1fix_targets.csv
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
logger = logging.getLogger(__name__)

_DEFAULT_CATALOG = Path("data/output/orbdet/mass_catalog_b1fix_targets.csv")

# Presentation vocabulary used in the paper's Table 8.
_PAPER_STATUS = {
    "measured": "measured",
    "not_identifiable": "bound",
    "non_physical": "non-physical",
    "unknown": "unconstrained",
}


def _f(row: dict, key: str) -> float | None:
    """Parse a possibly-empty CSV cell as float, else None."""
    v = row.get(key, "")
    if v is None or str(v).strip() == "":
        return None
    try:
        return float(v)
    except ValueError:
        return None


def classify(row: dict, *, min_snr_jack: float, min_n_jack: int, max_leverage: float) -> str:
    """Return the engine status for one catalogue row from the published rule."""
    m = _f(row, "mass_fit_kg")
    s_jack = _f(row, "sigma_jack_kg")
    s_boot = _f(row, "sigma_boot_kg")
    snr_jack = _f(row, "snr_jack")
    snr_jack_excl = _f(row, "snr_jack_excl_top1")
    snr_boot = _f(row, "snr_boot")
    lev = _f(row, "jack_leverage_top1")
    n = _f(row, "n_targets")
    n = int(n) if n is not None else 0

    if m is None or m <= 0:
        return "non_physical"
    if s_jack is None:
        return "unknown"
    use_boot = s_boot is not None and lev is not None and lev > max_leverage
    if use_boot:
        ok = snr_boot is not None and snr_boot >= min_snr_jack and n >= min_n_jack
        return "measured" if ok else "not_identifiable"
    ok = (
        snr_jack is not None
        and snr_jack >= min_snr_jack
        and (snr_jack_excl is None or snr_jack_excl >= min_snr_jack)
        and n >= min_n_jack
    )
    return "measured" if ok else "not_identifiable"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", type=Path, default=_DEFAULT_CATALOG)
    parser.add_argument("--min-snr-jack", type=float, default=3.0)
    parser.add_argument("--min-n-jack", type=int, default=10)
    parser.add_argument("--max-leverage", type=float, default=0.5)
    args = parser.parse_args()

    rows = list(csv.DictReader(args.catalog.open()))
    n_mismatch = 0
    for r in rows:
        derived = classify(
            r,
            min_snr_jack=args.min_snr_jack,
            min_n_jack=args.min_n_jack,
            max_leverage=args.max_leverage,
        )
        stored = r.get("mass_status", "")
        flag = "" if derived == stored else "  <-- MISMATCH"
        if flag:
            n_mismatch += 1
        logger.info(
            "(%s) %-11s stored=%-16s derived=%-16s paper=%-12s%s",
            r.get("perturber", "?"),
            r.get("name", "?"),
            stored,
            derived,
            _PAPER_STATUS.get(derived, "?"),
            flag,
        )

    logger.info(
        "Status-rule reproduction: %d/%d rows match published criteria "
        "(min_snr_jack=%.1f, min_n_jack=%d, max_leverage=%.2f)",
        len(rows) - n_mismatch,
        len(rows),
        args.min_snr_jack,
        args.min_n_jack,
        args.max_leverage,
    )
    logger.info("%s", "PASS" if n_mismatch == 0 else "FAIL")
    return 0 if n_mismatch == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
