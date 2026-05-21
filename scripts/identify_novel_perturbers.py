"""Identify which detected perturbers are GENUINELY novel mass candidates.

Cross-checks the perturbers in ``data/output/deflection_detections.csv``
against the literature encounter catalogs already loaded in the pipeline:

  - Fienga et al. (2003), VizieR J/A+A/406/751
  - Galád & Gray (2002), parsed from A&A 391, 1115

Perturbers that appear in either catalog have a published mass (or at least
were considered for mass determination). Perturbers that DON'T appear in
either are candidates for genuinely new mass measurements.

Output: ``data/output/novel_vs_known_perturbers.csv``

Usage
-----
    docker compose run --rm pipeline python -m scripts.identify_novel_perturbers
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import polars as pl

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)


def _load_literature_perturbers(raw_dir: Path) -> tuple[set[int], set[int]]:
    """Load the set of MPC numbers that appear as perturbers in Fienga / Galad."""
    fienga_path = raw_dir / "fienga_2003_encounters.parquet"
    galad_path = raw_dir / "galad_2002_encounters.parquet"

    fienga: set[int] = set()
    if fienga_path.exists():
        df = pl.read_parquet(fienga_path)
        for x in df["Perturber"].to_list():
            if x is not None:
                fienga.add(int(x))
        logger.info("Fienga 2003: %d unique perturbers", len(fienga))
    else:
        logger.warning("Fienga 2003 parquet not found at %s", fienga_path)

    galad: set[int] = set()
    if galad_path.exists():
        df = pl.read_parquet(galad_path)
        for x in df["perturber_number"].to_list():
            if x is not None:
                galad.add(int(x))
        logger.info("Galad 2002: %d unique perturbers", len(galad))
    else:
        logger.warning("Galad 2002 parquet not found at %s", galad_path)

    return fienga, galad


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument(
        "--detections",
        type=Path,
        default=Path("data/output/deflection_detections.csv"),
    )
    p.add_argument(
        "--raw",
        type=Path,
        default=Path("data/raw"),
    )
    p.add_argument(
        "--output",
        type=Path,
        default=Path("data/output/novel_vs_known_perturbers.csv"),
    )
    args = p.parse_args()

    if not args.detections.exists():
        logger.error("Detections CSV not found: %s", args.detections)
        return 1

    fienga, galad = _load_literature_perturbers(args.raw)
    literature_union = fienga | galad
    logger.info("Total literature perturbers (Fienga ∪ Galad): %d", len(literature_union))

    det = pl.read_csv(args.detections)
    detected = det.filter(pl.col("detection") == "yes")
    logger.info(
        "Detected perturbers in our pipeline (≥3σ): %d unique",
        detected["perturber_number"].n_unique(),
    )

    perturbers = detected.select(["perturber_number", "perturber_name"]).unique(
        subset=["perturber_number"]
    )

    rows: list[dict] = []
    for r in perturbers.iter_rows(named=True):
        num = int(r["perturber_number"])
        in_fienga = num in fienga
        in_galad = num in galad
        rows.append(
            {
                "perturber_number": num,
                "perturber_name": r["perturber_name"],
                "in_fienga_2003": in_fienga,
                "in_galad_2002": in_galad,
                "in_literature": in_fienga or in_galad,
                "genuinely_novel": not (in_fienga or in_galad),
            }
        )

    out = pl.DataFrame(rows).sort("perturber_number")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.write_csv(args.output)
    logger.info("Wrote %d rows to %s", out.height, args.output)

    n_total = out.height
    n_novel = int(out["genuinely_novel"].sum())
    n_known = n_total - n_novel
    logger.info("")
    logger.info("Summary:")
    logger.info("  Total detected perturbers:    %d", n_total)
    logger.info("  Already in literature:        %d (Fienga/Galad)", n_known)
    logger.info("  GENUINELY NOVEL:              %d", n_novel)
    logger.info("")
    logger.info("Already in literature (calibration set):")
    for r in out.filter(pl.col("in_literature")).iter_rows(named=True):
        srcs = []
        if r["in_fienga_2003"]:
            srcs.append("Fienga")
        if r["in_galad_2002"]:
            srcs.append("Galád")
        logger.info(
            "  (%d) %-20s  [%s]", r["perturber_number"], r["perturber_name"], "+".join(srcs)
        )
    logger.info("")
    logger.info("Genuinely novel (potential new mass determinations):")
    for r in out.filter(pl.col("genuinely_novel")).iter_rows(named=True):
        logger.info("  (%d) %s", r["perturber_number"], r["perturber_name"])

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
