"""Three-way validation against JPL Horizons for encounter-pair distances.

For each pair matched against the Fienga 2003 or Galád & Gray 2002 catalogs
(i.e. each row in ``data/output/fienga_2003_matches.csv`` and
``galad_2002_matches.csv``), query JPL Horizons for the heliocentric state
vectors of both asteroids in a dense window around the encounter epoch,
locate the true minimum-distance epoch, and report:

    pair  |  ours  |  literature  |  JPL  |  Δ(ours-JPL)  |  Δ(lit-JPL)

JPL Horizons is the highest-precision reference available (DE440 ephemeris +
numerical integration with full perturbation model), so its result is the
ground truth.  Both our pipeline output and the literature catalog are
compared against it.

Rate-limited: ``cfg.sources.jpl_horizons.rate_limit_seconds`` between requests.

Usage
-----
    docker compose run --rm pipeline python -m scripts.validate_jpl_horizons
"""

from __future__ import annotations

import argparse
import logging
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl
from astropy.time import Time

from src.utils.config import load_config

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)


def _horizons_vectors(asteroid_id: int, jd_center: float, half_window_days: float = 0.5):
    """Fetch heliocentric (location=@10, the Sun centre) state vectors around *jd_center*.

    Returns (jd_array, xyz_array) as numpy arrays (xyz in AU).
    """
    from astroquery.jplhorizons import Horizons

    t_center = Time(jd_center, format="jd", scale="tdb")
    epochs = {
        "start": (t_center - half_window_days).utc.iso[:19],
        "stop": (t_center + half_window_days).utc.iso[:19],
        "step": "1h",
    }
    h = Horizons(id=str(asteroid_id), location="@10", epochs=epochs, id_type="smallbody")
    tbl = h.vectors(refplane="ecliptic")
    jd = np.array(tbl["datetime_jd"], dtype=float)
    xyz = np.column_stack(
        [
            np.array(tbl["x"], dtype=float),
            np.array(tbl["y"], dtype=float),
            np.array(tbl["z"], dtype=float),
        ]
    )
    return jd, xyz


def _jpl_min_distance(
    a: int, b: int, jd_center: float, half_window_days: float = 1.0
) -> tuple[float, float]:
    """Return the minimum (distance_AU, jd_TDB) between *a* and *b* near *jd_center*."""
    jd_a, xyz_a = _horizons_vectors(a, jd_center, half_window_days)
    jd_b, xyz_b = _horizons_vectors(b, jd_center, half_window_days)
    # The two queries use the same epochs grid (server-side), so jd_a == jd_b.
    if not np.allclose(jd_a, jd_b):
        # Trim to common epochs by intersection — should never fire in practice.
        common = np.intersect1d(jd_a, jd_b)
        idx_a = np.array([np.where(jd_a == t)[0][0] for t in common])
        idx_b = np.array([np.where(jd_b == t)[0][0] for t in common])
        xyz_a = xyz_a[idx_a]
        xyz_b = xyz_b[idx_b]
        jd = common
    else:
        jd = jd_a
    d = np.linalg.norm(xyz_a - xyz_b, axis=1)
    k = int(np.argmin(d))
    return float(d[k]), float(jd[k])


def _load_matches(cfg) -> pl.DataFrame:
    """Load and unify Fienga + Galád match CSVs."""
    out_dir = Path(cfg.paths.output)
    fienga_csv = out_dir / "fienga_2003_matches.csv"
    galad_csv = out_dir / "galad_2002_matches.csv"

    frames: list[pl.DataFrame] = []
    if fienga_csv.exists():
        df = pl.read_csv(fienga_csv).rename(
            {"fienga_impact_au": "lit_dist_au", "fienga_date": "lit_date"}
        )
        df = df.with_columns(pl.lit("Fienga 2003").alias("source"))
        keep = [
            "source",
            "perturber",
            "target",
            "lit_date",
            "lit_dist_au",
            "our_date",
            "our_dist_au",
        ]
        df = df.select([c for c in keep if c in df.columns])
        frames.append(df)
    if galad_csv.exists():
        df = pl.read_csv(galad_csv).rename({"galad_r_au": "lit_dist_au", "galad_date": "lit_date"})
        df = df.with_columns(pl.lit("Galád 2002").alias("source"))
        keep = [
            "source",
            "perturber",
            "target",
            "lit_date",
            "lit_dist_au",
            "our_date",
            "our_dist_au",
        ]
        df = df.select([c for c in keep if c in df.columns])
        frames.append(df)
    if not frames:
        raise FileNotFoundError(
            f"No match CSV found under {out_dir}. Run validate_fienga_2003 and validate_galad_2002 first."
        )
    return pl.concat(frames, how="diagonal")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--config", default="config.yaml")
    p.add_argument(
        "--half-window-days",
        type=float,
        default=1.0,
        help="Half-width (days) of the dense window queried from Horizons around each "
        "encounter epoch. Default: 1.0 (i.e. 2-day window).",
    )
    args = p.parse_args()

    cfg = load_config(args.config)
    rate_limit = float(cfg.sources.jpl_horizons.rate_limit_seconds)

    matches = _load_matches(cfg)
    logger.info("Loaded %d match rows for JPL cross-check", len(matches))

    rows: list[dict] = []
    for i, row in enumerate(matches.iter_rows(named=True), 1):
        a, b = int(row["perturber"]), int(row["target"])
        our_dist = float(row["our_dist_au"])
        lit_dist = float(row["lit_dist_au"])
        # Use the date *we* found as the centre — it is day-precision unlike
        # Fienga's monthly approximation, so the dense JPL window catches the
        # true minimum without needing a 30-day-wide query.
        jd_center = float(Time(row["our_date"], scale="utc").tdb.jd)

        try:
            jpl_dist, jpl_jd = _jpl_min_distance(a, b, jd_center, args.half_window_days)
        except Exception as exc:
            logger.warning("  ✗ Horizons query failed for (%d, %d): %s", a, b, exc)
            jpl_dist, jpl_jd = float("nan"), float("nan")

        jpl_date = (
            Time(jpl_jd, format="jd", scale="tdb").utc.iso[:10] if np.isfinite(jpl_jd) else "—"
        )
        d_ours = our_dist - jpl_dist
        d_lit = lit_dist - jpl_dist
        rows.append(
            {
                "source": row["source"],
                "perturber": a,
                "target": b,
                "our_date": row["our_date"],
                "our_dist_au": our_dist,
                "lit_date": row["lit_date"],
                "lit_dist_au": lit_dist,
                "jpl_date": jpl_date,
                "jpl_dist_au": jpl_dist,
                "delta_ours_minus_jpl_au": d_ours,
                "delta_lit_minus_jpl_au": d_lit,
            }
        )
        logger.info(
            "[%d/%d] %-12s (%d, %d)  ours=%.5f  lit=%.5f  JPL=%.5f  Δours=%+.5f  Δlit=%+.5f",
            i,
            len(matches),
            row["source"],
            a,
            b,
            our_dist,
            lit_dist,
            jpl_dist,
            d_ours,
            d_lit,
        )
        if i < len(matches):
            time.sleep(rate_limit)

    report = pl.DataFrame(rows)
    out_path = Path(cfg.paths.output) / "jpl_horizons_validation.csv"
    report.write_csv(out_path)
    logger.info("Wrote %d rows to %s", len(report), out_path)

    finite = report.filter(pl.col("jpl_dist_au").is_finite())
    if len(finite) > 0:
        mae_ours_raw = finite["delta_ours_minus_jpl_au"].abs().mean()
        mae_lit_raw = finite["delta_lit_minus_jpl_au"].abs().mean()
        mae_ours = float(mae_ours_raw) if isinstance(mae_ours_raw, (int, float)) else 0.0
        mae_lit = float(mae_lit_raw) if isinstance(mae_lit_raw, (int, float)) else 0.0
        logger.info("Summary over %d successful queries:", len(finite))
        logger.info("  MAE(ours − JPL)        = %.5f AU", mae_ours)
        logger.info("  MAE(literature − JPL)  = %.5f AU", mae_lit)
    return 0


if __name__ == "__main__":
    sys.exit(main())
