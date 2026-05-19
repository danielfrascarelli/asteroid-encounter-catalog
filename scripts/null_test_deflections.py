"""NULL TEST: verify that the encounter-coincident shifts are not coincidence.

A real perturbation should produce a shift coincident with the *actual* encounter
date. If our detection method is biased (e.g. by linear orbital drift in
Horizons's target orbit), we'd see similar shifts at *arbitrary* dates.

This script repeats the detection but using *shuffled* encounter dates: each
candidate is assigned the date of a different candidate's encounter. The Gaia
observations of the original target are used, but split into "before/after"
using the wrong date.

If the method is real, the |t| distribution of the shuffled run should be
near the χ² null (most |t| < 3, few outliers). If the original |t|
distribution is much higher than the shuffled, the encounter coincidence is
real.

Output: ``data/output/deflection_null_test.csv`` plus printed comparison.

Usage
-----
    docker compose run --rm pipeline python -m scripts.null_test_deflections
"""

from __future__ import annotations

import argparse
import logging
import math
import random
import time
from pathlib import Path

import numpy as np
import polars as pl
from astropy.time import Time
from astroquery.jplhorizons import Horizons
from astroquery.utils.tap.core import TapPlus

from src.utils.config import load_config

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)


_HALF_WINDOW_DAYS = 180.0
_BLACKOUT_DAYS = 7.0
_J2010_TCB_JD = 2455197.5
_GAIA_OBSERVER = "500@-139479"


def fetch_gaia(archive_url: str, target: int, d_min: float, d_max: float) -> pl.DataFrame:
    adql = (
        "SELECT number_mp, epoch, ra, dec "
        "FROM gaiadr3.sso_observation "
        f"WHERE number_mp = {target} "
        f"AND epoch BETWEEN {d_min:.6f} AND {d_max:.6f} "
        "ORDER BY epoch"
    )
    tap = TapPlus(url=archive_url)
    job = tap.launch_job_async(adql)
    tbl = job.get_results()
    df = pl.from_pandas(tbl.to_pandas())
    df = df.rename({c: c.lower() for c in df.columns if c != c.lower()})
    return df


def horizons_apparent_radec(
    target: int, jd_tdb: np.ndarray, rate_limit_s: float
) -> tuple[np.ndarray, np.ndarray]:
    ra_out = np.empty(len(jd_tdb), dtype=float)
    dec_out = np.empty(len(jd_tdb), dtype=float)
    chunk = 50
    for start in range(0, len(jd_tdb), chunk):
        idx = np.arange(start, min(start + chunk, len(jd_tdb)))
        epochs = jd_tdb[idx].tolist()
        h = Horizons(
            id=str(target),
            location=_GAIA_OBSERVER,
            epochs=epochs,
            id_type="smallbody",
        )
        eph = h.ephemerides()
        ra_out[idx] = np.array(eph["RA"], dtype=float)
        dec_out[idx] = np.array(eph["DEC"], dtype=float)
        time.sleep(rate_limit_s)
    return ra_out, dec_out


def _t_welch(vals: np.ndarray, mask_before: np.ndarray, mask_after: np.ndarray) -> float:
    nb = int(mask_before.sum())
    na = int(mask_after.sum())
    if nb < 2 or na < 2:
        return float("nan")
    mb = float(np.mean(vals[mask_before]))
    ma = float(np.mean(vals[mask_after]))
    sb = float(np.std(vals[mask_before]))
    sa = float(np.std(vals[mask_after]))
    se = math.sqrt((sb * sb) / nb + (sa * sa) / na)
    return (ma - mb) / se if se > 0 else float("nan")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--config", default="config.yaml")
    p.add_argument(
        "--candidates",
        type=Path,
        default=Path("data/output/publishable_mass_candidates.csv"),
    )
    p.add_argument(
        "--top",
        type=int,
        default=10,
        help="Top-N candidates to null-test (each will use another candidate's date).",
    )
    p.add_argument("--seed", type=int, default=42)
    p.add_argument(
        "--output",
        type=Path,
        default=Path("data/output/deflection_null_test.csv"),
    )
    args = p.parse_args()

    cfg = load_config(args.config)
    archive_url = cfg.sources.gaia_sso.archive_url
    rate_limit = float(cfg.sources.jpl_horizons.rate_limit_seconds)

    cand = pl.read_csv(args.candidates).head(args.top)
    n = cand.height
    rng = random.Random(args.seed)

    # Shuffle dates: build a derangement so no candidate keeps its own date.
    indices = list(range(n))
    for _ in range(100):
        shuffled = indices.copy()
        rng.shuffle(shuffled)
        if all(shuffled[i] != i for i in range(n)):
            break
    else:
        logger.error("Could not produce a derangement")
        return 1

    dates = cand["date_utc"].to_list()
    fake_dates = [dates[shuffled[i]] for i in range(n)]

    rows: list[dict] = []
    for i, row in enumerate(cand.iter_rows(named=True)):
        perturber = int(row["perturber_number"])
        target_no = row.get("target_number")
        target_des = row["target_designation"]
        real_date = row["date_utc"]
        fake_date = fake_dates[i]

        if target_no is None or (isinstance(target_no, float) and math.isnan(target_no)):
            continue
        target_no = int(target_no)

        logger.info(
            "[%d/%d] (%d) + (%d) %s — real %s, FAKE %s",
            i + 1, n, perturber, target_no, target_des, real_date, fake_date,
        )

        # Use the fake date to set the window
        fake_days = float(Time(fake_date, scale="utc").tcb.jd) - _J2010_TCB_JD
        d_min = fake_days - _HALF_WINDOW_DAYS
        d_max = fake_days + _HALF_WINDOW_DAYS

        try:
            obs = fetch_gaia(archive_url, target_no, d_min, d_max)
        except Exception as exc:  # noqa: BLE001
            logger.warning("  Gaia query failed: %s", exc)
            continue
        if obs.height == 0:
            logger.info("  no Gaia transits at fake date")
            rows.append({
                "perturber_number": perturber, "target_number": target_no,
                "target_designation": target_des, "real_date": real_date,
                "fake_date": fake_date, "n_transits": 0,
                "t_fake_ra": float("nan"), "t_fake_dec": float("nan"),
                "detected_fake": False,
            })
            continue

        epochs_days = obs["epoch"].to_numpy()
        jd_tcb = epochs_days + _J2010_TCB_JD
        jd_tdb = Time(jd_tcb, format="jd", scale="tcb").tdb.jd.astype(float)

        try:
            ra_pred, dec_pred = horizons_apparent_radec(target_no, jd_tdb, rate_limit)
        except Exception as exc:  # noqa: BLE001
            logger.warning("  Horizons query failed: %s", exc)
            continue

        ra_obs = obs["ra"].to_numpy().astype(float)
        dec_obs = obs["dec"].to_numpy().astype(float)
        deg = np.pi / 180.0
        dra_mas = ((ra_obs - ra_pred + 540.0) % 360.0 - 180.0) * np.cos(dec_pred * deg) * 3_600_000.0
        ddec_mas = (dec_obs - dec_pred) * 3_600_000.0

        fake_jd_tdb = float(Time(fake_date, scale="utc").tdb.jd)
        days_from_fake = jd_tdb - fake_jd_tdb
        before = days_from_fake < -_BLACKOUT_DAYS
        after = days_from_fake > _BLACKOUT_DAYS

        t_ra = _t_welch(dra_mas, before, after)
        t_dec = _t_welch(ddec_mas, before, after)
        detected_fake = (
            (math.isfinite(t_ra) and abs(t_ra) >= 3.0)
            or (math.isfinite(t_dec) and abs(t_dec) >= 3.0)
        )

        rows.append({
            "perturber_number": perturber,
            "target_number": target_no,
            "target_designation": target_des,
            "real_date": real_date,
            "fake_date": fake_date,
            "n_transits": obs.height,
            "t_fake_ra": t_ra,
            "t_fake_dec": t_dec,
            "detected_fake": detected_fake,
        })
        logger.info("  t_fake_RA = %+.2fσ  t_fake_Dec = %+.2fσ  → %s",
                    t_ra, t_dec, "DETECTED" if detected_fake else "no")

    out = pl.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.write_csv(args.output)
    logger.info("Wrote %d rows to %s", out.height, args.output)

    # Compare to real detections
    real = pl.read_csv("data/output/deflection_detections.csv").head(args.top)
    real_detected = int((real["detection"] == "yes").sum())
    fake_detected = int(out["detected_fake"].sum())
    logger.info("")
    logger.info("NULL TEST RESULT:")
    logger.info("  Real encounter dates:  %d / %d detected (≥3σ)", real_detected, args.top)
    logger.info("  Fake encounter dates:  %d / %d detected (≥3σ)", fake_detected, out.height)
    if fake_detected < real_detected:
        logger.info("  → Real detections exceed null — encounter coincidence is real ✓")
    else:
        logger.warning("  → No excess over null — detections may be systematic")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
