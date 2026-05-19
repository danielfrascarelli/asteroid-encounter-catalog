"""Run the deflection-detection demo over the top-N publishable candidates.

Generalises ``demo_ate_clean.py`` to any subset of the candidate catalog.
For each (perturber, target, date) tuple in
``data/output/publishable_mass_candidates.csv``:

  1. Pull Gaia DR3 transits of the target in ±180 days of the encounter.
  2. Query JPL Horizons for the apparent RA/Dec of the target as observed
     from Gaia (location code 500@-139479).  Horizons handles light-time +
     aberration internally.
  3. Compute residuals (Gaia_obs − Horizons_pred) in mas.
  4. Split into pre- vs post-encounter (with ±7-day blackout).
  5. Welch t-statistic on the shift of the mean residual across the
     encounter date, per axis.

Horizons does NOT include the perturber asteroid in its N-body model
(except for the big-4 Ceres/Pallas/Vesta/Hygiea), so a non-zero residual
shift coincident with the encounter date is the perturbation signature.

Output: ``data/output/deflection_detections.csv``

Usage
-----
    docker compose run --rm pipeline python -m scripts.detect_deflections
    docker compose run --rm pipeline python -m scripts.detect_deflections --top 5
"""

from __future__ import annotations

import argparse
import logging
import math
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
_DEFAULT_INPUT = Path("data/output/publishable_mass_candidates.csv")
_DEFAULT_OUTPUT = Path("data/output/deflection_detections.csv")
_PER_CANDIDATE_RESIDUALS_DIR = Path("data/output/deflection_residuals")


def fetch_gaia(archive_url: str, target: int, d_min: float, d_max: float) -> pl.DataFrame:
    adql = (
        "SELECT number_mp, epoch, ra, dec, g_mag "
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


def analyze_one_candidate(
    archive_url: str,
    rate_limit_s: float,
    perturber_number: int,
    perturber_name: str,
    target_number: int,
    target_designation: str,
    date_utc: str,
    expected_muas: float,
    residual_dir: Path,
) -> dict:
    """Run the deflection demo on one (perturber, target, date) tuple."""
    enc_jd_tdb = float(Time(date_utc, scale="utc").tdb.jd)
    enc_days = float(Time(date_utc, scale="utc").tcb.jd) - _J2010_TCB_JD

    try:
        obs = fetch_gaia(
            archive_url,
            target_number,
            enc_days - _HALF_WINDOW_DAYS,
            enc_days + _HALF_WINDOW_DAYS,
        )
    except Exception as exc:  # noqa: BLE001
        logger.warning("  Gaia query failed: %s", exc)
        return _empty_row(
            perturber_number, perturber_name, target_number, target_designation,
            date_utc, expected_muas, note=f"gaia query failed: {exc}",
        )

    if obs.height == 0:
        logger.warning("  no Gaia transits in window")
        return _empty_row(
            perturber_number, perturber_name, target_number, target_designation,
            date_utc, expected_muas, note="no Gaia transits",
        )

    epochs_days = obs["epoch"].to_numpy()
    jd_tcb = epochs_days + _J2010_TCB_JD
    jd_tdb = Time(jd_tcb, format="jd", scale="tcb").tdb.jd.astype(float)

    try:
        ra_pred, dec_pred = horizons_apparent_radec(target_number, jd_tdb, rate_limit_s)
    except Exception as exc:  # noqa: BLE001
        logger.warning("  Horizons query failed: %s", exc)
        return _empty_row(
            perturber_number, perturber_name, target_number, target_designation,
            date_utc, expected_muas, note=f"horizons failed: {exc}",
        )

    ra_obs = obs["ra"].to_numpy().astype(float)
    dec_obs = obs["dec"].to_numpy().astype(float)

    deg = np.pi / 180.0
    dra_mas = ((ra_obs - ra_pred + 540.0) % 360.0 - 180.0) * np.cos(dec_pred * deg) * 3_600_000.0
    ddec_mas = (dec_obs - dec_pred) * 3_600_000.0
    sep_mas = np.sqrt(dra_mas**2 + ddec_mas**2)
    days_from_enc = jd_tdb - enc_jd_tdb

    # Save per-candidate residuals
    residual_dir.mkdir(parents=True, exist_ok=True)
    safe = f"{perturber_number:06d}_{target_number:06d}"
    pl.DataFrame(
        {
            "jd_tdb": jd_tdb,
            "days_from_encounter": days_from_enc,
            "ra_obs_deg": ra_obs,
            "dec_obs_deg": dec_obs,
            "ra_pred_deg": ra_pred,
            "dec_pred_deg": dec_pred,
            "dra_mas": dra_mas,
            "ddec_mas": ddec_mas,
            "sep_mas": sep_mas,
            "g_mag": obs["g_mag"].to_numpy().astype(float),
        }
    ).write_csv(residual_dir / f"{safe}.csv")

    before_mask = days_from_enc < -_BLACKOUT_DAYS
    after_mask = days_from_enc > _BLACKOUT_DAYS
    n_before = int(before_mask.sum())
    n_after = int(after_mask.sum())

    def _axis_stats(mask: np.ndarray, vals: np.ndarray) -> tuple[float, float]:
        if mask.sum() == 0:
            return float("nan"), float("nan")
        return float(np.mean(vals[mask])), float(np.std(vals[mask]))

    mu_dra_b, sd_dra_b = _axis_stats(before_mask, dra_mas)
    mu_dra_a, sd_dra_a = _axis_stats(after_mask, dra_mas)
    mu_ddec_b, sd_ddec_b = _axis_stats(before_mask, ddec_mas)
    mu_ddec_a, sd_ddec_a = _axis_stats(after_mask, ddec_mas)

    def _t_welch(mb, ma, sb, sa, nb, na) -> float:
        if nb < 2 or na < 2:
            return float("nan")
        se = math.sqrt((sb * sb) / nb + (sa * sa) / na)
        return (ma - mb) / se if se > 0 else float("nan")

    t_dra = _t_welch(mu_dra_b, mu_dra_a, sd_dra_b, sd_dra_a, n_before, n_after)
    t_ddec = _t_welch(mu_ddec_b, mu_ddec_a, sd_ddec_b, sd_ddec_a, n_before, n_after)

    return {
        "perturber_number": perturber_number,
        "perturber_name": perturber_name,
        "target_number": target_number,
        "target_designation": target_designation,
        "date_utc": date_utc,
        "expected_muas": expected_muas,
        "n_obs_before": n_before,
        "n_obs_after": n_after,
        "median_sep_mas": float(np.median(sep_mas)),
        "mu_dra_before_mas": mu_dra_b,
        "mu_dra_after_mas": mu_dra_a,
        "shift_dra_mas": mu_dra_a - mu_dra_b,
        "sd_dra_before": sd_dra_b,
        "sd_dra_after": sd_dra_a,
        "t_dra": t_dra,
        "mu_ddec_before_mas": mu_ddec_b,
        "mu_ddec_after_mas": mu_ddec_a,
        "shift_ddec_mas": mu_ddec_a - mu_ddec_b,
        "sd_ddec_before": sd_ddec_b,
        "sd_ddec_after": sd_ddec_a,
        "t_ddec": t_ddec,
        "detection": (
            "yes" if (math.isfinite(t_dra) and abs(t_dra) >= 3.0)
            or (math.isfinite(t_ddec) and abs(t_ddec) >= 3.0)
            else "no"
        ),
        "note": "",
    }


def _empty_row(perturber_number, perturber_name, target_number, target_designation,
               date_utc, expected_muas, note: str) -> dict:
    return {
        "perturber_number": perturber_number,
        "perturber_name": perturber_name,
        "target_number": target_number,
        "target_designation": target_designation,
        "date_utc": date_utc,
        "expected_muas": expected_muas,
        "n_obs_before": 0,
        "n_obs_after": 0,
        "median_sep_mas": float("nan"),
        "mu_dra_before_mas": float("nan"),
        "mu_dra_after_mas": float("nan"),
        "shift_dra_mas": float("nan"),
        "sd_dra_before": float("nan"),
        "sd_dra_after": float("nan"),
        "t_dra": float("nan"),
        "mu_ddec_before_mas": float("nan"),
        "mu_ddec_after_mas": float("nan"),
        "shift_ddec_mas": float("nan"),
        "sd_ddec_before": float("nan"),
        "sd_ddec_after": float("nan"),
        "t_ddec": float("nan"),
        "detection": "no",
        "note": note,
    }


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--config", default="config.yaml")
    p.add_argument(
        "--input",
        type=Path,
        default=_DEFAULT_INPUT,
        help="Candidates CSV (must have target_number, date_utc, deflection_muas).",
    )
    p.add_argument(
        "--output",
        type=Path,
        default=_DEFAULT_OUTPUT,
    )
    p.add_argument("--top", type=int, default=10, help="Process only the top-N candidates.")
    p.add_argument(
        "--residuals-dir",
        type=Path,
        default=_PER_CANDIDATE_RESIDUALS_DIR,
    )
    args = p.parse_args()

    cfg = load_config(args.config)
    archive_url = cfg.sources.gaia_sso.archive_url
    rate_limit = float(cfg.sources.jpl_horizons.rate_limit_seconds)

    if not args.input.exists():
        logger.error("Candidates CSV not found: %s", args.input)
        return 1

    df = pl.read_csv(args.input)
    df = df.head(args.top)
    logger.info("Processing top %d candidates from %s", df.height, args.input)

    rows: list[dict] = []
    for i, row in enumerate(df.iter_rows(named=True), start=1):
        perturber = int(row["perturber_number"])
        perturber_name = row["perturber_name"]
        target_no = row.get("target_number")
        target_des = row["target_designation"]
        date_utc = row["date_utc"]
        expected_muas = float(row.get("deflection_muas", float("nan")))

        if target_no is None or (isinstance(target_no, float) and math.isnan(target_no)):
            logger.warning(
                "[%d/%d] (%d) %s + %s: target has no MPC number — skipping",
                i, df.height, perturber, perturber_name, target_des,
            )
            rows.append(_empty_row(
                perturber, perturber_name, None, target_des,
                date_utc, expected_muas, note="no target_number",
            ))
            continue

        target_no = int(target_no)
        logger.info(
            "[%d/%d] (%d) %s + (%d) %s on %s  (expected δ = %.0f μas)…",
            i, df.height, perturber, perturber_name, target_no, target_des,
            date_utc, expected_muas,
        )

        result = analyze_one_candidate(
            archive_url=archive_url,
            rate_limit_s=rate_limit,
            perturber_number=perturber,
            perturber_name=perturber_name,
            target_number=target_no,
            target_designation=target_des,
            date_utc=date_utc,
            expected_muas=expected_muas,
            residual_dir=args.residuals_dir,
        )
        rows.append(result)

        if result["detection"] == "yes":
            logger.info(
                "  → DETECTION  shift_RA=%+.1f mas (t=%+.2f), shift_Dec=%+.1f mas (t=%+.2f), σ_RA(before)=%.0f mas",
                result["shift_dra_mas"], result["t_dra"],
                result["shift_ddec_mas"], result["t_ddec"],
                result["sd_dra_before"],
            )
        else:
            logger.info(
                "  → no detection  shift_RA=%+.1f mas (t=%+.2f), shift_Dec=%+.1f mas (t=%+.2f)",
                result["shift_dra_mas"], result["t_dra"],
                result["shift_ddec_mas"], result["t_ddec"],
            )

    out = pl.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.write_csv(args.output)
    logger.info("Wrote %d rows to %s", out.height, args.output)

    # Summary
    n_detected = int((out["detection"] == "yes").sum())
    logger.info("Summary: %d / %d candidates with |t| ≥ 3σ shift on either axis",
                n_detected, out.height)

    finite = out.filter(pl.col("t_dra").is_finite())
    if finite.height > 0:
        logger.info("")
        logger.info("Per-candidate summary (sorted by |t_dra|):")
        for r in finite.sort(pl.col("t_dra").abs(), descending=True).iter_rows(named=True):
            tag = "✓" if r["detection"] == "yes" else " "
            logger.info(
                " %s (%d) %-22s + (%s) %-15s  shift=%+7.1f mas  t=%+5.2fσ  (expected δ=%.0f μas)",
                tag,
                r["perturber_number"],
                (r["perturber_name"] or "")[:22],
                str(r["target_number"]) if r["target_number"] is not None else "—",
                (r["target_designation"] or "")[:15],
                r["shift_dra_mas"],
                r["t_dra"],
                r["expected_muas"],
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
