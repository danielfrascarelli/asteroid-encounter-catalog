"""Check Gaia DR3 SSO observations for mass-determination candidate targets.

For each Category B candidate in ``data/output/mass_candidates.csv`` (or the
filtered set in ``relevant_novel_encounters.csv``), queries the Gaia archive
TAP service for transits of the *target* asteroid in a window bracketing the
encounter date.

Mass determination from a close encounter requires astrometric measurements
of the deflected body BOTH before and after the encounter:

    pre-encounter (constrains baseline orbit)
        ↓
    ENCOUNTER (the perturbation we want to measure)
        ↓
    post-encounter (detects the deflection)

The default window is ±180 days around the encounter, excluding a ±7-day
"blackout" centred on the encounter itself (during which the encounter
geometry is still resolving). A candidate is flagged ``viable_obs`` when both
``n_obs_before >= 3`` and ``n_obs_after >= 3`` — three transits on each side
is the minimum to constrain a 2D astrometric fit.

This script makes ONE batched TAP query covering every candidate target, then
slices the result per-target locally. Fast (seconds, not hours) because it
fetches only the rows that matter.

Output: ``data/output/gaia_observations_check.csv``

Usage
-----
    docker compose run --rm pipeline python -m scripts.check_gaia_observations
    docker compose run --rm pipeline python -m scripts.check_gaia_observations \\
        --input data/output/mass_candidates.csv \\
        --half-window-days 180 --blackout-days 7
"""

from __future__ import annotations

import argparse
import logging
import math
from pathlib import Path

import polars as pl
from astropy.time import Time
from astroquery.utils.tap.core import TapPlus

from src.utils.config import load_config

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)


_DEFAULT_INPUT = Path("data/output/mass_candidates.csv")
_DEFAULT_OUTPUT = Path("data/output/gaia_observations_check.csv")
_MIN_OBS_EACH_SIDE = 3  # classical mass-det minimum (3 transits before + 3 after)

# Gaia DR3 SSO `epoch` column is "days since J2010.0 TCB" (Julian Date - 2455197.5),
# despite the column metadata claiming it is a JD.  Empirically: Ceres observations
# in 2015-12 have epoch ≈ 2167, which only makes sense as days since J2010.
_J2010_TCB_JD = 2455197.5


# ---------------------------------------------------------------------------
# Time helpers
# ---------------------------------------------------------------------------


def _date_utc_to_days_since_j2010(date_utc: str) -> float:
    """Convert ISO UTC string to *days since J2010.0 TCB*.

    Gaia DR3 stores ``sso_observation.epoch`` as ``JD_TCB − 2455197.5``.
    The TCB↔UTC offset (~50 s over the DR3 window) is negligible at the
    ±180-day timescales used here.
    """
    return float(Time(date_utc, format="iso", scale="utc").tcb.jd) - _J2010_TCB_JD


def _days_since_j2010_to_iso(days: float) -> str:
    """Convert days-since-J2010-TCB back to ISO UTC calendar date."""
    return Time(days + _J2010_TCB_JD, format="jd", scale="tcb").utc.iso[:10]


# ---------------------------------------------------------------------------
# Gaia TAP query
# ---------------------------------------------------------------------------


def fetch_gaia_observations(
    archive_url: str,
    target_numbers: list[int],
    days_min: float,
    days_max: float,
) -> pl.DataFrame:
    """Fetch transits of every *target_numbers* asteroid in [days_min, days_max].

    Time bounds are in *days since J2010.0 TCB* (the ``epoch`` column's units).
    Returns columns ``number_mp``, ``epoch``, ``epoch_utc``, ``ra``, ``dec``, ``g_mag``.
    """
    schema = {
        "number_mp": pl.Int64,
        "epoch": pl.Float64,
        "epoch_utc": pl.Float64,
        "ra": pl.Float64,
        "dec": pl.Float64,
        "g_mag": pl.Float64,
    }
    if not target_numbers:
        return pl.DataFrame(schema=schema)

    in_list = ", ".join(str(int(n)) for n in sorted(set(target_numbers)))
    adql = (
        "SELECT number_mp, epoch, epoch_utc, ra, dec, g_mag "
        "FROM gaiadr3.sso_observation "
        f"WHERE number_mp IN ({in_list}) "
        f"AND epoch BETWEEN {days_min:.6f} AND {days_max:.6f}"
    )
    logger.info(
        "Querying Gaia TAP for %d targets, epoch (d since J2010 TCB) ∈ [%.1f, %.1f]…",
        len(set(target_numbers)),
        days_min,
        days_max,
    )
    tap = TapPlus(url=archive_url)
    # Use the async endpoint so we are not capped at the 2000-row sync limit.
    job = tap.launch_job_async(adql)
    table = job.get_results()
    if len(table) == 0:
        logger.warning("Gaia TAP returned 0 rows")
        return pl.DataFrame(schema=schema)
    df = pl.from_pandas(table.to_pandas())
    df = df.rename({c: c.lower() for c in df.columns if c != c.lower()})
    logger.info("Received %d transit rows from Gaia", df.height)
    return df


# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------


def _is_missing_number(value: object) -> bool:
    if value is None:
        return True
    if isinstance(value, float) and math.isnan(value):
        return True
    return False


def analyze_candidates(
    candidates: pl.DataFrame,
    obs: pl.DataFrame,
    half_window_days: float,
    blackout_days: float,
) -> pl.DataFrame:
    """For each candidate, count Gaia transits before/after the encounter.

    A transit is "before" if its epoch is in [t_enc - half_window, t_enc - blackout],
    "after" if in [t_enc + blackout, t_enc + half_window].
    """
    obs_by_target: dict[int, pl.DataFrame] = {}
    if obs.height > 0:
        for k, g in obs.partition_by("number_mp", as_dict=True).items():
            # polars >= 0.20 returns tuple keys; older versions returned scalars
            number_mp = k[0] if isinstance(k, tuple) else k
            obs_by_target[int(number_mp)] = g

    rows: list[dict] = []
    for row in candidates.iter_rows(named=True):
        rank = row.get("rank")
        perturber = int(row["perturber_number"])
        perturber_name = row["perturber_name"]
        target_no_raw = row.get("target_number")
        target_designation = row["target_designation"]
        date_utc = row["date_utc"]
        dist_au = float(row["dist_au"])
        deflection_muas = float(row.get("deflection_muas", float("nan")))

        if _is_missing_number(target_no_raw):
            rows.append(
                {
                    "rank": rank,
                    "perturber_number": perturber,
                    "perturber_name": perturber_name,
                    "target_number": None,
                    "target_designation": target_designation,
                    "date_utc": date_utc,
                    "dist_au": dist_au,
                    "deflection_muas": deflection_muas,
                    "n_obs_before": 0,
                    "n_obs_after": 0,
                    "n_obs_total": 0,
                    "first_obs_date": None,
                    "last_obs_date": None,
                    "median_g_mag": None,
                    "viable_obs": False,
                    "note": "target has no MPC number",
                }
            )
            continue

        target_no = int(target_no_raw)
        target_obs = obs_by_target.get(target_no, pl.DataFrame(schema=obs.schema))

        d_center = _date_utc_to_days_since_j2010(date_utc)
        d_lo = d_center - half_window_days
        d_hi = d_center + half_window_days
        d_blackout_lo = d_center - blackout_days
        d_blackout_hi = d_center + blackout_days

        before = target_obs.filter((pl.col("epoch") >= d_lo) & (pl.col("epoch") < d_blackout_lo))
        after = target_obs.filter((pl.col("epoch") > d_blackout_hi) & (pl.col("epoch") <= d_hi))
        in_window = target_obs.filter((pl.col("epoch") >= d_lo) & (pl.col("epoch") <= d_hi))

        n_before = before.height
        n_after = after.height
        n_total = in_window.height

        first_date = (
            _days_since_j2010_to_iso(float(in_window["epoch"].min())) if n_total > 0 else None
        )
        last_date = (
            _days_since_j2010_to_iso(float(in_window["epoch"].max())) if n_total > 0 else None
        )
        median_g = float(in_window["g_mag"].median()) if n_total > 0 else None

        viable_obs = n_before >= _MIN_OBS_EACH_SIDE and n_after >= _MIN_OBS_EACH_SIDE

        rows.append(
            {
                "rank": rank,
                "perturber_number": perturber,
                "perturber_name": perturber_name,
                "target_number": target_no,
                "target_designation": target_designation,
                "date_utc": date_utc,
                "dist_au": dist_au,
                "deflection_muas": deflection_muas,
                "n_obs_before": n_before,
                "n_obs_after": n_after,
                "n_obs_total": n_total,
                "first_obs_date": first_date,
                "last_obs_date": last_date,
                "median_g_mag": median_g,
                "viable_obs": viable_obs,
                "note": (
                    ""
                    if viable_obs
                    else (
                        "no Gaia observations"
                        if n_total == 0
                        else "insufficient bracketing transits"
                    )
                ),
            }
        )

    return pl.DataFrame(rows)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--config", default="config.yaml", help="Path to config YAML.")
    p.add_argument(
        "--input",
        type=Path,
        default=_DEFAULT_INPUT,
        help="Path to candidates CSV (must have target_number, date_utc columns).",
    )
    p.add_argument(
        "--output",
        type=Path,
        default=_DEFAULT_OUTPUT,
        help="Output CSV path.",
    )
    p.add_argument(
        "--half-window-days",
        type=float,
        default=180.0,
        help="Half-width (days) of the observation window around the encounter.",
    )
    p.add_argument(
        "--blackout-days",
        type=float,
        default=7.0,
        help="Half-width (days) of the blackout interval excluded around the encounter.",
    )
    args = p.parse_args()

    if not args.input.exists():
        logger.error("Input CSV not found: %s", args.input)
        return 1

    cfg = load_config(args.config)
    archive_url = cfg.sources.gaia_sso.archive_url

    candidates = pl.read_csv(args.input)
    logger.info("Loaded %d candidates from %s", candidates.height, args.input)

    # Collect unique target numbers (skip missing)
    target_numbers: list[int] = []
    d_min_list: list[float] = []
    d_max_list: list[float] = []
    for row in candidates.iter_rows(named=True):
        tn = row.get("target_number")
        if _is_missing_number(tn):
            continue
        target_numbers.append(int(tn))
        d_center = _date_utc_to_days_since_j2010(row["date_utc"])
        d_min_list.append(d_center - args.half_window_days)
        d_max_list.append(d_center + args.half_window_days)

    if not target_numbers:
        logger.warning("No candidate has a usable target_number; nothing to query.")
        return 0

    # Single batched query covering the union of all per-target windows.
    d_min = min(d_min_list)
    d_max = max(d_max_list)

    obs = fetch_gaia_observations(archive_url, target_numbers, d_min, d_max)

    if obs.height > 0:
        unique_targets_returned = obs["number_mp"].n_unique()
        logger.info(
            "Returned data for %d of %d target asteroids",
            unique_targets_returned,
            len(set(target_numbers)),
        )

    result = analyze_candidates(
        candidates,
        obs,
        half_window_days=args.half_window_days,
        blackout_days=args.blackout_days,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.write_csv(args.output)
    logger.info("Wrote %d rows to %s", result.height, args.output)

    # Also write the publishable subset: viable_obs == True, sorted by deflection.
    publishable_path = args.output.parent / "mass_followup_candidates.csv"
    publishable = result.filter(pl.col("viable_obs")).sort("deflection_muas", descending=True)
    publishable.write_csv(publishable_path)
    logger.info(
        "Wrote %d publishable candidates (viable_obs == True) → %s",
        publishable.height,
        publishable_path,
    )

    # Summary
    n_viable = int(result["viable_obs"].sum())
    n_with_any_obs = int(result.filter(pl.col("n_obs_total") > 0).height)
    logger.info(
        "Summary: %d / %d candidates have ≥%d transits on each side  (viable_obs)",
        n_viable,
        result.height,
        _MIN_OBS_EACH_SIDE,
    )
    logger.info(
        "         %d / %d candidates have any Gaia transit in the window",
        n_with_any_obs,
        result.height,
    )

    # Pretty top-10 table
    logger.info("Top 10 candidates by deflection_muas, with Gaia coverage:")
    header = (
        f"{'rk':>3}  {'perturber':<22}  {'target':<18}  {'date':<10}  "
        f"{'δ_μas':>8}  {'before':>6}  {'after':>6}  {'mag':>5}  {'viable':>6}"
    )
    logger.info(header)
    logger.info("-" * len(header))
    top = result.sort("deflection_muas", descending=True).head(10)
    for r in top.iter_rows(named=True):
        pert_label = f"({r['perturber_number']}) {r['perturber_name']}"
        mag_str = f"{r['median_g_mag']:.1f}" if r["median_g_mag"] is not None else "  —"
        logger.info(
            "%3s  %-22s  %-18s  %-10s  %8.0f  %6d  %6d  %5s  %6s",
            r["rank"] if r["rank"] is not None else "—",
            pert_label[:22],
            (r["target_designation"] or "")[:18],
            r["date_utc"][:10],
            r["deflection_muas"] if r["deflection_muas"] is not None else float("nan"),
            r["n_obs_before"],
            r["n_obs_after"],
            mag_str,
            "yes" if r["viable_obs"] else "no",
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
