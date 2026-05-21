"""Audit Gaia DR3 transit coverage for every novel-encounter candidate.

For each encounter in data/output/relevant_novel_encounters.csv this script:

1. Sends ONE batch ADQL query to the Gaia TAP service to fetch transit epochs
   for all 379 target asteroids at once.
2. For each encounter computes:
   - n_pre_transits   : Gaia transits in [enc_date − 365 d, enc_date − 30 d]
   - n_post_transits  : Gaia transits in [enc_date + 5 d,   enc_date + 180 d]
   - nearest_pre_days : days from encounter to closest pre-encounter transit
   - nearest_post_days: days from encounter to closest post-encounter transit
   - n_window_transits: all transits within ±180 d of encounter
   - viable_coverage  : n_pre ≥ 5 AND n_post ≥ 3
3. Writes data/output/gaia_coverage_audit.csv ranked by deflection_score.

Usage
-----
    docker compose run --rm pipeline python -m scripts.audit_gaia_coverage

The script fetches data from the Gaia archive (requires network access).
Typical run time: 3–10 min for the full 379-target batch query.
"""

from __future__ import annotations

import logging
import math
from pathlib import Path

import numpy as np
import polars as pl
import yaml
from astroquery.utils.tap.core import TapPlus

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Paths & constants
# ---------------------------------------------------------------------------

_ENCOUNTERS_PATH = Path("data/output/relevant_novel_encounters.csv")
_OUTPUT_PATH = Path("data/output/gaia_coverage_audit.csv")
_CONFIG_PATH = Path("config.yaml")

# Gaia time reference: epoch column = BJD_TCB − J2010.0
_J2010_TCB_JD: float = 2455197.5
_GAIA_START_JD_TCB: float = 2456863.5  # 2014-07-25
_GAIA_END_JD_TCB: float = 2457910.5  # 2017-05-28

# Coverage windows (days relative to encounter)
_PRE_MIN_DAYS: float = 30.0  # at least 30 d before (for clean pre-arc)
_PRE_MAX_DAYS: float = 365.0  # up to 1 year before
_POST_MIN_DAYS: float = 5.0  # skip first 5 d (too close to encounter numerically)
_POST_MAX_DAYS: float = 180.0  # up to 6 months after
_WIN_HALF_DAYS: float = 180.0  # ±window for n_window_transits

# Minimum transit counts to call coverage "viable" for LOO mass fit
_MIN_PRE_TRANSITS: int = 5
_MIN_POST_TRANSITS: int = 3

# Gaia TAP chunk size (IN clauses have practical limits)
_CHUNK_SIZE: int = 200


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _date_to_gaia_epoch(date_utc: str) -> float:
    """Convert ISO date string to Gaia epoch units (BJD_TCB − J2010.0)."""
    from astropy.time import Time

    jd_utc = float(Time(date_utc, format="iso", scale="utc").jd)
    return jd_utc - _J2010_TCB_JD


def _batch_query(
    archive_url: str,
    numbers: list[int],
) -> pl.DataFrame:
    """Fetch epoch + number_mp for a batch of asteroid numbers from Gaia TAP."""
    d_min = _GAIA_START_JD_TCB - _J2010_TCB_JD
    d_max = _GAIA_END_JD_TCB - _J2010_TCB_JD
    nums_str = ",".join(str(n) for n in numbers)
    adql = (
        "SELECT number_mp, epoch "
        "FROM gaiadr3.sso_observation "
        f"WHERE number_mp IN ({nums_str}) "
        f"AND epoch BETWEEN {d_min:.6f} AND {d_max:.6f} "
        "ORDER BY number_mp, epoch"
    )
    logger.info("Querying Gaia TAP for %d asteroids…", len(numbers))
    tap = TapPlus(url=archive_url)
    job = tap.launch_job_async(adql)
    raw = job.get_results().to_pandas()
    df = pl.from_pandas(raw)
    df = df.rename({c: c.lower() for c in df.columns if c != c.lower()})
    logger.info("  → %d rows returned", len(df))
    return df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    cfg = yaml.safe_load(_CONFIG_PATH.read_text())
    archive_url: str = cfg["sources"]["gaia_sso"]["archive_url"]

    enc = pl.read_csv(_ENCOUNTERS_PATH)
    logger.info("Loaded %d encounters from %s", len(enc), _ENCOUNTERS_PATH)

    target_numbers: list[int] = enc["number_2"].to_list()

    # ── Batch-fetch Gaia epochs for all targets ────────────────────────────
    # Split into chunks to stay within TAP IN-clause limits
    all_frames: list[pl.DataFrame] = []
    n_chunks = math.ceil(len(target_numbers) / _CHUNK_SIZE)
    for i in range(n_chunks):
        chunk = target_numbers[i * _CHUNK_SIZE : (i + 1) * _CHUNK_SIZE]
        logger.info("Chunk %d/%d: %d numbers", i + 1, n_chunks, len(chunk))
        df_chunk = _batch_query(archive_url, chunk)
        all_frames.append(df_chunk)

    if all_frames:
        gaia_all = pl.concat(all_frames)
    else:
        gaia_all = pl.DataFrame(
            {"number_mp": pl.Series([], dtype=pl.Int64), "epoch": pl.Series([], dtype=pl.Float64)}
        )

    logger.info(
        "Total Gaia transits fetched: %d across %d unique target asteroids",
        len(gaia_all),
        gaia_all["number_mp"].n_unique() if len(gaia_all) else 0,
    )

    # Build per-asteroid lookup: number_mp → sorted epoch array
    gaia_by_number: dict[int, np.ndarray] = {}
    for num, grp in gaia_all.group_by("number_mp"):
        gaia_by_number[int(num[0])] = grp["epoch"].sort().to_numpy()

    # ── Analyse coverage per encounter ────────────────────────────────────
    results: list[dict] = []

    for row in enc.iter_rows(named=True):
        perturber_num = int(row["number_1"])
        perturber_name = str(row["designation_1"])
        target_num = int(row["number_2"])
        target_desg = str(row["designation_2"])
        date_utc = str(row["date_utc"])
        dist_au = float(row["dist_au"])
        deflection_score = float(row.get("deflection_score", 0.0))

        enc_epoch = _date_to_gaia_epoch(date_utc)
        epochs = gaia_by_number.get(target_num)

        if epochs is None or len(epochs) == 0:
            results.append(
                {
                    "perturber_number": perturber_num,
                    "perturber_name": perturber_name,
                    "target_number": target_num,
                    "target_designation": target_desg,
                    "date_utc": date_utc,
                    "dist_au": dist_au,
                    "deflection_score": deflection_score,
                    "has_gaia_data": False,
                    "n_total_transits": 0,
                    "n_pre_transits": 0,
                    "n_post_transits": 0,
                    "n_window_transits": 0,
                    "nearest_pre_days": None,
                    "nearest_post_days": None,
                    "viable_coverage": False,
                    "note": "target has no Gaia DR3 transits",
                }
            )
            continue

        delta = epochs - enc_epoch  # positive = after encounter

        pre_mask = (delta >= -_PRE_MAX_DAYS) & (delta <= -_PRE_MIN_DAYS)
        post_mask = (delta >= _POST_MIN_DAYS) & (delta <= _POST_MAX_DAYS)
        win_mask = np.abs(delta) <= _WIN_HALF_DAYS

        n_pre = int(pre_mask.sum())
        n_post = int(post_mask.sum())
        n_win = int(win_mask.sum())

        nearest_pre: float | None = None
        nearest_post: float | None = None

        if n_pre > 0:
            # Most recent pre-encounter transit (closest to encounter)
            nearest_pre = float(np.max(delta[pre_mask]))  # least-negative = closest
            nearest_pre = abs(nearest_pre)

        if n_post > 0:
            nearest_post = float(np.min(delta[post_mask]))  # first post-encounter transit

        viable = (n_pre >= _MIN_PRE_TRANSITS) and (n_post >= _MIN_POST_TRANSITS)

        note_parts: list[str] = []
        if n_pre < _MIN_PRE_TRANSITS:
            note_parts.append(f"pre: {n_pre}<{_MIN_PRE_TRANSITS}")
        if n_post < _MIN_POST_TRANSITS:
            note_parts.append(f"post: {n_post}<{_MIN_POST_TRANSITS}")
        note = "; ".join(note_parts) if note_parts else "ok"

        results.append(
            {
                "perturber_number": perturber_num,
                "perturber_name": perturber_name,
                "target_number": target_num,
                "target_designation": target_desg,
                "date_utc": date_utc,
                "dist_au": dist_au,
                "deflection_score": deflection_score,
                "has_gaia_data": True,
                "n_total_transits": int(len(epochs)),
                "n_pre_transits": n_pre,
                "n_post_transits": n_post,
                "n_window_transits": n_win,
                "nearest_pre_days": nearest_pre,
                "nearest_post_days": nearest_post,
                "viable_coverage": viable,
                "note": note,
            }
        )

    result_df = (
        pl.DataFrame(results)
        .sort("deflection_score", descending=True)
        .sort("viable_coverage", descending=True)
    )

    _OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    result_df.write_csv(_OUTPUT_PATH)

    n_viable = int(result_df["viable_coverage"].sum())
    n_with_data = int(result_df["has_gaia_data"].sum())
    logger.info("─" * 60)
    logger.info("Wrote %d rows → %s", len(result_df), _OUTPUT_PATH)
    logger.info("Targets with Gaia data: %d / %d", n_with_data, len(result_df))
    logger.info("Viable for LOO mass fit: %d / %d", n_viable, len(result_df))

    if n_viable > 0:
        viable_df = result_df.filter(pl.col("viable_coverage"))
        logger.info("\nTop viable candidates:")
        for r in viable_df.head(10).iter_rows(named=True):
            logger.info(
                "  (%d) %-14s → %-20s  %s  %.5f AU  pre=%d  post=%d  gap=%.0f d",
                r["perturber_number"],
                r["perturber_name"],
                r["target_designation"],
                r["date_utc"],
                r["dist_au"],
                r["n_pre_transits"],
                r["n_post_transits"],
                r["nearest_post_days"] or 0,
            )


if __name__ == "__main__":
    main()
