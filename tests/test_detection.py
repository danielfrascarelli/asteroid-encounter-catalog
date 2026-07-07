"""Tests for src/detect/ — prefilter, KD-tree scan, refinement, pipeline.

All tests use synthetic orbital elements and require no real data files.

Synthetic setup
---------------
Three asteroids, epoch JD 2457000.0 TDB (2014-12-09):

  idx  a_au   e     i_deg   Omega  omega   M_deg
   0   2.500  0.0   0.0     0.0    0.0     0.000
   1   2.500  0.0   0.0     0.0    0.0     0.010   ← very close to idx=0
   2   3.500  0.0   0.0     0.0    0.0     0.000   ← different orbit

At t=epoch:
  • separation(0, 1) ≈ 2.5 × 0.010° × π/180 ≈ 4.36 × 10⁻⁴ AU  → encounter
  • separation(0, 2) ≈ 1.0 AU                                    → no encounter

Prefilter (semimajor_diff_max_au=0.5, inclination_diff_max_deg=30°):
  • (0,1) compatible  — |Δa|=0,   |Δi|=0°
  • (0,2) filtered    — |Δa|=1.0 > 0.5
  • (1,2) filtered    — |Δa|=1.0 > 0.5
"""

from __future__ import annotations

import numpy as np
import polars as pl
import pytest

from src.detect.kdtree_scan import scan_time_grid
from src.detect.pipeline import (
    detect_encounters,
    effective_prefilter_mode,
    load_scan_checkpoint,
    write_scan_checkpoint,
)
from src.detect.prefilter import compatible_pairs
from src.detect.refine import _quadratic_min, refine_candidates
from src.propagate.grid import make_time_grid

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

_EPOCH = 2457000.0  # JD TDB
_THRESHOLD = 0.05  # AU — synthetic test value, not production config

_DETECT_KWARGS: dict = dict(
    semimajor_diff_max_au=0.5,
    inclination_diff_max_deg=30.0,
    leaf_size=30,
    fine_step_seconds=60.0,
    window_hours=6.0,
    prefilter_enabled=True,
    refinement_enabled=True,
    n_workers=1,
    chunk_size_days=30.0,
)


def _make_elements(
    *,
    a: list[float],
    e: list[float] | None = None,
    i_deg: list[float] | None = None,
    omega_big_deg: list[float] | None = None,
    omega_deg: list[float] | None = None,
    m_deg: list[float] | None = None,
    epoch_jd: list[float] | None = None,
    numbers: list[int] | None = None,
    designations: list[str] | None = None,
) -> pl.DataFrame:
    n = len(a)
    return pl.DataFrame(
        {
            "number": pl.Series(numbers if numbers is not None else list(range(n)), dtype=pl.Int32),
            "designation": (
                designations if designations is not None else [f"A{k}" for k in range(n)]
            ),
            "a_au": a,
            "e": e if e is not None else [0.0] * n,
            "i_deg": i_deg if i_deg is not None else [0.0] * n,
            "Omega_deg": omega_big_deg if omega_big_deg is not None else [0.0] * n,
            "omega_deg": omega_deg if omega_deg is not None else [0.0] * n,
            "M_deg": m_deg if m_deg is not None else [0.0] * n,
            "epoch_jd": epoch_jd if epoch_jd is not None else [_EPOCH] * n,
        },
        schema={
            "number": pl.Int32,
            "designation": pl.Utf8,
            "a_au": pl.Float64,
            "e": pl.Float64,
            "i_deg": pl.Float64,
            "Omega_deg": pl.Float64,
            "omega_deg": pl.Float64,
            "M_deg": pl.Float64,
            "epoch_jd": pl.Float64,
        },
    )


@pytest.fixture()
def three_asteroids() -> pl.DataFrame:
    """Standard three-asteroid test set: (0,1) close, (0,2)/(1,2) far."""
    return _make_elements(
        a=[2.5, 2.5, 3.5],
        m_deg=[0.0, 0.01, 0.0],
        numbers=[1, 2, 3],
        designations=["Ceres", "Near", "Far"],
    )


@pytest.fixture()
def single_step_grid() -> np.ndarray:
    """One-point time grid at the epoch."""
    return np.array([_EPOCH])


# ===========================================================================
# prefilter.py
# ===========================================================================


def test_prefilter_compatible_pair(three_asteroids: pl.DataFrame) -> None:
    pairs = compatible_pairs(three_asteroids, semimajor_diff_max_au=0.5)
    # Only (0,1) should survive
    assert len(pairs) == 1
    assert tuple(pairs[0]) == (0, 1)


def test_prefilter_filters_by_semimajor(three_asteroids: pl.DataFrame) -> None:
    pairs = compatible_pairs(three_asteroids, semimajor_diff_max_au=0.5)
    pair_tuples = {tuple(p) for p in pairs}
    assert (0, 2) not in pair_tuples
    assert (1, 2) not in pair_tuples


def test_prefilter_filters_by_inclination() -> None:
    elems = _make_elements(
        a=[2.5, 2.5],
        i_deg=[0.0, 35.0],
    )
    pairs = compatible_pairs(elems, inclination_diff_max_deg=30.0)
    assert len(pairs) == 0


def test_prefilter_passes_within_inclination_limit() -> None:
    elems = _make_elements(
        a=[2.5, 2.5],
        i_deg=[0.0, 29.0],
    )
    pairs = compatible_pairs(elems, inclination_diff_max_deg=30.0)
    assert len(pairs) == 1


def test_prefilter_dtype_is_int32(three_asteroids: pl.DataFrame) -> None:
    pairs = compatible_pairs(three_asteroids)
    assert pairs.dtype == np.int32


def test_prefilter_no_pairs_for_single_asteroid() -> None:
    elems = _make_elements(a=[2.5])
    pairs = compatible_pairs(elems)
    assert pairs.shape == (0, 2)


def test_effective_prefilter_mode_disabled() -> None:
    # Config flag off → "disabled" regardless of N.
    assert effective_prefilter_mode(10, enabled=False) == "disabled"
    assert effective_prefilter_mode(1_000_000, enabled=False) == "disabled"


def test_effective_prefilter_mode_applied_small_n() -> None:
    # Enabled and N ≤ max_n → pair-list built.
    assert effective_prefilter_mode(5_000, enabled=True, max_n=5_000) == "applied"
    assert effective_prefilter_mode(100, enabled=True) == "applied"


def test_effective_prefilter_mode_skipped_large_n() -> None:
    # Enabled but N > max_n → pair-list skipped, KD-tree only (the frozen-run case).
    assert effective_prefilter_mode(5_001, enabled=True, max_n=5_000) == "skipped_large_n"
    assert effective_prefilter_mode(150_000, enabled=True) == "skipped_large_n"


def test_prefilter_no_pairs_for_empty_elements() -> None:
    elems = pl.DataFrame(
        schema={
            "number": pl.Int32,
            "designation": pl.Utf8,
            "a_au": pl.Float64,
            "e": pl.Float64,
            "i_deg": pl.Float64,
            "Omega_deg": pl.Float64,
            "omega_deg": pl.Float64,
            "M_deg": pl.Float64,
            "epoch_jd": pl.Float64,
        }
    )
    pairs = compatible_pairs(elems)
    assert pairs.shape == (0, 2)


# ===========================================================================
# refine.py — _quadratic_min
# ===========================================================================


def test_quadratic_min_known_parabola() -> None:
    # d(t) = (t - 5)²  →  minimum at t=5, d=0
    t_min, d_min = _quadratic_min(4.0, 5.0, 6.0, 1.0, 0.0, 1.0)
    assert abs(t_min - 5.0) < 1e-10
    assert abs(d_min) < 1e-10


def test_quadratic_min_shifted_vertex() -> None:
    # d(t) = (t - 5.3)²,  t = [4, 5, 6]
    # d0 = 1.69,  d1 = 0.09,  d2 = 0.49
    d0, d1, d2 = (4.0 - 5.3) ** 2, (5.0 - 5.3) ** 2, (6.0 - 5.3) ** 2
    t_min, d_min = _quadratic_min(4.0, 5.0, 6.0, d0, d1, d2)
    assert abs(t_min - 5.3) < 1e-6
    assert abs(d_min - 0.0) < 1e-10


def test_quadratic_min_downward_parabola_returns_argmin() -> None:
    # Downward parabola: d(t) = -(t-5)² + 10 — minimum on boundary
    t_min, d_min = _quadratic_min(4.0, 5.0, 6.0, 9.0, 10.0, 9.0)
    # denom = 9 - 20 + 9 = -2 < 0  → fall back to argmin (d=9 at t=4 or t=6)
    assert d_min == 9.0
    assert t_min in (4.0, 6.0)


def test_quadratic_min_clamped_to_interval() -> None:
    # Vertex would fall outside [t0, t2]; must be clamped
    t_min, d_min = _quadratic_min(0.0, 1.0, 2.0, 0.1, 1.0, 2.0)
    assert 0.0 <= t_min <= 2.0


# ===========================================================================
# kdtree_scan.py
# ===========================================================================


def test_kdtree_scan_finds_close_pair(
    three_asteroids: pl.DataFrame,
    single_step_grid: np.ndarray,
) -> None:
    pairs = compatible_pairs(three_asteroids)
    results = scan_time_grid(three_asteroids, single_step_grid, pairs, _THRESHOLD)
    # (0,1) should be found
    assert len(results) == 1
    idx_i, idx_j, _t, dist = results[0]
    assert (idx_i, idx_j) == (0, 1)
    assert dist < _THRESHOLD


def test_kdtree_scan_distance_is_accurate(
    three_asteroids: pl.DataFrame,
    single_step_grid: np.ndarray,
) -> None:
    """Reported distance should match the analytical arc-length estimate."""
    pairs = compatible_pairs(three_asteroids)
    results = scan_time_grid(three_asteroids, single_step_grid, pairs, _THRESHOLD)
    _, _, _, dist = results[0]
    # arc ≈ 2.5 AU × 0.01° × π/180
    expected = 2.5 * 0.01 * np.pi / 180.0
    assert abs(dist - expected) < 1e-6


def test_kdtree_scan_respects_pairs_filter(
    three_asteroids: pl.DataFrame,
    single_step_grid: np.ndarray,
) -> None:
    """(0,1) is close but not in the pairs list → must NOT appear in results."""
    empty_pairs = np.empty((0, 2), dtype=np.int32)
    results = scan_time_grid(three_asteroids, single_step_grid, empty_pairs, _THRESHOLD)
    assert results == []


def test_kdtree_scan_no_encounters() -> None:
    """All pairs too far apart — scan returns empty."""
    elems = _make_elements(a=[2.5, 3.5], m_deg=[0.0, 0.0])
    pairs = np.array([[0, 1]], dtype=np.int32)
    grid = np.array([_EPOCH])
    results = scan_time_grid(elems, grid, pairs, threshold_au=0.001)
    assert results == []


def test_kdtree_scan_empty_pairs_returns_empty(
    three_asteroids: pl.DataFrame,
    single_step_grid: np.ndarray,
) -> None:
    results = scan_time_grid(
        three_asteroids, single_step_grid, np.empty((0, 2), dtype=np.int32), 1.0
    )
    assert results == []


# ===========================================================================
# refine.py — refine_candidates
# ===========================================================================


def test_refine_finds_minimum(
    three_asteroids: pl.DataFrame,
) -> None:
    candidates = [(0, 1, _EPOCH, 5e-4)]
    result = refine_candidates(
        three_asteroids,
        candidates,
        threshold_au=_THRESHOLD,
        fine_step_seconds=60.0,
        window_hours=1.0,
    )
    assert len(result) == 1
    assert result["dist_au"][0] < _THRESHOLD
    assert result["rel_vel_au_day"][0] >= 0.0  # non-negative


def test_refine_excludes_above_threshold(
    three_asteroids: pl.DataFrame,
) -> None:
    """Candidate whose refined distance exceeds threshold is dropped."""
    candidates = [(0, 2, _EPOCH, 0.9)]  # (0,2) far apart — will refine to > threshold
    result = refine_candidates(
        three_asteroids,
        candidates,
        threshold_au=_THRESHOLD,
        fine_step_seconds=60.0,
        window_hours=1.0,
    )
    assert len(result) == 0


def test_refine_empty_candidates_returns_empty_df(
    three_asteroids: pl.DataFrame,
) -> None:
    result = refine_candidates(three_asteroids, [], threshold_au=_THRESHOLD)
    assert len(result) == 0
    assert set(result.columns) == {
        "number_1",
        "number_2",
        "designation_1",
        "designation_2",
        "jd_tdb",
        "dist_au",
        "rel_vel_au_day",
    }


def test_refine_output_schema(three_asteroids: pl.DataFrame) -> None:
    candidates = [(0, 1, _EPOCH, 5e-4)]
    result = refine_candidates(three_asteroids, candidates, threshold_au=_THRESHOLD)
    assert result.schema == {
        "number_1": pl.Int32,
        "number_2": pl.Int32,
        "designation_1": pl.Utf8,
        "designation_2": pl.Utf8,
        "jd_tdb": pl.Float64,
        "dist_au": pl.Float64,
        "rel_vel_au_day": pl.Float64,
    }


# ===========================================================================
# refine.py — B1 regression: minimum beyond the fine window (tribunal 2026-07-04)
#
# The coarse scan samples every coarse_step (12 h in production); the true
# minimum can lie up to coarse_step/2 from the nearest coarse sample.  Before
# the fix, a ±2 h fine window clipped the epoch at the window edge without
# interpolating and biased the distance high for ~60 % of the frozen catalog.
# ===========================================================================


def _crossing_pair(delta_deg: float) -> pl.DataFrame:
    """Two circular orbits (a=2.5 AU) with intersecting planes (±15° tilt).

    Both bodies pass the mutual intersection point; body 0 arrives
    ``delta_deg`` of mean anomaly late, producing a close approach with
    realistic relative velocity |v_rel| = 2·a·n·sin(15°) ≈ 5.6×10⁻³ AU/d
    (~9.8 km/s) and miss distance ≈ a·cos(15°)·Δ (rad).
    """
    return _make_elements(
        a=[2.5, 2.5],
        i_deg=[15.0, 15.0],
        omega_big_deg=[0.0, 180.0],
        m_deg=[180.0 - delta_deg, 0.0],
    )


def _dense_scan_minimum(
    elems: pl.DataFrame, t_lo: float, t_hi: float, step_seconds: float = 5.0
) -> tuple[float, float]:
    """Brute-force ground-truth minimum via dense sampling, independent of refine.py."""
    from src.propagate.kepler import kepler_to_cartesian

    t = np.arange(t_lo, t_hi, step_seconds / 86400.0)

    def _pos(k: int) -> np.ndarray:
        row = elems.row(k, named=True)
        return kepler_to_cartesian(
            a_au=row["a_au"],
            e=row["e"],
            i_rad=np.radians(row["i_deg"]),
            Omega_rad=np.radians(row["Omega_deg"]),
            omega_rad=np.radians(row["omega_deg"]),
            M0_rad=np.radians(row["M_deg"]),
            epoch_jd=row["epoch_jd"],
            t_jd=t,
        )

    d = np.linalg.norm(_pos(0) - _pos(1), axis=1)
    k = int(np.argmin(d))
    return float(t[k]), float(d[k])


@pytest.mark.parametrize("window_hours", [6.0, 2.0])
def test_refine_recovers_minimum_beyond_window(window_hours: float) -> None:
    """B1 gate: a minimum ~5 h past the coarse sample (beyond the old ±2 h
    window) must be recovered with |Δt| ≤ fine step and |Δd| ≤ 1 μAU.

    With window_hours=6 the minimum falls inside the first window; with the
    old window_hours=2 it is only reachable through the edge re-centring loop.
    """
    elems = _crossing_pair(delta_deg=0.0237)
    # Ground truth: dense 5 s scan around the plane-crossing epoch.
    t_true, d_true = _dense_scan_minimum(elems, _EPOCH - 0.2, _EPOCH + 0.3)
    assert d_true < 0.0015  # sanity: geometry produces a sub-threshold approach

    # Simulate the coarse sample sitting 5 h before the true minimum.
    t_coarse = t_true - 5.0 / 24.0
    fine_step_seconds = 60.0
    result = refine_candidates(
        elems,
        [(0, 1, t_coarse, 1.0e-2)],
        threshold_au=0.05,
        fine_step_seconds=fine_step_seconds,
        window_hours=window_hours,
    )
    assert len(result) == 1
    dt_days = abs(result["jd_tdb"][0] - t_true)
    dd_au = abs(result["dist_au"][0] - d_true)
    assert dt_days <= fine_step_seconds / 86400.0, f"epoch clipped: Δt = {dt_days * 24:.3f} h"
    assert dd_au <= 1.0e-6, f"distance biased: Δd = {dd_au:.3e} AU"


def test_refine_edge_bias_would_be_caught() -> None:
    """Documents the magnitude of the pre-fix bias: the distance at the old
    ±2 h window edge is > 100 μAU above the true minimum for this geometry, so
    the 1 μAU tolerance above genuinely discriminates the buggy behaviour."""
    elems = _crossing_pair(delta_deg=0.0237)
    t_true, d_true = _dense_scan_minimum(elems, _EPOCH - 0.2, _EPOCH + 0.3)
    t_coarse = t_true - 5.0 / 24.0
    # Distance at the old window edge (t_coarse + 2 h, i.e. 3 h before t_true)
    _, d_at_edge = _dense_scan_minimum(
        elems, t_coarse + 2.0 / 24.0, t_coarse + 2.0 / 24.0 + 1.0 / 86400.0, step_seconds=1.0
    )
    assert d_at_edge - d_true > 1.0e-4


def test_pipeline_rejects_window_narrower_than_half_step(
    three_asteroids: pl.DataFrame,
) -> None:
    """B1 gate: config validation must fail when window_hours < coarse_step/2."""
    grid = make_time_grid(_EPOCH, _EPOCH + 2.0, step_hours=12.0)
    with pytest.raises(ValueError, match="window_hours"):
        detect_encounters(
            three_asteroids,
            grid,
            threshold_au=_THRESHOLD,
            **{**_DETECT_KWARGS, "window_hours": 2.0},
        )


# ===========================================================================
# pipeline.py — scan checkpoint / resume (insurance against refinement crash)
# ===========================================================================


def test_scan_checkpoint_roundtrip(tmp_path) -> None:
    """write_scan_checkpoint → load_scan_checkpoint preserves candidate tuples."""
    cands = [(0, 1, 2457000.0, 5e-4), (3, 7, 2457001.5, 1e-3)]
    p = tmp_path / "scan.parquet"
    write_scan_checkpoint(cands, p)
    loaded = load_scan_checkpoint(p)
    assert loaded == cands


def test_resume_from_scan_matches_full_run(three_asteroids: pl.DataFrame, tmp_path) -> None:
    """A run that writes a scan checkpoint and a resume-from-scan run must yield
    the identical catalog — the insurance path is behaviourally transparent."""
    grid = make_time_grid(_EPOCH, _EPOCH + 0.5, step_hours=6.0)
    ckpt = tmp_path / "scan_candidates.parquet"

    full = detect_encounters(
        three_asteroids,
        grid,
        threshold_au=_THRESHOLD,
        **{**_DETECT_KWARGS, "scan_checkpoint_path": str(ckpt)},
    )
    assert ckpt.exists()  # checkpoint written before refinement

    resumed = detect_encounters(
        three_asteroids,
        grid,
        threshold_au=_THRESHOLD,
        **{**_DETECT_KWARGS, "resume_from_scan": str(ckpt)},
    )
    assert full.to_dicts() == resumed.to_dicts()


# ===========================================================================
# ooc.py — out-of-core detection (disk shards + DuckDB dedup + batched refine)
# ===========================================================================


def test_ooc_dedup_matches_in_memory(tmp_path) -> None:
    """dedup_shards keeps the same min-(dist,t) row per pair as _merge_into."""
    from src.detect.ooc import _SHARD_SCHEMA, dedup_shards
    from src.detect.parallel import _merge_candidates

    shards = [
        [(0, 1, 10.0, 0.03), (0, 1, 9.0, 0.02), (2, 5, 1.0, 0.04)],
        [(0, 1, 8.0, 0.02), (2, 5, 3.0, 0.01)],  # (0,1): tie dist 0.02 → earlier t=8
    ]
    d = tmp_path / "shards"
    d.mkdir()
    for i, ch in enumerate(shards):
        ii, jj, t, dd = zip(*ch)
        pl.DataFrame(
            {"idx_i": ii, "idx_j": jj, "t_coarse_jd": t, "d_coarse_au": dd},
            schema=_SHARD_SCHEMA,
        ).write_parquet(d / f"chunk_{i:05d}.parquet")

    out = tmp_path / "dedup.parquet"
    n = dedup_shards(d, out)
    got = {
        (r["idx_i"], r["idx_j"]): (r["t_coarse_jd"], r["d_coarse_au"])
        for r in pl.read_parquet(out).to_dicts()
    }
    ref = {(a, b): (t, dd) for a, b, t, dd in _merge_candidates(shards)}
    assert n == len(ref)
    assert got == ref  # (0,1)→(8.0,0.02), (2,5)→(3.0,0.01)


def test_ooc_shard_sink_named_by_chunk_id(tmp_path) -> None:
    """Shards are named by chunk id and re-discovered for resume."""
    from src.detect.ooc import existing_shard_ids, make_shard_sink

    d = tmp_path / "shards"
    sink = make_shard_sink(d)
    sink(0, [(0, 1, 10.0, 0.02)])
    sink(20, [(2, 5, 11.0, 0.03)])
    sink(40, [])  # empty chunk still writes a shard so its id counts as done
    assert existing_shard_ids(d) == {0, 20, 40}
    assert (d / "chunk_000000.parquet").exists()
    assert (d / "chunk_000040.parquet").exists()


def test_ooc_scan_skips_done_chunks(three_asteroids: pl.DataFrame, tmp_path) -> None:
    """scan_parallel(skip_chunk_ids=...) drops chunks already scanned (resume)."""
    from src.detect.ooc import existing_shard_ids, make_shard_sink
    from src.detect.parallel import scan_parallel

    grid = make_time_grid(_EPOCH, _EPOCH + 2.0, step_hours=12.0)  # 5 steps
    d = tmp_path / "shards"
    sink = make_shard_sink(d)
    # First pass: scan everything into shards.
    scan_parallel(three_asteroids, grid, None, _THRESHOLD, 30, 2, 1.0, on_chunk=sink)
    done = existing_shard_ids(d)
    assert done  # some shards written
    # Second pass skipping all done ids: nothing new dispatched.
    before = len(list(d.glob("chunk_*.parquet")))
    scan_parallel(
        three_asteroids, grid, None, _THRESHOLD, 30, 2, 1.0, on_chunk=sink, skip_chunk_ids=done
    )
    after = len(list(d.glob("chunk_*.parquet")))
    assert after == before  # no re-scan


def test_ooc_detect_matches_in_memory(three_asteroids: pl.DataFrame, tmp_path) -> None:
    """Full out-of-core detection recovers the same encounter as the in-memory
    path (same pair, refined distance within 1 nAU)."""
    from src.detect.ooc import detect_encounters_ooc

    grid = make_time_grid(_EPOCH, _EPOCH + 0.5, step_hours=6.0)
    ref = detect_encounters(three_asteroids, grid, threshold_au=_THRESHOLD, **_DETECT_KWARGS)

    out = tmp_path / "ooc_catalog.parquet"
    n = detect_encounters_ooc(
        three_asteroids,
        grid,
        threshold_au=_THRESHOLD,
        leaf_size=30,
        fine_step_seconds=60.0,
        window_hours=6.0,
        n_workers=2,
        chunk_size_days=30.0,
        out_path=out,
        spill_dir=tmp_path / "spill",
    )
    ooc = pl.read_parquet(out)
    assert n == len(ref) == len(ooc)
    ref_pairs = set(zip(ref["number_1"].to_list(), ref["number_2"].to_list()))
    ooc_pairs = set(zip(ooc["number_1"].to_list(), ooc["number_2"].to_list()))
    assert ref_pairs == ooc_pairs
    assert abs(ref["dist_au"].min() - ooc["dist_au"].min()) < 1e-9


# ===========================================================================
# pipeline.py — end-to-end
# ===========================================================================


def test_pipeline_end_to_end(
    three_asteroids: pl.DataFrame,
    single_step_grid: np.ndarray,
) -> None:
    """Full pipeline returns exactly one encounter for the synthetic set."""
    result = detect_encounters(
        three_asteroids,
        single_step_grid,
        threshold_au=_THRESHOLD,
        **_DETECT_KWARGS,
    )
    assert len(result) == 1
    row = result.row(0, named=True)
    assert row["number_1"] == 1  # asteroid idx=0 has number=1
    assert row["number_2"] == 2  # asteroid idx=1 has number=2
    assert row["dist_au"] < _THRESHOLD


def test_pipeline_no_compatible_pairs_returns_empty() -> None:
    """When all pairs are filtered out the catalog is empty."""
    elems = _make_elements(a=[2.5, 4.0], m_deg=[0.0, 0.0])
    grid = np.array([_EPOCH])
    result = detect_encounters(elems, grid, threshold_au=_THRESHOLD, **_DETECT_KWARGS)
    assert len(result) == 0


def test_pipeline_output_sorted_by_dist(three_asteroids: pl.DataFrame) -> None:
    """Output must be sorted by dist_au ascending."""
    grid = np.array([_EPOCH])
    result = detect_encounters(three_asteroids, grid, threshold_au=_THRESHOLD, **_DETECT_KWARGS)
    dists = result["dist_au"].to_list()
    assert dists == sorted(dists)


def test_pipeline_output_schema(
    three_asteroids: pl.DataFrame, single_step_grid: np.ndarray
) -> None:
    result = detect_encounters(
        three_asteroids, single_step_grid, threshold_au=_THRESHOLD, **_DETECT_KWARGS
    )
    assert result.schema == {
        "number_1": pl.Int32,
        "number_2": pl.Int32,
        "designation_1": pl.Utf8,
        "designation_2": pl.Utf8,
        "jd_tdb": pl.Float64,
        "dist_au": pl.Float64,
        "rel_vel_au_day": pl.Float64,
    }


def test_pipeline_no_duplicate_pairs(three_asteroids: pl.DataFrame) -> None:
    """Each (number_1, number_2) pair appears at most once."""
    grid = make_time_grid(_EPOCH, _EPOCH + 0.5, step_hours=6.0)
    result = detect_encounters(three_asteroids, grid, threshold_au=_THRESHOLD, **_DETECT_KWARGS)
    pairs = list(zip(result["number_1"].to_list(), result["number_2"].to_list()))
    assert len(pairs) == len(set(pairs))


def test_pipeline_prefilter_disabled(
    three_asteroids: pl.DataFrame, single_step_grid: np.ndarray
) -> None:
    """With prefilter disabled the encounter is still detected."""
    result = detect_encounters(
        three_asteroids,
        single_step_grid,
        threshold_au=_THRESHOLD,
        **{**_DETECT_KWARGS, "prefilter_enabled": False},
    )
    assert len(result) >= 1
    assert result["dist_au"][0] < _THRESHOLD


def test_pipeline_refinement_disabled(
    three_asteroids: pl.DataFrame, single_step_grid: np.ndarray
) -> None:
    """With refinement disabled the coarse epoch and distance are returned."""
    result = detect_encounters(
        three_asteroids,
        single_step_grid,
        threshold_au=_THRESHOLD,
        **{**_DETECT_KWARGS, "refinement_enabled": False},
    )
    assert len(result) == 1
    assert result["dist_au"][0] < _THRESHOLD
    # rel_vel is NaN when refinement is skipped
    assert np.isnan(result["rel_vel_au_day"][0])


def test_pipeline_designation_preserved(
    three_asteroids: pl.DataFrame, single_step_grid: np.ndarray
) -> None:
    result = detect_encounters(
        three_asteroids, single_step_grid, threshold_au=_THRESHOLD, **_DETECT_KWARGS
    )
    assert result["designation_1"][0] == "Ceres"
    assert result["designation_2"][0] == "Near"
