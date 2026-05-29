"""Tests for null-perturber sampling used by the Stage 3 specificity test."""

from __future__ import annotations

import polars as pl

from src.mass.null_perturbers import asteroids_close_to_target, sample_null_perturbers


def _toy_mpcorb() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "number": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
            "a_au": [2.50, 2.55, 2.60, 2.70, 3.50, 2.45, 2.65, 2.40, 2.80, 2.20],
            "H": [9.0, 9.1, 9.2, 9.3, 5.0, 9.0, 9.2, 9.4, 9.5, 7.0],
        }
    )


def _toy_encounters(target_number: int, close: list[int], dist: float = 0.04) -> pl.DataFrame:
    rows = []
    for other in close:
        rows.append(
            {
                "number_1": target_number,
                "number_2": int(other),
                "dist_au": dist,
            }
        )
    rows.append(
        {
            "number_1": target_number,
            "number_2": 99,
            "dist_au": 0.20,  # too far to count as a real encounter
        }
    )
    return pl.DataFrame(rows)


def test_asteroids_close_to_target_filters_by_distance() -> None:
    catalog = _toy_encounters(target_number=42, close=[1, 5])
    out = asteroids_close_to_target(catalog, target_number=42, min_separation_au=0.1)
    assert out == {1, 5}


def test_asteroids_close_to_target_handles_both_columns() -> None:
    # Pair stored with target on the right side as well.
    catalog = pl.DataFrame(
        {
            "number_1": [7, 42],
            "number_2": [42, 3],
            "dist_au": [0.02, 0.03],
        }
    )
    out = asteroids_close_to_target(catalog, target_number=42, min_separation_au=0.1)
    assert out == {3, 7}


def test_sample_null_perturbers_respects_a_window() -> None:
    mpcorb = _toy_mpcorb()
    encounters = pl.DataFrame(
        {"number_1": [100], "number_2": [999], "dist_au": [0.5]}  # no real overlap
    )
    chosen = sample_null_perturbers(
        target_number=100,
        real_perturber_number=2,
        real_perturber_a_au=2.55,
        mpcorb=mpcorb,
        encounters=encounters,
        n_nulls=10,
        a_window_au=0.2,
        min_separation_au=0.1,
        h_window_mag=None,
        seed=0,
    )
    # a in [2.35, 2.75]: numbers 1, 3, 4, 6, 7, 8 (excluding 2 itself)
    assert set(chosen).issubset({1, 3, 4, 6, 7, 8})
    assert 2 not in chosen
    assert 5 not in chosen  # a = 3.5, outside window
    assert 9 not in chosen  # a = 2.80, outside window


def test_sample_null_perturbers_excludes_real_close_pairs() -> None:
    mpcorb = _toy_mpcorb()
    encounters = _toy_encounters(target_number=100, close=[1, 3, 7])
    chosen = sample_null_perturbers(
        target_number=100,
        real_perturber_number=2,
        real_perturber_a_au=2.55,
        mpcorb=mpcorb,
        encounters=encounters,
        n_nulls=10,
        a_window_au=0.5,
        min_separation_au=0.1,
        h_window_mag=None,
        seed=0,
    )
    # blocked: 1, 3, 7 (real encounters). a in [2.05, 3.05]: 4, 6, 8, 9, 10
    assert not (set(chosen) & {1, 3, 7})


def test_sample_null_perturbers_respects_h_window() -> None:
    mpcorb = _toy_mpcorb()
    encounters = pl.DataFrame({"number_1": [100], "number_2": [999], "dist_au": [0.5]})
    # real perturber 2 has H = 9.1; with h_window=0.2, only H in [8.9, 9.3] allowed.
    chosen = sample_null_perturbers(
        target_number=100,
        real_perturber_number=2,
        real_perturber_a_au=2.55,
        mpcorb=mpcorb,
        encounters=encounters,
        n_nulls=10,
        a_window_au=1.0,
        min_separation_au=0.1,
        h_window_mag=0.2,
        seed=0,
    )
    # H in [8.9, 9.3]: 1 (9.0), 3 (9.2), 4 (9.3), 6 (9.0), 7 (9.2).
    # Excluded: 5 (H=5), 8 (9.4), 9 (9.5), 10 (7). Self: 2.
    assert set(chosen) <= {1, 3, 4, 6, 7}


def test_sample_null_perturbers_deterministic_by_seed() -> None:
    mpcorb = _toy_mpcorb()
    encounters = pl.DataFrame({"number_1": [100], "number_2": [999], "dist_au": [0.5]})
    a = sample_null_perturbers(
        target_number=100,
        real_perturber_number=2,
        real_perturber_a_au=2.55,
        mpcorb=mpcorb,
        encounters=encounters,
        n_nulls=3,
        a_window_au=0.5,
        min_separation_au=0.1,
        h_window_mag=None,
        seed=7,
    )
    b = sample_null_perturbers(
        target_number=100,
        real_perturber_number=2,
        real_perturber_a_au=2.55,
        mpcorb=mpcorb,
        encounters=encounters,
        n_nulls=3,
        a_window_au=0.5,
        min_separation_au=0.1,
        h_window_mag=None,
        seed=7,
    )
    assert a == b


def test_sample_null_perturbers_returns_all_when_pool_small() -> None:
    mpcorb = _toy_mpcorb().head(3)  # only 3 rows
    encounters = pl.DataFrame({"number_1": [100], "number_2": [999], "dist_au": [0.5]})
    chosen = sample_null_perturbers(
        target_number=100,
        real_perturber_number=2,
        real_perturber_a_au=2.55,
        mpcorb=mpcorb,
        encounters=encounters,
        n_nulls=10,
        a_window_au=0.5,
        min_separation_au=0.1,
        h_window_mag=None,
        seed=0,
    )
    # Eligible after dropping self (2): 1 and 3
    assert sorted(chosen) == [1, 3]
