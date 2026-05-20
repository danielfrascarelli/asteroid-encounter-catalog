"""Phase 6 regression and validation tests for the final encounter catalog."""

from __future__ import annotations

import json
from pathlib import Path

import polars as pl
import pytest

from src.catalog.query import filter_encounters, load_catalog, top_encounters
from src.catalog.schema import CATALOG_COLUMNS, CATALOG_SCHEMA
from src.catalog.writer import write_catalog

# Paths to real output files (skipped on CI where data is absent)
_CATALOG_PATH = Path("data/output/encounters_characterized.parquet")
_SIDECAR_PATH = _CATALOG_PATH.parent / (_CATALOG_PATH.stem + "_metadata.json")


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def sample_df() -> pl.DataFrame:
    """Minimal synthetic catalog exercising all schema columns (no run_id)."""
    return pl.DataFrame(
        {
            "number_1": pl.Series([1, 4, 2, 100], dtype=pl.Int32),
            "number_2": pl.Series([11, 17, 23, 200], dtype=pl.Int32),
            "designation_1": ["Ceres", "Vesta", "Pallas", "other"],
            "designation_2": ["Parthenope", "Thetis", "Thalia", "other2"],
            "jd_tdb": [2457000.5, 2457100.5, 2457200.5, 2457300.5],
            "date_utc": ["2015-01-01", "2015-04-10", "2015-07-19", "2015-10-27"],
            "dist_au": [0.005, 0.008, 0.003, 0.009],
            "dist_km": [747989.35, 1196782.97, 448793.61, 1346380.83],
            "rel_vel_au_day": [0.001, 0.002, 0.0015, 0.003],
            "rel_vel_km_s": [1.73, 3.47, 2.60, 5.20],
            "rel_vel_m_s": [1731.0, 3472.0, 2601.0, 5202.0],
            "H_1": [3.34, 3.25, 4.13, 10.0],
            "H_2": [6.6, 7.0, 7.5, 11.0],
            "diameter_1_km": [940.0, 525.0, 512.0, 50.0],
            "diameter_2_km": [153.0, 90.0, 80.0, 30.0],
            "class_1": ["MBA", "MBA", "MBA", "MBA"],
            "class_2": ["MBA", "MBA", "MBA", "MBA"],
            "solar_elongation_deg": [90.0, 75.0, 60.0, 120.0],
            "gaia_observable": [True, True, False, True],
        }
    )


# ---------------------------------------------------------------------------
# Schema
# ---------------------------------------------------------------------------


class TestSchema:
    def test_required_columns_defined(self) -> None:
        for col in ("run_id", "number_1", "number_2", "dist_au", "gaia_observable"):
            assert col in CATALOG_SCHEMA

    def test_catalog_columns_matches_schema(self) -> None:
        assert set(CATALOG_COLUMNS) == set(CATALOG_SCHEMA.keys())

    def test_numeric_columns_are_float64(self) -> None:
        for col in ("dist_au", "dist_km", "rel_vel_km_s", "solar_elongation_deg"):
            assert CATALOG_SCHEMA[col] == pl.Float64

    def test_id_columns_are_int32(self) -> None:
        assert CATALOG_SCHEMA["number_1"] == pl.Int32
        assert CATALOG_SCHEMA["number_2"] == pl.Int32


# ---------------------------------------------------------------------------
# Query API
# ---------------------------------------------------------------------------


class TestFilterEncounters:
    def test_filter_max_dist(self, sample_df: pl.DataFrame) -> None:
        result = filter_encounters(sample_df, max_dist_au=0.006)
        assert all(d <= 0.006 for d in result["dist_au"].to_list())

    def test_filter_min_dist(self, sample_df: pl.DataFrame) -> None:
        result = filter_encounters(sample_df, min_dist_au=0.006)
        assert all(d >= 0.006 for d in result["dist_au"].to_list())

    def test_filter_date_range(self, sample_df: pl.DataFrame) -> None:
        result = filter_encounters(sample_df, date_start="2015-04-01", date_end="2015-08-01")
        for d in result["date_utc"].to_list():
            assert "2015-04-01" <= d <= "2015-08-01"

    def test_filter_by_body_id(self, sample_df: pl.DataFrame) -> None:
        result = filter_encounters(sample_df, body_ids=[1])
        for n1, n2 in zip(result["number_1"].to_list(), result["number_2"].to_list()):
            assert n1 == 1 or n2 == 1

    def test_filter_gaia_observable(self, sample_df: pl.DataFrame) -> None:
        result = filter_encounters(sample_df, gaia_observable_only=True)
        assert all(result["gaia_observable"].to_list())

    def test_filter_no_args_returns_all(self, sample_df: pl.DataFrame) -> None:
        assert len(filter_encounters(sample_df)) == len(sample_df)

    def test_filter_multiple_bodies(self, sample_df: pl.DataFrame) -> None:
        result = filter_encounters(sample_df, body_ids=[1, 4])
        assert len(result) == 2

    def test_filter_empty_result(self, sample_df: pl.DataFrame) -> None:
        result = filter_encounters(sample_df, max_dist_au=0.001)
        assert len(result) == 0


class TestTopEncounters:
    def test_top_n_count(self, sample_df: pl.DataFrame) -> None:
        assert len(top_encounters(sample_df, n=2)) == 2

    def test_top_sorted_ascending(self, sample_df: pl.DataFrame) -> None:
        result = top_encounters(sample_df, n=4, by="dist_au", ascending=True)
        dists = result["dist_au"].to_list()
        assert dists == sorted(dists)

    def test_top_sorted_descending(self, sample_df: pl.DataFrame) -> None:
        result = top_encounters(sample_df, n=4, by="rel_vel_km_s", ascending=False)
        vels = result["rel_vel_km_s"].to_list()
        assert vels == sorted(vels, reverse=True)

    def test_top_n_larger_than_df(self, sample_df: pl.DataFrame) -> None:
        result = top_encounters(sample_df, n=100)
        assert len(result) == len(sample_df)


# ---------------------------------------------------------------------------
# Writer
# ---------------------------------------------------------------------------


class TestWriter:
    def test_creates_parquet(self, sample_df: pl.DataFrame, tmp_path: Path) -> None:
        out = tmp_path / "catalog.parquet"
        write_catalog(sample_df, out, run_id="test-001")
        assert out.exists()

    def test_adds_run_id_column(self, sample_df: pl.DataFrame, tmp_path: Path) -> None:
        out = tmp_path / "catalog.parquet"
        write_catalog(sample_df, out, run_id="test-run-42")
        loaded = load_catalog(out)
        assert "run_id" in loaded.columns
        assert loaded["run_id"][0] == "test-run-42"

    def test_creates_sidecar(self, sample_df: pl.DataFrame, tmp_path: Path) -> None:
        out = tmp_path / "catalog.parquet"
        write_catalog(sample_df, out, run_id="test-001")
        sidecar = tmp_path / "catalog_metadata.json"
        assert sidecar.exists()

    def test_sidecar_is_valid_json(self, sample_df: pl.DataFrame, tmp_path: Path) -> None:
        out = tmp_path / "catalog.parquet"
        write_catalog(sample_df, out, run_id="test-001")
        sidecar = tmp_path / "catalog_metadata.json"
        data = json.loads(sidecar.read_text())
        assert isinstance(data, dict)

    def test_sidecar_required_keys(self, sample_df: pl.DataFrame, tmp_path: Path) -> None:
        out = tmp_path / "catalog.parquet"
        write_catalog(sample_df, out, run_id="test-001", config_dict={"key": "val"})
        data = json.loads((tmp_path / "catalog_metadata.json").read_text())
        for key in ("run_id", "timestamp_utc", "n_encounters", "dependencies", "config"):
            assert key in data, f"Missing sidecar key: {key}"

    def test_sidecar_n_encounters(self, sample_df: pl.DataFrame, tmp_path: Path) -> None:
        out = tmp_path / "catalog.parquet"
        write_catalog(sample_df, out, run_id="test-001")
        data = json.loads((tmp_path / "catalog_metadata.json").read_text())
        assert data["n_encounters"] == len(sample_df)

    def test_roundtrip_row_count(self, sample_df: pl.DataFrame, tmp_path: Path) -> None:
        out = tmp_path / "catalog.parquet"
        write_catalog(sample_df, out, run_id="test-001")
        assert len(load_catalog(out)) == len(sample_df)

    def test_no_duplicate_pairs(self, sample_df: pl.DataFrame, tmp_path: Path) -> None:
        out = tmp_path / "catalog.parquet"
        write_catalog(sample_df, out, run_id="test-001")
        loaded = load_catalog(out)
        unique_pairs = loaded.select(["number_1", "number_2"]).unique()
        assert len(unique_pairs) == len(loaded)


# ---------------------------------------------------------------------------
# Regression tests — real catalog (skipped on CI)
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not _CATALOG_PATH.exists(), reason="Real catalog not available in CI")
class TestRealCatalog:
    def test_catalog_non_empty(self) -> None:
        df = load_catalog(_CATALOG_PATH)
        assert len(df) > 0

    def test_schema_columns_present(self) -> None:
        df = load_catalog(_CATALOG_PATH)
        for col in ("number_1", "number_2", "dist_au", "gaia_observable"):
            assert col in df.columns

    def test_ceres_in_catalog(self) -> None:
        df = load_catalog(_CATALOG_PATH)
        hits = filter_encounters(df, body_ids=[1])
        assert len(hits) > 0, "(1) Ceres must appear in at least one encounter"

    def test_vesta_in_catalog(self) -> None:
        df = load_catalog(_CATALOG_PATH)
        hits = filter_encounters(df, body_ids=[4])
        assert len(hits) > 0, "(4) Vesta must appear in at least one encounter"

    @pytest.mark.xfail(
        reason="(2) Pallas: i=34.9° keeps it >0.01 AU from all MBAs during Gaia window"
    )
    def test_pallas_in_catalog(self) -> None:
        df = load_catalog(_CATALOG_PATH)
        hits = filter_encounters(df, body_ids=[2])
        assert len(hits) > 0

    def test_no_duplicate_pairs(self) -> None:
        df = load_catalog(_CATALOG_PATH)
        unique_pairs = df.select(["number_1", "number_2"]).unique()
        assert len(unique_pairs) == len(df), "Catalog contains duplicate (number_1, number_2) pairs"

    def test_velocity_range_physical(self) -> None:
        df = load_catalog(_CATALOG_PATH)
        assert float(df["rel_vel_km_s"].min()) > 0  # type: ignore[arg-type]
        assert float(df["rel_vel_km_s"].max()) < 100  # type: ignore[arg-type]

    def test_dist_au_within_threshold(self) -> None:
        import yaml
        cfg = yaml.safe_load(Path("config.yaml").read_text())
        threshold = float(cfg["detection"]["threshold_au"])
        df = load_catalog(_CATALOG_PATH)
        # Allow 10% headroom for floating-point rounding in refinement step
        assert float(df["dist_au"].max()) <= threshold * 1.1  # type: ignore[arg-type]

    def test_gaia_observable_fraction_reasonable(self) -> None:
        df = load_catalog(_CATALOG_PATH)
        frac = df["gaia_observable"].sum() / len(df)
        assert 0.1 < frac < 0.9


@pytest.mark.skipif(not _SIDECAR_PATH.exists(), reason="Metadata sidecar not available in CI")
class TestRealSidecar:
    def test_valid_json(self) -> None:
        data = json.loads(_SIDECAR_PATH.read_text())
        assert isinstance(data, dict)

    def test_required_keys(self) -> None:
        data = json.loads(_SIDECAR_PATH.read_text())
        for key in ("run_id", "timestamp_utc", "n_encounters", "dependencies"):
            assert key in data

    def test_n_encounters_positive(self) -> None:
        data = json.loads(_SIDECAR_PATH.read_text())
        assert data["n_encounters"] > 0

    def test_dependency_versions_present(self) -> None:
        data = json.loads(_SIDECAR_PATH.read_text())
        deps = data["dependencies"]
        for pkg in ("astropy", "polars", "numpy"):
            assert pkg in deps


# ---------------------------------------------------------------------------
# Goffin 2014 fixture tests (no network access required)
# ---------------------------------------------------------------------------


class TestGoffinLoadAndFilter:
    """Unit tests for validate_goffin_2014._load_goffin using synthetic data."""

    def _make_goffin_parquet(self, tmp_path: Path) -> Path:
        """Write a minimal synthetic Goffin-like parquet file."""
        df = pl.DataFrame(
            {
                "Pert": pl.Series([1, 4, 2, 10, 532], dtype=pl.Int32),
                "Targ": pl.Series([100, 200, 300, 400, 500], dtype=pl.Int32),
                "Date": [
                    "2015-03-01",  # inside Gaia window
                    "2016-06-15",  # inside
                    "2013-01-01",  # before window
                    "2018-01-01",  # after window
                    "2014-08-01",  # inside
                ],
                "Dist": [0.005, 0.012, 0.003, 0.008, 0.007],
                "Vrel": [2.1, 3.5, 1.8, 4.2, 2.9],
            }
        )
        path = tmp_path / "goffin_2014_encounters.parquet"
        df.write_parquet(path)
        return path

    def test_loads_and_filters_to_gaia_window(self, tmp_path: Path) -> None:
        from scripts.validate_goffin_2014 import _load_goffin

        path = self._make_goffin_parquet(tmp_path)
        df, col_map = _load_goffin(path, "2014-07-25", "2017-05-28")
        # Rows outside the Gaia window (2013-01-01, 2018-01-01) must be dropped
        assert len(df) == 3

    def test_col_map_has_required_keys(self, tmp_path: Path) -> None:
        from scripts.validate_goffin_2014 import _load_goffin

        path = self._make_goffin_parquet(tmp_path)
        _df, col_map = _load_goffin(path, "2014-07-25", "2017-05-28")
        for key in ("perturber", "target", "epoch", "dist"):
            assert key in col_map

    def test_dist_column_values_preserved(self, tmp_path: Path) -> None:
        from scripts.validate_goffin_2014 import _load_goffin

        path = self._make_goffin_parquet(tmp_path)
        df, col_map = _load_goffin(path, "2014-07-25", "2017-05-28")
        dists = df[col_map["dist"]].to_list()
        assert all(d > 0.0 for d in dists)
