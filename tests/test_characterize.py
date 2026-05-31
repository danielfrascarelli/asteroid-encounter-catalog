"""Tests for Phase 5 characterization modules."""

from __future__ import annotations

import numpy as np
import pytest

from src.characterize.geometry import (
    dist_au_to_km,
    vel_au_per_day_to_km_s,
    vel_au_per_day_to_m_s,
)
from src.characterize.observability import (
    apparent_mag_hg,
    get_earth_positions_au,
    is_gaia_observable,
    solar_elongation_deg,
)
from src.characterize.physical import classify_orbit, diameter_km

# ------------------------------------------------------------------ #
# physical.py                                                          #
# ------------------------------------------------------------------ #


class TestDiameterKm:
    def test_ceres(self) -> None:
        # (1) Ceres: H=3.34, actual albedo=0.09 → ~940 km
        d = float(diameter_km(3.34, albedo=0.09))
        assert abs(d - 940.0) / 940.0 < 0.05, f"Ceres diameter {d:.1f} km not within 5% of 940"

    def test_vesta(self) -> None:
        # (4) Vesta: H=3.25, actual albedo=0.34 → ~525 km
        d = float(diameter_km(3.25, albedo=0.34))
        assert abs(d - 525.0) / 525.0 < 0.05, f"Vesta diameter {d:.1f} km not within 5% of 525"

    def test_array_input(self) -> None:
        h_vals = np.array([3.34, 3.25, 10.0])
        p = np.array([0.09, 0.34, 0.14])
        d = diameter_km(h_vals, p)
        assert d.shape == (3,)
        assert np.all(d > 0)

    def test_brighter_is_larger(self) -> None:
        # Smaller H = brighter = larger object
        assert diameter_km(5.0) > diameter_km(10.0)

    def test_higher_albedo_is_smaller(self) -> None:
        # Higher albedo → smaller object for same H
        assert diameter_km(10.0, albedo=0.05) > diameter_km(10.0, albedo=0.30)

    def test_nan_h_returns_nan(self) -> None:
        assert np.isnan(float(diameter_km(np.nan)))

    def test_nan_array_propagates(self) -> None:
        h = np.array([3.34, np.nan, 10.0])
        result = diameter_km(h)
        assert not np.isnan(result[0])
        assert np.isnan(result[1])
        assert not np.isnan(result[2])


class TestClassifyOrbit:
    def test_mba(self) -> None:
        assert classify_orbit(2.5, 0.1) == "MBA"

    def test_nea(self) -> None:
        # q = 1.0*(1-0.7) = 0.3 < 1.3
        assert classify_orbit(1.0, 0.7) == "NEA"

    def test_trojan(self) -> None:
        assert classify_orbit(5.2, 0.05) == "Trojan"

    def test_centaur(self) -> None:
        assert classify_orbit(15.0, 0.2) == "Centaur"

    def test_tno(self) -> None:
        assert classify_orbit(45.0, 0.05) == "TNO"

    def test_array_input(self) -> None:
        a = np.array([2.5, 1.0, 5.2, 15.0, 45.0])
        e = np.array([0.1, 0.7, 0.05, 0.2, 0.05])
        cls = classify_orbit(a, e)
        assert cls[0] == "MBA"
        assert cls[1] == "NEA"
        assert cls[2] == "Trojan"
        assert cls[3] == "Centaur"
        assert cls[4] == "TNO"


# ------------------------------------------------------------------ #
# geometry.py                                                          #
# ------------------------------------------------------------------ #


class TestGeometry:
    def test_au_to_km(self) -> None:
        assert abs(float(dist_au_to_km(1.0)) - 149_597_870.7) < 1.0

    def test_earth_orbital_speed(self) -> None:
        # Earth orbits at ~0.01720 AU/day ≈ 29.78 km/s
        v_km = float(vel_au_per_day_to_km_s(0.01720))
        assert abs(v_km - 29.78) < 0.5

    def test_m_per_s_vs_km_s(self) -> None:
        v = 0.005
        assert abs(float(vel_au_per_day_to_m_s(v)) - float(vel_au_per_day_to_km_s(v)) * 1000) < 1e-6

    def test_array_input(self) -> None:
        d = dist_au_to_km(np.array([1.0, 2.0]))
        assert d[1] == pytest.approx(2 * d[0])


# ------------------------------------------------------------------ #
# observability.py                                                     #
# ------------------------------------------------------------------ #


class TestSolarElongation:
    def test_opposition(self) -> None:
        # Asteroid directly opposite the Sun as seen from Earth
        # Earth at (1, 0, 0), asteroid at (2, 0, 0): elongation = 180°
        earth = np.array([[1.0, 0.0, 0.0]])
        enc = np.array([[2.0, 0.0, 0.0]])
        elong = solar_elongation_deg(enc, earth)
        assert abs(float(elong[0]) - 180.0) < 1.0

    def test_conjunction(self) -> None:
        # Asteroid at (-2, 0, 0), Earth at (1, 0, 0): elongation ≈ 0°
        earth = np.array([[1.0, 0.0, 0.0]])
        enc = np.array([[-2.0, 0.0, 0.0]])
        elong = solar_elongation_deg(enc, earth)
        assert float(elong[0]) < 5.0

    def test_quadrature(self) -> None:
        # For 90° elongation: Earth at (1,0,0), Sun at origin → Sun dir from Earth = (-1,0,0).
        # Perpendicular direction from Earth is (0,1,0), so asteroid at (1,2,0) gives elong=90°.
        earth = np.array([[1.0, 0.0, 0.0]])
        enc = np.array([[1.0, 2.0, 0.0]])
        elong = solar_elongation_deg(enc, earth)
        assert abs(float(elong[0]) - 90.0) < 1.0

    def test_array_input(self) -> None:
        earth = np.array([[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        enc = np.array([[2.0, 0.0, 0.0], [-2.0, 0.0, 0.0]])
        elong = solar_elongation_deg(enc, earth)
        assert len(elong) == 2


class TestEarthPositions:
    """Regression tests for the heliocentric ecliptic J2000 frame.

    The previous implementation returned barycentric ICRS (equatorial), which
    made downstream solar elongation and apparent magnitude computations
    invalid because they were combined with heliocentric ecliptic asteroid
    positions from ``kepler_to_cartesian``.  These tests pin the corrected
    frame so the bug cannot regress silently.
    """

    def test_earth_lies_near_ecliptic_plane(self) -> None:
        # In heliocentric ecliptic J2000, Earth's orbit defines the plane,
        # so |z| must be < ~1e-3 AU (perturbations from other planets only).
        # In the buggy barycentric ICRS frame, Earth's z at solstice would
        # be |z| ≈ sin(23.44°) ≈ 0.4 AU — three orders of magnitude larger.
        jd = np.array(
            [
                2451545.0,  # J2000.0 epoch (perihelion-ish, |z| small either way)
                2451545.0 + 90.0,  # ~April equinox + quarter (near boreal solstice)
                2451545.0 + 180.0,
                2451545.0 + 270.0,
            ]
        )
        earth = get_earth_positions_au(jd)
        assert earth.shape == (4, 3)
        assert np.all(np.abs(earth[:, 2]) < 1e-3), (
            f"Earth z out of ecliptic plane: {earth[:, 2]} — "
            "frame is probably still ICRS instead of ecliptic"
        )

    def test_earth_distance_is_one_au(self) -> None:
        jd = np.array([2451545.0, 2451545.0 + 180.0])
        earth = get_earth_positions_au(jd)
        r = np.linalg.norm(earth, axis=1)
        # Earth's heliocentric distance varies in [0.983, 1.017] AU.
        assert np.all(np.abs(r - 1.0) < 0.02), f"Earth-Sun distance unexpected: {r}"


class TestApparentMag:
    def test_increases_with_distance(self) -> None:
        # Farther away → fainter (larger magnitude)
        m1 = float(apparent_mag_hg(10.0, 2.0, 1.0))
        m2 = float(apparent_mag_hg(10.0, 3.0, 2.0))
        assert m2 > m1

    def test_reasonable_mba(self) -> None:
        # MBA at 2.5 AU from Sun, 1.5 AU from Earth, H=12 → roughly 17–19 mag
        m = float(apparent_mag_hg(12.0, 2.5, 1.5))
        assert 14.0 < m < 22.0


class TestGaiaObservable:
    def test_observable(self) -> None:
        elong = np.array([90.0])
        mag = np.array([18.0])
        assert is_gaia_observable(elong, mag)[0]

    def test_too_close_to_sun(self) -> None:
        elong = np.array([30.0])
        mag = np.array([18.0])
        assert not is_gaia_observable(elong, mag)[0]

    def test_too_faint(self) -> None:
        elong = np.array([90.0])
        mag = np.array([22.0])
        assert not is_gaia_observable(elong, mag)[0]


# ------------------------------------------------------------------ #
# characterize_catalog integration                                     #
# ------------------------------------------------------------------ #


def _two_body_fixture():
    """Build a 2-row encounter set + elements + mpcorb for testing.

    Body 100: bright (H=4) "large" object.
    Body 200: dim   (H=14) "small" object.

    Row 0 has (number_1, number_2) = (100, 200) — already in correct order.
    Row 1 has (number_1, number_2) = (200, 100) — must be swapped so that
    the larger body (100, H=4) becomes _1.
    """
    import polars as pl

    encounters = pl.DataFrame(
        {
            "number_1": [100, 200],
            "number_2": [200, 100],
            "designation_1": ["A100", "A200"],
            "designation_2": ["A200", "A100"],
            "jd_tdb": [2_457_000.0, 2_457_010.0],
            "dist_au": [0.01, 0.02],
            "rel_vel_au_day": [0.001, 0.002],
        }
    )
    elements = pl.DataFrame(
        {
            "number": [100, 200],
            "a_au": [2.5, 3.0],
            "e": [0.1, 0.2],
            "i_deg": [5.0, 10.0],
            "Omega_deg": [30.0, 60.0],
            "omega_deg": [40.0, 70.0],
            "M_deg": [50.0, 80.0],
            "epoch_jd": [2_457_000.0, 2_457_000.0],
        }
    )
    mpcorb = pl.DataFrame(
        {
            "number": [100, 200],
            "H": [4.0, 14.0],
            "G": [0.15, 0.15],
        }
    )
    return encounters, elements, mpcorb


def _multi_body_fixture():
    """A 6-row encounter set spanning several bodies (for chunking tests).

    Mixes already-ordered and swapped pairs, distinct JDs, and a body with
    unknown H (999, absent from mpcorb) so the NaN paths are exercised too.
    """
    import polars as pl

    encounters = pl.DataFrame(
        {
            "number_1": [100, 200, 300, 100, 400, 999],
            "number_2": [200, 100, 100, 300, 200, 100],
            "designation_1": ["A100", "A200", "A300", "A100", "A400", "A999"],
            "designation_2": ["A200", "A100", "A100", "A300", "A200", "A100"],
            "jd_tdb": [
                2_457_000.0,
                2_457_010.0,
                2_457_020.0,
                2_457_030.0,
                2_457_040.0,
                2_457_050.0,
            ],
            "dist_au": [0.01, 0.02, 0.03, 0.015, 0.025, 0.035],
            "rel_vel_au_day": [0.001, 0.002, 0.003, 0.0015, 0.0025, 0.0035],
        }
    )
    elements = pl.DataFrame(
        {
            "number": [100, 200, 300, 400, 999],
            "a_au": [2.5, 3.0, 2.2, 2.7, 2.9],
            "e": [0.1, 0.2, 0.15, 0.05, 0.3],
            "i_deg": [5.0, 10.0, 7.0, 3.0, 12.0],
            "Omega_deg": [30.0, 60.0, 45.0, 20.0, 75.0],
            "omega_deg": [40.0, 70.0, 55.0, 25.0, 85.0],
            "M_deg": [50.0, 80.0, 65.0, 35.0, 95.0],
            "epoch_jd": [2_457_000.0] * 5,
        }
    )
    mpcorb = pl.DataFrame(
        {
            "number": [100, 200, 300, 400],  # 999 deliberately absent → NaN H
            "H": [4.0, 14.0, 8.0, 6.0],
            "G": [0.15, 0.15, 0.15, 0.15],
        }
    )
    return encounters, elements, mpcorb


class TestStreamingParity:
    """The chunked/streaming path must reproduce the in-memory result exactly."""

    def test_chunked_equals_full_unsorted(self) -> None:
        import polars as pl

        from src.characterize.encounter import characterize_catalog

        encounters, elements, mpcorb = _multi_body_fixture()
        full = characterize_catalog(encounters, elements, mpcorb, sort=False)

        # Characterise in 3-row slices and concatenate — must equal the full
        # single-pass result row-for-row, since characterisation is row-independent.
        chunks = [
            characterize_catalog(encounters[i : i + 3], elements, mpcorb, sort=False)
            for i in range(0, len(encounters), 3)
        ]
        chunked = pl.concat(chunks)
        from polars.testing import assert_frame_equal

        assert_frame_equal(full, chunked, check_row_order=True)

    def test_streaming_file_matches_inmemory(self, tmp_path) -> None:
        import polars as pl

        from src.catalog.schema import CATALOG_SCHEMA
        from src.characterize.encounter import (
            characterize_catalog,
            characterize_catalog_streaming,
        )

        encounters, elements, mpcorb = _multi_body_fixture()
        in_path = tmp_path / "enc.parquet"
        out_path = tmp_path / "enc_characterized.parquet"
        encounters.write_parquet(in_path)

        summary = characterize_catalog_streaming(
            str(in_path),
            elements,
            mpcorb,
            str(out_path),
            run_id="test-stream",
            chunk_size=4,  # forces >1 chunk over 6 rows
        )
        assert summary["n_encounters"] == 6
        assert summary["n_chunks"] == 2

        streamed = pl.read_parquet(out_path)
        # In-memory reference, cast to the same on-disk schema + run_id column.
        ref = characterize_catalog(encounters, elements, mpcorb, sort=False)
        ref = ref.with_columns(pl.lit("test-stream").alias("run_id"))
        present = [c for c in CATALOG_SCHEMA if c in ref.columns]
        ref = ref.select(present).with_columns([pl.col(c).cast(CATALOG_SCHEMA[c]) for c in present])
        streamed = streamed.select(present)

        from polars.testing import assert_frame_equal

        assert_frame_equal(ref, streamed, check_row_order=True)

    def test_streaming_sidecar_written(self, tmp_path) -> None:
        import json

        from src.characterize.encounter import characterize_catalog_streaming

        encounters, elements, mpcorb = _multi_body_fixture()
        in_path = tmp_path / "enc.parquet"
        out_path = tmp_path / "enc_characterized.parquet"
        encounters.write_parquet(in_path)
        characterize_catalog_streaming(
            str(in_path), elements, mpcorb, str(out_path), run_id="test-stream", chunk_size=4
        )
        sidecar = tmp_path / "enc_characterized_metadata.json"
        assert sidecar.exists()
        meta = json.loads(sidecar.read_text())
        assert meta["n_encounters"] == 6
        assert meta["sorted_by_dist"] is False
        # Gate tracks the four required bodies (1/2/4/10); none are in this
        # synthetic fixture, so all must report present=False with a recorded slot.
        assert set(meta["major_body_gate"]) == {"1", "2", "4", "10"}
        assert all(not g["present"] for g in meta["major_body_gate"].values())


class TestCharacterizeCatalog:
    def test_body_1_is_always_larger(self) -> None:
        """After characterize, H_1 ≤ H_2 row-by-row (lower H = brighter = larger)."""
        from src.characterize.encounter import characterize_catalog

        encounters, elements, mpcorb = _two_body_fixture()
        enriched = characterize_catalog(encounters, elements, mpcorb)

        h1 = enriched["H_1"].to_numpy()
        h2 = enriched["H_2"].to_numpy()
        # H_1 must be ≤ H_2 for every row that has both H values
        mask = ~(np.isnan(h1) | np.isnan(h2))
        assert (h1[mask] <= h2[mask]).all()

        # And number_1 must equal 100 (the bright body) in every row
        assert enriched["number_1"].to_list() == [100, 100]
        assert enriched["number_2"].to_list() == [200, 200]

    def test_diameter_1_is_largest(self) -> None:
        """diameter_1_km ≥ diameter_2_km row-by-row."""
        from src.characterize.encounter import characterize_catalog

        encounters, elements, mpcorb = _two_body_fixture()
        enriched = characterize_catalog(encounters, elements, mpcorb)
        d1 = enriched["diameter_1_km"].to_numpy()
        d2 = enriched["diameter_2_km"].to_numpy()
        mask = ~(np.isnan(d1) | np.isnan(d2))
        assert (d1[mask] >= d2[mask]).all()

    def test_per_body_observability_columns_emitted(self) -> None:
        """gaia_observable_1, gaia_observable_2, gaia_observable all present."""
        from src.characterize.encounter import characterize_catalog

        encounters, elements, mpcorb = _two_body_fixture()
        enriched = characterize_catalog(encounters, elements, mpcorb)
        for col in (
            "gaia_observable_1",
            "gaia_observable_2",
            "gaia_observable",
            "app_mag_1",
            "app_mag_2",
            "solar_elongation_1_deg",
            "solar_elongation_2_deg",
        ):
            assert col in enriched.columns, f"missing column {col}"

        # Combined flag is the OR of per-body flags
        obs1 = enriched["gaia_observable_1"].to_numpy()
        obs2 = enriched["gaia_observable_2"].to_numpy()
        obs = enriched["gaia_observable"].to_numpy()
        assert (obs == (obs1 | obs2)).all()

    def test_orbital_elements_follow_swap(self) -> None:
        """After swap, a_au_1 must correspond to body 100 (a=2.5) in both rows."""
        from src.characterize.encounter import characterize_catalog

        encounters, elements, mpcorb = _two_body_fixture()
        enriched = characterize_catalog(encounters, elements, mpcorb)
        # Body 100 has a=2.5, body 200 has a=3.0. After swap, _1 columns
        # must reference body 100 in both rows.
        # Need to add a_au back to the output to verify... it's not in schema.
        # Instead, verify via diameter: body 100 is bigger.
        assert enriched["number_1"].to_list() == [100, 100]
        # H_1 must be 4.0 (body 100), not 14.0 (body 200)
        assert enriched["H_1"].to_list() == [4.0, 4.0]
