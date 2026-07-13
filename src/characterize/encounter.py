"""Encounter characterization — enriches the raw detection catalog.

Takes the output of the detection pipeline and adds physical properties,
unit conversions, orbit classification, and Gaia observability estimates.

Public entry point: :func:`characterize_catalog`.
"""

from __future__ import annotations

import logging

import numpy as np
import polars as pl
from astropy.time import Time

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
from src.characterize.physical import (
    classify_orbit,
    deflection_dv_m_s,
    diameter_km_with_source,
)
from src.propagate.kepler import kepler_to_cartesian

logger = logging.getLogger(__name__)

_DEG = np.pi / 180.0

# Columns from elements DataFrame needed for position propagation
_ELEM_COLS = ["number", "a_au", "e", "i_deg", "Omega_deg", "omega_deg", "M_deg", "epoch_jd"]


def _join_elements(
    encounters: pl.DataFrame,
    elements: pl.DataFrame,
    suffix: str,
) -> pl.DataFrame:
    """Left-join orbital elements onto encounters for one body (1 or 2).

    Parameters
    ----------
    encounters:
        Encounters DataFrame with ``number_{suffix}`` column.
    elements:
        Orbital elements DataFrame (must include ``number`` and ``_ELEM_COLS``).
    suffix:
        ``"1"`` or ``"2"``.

    Returns
    -------
    pl.DataFrame
        encounters with added columns ``a_au_{suffix}``, ``e_{suffix}``, etc.
    """
    elem_cols = [c for c in _ELEM_COLS if c != "number"]
    renamed = elements.select(["number"] + elem_cols).rename(
        {c: f"{c}_{suffix}" for c in elem_cols}
    )
    return encounters.join(renamed, left_on=f"number_{suffix}", right_on="number", how="left")


def _apply_swap(df: pl.DataFrame, swap: np.ndarray) -> pl.DataFrame:
    """Return df with _1/_2 paired columns swapped where ``swap`` is True.

    Touches every paired column that this module joins or reads: number,
    designation, H, G, and the orbital elements (a_au, e, i_deg, Omega_deg,
    omega_deg, M_deg, epoch_jd). Columns that are unaffected by the swap
    (jd_tdb, dist_au, rel_vel_au_day) are passed through unchanged.
    """
    paired_bases = ["number", "designation", "H", "G"] + [c for c in _ELEM_COLS if c != "number"]
    new_cols = {}
    for base in paired_bases:
        c1, c2 = f"{base}_1", f"{base}_2"
        if c1 not in df.columns or c2 not in df.columns:
            continue
        v1 = df[c1].to_numpy()
        v2 = df[c2].to_numpy()
        new_v1 = np.where(swap, v2, v1)
        new_v2 = np.where(swap, v1, v2)
        new_cols[c1] = pl.Series(c1, new_v1)
        new_cols[c2] = pl.Series(c2, new_v2)
    return df.with_columns(list(new_cols.values()))


def _positions_from_elements(df: pl.DataFrame, suffix: str) -> np.ndarray:
    """Vectorised Kepler propagation for body *suffix* at each encounter JD."""
    return kepler_to_cartesian(
        a_au=df[f"a_au_{suffix}"].to_numpy(),
        e=df[f"e_{suffix}"].to_numpy(),
        i_rad=df[f"i_deg_{suffix}"].to_numpy() * _DEG,
        Omega_rad=df[f"Omega_deg_{suffix}"].to_numpy() * _DEG,
        omega_rad=df[f"omega_deg_{suffix}"].to_numpy() * _DEG,
        M0_rad=df[f"M_deg_{suffix}"].to_numpy() * _DEG,
        epoch_jd=df[f"epoch_jd_{suffix}"].to_numpy(),
        t_jd=df["jd_tdb"].to_numpy(),
    )


def characterize_catalog(
    encounters: pl.DataFrame,
    elements: pl.DataFrame,
    mpcorb: pl.DataFrame,
    albedo: float = 0.14,
    *,
    physical: pl.DataFrame | None = None,
    sort: bool = True,
) -> pl.DataFrame:
    """Enrich the detection catalog with physical and observational properties.

    Parameters
    ----------
    encounters:
        Output of :func:`src.detect.pipeline.detect_encounters`.  Required
        columns: ``number_1``, ``number_2``, ``designation_1``,
        ``designation_2``, ``jd_tdb``, ``dist_au``, ``rel_vel_au_day``.
    elements:
        Orbital elements (gaia_orbits + Horizons supplement).  Must include
        all columns in ``_ELEM_COLS``.
    mpcorb:
        Parsed MPCORB DataFrame with columns ``number``, ``H``, ``G``.
    albedo:
        Default geometric albedo for diameter estimation when an
        asteroid-specific value is unavailable.
    physical:
        Optional measured-physical-data table (from
        ``scripts.ingest.download_sbdb_physical``) with columns ``number``
        (Int32), ``diameter_km`` (Float64, NaN/null where unmeasured) and
        ``albedo`` (Float64).  When provided, diameters follow the B3 priority
        chain (measured D > D(H, measured albedo) > D(H, zone albedo) >
        D(H, default)) and the provenance is emitted in
        ``diameter_source_1/2``.  When None, diameters fall back to the zone
        albedo (source ``zone_albedo``/``default_albedo``).
    sort:
        When True (default) the result is sorted by ``dist_au`` ascending.
        Set to False for the streaming path (:func:`characterize_catalog_streaming`),
        where a global sort across chunks is neither possible nor needed — each
        encounter is characterised independently, so chunked output preserves
        the input row order.

    Returns
    -------
    pl.DataFrame
        Enriched catalog with columns:

        ``number_1``, ``number_2``, ``designation_1``, ``designation_2``,
        ``jd_tdb``, ``date_utc``,
        ``dist_au``, ``dist_km``,
        ``rel_vel_au_day``, ``rel_vel_km_s``, ``rel_vel_m_s``,
        ``H_1``, ``H_2``,
        ``diameter_1_km``, ``diameter_2_km``,
        ``class_1``, ``class_2``,
        ``solar_elongation_deg`` (midpoint, for backwards compatibility),
        ``solar_elongation_1_deg``, ``solar_elongation_2_deg``,
        ``app_mag_1``, ``app_mag_2``,
        ``gaia_observable_1``, ``gaia_observable_2``,
        ``gaia_observable`` (= ``gaia_observable_1 | gaia_observable_2``).

    Notes
    -----
    Bodies are reordered so that ``_1`` is always the **larger** body
    (the putative perturber): the pair is swapped when ``H_2 < H_1`` —
    or when ``H_1`` is NaN and ``H_2`` is not — so any downstream code
    that treats ``number_1`` / ``diameter_1_km`` as the perturber is
    correct by construction. Ties (equal H, or both NaN) keep the
    detector's index ordering.
    """
    n = len(encounters)
    logger.info("Characterizing %d encounters…", n)

    # ------------------------------------------------------------------ #
    # 1. Join orbital elements for both bodies                            #
    # ------------------------------------------------------------------ #
    df = _join_elements(encounters, elements, "1")
    df = _join_elements(df, elements, "2")

    # ------------------------------------------------------------------ #
    # 2. Join H and G from MPCORB                                         #
    # ------------------------------------------------------------------ #
    h_cols = mpcorb.select(["number", "H", "G"])
    df = df.join(
        h_cols.rename({"H": "H_1", "G": "G_1"}), left_on="number_1", right_on="number", how="left"
    )
    df = df.join(
        h_cols.rename({"H": "H_2", "G": "G_2"}), left_on="number_2", right_on="number", how="left"
    )

    # ------------------------------------------------------------------ #
    # 3. Unit conversions                                                  #
    # ------------------------------------------------------------------ #
    dist_km = dist_au_to_km(df["dist_au"].to_numpy())
    vel_km_s = vel_au_per_day_to_km_s(df["rel_vel_au_day"].to_numpy())
    vel_m_s = vel_au_per_day_to_m_s(df["rel_vel_au_day"].to_numpy())

    date_utc = [Time(jd, format="jd", scale="tdb").utc.iso[:10] for jd in df["jd_tdb"].to_list()]

    # ------------------------------------------------------------------ #
    # 4. Diameters                                                         #
    # ------------------------------------------------------------------ #
    h1 = df["H_1"].to_numpy()
    h2 = df["H_2"].to_numpy()

    # ------------------------------------------------------------------ #
    # 4b. Perturber/target ordering                                        #
    # ------------------------------------------------------------------ #
    # The detector emits pairs by index order (i < j), not by mass. Anything
    # downstream that scores or filters on `diameter_1_km` (mass-fit
    # candidate selection, deflection-score ranking) implicitly assumes
    # `_1` is the perturber. Enforce that here: swap whenever body 2 is
    # larger (lower H) than body 1, or when body 1's H is unknown but
    # body 2's is known.
    swap = (h2 < h1) | (np.isnan(h1) & ~np.isnan(h2))
    if swap.any():
        logger.info("Reordering %d/%d pairs so body_1 is the larger body", swap.sum(), n)
        df = _apply_swap(df, swap)
        h1 = df["H_1"].to_numpy()
        h2 = df["H_2"].to_numpy()

    # Measured diameters/albedos (B3): join post-swap so `_1`/`_2` line up with
    # the final perturber/target ordering, then apply the priority chain.
    a1_diam = df["a_au_1"].to_numpy()
    a2_diam = df["a_au_2"].to_numpy()
    dmeas_1 = ameas_1 = dmeas_2 = ameas_2 = None
    if physical is not None:
        phys = physical.select(
            pl.col("number").cast(pl.Int32),
            pl.col("diameter_km"),
            pl.col("albedo"),
        )
        df = df.join(
            phys.rename({"diameter_km": "_diam_meas_1", "albedo": "_albedo_meas_1"}),
            left_on="number_1",
            right_on="number",
            how="left",
        ).join(
            phys.rename({"diameter_km": "_diam_meas_2", "albedo": "_albedo_meas_2"}),
            left_on="number_2",
            right_on="number",
            how="left",
        )
        dmeas_1 = df["_diam_meas_1"].to_numpy()
        ameas_1 = df["_albedo_meas_1"].to_numpy()
        dmeas_2 = df["_diam_meas_2"].to_numpy()
        ameas_2 = df["_albedo_meas_2"].to_numpy()

    diam_1, diam_src_1 = diameter_km_with_source(
        h1, a1_diam, dmeas_1, ameas_1, default_albedo=albedo
    )
    diam_2, diam_src_2 = diameter_km_with_source(
        h2, a2_diam, dmeas_2, ameas_2, default_albedo=albedo
    )

    # Señal de deflexión por par (M7): kick Δv del perturbador (cuerpo 1, el
    # grande tras el reordenamiento) sobre el cuerpo 2. Métrica de ranking de
    # utilidad para determinación de masas.
    deflection_dv = deflection_dv_m_s(diam_1, a1_diam, dist_km * 1e3, vel_m_s)

    # ------------------------------------------------------------------ #
    # 5. Orbit classification                                              #
    # ------------------------------------------------------------------ #
    a1 = df["a_au_1"].to_numpy()
    e1 = df["e_1"].to_numpy()
    a2 = df["a_au_2"].to_numpy()
    e2 = df["e_2"].to_numpy()
    cls1 = classify_orbit(a1, e1).astype(str)
    cls2 = classify_orbit(a2, e2).astype(str)

    # ------------------------------------------------------------------ #
    # 6. Heliocentric positions at encounter epoch (vectorised Kepler)    #
    # ------------------------------------------------------------------ #
    logger.info("Propagating asteroid positions for observability…")
    pos1 = _positions_from_elements(df, "1")  # (N, 3)
    pos2 = _positions_from_elements(df, "2")  # (N, 3)
    enc_pos_mid = 0.5 * (pos1 + pos2)  # midpoint, kept for backwards compatibility

    r_helio_1 = np.linalg.norm(pos1, axis=1)
    r_helio_2 = np.linalg.norm(pos2, axis=1)

    # ------------------------------------------------------------------ #
    # 7. Earth positions and solar elongation                             #
    # ------------------------------------------------------------------ #
    # get_earth_positions_au returns heliocentric ecliptic J2000 to match
    # kepler_to_cartesian (which produced pos1, pos2 above) — both vectors
    # live in the same frame, so subtractions and dot products are valid.
    logger.info("Fetching Earth positions from astropy ephemeris…")
    jd_arr = df["jd_tdb"].to_numpy()
    earth_pos = get_earth_positions_au(jd_arr)  # (N, 3)

    elong_mid = solar_elongation_deg(enc_pos_mid, earth_pos)
    elong_1 = solar_elongation_deg(pos1, earth_pos)
    elong_2 = solar_elongation_deg(pos2, earth_pos)

    # ------------------------------------------------------------------ #
    # 8. Per-body apparent magnitude and Gaia observability               #
    # ------------------------------------------------------------------ #
    # Each body has its own Earth distance, heliocentric distance, H and G.
    # Observability is reported per body; the legacy combined flag is the
    # logical OR (either body in the encounter was detectable by Gaia).
    g1 = np.where(np.isnan(df["G_1"].to_numpy()), 0.15, df["G_1"].to_numpy())
    g2 = np.where(np.isnan(df["G_2"].to_numpy()), 0.15, df["G_2"].to_numpy())
    delta_1 = np.linalg.norm(pos1 - earth_pos, axis=1)
    delta_2 = np.linalg.norm(pos2 - earth_pos, axis=1)
    app_mag_1 = np.where(
        np.isnan(h1),
        np.nan,
        apparent_mag_hg(np.where(np.isnan(h1), 20.0, h1), r_helio_1, delta_1, G=g1),
    )
    app_mag_2 = np.where(
        np.isnan(h2),
        np.nan,
        apparent_mag_hg(np.where(np.isnan(h2), 20.0, h2), r_helio_2, delta_2, G=g2),
    )
    observable_1 = is_gaia_observable(elong_1, app_mag_1)
    observable_2 = is_gaia_observable(elong_2, app_mag_2)
    observable = observable_1 | observable_2

    # ------------------------------------------------------------------ #
    # 9. Assemble output DataFrame                                         #
    # ------------------------------------------------------------------ #
    logger.info("Assembling enriched catalog…")
    enriched = pl.DataFrame(
        {
            "number_1": df["number_1"],
            "number_2": df["number_2"],
            "designation_1": df["designation_1"],
            "designation_2": df["designation_2"],
            "jd_tdb": df["jd_tdb"],
            "date_utc": date_utc,
            "dist_au": df["dist_au"],
            "dist_km": dist_km,
            "rel_vel_au_day": df["rel_vel_au_day"],
            "rel_vel_km_s": vel_km_s,
            "rel_vel_m_s": vel_m_s,
            "H_1": h1,
            "H_2": h2,
            "diameter_1_km": diam_1,
            "diameter_2_km": diam_2,
            "diameter_source_1": diam_src_1.astype(str),
            "diameter_source_2": diam_src_2.astype(str),
            "deflection_dv_m_s": deflection_dv,
            "class_1": cls1,
            "class_2": cls2,
            "solar_elongation_deg": elong_mid,
            "solar_elongation_1_deg": elong_1,
            "solar_elongation_2_deg": elong_2,
            "app_mag_1": app_mag_1,
            "app_mag_2": app_mag_2,
            "gaia_observable_1": observable_1,
            "gaia_observable_2": observable_2,
            "gaia_observable": observable,
        }
    )
    return enriched.sort("dist_au") if sort else enriched


# Detection columns read from the input catalog — the only ones
# characterize_catalog needs. Reading just these keeps memory bounded even when
# the source carries the wide hybrid schema (refinement_method, *_kepler, …).
_DETECTION_COLS = [
    "number_1",
    "number_2",
    "designation_1",
    "designation_2",
    "jd_tdb",
    "dist_au",
    "rel_vel_au_day",
]


def characterize_catalog_streaming(
    input_path: str,
    elements: pl.DataFrame,
    mpcorb: pl.DataFrame,
    out_path: str,
    run_id: str,
    *,
    albedo: float = 0.14,
    physical: pl.DataFrame | None = None,
    chunk_size: int = 1_000_000,
    mpcorb_path: object | None = None,
    config_dict: dict | None = None,
) -> dict:
    """Characterise an arbitrarily large catalog by streaming it in chunks.

    The in-memory :func:`characterize_catalog` materialises ~10 N-length arrays
    (two Kepler position fields, Earth positions, …) and OOMs on the 72 M-row
    frozen catalog (~31 GB peak). Characterisation is **row-independent** (the
    only cross-row step is the final sort), so this reads the input parquet in
    ``chunk_size`` batches, characterises each (``sort=False``), and appends to a
    single output parquet via a streaming writer — peak RAM is set by one chunk,
    not the whole catalog.

    Only :data:`_DETECTION_COLS` are read from the input, so a wide hybrid-schema
    source is handled without loading its extra columns. The output is **not**
    globally sorted by ``dist_au`` (unlike the in-memory path); it preserves the
    input row order. A provenance sidecar (``<stem>_metadata.json``) is written
    next to the output, mirroring :func:`src.catalog.writer.write_catalog`.

    Parameters
    ----------
    input_path:
        Detection catalog parquet (Kepler or hybrid). Read in batches.
    elements, mpcorb, albedo:
        As for :func:`characterize_catalog`.
    out_path:
        Destination parquet for the enriched full catalog.
    run_id:
        Run identifier written to every row and the sidecar.
    chunk_size:
        Rows per batch. 1 M ≈ <2 GB peak per chunk.
    mpcorb_path, config_dict:
        Optional provenance for the sidecar.

    Returns
    -------
    dict
        Summary: ``n_encounters``, ``n_gaia_observable``, ``n_chunks``, and the
        major-body ``gate`` (presence + closest approach per required body).
    """
    import json
    from datetime import UTC, datetime

    import pyarrow.parquet as pq

    from src.catalog.schema import CATALOG_SCHEMA
    from src.catalog.writer import _dep_versions, _hash_file

    required_bodies = (1, 2, 4, 10)

    pf = pq.ParquetFile(input_path)
    writer: pq.ParquetWriter | None = None
    n_total = 0
    n_obs = 0
    n_chunks = 0
    gate: dict[int, dict[str, float]] = {
        b: {"n_encounters": 0, "closest_au": float("inf")} for b in required_bodies
    }

    logger.info(
        "Streaming characterisation: %s → %s (chunk_size=%d)", input_path, out_path, chunk_size
    )
    try:
        for batch in pf.iter_batches(batch_size=chunk_size, columns=_DETECTION_COLS):
            chunk = pl.from_arrow(batch)
            if isinstance(chunk, pl.Series):  # single-column safety (never here)
                chunk = chunk.to_frame()
            enriched = characterize_catalog(
                chunk, elements, mpcorb, albedo=albedo, physical=physical, sort=False
            )
            enriched = enriched.with_columns(pl.lit(run_id).alias("run_id"))
            present = [c for c in CATALOG_SCHEMA if c in enriched.columns]
            enriched = enriched.select(present).with_columns(
                [pl.col(c).cast(CATALOG_SCHEMA[c]) for c in present]
            )

            n_total += len(enriched)
            n_obs += int(enriched["gaia_observable"].sum())
            n_chunks += 1
            for b in required_bodies:
                hits = enriched.filter((pl.col("number_1") == b) | (pl.col("number_2") == b))
                if len(hits):
                    gate[b]["n_encounters"] += len(hits)
                    closest = float(hits["dist_au"].min())  # type: ignore[arg-type]
                    gate[b]["closest_au"] = min(gate[b]["closest_au"], closest)

            table = enriched.to_arrow()
            if writer is None:
                writer = pq.ParquetWriter(out_path, table.schema, compression="zstd")
            else:
                table = table.cast(writer.schema)  # guard against per-chunk type drift
            writer.write_table(table)
            logger.info("  chunk %d: %d rows (total %d)", n_chunks, len(enriched), n_total)
    finally:
        if writer is not None:
            writer.close()

    gate_out = {
        str(b): {
            "present": gate[b]["n_encounters"] > 0,
            "n_encounters": gate[b]["n_encounters"],
            "closest_au": (
                None if gate[b]["closest_au"] == float("inf") else gate[b]["closest_au"]
            ),
        }
        for b in required_bodies
    }

    meta: dict = {
        "run_id": run_id,
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "n_encounters": n_total,
        "n_gaia_observable": n_obs,
        "n_chunks": n_chunks,
        "chunk_size": chunk_size,
        "input_catalog": str(input_path),
        "catalog_path": str(out_path),
        "sorted_by_dist": False,
        "schema_columns": list(CATALOG_SCHEMA),
        "major_body_gate": gate_out,
        "dependencies": _dep_versions(),
    }
    from pathlib import Path as _Path

    mp = _Path(str(mpcorb_path)) if mpcorb_path is not None else None
    if mp is not None and mp.exists():
        meta["mpcorb"] = {
            "path": str(mp),
            "sha256": _hash_file(mp),
            "size_bytes": mp.stat().st_size,
        }
    if config_dict is not None:
        meta["config"] = config_dict

    sidecar = _Path(out_path).parent / (_Path(out_path).stem + "_metadata.json")
    sidecar.write_text(json.dumps(meta, indent=2, default=str))
    logger.info(
        "Streaming characterisation done: %d rows in %d chunks → %s (sidecar %s)",
        n_total,
        n_chunks,
        out_path,
        sidecar,
    )
    return {
        "n_encounters": n_total,
        "n_gaia_observable": n_obs,
        "n_chunks": n_chunks,
        "gate": gate_out,
    }
