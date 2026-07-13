"""Descarga diámetros y albedos medidos de JPL SBDB para asteroides numerados.

Fuente: SBDB Query API (``ssd-api.jpl.nasa.gov/sbdb_query.api``), que compila
mediciones IRAS/AKARI/NEOWISE/ocultaciones. Motivación: hallazgo B3 del tribunal
(2026-07-04) — la caracterización aplicaba albedo fijo 0.14 a todos los cuerpos
(Ceres 763 km en vez de 939). Este dataset alimenta la cadena de prioridades de
:func:`src.characterize.physical.diameter_km_with_source`.

Salida
------
``data/raw/sbdb_physical.parquet`` con columnas ``number`` (Int32),
``diameter_km`` (Float64, null si no medido), ``diameter_sigma_km`` (Float64),
``albedo`` (Float64, null si no medido), más un sidecar JSON con URL, fecha,
SHA-256 y conteos.

Uso
---
    docker compose run --rm pipeline python -m scripts.ingest.download_sbdb_physical
"""

from __future__ import annotations

import hashlib
import json
import logging
import urllib.parse
import urllib.request
from datetime import UTC, datetime
from pathlib import Path

import polars as pl

logger = logging.getLogger(__name__)

_API_URL = "https://ssd-api.jpl.nasa.gov/sbdb_query.api"
_OUT_PATH = Path("data/raw/sbdb_physical.parquet")

# pdes = designación primaria (para numerados, el número como string).
_FIELDS = "pdes,diameter,diameter_sigma,albedo"


def _query_sbdb() -> dict:
    """Una sola query bulk: asteroides numerados, campos físicos."""
    params = {
        "fields": _FIELDS,
        "sb-ns": "n",  # numbered only
        "sb-kind": "a",  # asteroids only
        "full-prec": "false",
    }
    url = f"{_API_URL}?{urllib.parse.urlencode(params)}"
    logger.info("Query SBDB: %s", url)
    req = urllib.request.Request(url, headers={"User-Agent": "gaia-encounters-pipeline"})
    with urllib.request.urlopen(req, timeout=600) as resp:
        return json.loads(resp.read().decode())


def _to_frame(payload: dict) -> pl.DataFrame:
    fields = payload["fields"]
    idx = {name: k for k, name in enumerate(fields)}
    rows = payload["data"]

    def _col(name: str) -> list:
        j = idx[name]
        return [r[j] for r in rows]

    def _float(values: list) -> list[float | None]:
        out: list[float | None] = []
        for v in values:
            try:
                out.append(float(v) if v not in (None, "") else None)
            except (TypeError, ValueError):
                out.append(None)
        return out

    numbers: list[int | None] = []
    for v in _col("pdes"):
        try:
            numbers.append(int(v))
        except (TypeError, ValueError):
            numbers.append(None)

    df = pl.DataFrame(
        {
            "number": pl.Series(numbers, dtype=pl.Int32),
            "diameter_km": pl.Series(_float(_col("diameter")), dtype=pl.Float64),
            "diameter_sigma_km": pl.Series(_float(_col("diameter_sigma")), dtype=pl.Float64),
            "albedo": pl.Series(_float(_col("albedo")), dtype=pl.Float64),
        }
    )
    return df.filter(pl.col("number").is_not_null()).sort("number")


def main() -> int:
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s"
    )
    payload = _query_sbdb()
    df = _to_frame(payload)

    n_total = len(df)
    n_diam = int(df["diameter_km"].is_not_null().sum())
    n_albedo = int(df["albedo"].is_not_null().sum())
    logger.info("SBDB: %d numerados, %d con diámetro, %d con albedo", n_total, n_diam, n_albedo)

    # Sanity gates (B3): Ceres y Nysa con valores medidos correctos.
    ceres = df.filter(pl.col("number") == 1)
    if len(ceres) and ceres["diameter_km"][0] is not None:
        d_ceres = ceres["diameter_km"][0]
        assert 900.0 < d_ceres < 1000.0, f"Ceres diameter {d_ceres} fuera de rango"
        logger.info("Gate Ceres OK: D = %.1f km", d_ceres)

    _OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    df.write_parquet(_OUT_PATH, compression="zstd")

    sha = hashlib.sha256(_OUT_PATH.read_bytes()).hexdigest()
    sidecar = {
        "source": _API_URL,
        "fields": _FIELDS,
        "downloaded_utc": datetime.now(UTC).isoformat(),
        "n_rows": n_total,
        "n_with_diameter": n_diam,
        "n_with_albedo": n_albedo,
        "sha256": sha,
        "signature": payload.get("signature"),
    }
    _OUT_PATH.with_suffix(".json").write_text(json.dumps(sidecar, indent=2))
    logger.info("Escrito %s (+ sidecar)", _OUT_PATH)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
