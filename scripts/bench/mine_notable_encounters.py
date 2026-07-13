"""Mina el catálogo de encuentros caracterizado en busca de eventos notables.

Frente P1 del plan ``planning/PUBLISH_PUSH_PLAN.md``: recorrer el catálogo de
encuentros cercanos 3D (Gaia DR3) para aislar eventos únicos, publicables o no
catalogados en otro lado. El catálogo es de encuentros 3D reales (mínima
distancia física entre pares), no co-localizaciones aparentes en el plano del
cielo.

Categorías reportadas
---------------------
1. Encuentros grande-grande (ambos diámetros por encima de un umbral).
2. Extremos: mínima distancia absoluta, encuentros lentos (candidatos naturales
   a determinación de masa), y la combinación grande + lento + cercano.
3. Pares en la misma región dinámica (proximidad en a/e/i como proxy de familia),
   uniendo con un snapshot de elementos orbitales.
4. Perturbadores: separar lo ya cubierto por los 16 perturbadores estudiados de
   los cuerpos grandes fuera de esa lista (candidatos a masa nueva para F4).

Uso
---
Ejecutar dentro de Docker. Sólo ``./data`` y ``./scripts`` están montados en el
contenedor, así que el reporte se escribe a un directorio montado y luego se
mueve a ``docs/`` en el host::

    docker compose run --rm pipeline python -m scripts.bench.mine_notable_encounters \
        --out data/output/notable_encounters.md
    mv data/output/notable_encounters.md docs/notable_encounters.md

Notas
-----
La entrada es un parquet de ~5.7 GB con ~72 M filas; cada fila es un único evento
de máxima aproximación por par (ya deduplicado). Se usa ``polars.scan_parquet``
(lazy) con filtros empujados para no materializar el catálogo en memoria.
"""

from __future__ import annotations

import argparse
import logging
from dataclasses import dataclass
from pathlib import Path

import polars as pl

logger = logging.getLogger(__name__)

# --------------------------------------------------------------------------- #
# Constantes de dominio
# --------------------------------------------------------------------------- #

#: Los 16 perturbadores ya estudiados en la capa de masas del proyecto.
KNOWN_PERTURBERS: frozenset[int] = frozenset(
    {1, 2, 3, 4, 7, 10, 15, 16, 31, 52, 65, 87, 88, 107, 511, 704}
)

#: Ruta por defecto al catálogo de encuentros caracterizado (freeze b1fix del paper).
DEFAULT_CATALOG = "data/output/encounters_characterized_b1fix.parquet"

#: Snapshot de elementos orbitales (a, e, i) por número, para el proxy de familia.
DEFAULT_ELEMENTS = "data/raw/gaia_orbits.parquet"

#: Umbrales de diámetro (km) para la categoría grande-grande.
BIG_THRESHOLDS_KM: tuple[float, ...] = (50.0, 100.0)

#: Corte de velocidad relativa (km/s) para "encuentro lento".
SLOW_VREL_KM_S = 1.0

#: Corte de distancia (km) para "encuentro cercano" en la combinación de interés.
CLOSE_DIST_KM = 1.0e6

#: Cantidad de filas a mostrar por tabla.
TOP_N = 30

#: Columnas legibles reportadas en las tablas de eventos.
DISPLAY_COLS: tuple[str, ...] = (
    "number_1",
    "designation_1",
    "number_2",
    "designation_2",
    "date_utc",
    "dist_km",
    "rel_vel_km_s",
    "diameter_1_km",
    "diameter_2_km",
    "class_1",
    "class_2",
)


# --------------------------------------------------------------------------- #
# Utilidades de formato
# --------------------------------------------------------------------------- #


def _collect(lf: pl.LazyFrame) -> pl.DataFrame:
    """Colecta un ``LazyFrame`` con el motor de streaming (memoria acotada).

    El catálogo tiene ~72 M filas y el entorno dispone de poca RAM libre; el motor
    de streaming procesa ``sort``/``head`` (top-k) y joins sin materializar todo el
    frame en memoria.

    Parameters
    ----------
    lf : polars.LazyFrame
        Plan a ejecutar.

    Returns
    -------
    polars.DataFrame
        Resultado materializado.
    """
    return lf.collect(engine="streaming")


def _fmt(value: object) -> str:
    """Formatea un valor de celda para una tabla markdown, tolerando nulls.

    Parameters
    ----------
    value : object
        Valor de celda (posible ``None``).

    Returns
    -------
    str
        Representación legible; ``—`` para nulos.
    """
    if value is None:
        return "—"
    if isinstance(value, float):
        if value != value:  # NaN
            return "—"
        if abs(value) >= 1.0e5:
            return f"{value:.3e}"
        if abs(value) >= 100:
            return f"{value:.0f}"
        return f"{value:.3f}"
    return str(value)


def df_to_markdown(df: pl.DataFrame, columns: tuple[str, ...] = DISPLAY_COLS) -> str:
    """Convierte un ``DataFrame`` en una tabla markdown con columnas legibles.

    Parameters
    ----------
    df : polars.DataFrame
        Filas a renderizar.
    columns : tuple of str, optional
        Columnas a incluir, en orden.

    Returns
    -------
    str
        Tabla markdown. Cadena vacía con nota si el ``DataFrame`` no tiene filas.
    """
    cols = [c for c in columns if c in df.columns]
    if df.height == 0:
        return "_(sin filas para este corte)_\n"
    header = "| " + " | ".join(cols) + " |"
    sep = "| " + " | ".join("---" for _ in cols) + " |"
    lines = [header, sep]
    for row in df.select(cols).iter_rows(named=True):
        lines.append("| " + " | ".join(_fmt(row[c]) for c in cols) + " |")
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------- #
# Consultas
# --------------------------------------------------------------------------- #


def big_big_encounters(lf: pl.LazyFrame, threshold_km: float, top_n: int) -> pl.DataFrame:
    """Encuentros donde ambos cuerpos superan un umbral de diámetro.

    Parameters
    ----------
    lf : polars.LazyFrame
        Catálogo de encuentros.
    threshold_km : float
        Diámetro mínimo (km) de *ambos* cuerpos.
    top_n : int
        Número máximo de filas (ordenadas por distancia mínima ascendente).

    Returns
    -------
    polars.DataFrame
        Encuentros grande-grande ordenados por ``dist_km``.
    """
    return _collect(
        lf.filter(
            (pl.col("diameter_1_km") >= threshold_km) & (pl.col("diameter_2_km") >= threshold_km)
        )
        .sort("dist_km")
        .head(top_n)
    )


def closest_encounters(lf: pl.LazyFrame, top_n: int) -> pl.DataFrame:
    """Los encuentros de menor distancia absoluta del catálogo."""
    return _collect(lf.sort("dist_km").head(top_n))


def slowest_encounters(lf: pl.LazyFrame, top_n: int, min_diam_km: float = 5.0) -> pl.DataFrame:
    """Encuentros de menor velocidad relativa (candidatos a determinación de masa).

    Se aplica un piso de diámetro para descartar el ruido de cuerpos diminutos
    cuya masa es irrelevante como perturbador; la deflexión útil requiere que al
    menos uno de los dos tenga tamaño apreciable.

    Parameters
    ----------
    lf : polars.LazyFrame
        Catálogo de encuentros.
    top_n : int
        Número máximo de filas.
    min_diam_km : float, optional
        Diámetro mínimo (km) de al menos uno de los cuerpos.

    Returns
    -------
    polars.DataFrame
        Encuentros lentos ordenados por ``rel_vel_km_s`` ascendente.
    """
    return _collect(
        lf.filter(
            (pl.col("diameter_1_km") >= min_diam_km) | (pl.col("diameter_2_km") >= min_diam_km)
        )
        .sort("rel_vel_km_s")
        .head(top_n)
    )


def big_slow_close(
    lf: pl.LazyFrame,
    top_n: int,
    min_diam_km: float = 50.0,
    max_vrel_km_s: float = SLOW_VREL_KM_S,
    max_dist_km: float = CLOSE_DIST_KM,
) -> pl.DataFrame:
    """Combinación de máximo interés físico: grande + lento + cercano.

    Un cuerpo grande (perturbador masivo) que pasa lento y cerca de otro produce
    la máxima deflexión medible, el escenario ideal para determinación de masa.

    Parameters
    ----------
    lf : polars.LazyFrame
        Catálogo de encuentros.
    top_n : int
        Número máximo de filas.
    min_diam_km : float, optional
        Diámetro mínimo (km) de al menos un cuerpo (el perturbador).
    max_vrel_km_s : float, optional
        Velocidad relativa máxima (km/s).
    max_dist_km : float, optional
        Distancia máxima (km).

    Returns
    -------
    polars.DataFrame
        Eventos filtrados, ordenados por distancia mínima.
    """
    return _collect(
        lf.filter(
            ((pl.col("diameter_1_km") >= min_diam_km) | (pl.col("diameter_2_km") >= min_diam_km))
            & (pl.col("rel_vel_km_s") <= max_vrel_km_s)
            & (pl.col("dist_km") <= max_dist_km)
        )
        .sort("dist_km")
        .head(top_n)
    )


def same_dynamical_region(
    lf: pl.LazyFrame,
    elements_path: str,
    top_n: int,
    da_frac: float = 0.01,
    de: float = 0.02,
    di_deg: float = 1.0,
    max_dist_km: float = 5.0e6,
) -> pl.DataFrame:
    """Pares próximos en espacio orbital (a, e, i) como proxy de familia dinámica.

    No se dispone de asignación de familia; se une el catálogo con un snapshot de
    elementos orbitales por número y se marcan los pares cuyos elementos son
    mutuamente cercanos. Es un proxy grueso (no un clasificador de familias en el
    espacio de elementos propios), útil para señalar candidatos.

    Parameters
    ----------
    lf : polars.LazyFrame
        Catálogo de encuentros.
    elements_path : str
        Parquet de elementos orbitales con columnas ``number, a_au, e, i_deg``.
    top_n : int
        Número máximo de filas.
    da_frac : float, optional
        Diferencia fraccional máxima en semieje mayor.
    de : float, optional
        Diferencia absoluta máxima en excentricidad.
    di_deg : float, optional
        Diferencia absoluta máxima en inclinación (grados).
    max_dist_km : float, optional
        Prefiltro de distancia para acotar el volumen antes del join.

    Returns
    -------
    polars.DataFrame
        Pares co-familia candidatos, ordenados por distancia mínima. Incluye las
        columnas de elementos de ambos cuerpos.
    """
    elem = pl.scan_parquet(elements_path).select(
        pl.col("number").alias("num"),
        pl.col("a_au"),
        pl.col("e"),
        pl.col("i_deg"),
    )

    joined = (
        lf.filter(pl.col("dist_km") <= max_dist_km)
        .join(
            elem.rename({"num": "number_1", "a_au": "a1", "e": "e1", "i_deg": "i1"}),
            on="number_1",
            how="inner",
        )
        .join(
            elem.rename({"num": "number_2", "a_au": "a2", "e": "e2", "i_deg": "i2"}),
            on="number_2",
            how="inner",
        )
    )

    da = (pl.col("a1") - pl.col("a2")).abs() / ((pl.col("a1") + pl.col("a2")) / 2.0)
    return _collect(
        joined.filter(
            (da <= da_frac)
            & ((pl.col("e1") - pl.col("e2")).abs() <= de)
            & ((pl.col("i1") - pl.col("i2")).abs() <= di_deg)
        )
        .with_columns(
            pl.col("a1").round(4),
            pl.col("a2").round(4),
            pl.col("e1").round(4),
            pl.col("e2").round(4),
            pl.col("i1").round(3),
            pl.col("i2").round(3),
        )
        .sort("dist_km")
        .head(top_n)
    )


def known_perturber_encounters(lf: pl.LazyFrame, top_n: int) -> pl.DataFrame:
    """Encuentros donde alguno de los cuerpos es uno de los 16 perturbadores."""
    return _collect(
        lf.filter(
            pl.col("number_1").is_in(list(KNOWN_PERTURBERS))
            | pl.col("number_2").is_in(list(KNOWN_PERTURBERS))
        )
        .sort("dist_km")
        .head(top_n)
    )


def new_big_perturber_candidates(
    lf: pl.LazyFrame,
    top_n: int,
    min_diam_km: float = 100.0,
    max_vrel_km_s: float = 3.0,
    max_dist_km: float = 3.0e6,
) -> pl.DataFrame:
    """Cuerpos grandes fuera de los 16 con encuentros cercanos con objetivos pequeños.

    Alimenta F4 (masa nueva). Se cuentan, por cuerpo grande no catalogado, los
    encuentros cercanos y lentos con otros cuerpos: mayor número de buenos
    eventos ⇒ mejor candidato a determinación de masa.

    Parameters
    ----------
    lf : polars.LazyFrame
        Catálogo de encuentros.
    top_n : int
        Número máximo de candidatos.
    min_diam_km : float, optional
        Diámetro mínimo (km) para considerar "grande".
    max_vrel_km_s : float, optional
        Velocidad relativa máxima (km/s) para un evento "útil".
    max_dist_km : float, optional
        Distancia máxima (km) para un evento "útil".

    Returns
    -------
    polars.DataFrame
        Un renglón por cuerpo grande no catalogado, con el conteo de eventos
        útiles, la mejor (mínima) distancia y velocidad, ordenado por conteo.
    """
    known = list(KNOWN_PERTURBERS)
    useful = (pl.col("rel_vel_km_s") <= max_vrel_km_s) & (pl.col("dist_km") <= max_dist_km)
    # Normaliza cada evento a la perspectiva del cuerpo grande no catalogado.
    big1 = lf.filter(
        (pl.col("diameter_1_km") >= min_diam_km) & ~pl.col("number_1").is_in(known) & useful
    ).select(
        pl.col("number_1").alias("big_number"),
        pl.col("designation_1").alias("big_name"),
        pl.col("diameter_1_km").alias("big_diam_km"),
        pl.col("class_1").alias("big_class"),
        pl.col("dist_km"),
        pl.col("rel_vel_km_s"),
    )
    big2 = lf.filter(
        (pl.col("diameter_2_km") >= min_diam_km) & ~pl.col("number_2").is_in(known) & useful
    ).select(
        pl.col("number_2").alias("big_number"),
        pl.col("designation_2").alias("big_name"),
        pl.col("diameter_2_km").alias("big_diam_km"),
        pl.col("class_2").alias("big_class"),
        pl.col("dist_km"),
        pl.col("rel_vel_km_s"),
    )
    return _collect(
        pl.concat([big1, big2])
        .group_by("big_number")
        .agg(
            pl.col("big_name").first(),
            pl.col("big_diam_km").first().round(1),
            pl.col("big_class").first(),
            pl.len().alias("n_useful_events"),
            pl.col("dist_km").min().alias("best_dist_km"),
            pl.col("rel_vel_km_s").min().alias("min_vrel_km_s"),
        )
        .sort(["n_useful_events", "best_dist_km"], descending=[True, False])
        .head(top_n)
    )


# --------------------------------------------------------------------------- #
# Orquestación
# --------------------------------------------------------------------------- #


@dataclass
class MiningResults:
    """Contenedor de todos los ``DataFrame`` producidos por la minería."""

    big_big: dict[float, pl.DataFrame]
    closest: pl.DataFrame
    slowest: pl.DataFrame
    big_slow_close: pl.DataFrame
    family_proxy: pl.DataFrame | None
    known_perturbers: pl.DataFrame
    new_big_candidates: pl.DataFrame
    total_rows: int


def run_mining(
    catalog_path: str,
    elements_path: str,
    top_n: int = TOP_N,
) -> MiningResults:
    """Ejecuta todas las consultas de minería sobre el catálogo.

    Parameters
    ----------
    catalog_path : str
        Ruta al parquet de encuentros caracterizado.
    elements_path : str
        Ruta al parquet de elementos orbitales (para el proxy de familia).
    top_n : int, optional
        Filas máximas por tabla.

    Returns
    -------
    MiningResults
        Todos los resultados.
    """
    lf = pl.scan_parquet(catalog_path)
    total_rows = lf.select(pl.len()).collect().item()
    logger.info("Catálogo cargado: %s filas (lazy)", f"{total_rows:,}")

    big_big: dict[float, pl.DataFrame] = {}
    for thr in BIG_THRESHOLDS_KM:
        big_big[thr] = big_big_encounters(lf, thr, top_n)
        logger.info("Grande-grande D>=%g km: %d filas", thr, big_big[thr].height)

    closest = closest_encounters(lf, top_n)
    logger.info("Más cercanos: min dist_km = %.1f", closest["dist_km"].min())

    slowest = slowest_encounters(lf, top_n)
    logger.info("Más lentos: min v_rel = %.4f km/s", slowest["rel_vel_km_s"].min())

    bsc = big_slow_close(lf, top_n)
    logger.info("Grande+lento+cercano: %d filas", bsc.height)

    family_proxy: pl.DataFrame | None
    try:
        family_proxy = same_dynamical_region(lf, elements_path, top_n)
        logger.info("Proxy de familia: %d filas", family_proxy.height)
    except Exception as exc:  # noqa: BLE001 -- no bloquear el reporte
        logger.warning("Proxy de familia omitido (queda para P1b): %s", exc)
        family_proxy = None

    known = known_perturber_encounters(lf, top_n)
    logger.info("Encuentros de perturbadores conocidos: %d filas", known.height)

    new_big = new_big_perturber_candidates(lf, top_n)
    logger.info("Candidatos grandes fuera de los 16: %d filas", new_big.height)

    return MiningResults(
        big_big=big_big,
        closest=closest,
        slowest=slowest,
        big_slow_close=bsc,
        family_proxy=family_proxy,
        known_perturbers=known,
        new_big_candidates=new_big,
        total_rows=total_rows,
    )


def _in_known(df: pl.DataFrame) -> pl.Series:
    """Máscara de filas cuyo par toca uno de los 16 perturbadores."""
    known = list(KNOWN_PERTURBERS)
    return df["number_1"].is_in(known) | df["number_2"].is_in(known)


def build_report(results: MiningResults, catalog_path: str) -> str:
    """Construye el reporte markdown a partir de los resultados de minería.

    Parameters
    ----------
    results : MiningResults
        Salida de :func:`run_mining`.
    catalog_path : str
        Ruta del catálogo (para la cabecera de procedencia).

    Returns
    -------
    str
        Documento markdown completo.
    """
    today = "2026-07-02"
    parts: list[str] = []

    parts.append("# Encuentros notables en el catálogo Gaia DR3")
    parts.append("")
    parts.append("> **Estado:** 🟡 EN CURSO (frente P1 de `planning/PUBLISH_PUSH_PLAN.md`)")
    parts.append(f"> **Fecha:** {today}")
    parts.append(
        f"> **Fuente:** `{catalog_path}` — {results.total_rows:,} encuentros 3D "
        "reales (una fila por par, ya deduplicado por máxima aproximación)."
    )
    parts.append(">")
    parts.append(
        "> Este catálogo registra la **mínima distancia física en 3D** entre pares "
        "de asteroides durante la ventana Gaia DR3 (jul 2014 – may 2017), no "
        "co-localizaciones aparentes en el plano del cielo."
    )
    parts.append("")

    parts.append("## Metodología y limitaciones")
    parts.append("")
    parts.append(
        "- Los diámetros derivan de `H` con albedo por clase cuando no hay medida "
        "directa; deben leerse como estimaciones de orden de magnitud, no como "
        "diámetros medidos. Los cortes por tamaño son por tanto aproximados."
    )
    parts.append(
        "- La propagación de base del catálogo es Kepler de dos cuerpos; el "
        "refinamiento N-cuerpos se aplica sólo al subset de determinación de masas. "
        "Las distancias mínimas de esta minería pueden tener sesgo cerca de "
        "resonancias o encuentros planetarios. Cualquier evento seleccionado como "
        "candidato requiere revalidación N-cuerpos antes de publicar."
    )
    parts.append(
        "- El presupuesto de completitud del catálogo (censura ~0.70 %, recall "
        "prefiltro ~76 % en el tail adverso) implica que faltan algunos encuentros "
        "genuinos; las tablas de abajo son un piso, no un censo exhaustivo."
    )
    parts.append(
        "- El proxy de familia (§4) es proximidad en elementos osculantes "
        "`(a, e, i)`, **no** clasificación en elementos propios. Señala candidatos, "
        "no confirma pertenencia a familia."
    )
    parts.append("")

    # 1. Grande-grande
    parts.append("## 1. Encuentros grande-grande")
    parts.append("")
    parts.append(
        "Ambos cuerpos por encima del umbral de diámetro. Son los más raros: dos "
        "cuerpos masivos que se aproximan en 3D. Ordenados por distancia mínima."
    )
    for thr in BIG_THRESHOLDS_KM:
        df = results.big_big[thr]
        parts.append("")
        parts.append(
            f"### 1.{int(thr == 100.0) + 1} Ambos D ≳ {thr:g} km "
            f"({df.height} en el top {TOP_N})"
        )
        parts.append("")
        parts.append(df_to_markdown(df))
    parts.append("")

    # 2. Extremos
    parts.append("## 2. Encuentros extremos")
    parts.append("")
    parts.append("### 2a. Mínima distancia absoluta del catálogo")
    parts.append("")
    parts.append(df_to_markdown(results.closest))
    parts.append("")
    parts.append("### 2b. Encuentros más lentos (candidatos naturales a masa)")
    parts.append("")
    parts.append(
        "Velocidad relativa mínima (con al menos un cuerpo D ≳ 5 km). Una v_rel baja "
        "prolonga la interacción gravitatoria ⇒ deflexión mayor y más medible."
    )
    parts.append("")
    parts.append(df_to_markdown(results.slowest))
    parts.append("")
    parts.append("### 2c. Grande + lento + cercano (máximo interés físico)")
    parts.append("")
    parts.append(
        f"Al menos un cuerpo D ≳ 50 km, v_rel ≤ {SLOW_VREL_KM_S:g} km/s, "
        f"dist ≤ {CLOSE_DIST_KM:.0e} km."
    )
    parts.append("")
    parts.append(df_to_markdown(results.big_slow_close))
    parts.append("")

    # 3. Familia
    parts.append("## 3. Pares en la misma región dinámica (proxy de familia)")
    parts.append("")
    if results.family_proxy is None:
        parts.append(
            "_Omitido en esta corrida (queda para P1b): no se pudo unir con una "
            "fuente de elementos orbitales._"
        )
    else:
        parts.append(
            "Pares cuyos elementos osculantes `(a, e, i)` son mutuamente cercanos "
            "(Δa/a ≤ 1 %, Δe ≤ 0.02, Δi ≤ 1°) y que además tuvieron un encuentro "
            "físico cercano. Se añaden las columnas de elementos de ambos cuerpos."
        )
        parts.append("")
        fam_cols = DISPLAY_COLS + ("a1", "a2", "e1", "e2", "i1", "i2")
        parts.append(df_to_markdown(results.family_proxy, fam_cols))
    parts.append("")

    # 4. Perturbadores
    parts.append("## 4. Perturbadores")
    parts.append("")
    parts.append("### 4a. Encuentros que tocan uno de los 16 perturbadores estudiados")
    parts.append("")
    parts.append(
        "Lista de referencia (ya cubierta): "
        f"{', '.join(str(n) for n in sorted(KNOWN_PERTURBERS))}. "
        "Sirve para **separar lo ya trabajado** del descubrimiento."
    )
    parts.append("")
    parts.append(df_to_markdown(results.known_perturbers))
    parts.append("")
    parts.append("### 4b. Cuerpos grandes FUERA de los 16 (candidatos a masa nueva → F4)")
    parts.append("")
    parts.append(
        "Cuerpos con D ≳ 100 km, no incluidos en los 16, rankeados por número de "
        "encuentros *útiles* (v_rel ≤ 3 km/s y dist ≤ 3×10⁶ km). Más eventos buenos "
        "⇒ mejor candidato a determinación de masa."
    )
    parts.append("")
    ncols = (
        "big_number",
        "big_name",
        "big_diam_km",
        "big_class",
        "n_useful_events",
        "best_dist_km",
        "min_vrel_km_s",
    )
    parts.append(df_to_markdown(results.new_big_candidates, ncols))
    parts.append("")

    # 5. Candidatos a seguimiento
    parts.append("## 5. Candidatos a seguimiento")
    parts.append("")

    # (a) eventos notables / no catalogados
    bb100 = results.big_big[100.0]
    bb50_new = results.big_big[50.0].filter(~_in_known(results.big_big[50.0]))
    parts.append("### 5a. Eventos genuinamente notables / potencialmente no catalogados")
    parts.append("")
    parts.append(
        f"- **{bb100.height} encuentros grande-grande D ≳ 100 km** en el catálogo "
        "(sección 1.2). Cualquier par de este grupo que **no** toque a los 16 "
        "perturbadores conocidos es un evento de alto perfil sin cobertura previa "
        "obvia; ver la marca de perturbador en §4a."
    )
    if bb50_new.height:
        top = bb50_new.row(0, named=True)
        parts.append(
            f"- El encuentro grande-grande (D ≳ 50 km) más cercano **fuera de los "
            f"16** es {top['designation_1']} × {top['designation_2']} el "
            f"{top['date_utc']} a {_fmt(top['dist_km'])} km — candidato a revisión "
            f"N-cuerpos."
        )
    parts.append(
        "- Los encuentros más cercanos en términos absolutos (§2a) merecen "
        "revalidación N-cuerpos: a esas distancias la aproximación Kepler es más "
        "frágil, pero si sobreviven son los eventos geométricamente más notables "
        "del dataset."
    )
    parts.append("")

    # (b) mejores candidatos a masa
    parts.append("### 5b. Cuerpos grandes fuera de los 16 con más encuentros útiles (→ F4)")
    parts.append("")
    parts.append(
        "Ranking del §4b (top 10). Estos son los perturbadores no estudiados con más "
        "eventos de baja v_rel y corta distancia, es decir, con el mayor potencial "
        "de deflexión medible para una masa nueva:"
    )
    parts.append("")
    top10 = results.new_big_candidates.head(10)
    for row in top10.iter_rows(named=True):
        parts.append(
            f"- **{row['big_name']}** "
            f"(D≈{_fmt(row['big_diam_km'])} km, clase {_fmt(row['big_class'])}): "
            f"{row['n_useful_events']} eventos útiles; mejor par a "
            f"{_fmt(row['best_dist_km'])} km, v_rel mín {_fmt(row['min_vrel_km_s'])} km/s."
        )
    parts.append("")
    parts.append(
        "> Contrastar contra los candidatos F4 propuestos en "
        "`docs/mass_layer_f4_design.md` (24 Themis, 532 Herculina, 29 Amphitrite, "
        "354 Eleonora) y priorizar los que además aparezcan alto en este ranking."
    )
    parts.append("")

    parts.append("---")
    parts.append("")
    parts.append(f"_Generado por `scripts/bench/mine_notable_encounters.py` el {today}._")
    parts.append("")
    return "\n".join(parts)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parsea argumentos de línea de comandos."""
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--catalog", default=DEFAULT_CATALOG, help="parquet de encuentros")
    p.add_argument("--elements", default=DEFAULT_ELEMENTS, help="parquet de elementos orbitales")
    p.add_argument("--out", default="docs/notable_encounters.md", help="reporte markdown de salida")
    p.add_argument("--top-n", type=int, default=TOP_N, help="filas por tabla")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    """Punto de entrada CLI."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    args = parse_args(argv)
    logger.info("Minando %s", args.catalog)
    results = run_mining(args.catalog, args.elements, top_n=args.top_n)
    report = build_report(results, args.catalog)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(report, encoding="utf-8")
    logger.info("Reporte escrito en %s (%d bytes)", out_path, len(report))


if __name__ == "__main__":
    main()
