"""Ranking de perturbadores fuera de los 16 de la efeméride para el gate F4.

Cuenta, por cada número MPC que **no** está en los 16 grandes de ``sb441-n16.bsp``,
cuántos objetivos pequeños distintos tuvieron un encuentro < ``--max-dist-au`` con
él en el catálogo congelado (``encounters_catalog_hybrid_stageb.parquet``). El
conteo replica exactamente la selección de objetivos de
:func:`scripts.mass.orbdet_fit_realdata._read_targets_from_catalog`:

- un encuentro cuenta para el número que actúa como perturbador si el **otro**
  cuerpo (el objetivo) tiene número < ``--max-target-number``,
- se deduplica por objetivo (un objetivo cuenta una sola vez por perturbador),

de modo que el ``N`` reportado es el número máximo de objetivos que un fit F4
podría usar para ese candidato. Se listan los de mayor ``N`` — buenos candidatos
para el gate (≥ 1 perturbador fuera de los 16 ajustado con χ²_red ∈ [0.95, 1.05]).

Uso
---
    docker compose run --rm pipeline python -m scripts.mass.find_extra_perturber_candidates \\
        --catalog data/output/encounters_catalog_hybrid_stageb.parquet --top 30
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import polars as pl

from src.orbdet.dynamics_assist import BIG_ASTEROIDS

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
logger = logging.getLogger(__name__)

# Números MPC de los 16 perturbadores de la efeméride (a excluir de los candidatos).
_EPHEM_NUMBERS: frozenset[int] = frozenset(
    {1, 2, 3, 4, 7, 10, 15, 16, 31, 52, 65, 87, 88, 107, 511, 704}
)

# Nombres conocidos para anotar el ranking (puramente cosmético).
_KNOWN_NAMES: dict[int, str] = {
    13: "Egeria",
    19: "Fortuna",
    24: "Themis",
    29: "Amphitrite",
    45: "Eugenia",
    48: "Doris",
    354: "Eleonora",
    532: "Herculina",
}


def rank_candidates(
    catalog_path: Path,
    *,
    max_dist_au: float = 0.05,
    max_target_number: int = 100_000,
    max_perturber_number: int | None = None,
    top: int = 30,
) -> pl.DataFrame:
    """Devuelve los perturbadores no-16 con más objetivos < *max_dist_au*.

    Parameters
    ----------
    catalog_path:
        Ruta al catálogo de encuentros (parquet).
    max_dist_au:
        Distancia máxima de encuentro considerada (igual que el fit).
    max_target_number:
        Número MPC máximo de un objetivo para que cuente (igual que el fit).
    max_perturber_number:
        Si se da, sólo considera perturbadores candidatos con número ≤ este valor.
        Los cuerpos *físicamente* masivos (los que determinan masas por encuentros)
        son de número bajo; sin este corte el ranking crudo lo dominan asteroides
        pequeños de número alto que acumulan muchos encuentros de azar pero con masa
        despreciable. Un corte ~1000 aísla los candidatos grandes clásicos.
    top:
        Número de candidatos a devolver.

    Returns
    -------
    polars.DataFrame
        Columnas ``perturber`` (Int), ``name`` (str), ``n_targets`` (Int),
        ordenadas por ``n_targets`` descendente.
    """
    lf = pl.scan_parquet(catalog_path).select(["number_1", "number_2", "dist_au"])
    lf = lf.filter(pl.col("dist_au") < max_dist_au)

    # Dos vistas: cada cuerpo del par puede ser el perturbador; el otro es el objetivo.
    # (mismo criterio que _read_targets_from_catalog, aplicado a todos los números)
    a = lf.select(
        pl.col("number_1").alias("perturber"), pl.col("number_2").alias("target")
    )
    b = lf.select(
        pl.col("number_2").alias("perturber"), pl.col("number_1").alias("target")
    )
    pairs = pl.concat([a, b])

    ephem = list(_EPHEM_NUMBERS)
    pairs = pairs.filter(
        (~pl.col("perturber").is_in(ephem)) & (pl.col("target") < max_target_number)
    )
    if max_perturber_number is not None:
        pairs = pairs.filter(pl.col("perturber") <= max_perturber_number)
    # Deduplicar por (perturbador, objetivo): un objetivo cuenta una vez por perturbador.
    counts = (
        pairs.unique(subset=["perturber", "target"])
        .group_by("perturber")
        .agg(pl.len().alias("n_targets"))
        .sort("n_targets", descending=True)
        .head(top)
        .collect()
    )
    counts = counts.with_columns(
        pl.col("perturber")
        .map_elements(lambda n: _KNOWN_NAMES.get(int(n), ""), return_dtype=pl.String)
        .alias("name")
    ).select(["perturber", "name", "n_targets"])
    return counts


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument(
        "--catalog",
        type=Path,
        default=Path("data/output/encounters_catalog_hybrid_stageb.parquet"),
    )
    parser.add_argument("--max-dist-au", type=float, default=0.05)
    parser.add_argument("--max-target-number", type=int, default=100_000)
    parser.add_argument(
        "--max-perturber-number",
        type=int,
        default=None,
        help="sólo candidatos con número ≤ este valor; los cuerpos masivos (que "
        "determinan masas) son de número bajo. Recomendado ~1000 para aislar los "
        "candidatos grandes clásicos (Themis, Herculina, Amphitrite, …)",
    )
    parser.add_argument("--top", type=int, default=30)
    args = parser.parse_args()

    logger.info(
        "Excluyendo los 16 de la efeméride (%s)", ", ".join(sorted(BIG_ASTEROIDS))
    )
    df = rank_candidates(
        args.catalog,
        max_dist_au=args.max_dist_au,
        max_target_number=args.max_target_number,
        max_perturber_number=args.max_perturber_number,
        top=args.top,
    )
    print(f"\n=== candidatos F4 (perturbadores fuera de los 16), dist < {args.max_dist_au} AU ===")
    print(f"{'perturber':>10}  {'name':<14}  {'n_targets':>9}")
    for row in df.iter_rows(named=True):
        print(f"{row['perturber']:>10}  {row['name']:<14}  {row['n_targets']:>9}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
