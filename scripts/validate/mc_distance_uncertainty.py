"""MC de la incertidumbre de distancia por covarianza de elementos de entrada (B4).

El tribunal (B4) objeta que ninguna distancia del catálogo lleva σ, y que el término
dominante del presupuesto de completitud es la **incertidumbre de los elementos
orbitales de entrada** (no propagada). Este script la cuantifica para los pares de
las tablas de eventos extremos: baja la covarianza orbital de cada cuerpo de JPL
SBDB, muestrea órbitas del elipsoide de error, las propaga (2 cuerpos, el mismo
modelo del refinador Kepler del catálogo) a la época del encuentro, y reporta la
dispersión de la distancia mutua ``σ_dist``.

Como el generador y el catálogo comparten el propagador 2-cuerpos, ``σ_dist`` aísla
la contribución de la **incertidumbre de elementos** (no la del modelo de fuerzas),
que es justo el término que B4 pide acotar por objeto.

Fuente: SBDB API (``?cov=mat``) da la covarianza 6×6 en la base
``[e, q, tp, node, peri, i]`` en su época; se muestrea ahí y se propaga.

Uso
---
    docker compose run --rm pipeline python -m scripts.validate.mc_distance_uncertainty \\
        --pairs data/output/orbdet/extreme_pairs.json --n-mc 2000 --seed 42 \\
        --out data/output/orbdet/extreme_pairs_sigma.json
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import time
import urllib.request
from pathlib import Path

import numpy as np

from src.orbdet.kepler import KeplerElements, elements_to_state, mean_motion, propagate

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
logger = logging.getLogger(__name__)

_SBDB = "https://ssd-api.jpl.nasa.gov/sbdb.api?sstr={n}&cov=mat&full-prec=1"
_CACHE = Path("data/cache/sbdb_cov")
_COV_LABELS = ["e", "q", "tp", "node", "peri", "i"]  # orden que devuelve SBDB


def _fetch_cov(number: int) -> dict | None:
    """Baja (y cachea) la covarianza orbital de SBDB para *number*. None si falta."""
    _CACHE.mkdir(parents=True, exist_ok=True)
    cache = _CACHE / f"{number}.json"
    if cache.exists():
        d = json.loads(cache.read_text())
    else:
        url = _SBDB.format(n=number)
        for attempt in range(3):
            try:
                with urllib.request.urlopen(url, timeout=30) as r:
                    d = json.loads(r.read().decode())
                break
            except Exception as exc:  # noqa: BLE001
                if attempt == 2:
                    logger.warning("SBDB falló para %d: %s", number, exc)
                    return None
                time.sleep(2)
        cache.write_text(json.dumps(d))
    orb = d.get("orbit") or {}
    cov = orb.get("covariance")
    if not cov or "data" not in cov:
        return None
    return cov


def _cov_to_elements(vec: np.ndarray, epoch: float) -> KeplerElements | None:
    """Vector SBDB [e,q,tp,node,peri,i] (AU, día, grados) → KeplerElements en *epoch*."""
    e, q, tp, node, peri, i = vec
    if not (0.0 < e < 1.0) or q <= 0.0:
        return None
    a = q / (1.0 - e)
    if not (0.3 < a < 100.0):
        return None
    n = mean_motion(a)  # rad/día
    big_m = (n * (epoch - tp)) % (2.0 * math.pi)
    return KeplerElements(
        a=a,
        e=e,
        i=math.radians(i),
        Omega=math.radians(node),
        omega=math.radians(peri),
        M=big_m,
    )


def _sample_positions(cov: dict, enc_jd: float, n_mc: int, rng: np.random.Generator) -> np.ndarray:
    """N posiciones heliocéntricas (AU) del cuerpo en *enc_jd*, muestreadas de la covarianza."""
    epoch = float(cov["epoch"])
    labels = cov["labels"]
    # SBDB entrega los datos de la matriz como strings en notación científica.
    data = np.array([[float(x) for x in row] for row in cov["data"]], dtype=float)
    # Reordenar a _COV_LABELS por si SBDB cambia el orden.
    idx = [labels.index(name) for name in _COV_LABELS]
    C = data[np.ix_(idx, idx)]
    # cov["elements"] es una lista de dicts {name, value, ...}; SBDB usa 'om'/'w'
    # para node/peri en esa lista aunque la matriz los etiquete 'node'/'peri'.
    _ELEM_NAME = {"e": "e", "q": "q", "tp": "tp", "node": "om", "peri": "w", "i": "i"}
    elem_by_name = {e["name"]: float(e["value"]) for e in cov["elements"]}
    mean = np.array([elem_by_name[_ELEM_NAME[k]] for k in _COV_LABELS], dtype=float)
    # Muestreo multivariado; simetrizar C por seguridad numérica.
    C = 0.5 * (C + C.T)
    try:
        draws = rng.multivariate_normal(mean, C, size=n_mc)
    except np.linalg.LinAlgError:
        # covarianza no PSD: añadir jitter mínimo en la diagonal
        C = C + np.eye(6) * (np.trace(C) * 1e-12)
        draws = rng.multivariate_normal(mean, C, size=n_mc)
    pos = []
    for v in draws:
        el = _cov_to_elements(v, epoch)
        if el is None:
            continue
        el_enc = propagate(el, enc_jd - epoch)
        r, _v = elements_to_state(el_enc)
        pos.append(r)
    return np.asarray(pos, dtype=float)


def run(args: argparse.Namespace) -> int:
    pairs = json.loads(Path(args.pairs).read_text())
    rng = np.random.default_rng(args.seed)
    results = []
    for p in pairs:
        n1, n2, enc_jd, d_cat = p["n1"], p["n2"], float(p["jd_tdb"]), float(p["dist_au"])
        cov1, cov2 = _fetch_cov(n1), _fetch_cov(n2)
        if cov1 is None or cov2 is None:
            logger.warning("%d×%d: sin covarianza SBDB (%s)", n1, n2, p.get("grp"))
            results.append({**p, "sigma_dist_au": None, "note": "no_cov"})
            continue
        r1 = _sample_positions(cov1, enc_jd, args.n_mc, rng)
        r2 = _sample_positions(cov2, enc_jd, args.n_mc, rng)
        k = min(len(r1), len(r2))
        if k < 50:
            results.append({**p, "sigma_dist_au": None, "note": "too_few_valid"})
            continue
        dist = np.linalg.norm(r1[:k] - r2[:k], axis=1)  # AU, muestras independientes
        med = float(np.median(dist))
        sigma = float(dist.std(ddof=1))
        lo, hi = (float(x) for x in np.percentile(dist, [16, 84]))
        au_km = 1.495978707e8
        results.append(
            {
                **p,
                "sigma_dist_au": sigma,
                "median_dist_au": med,
                "p16_au": lo,
                "p84_au": hi,
                "sigma_dist_km": sigma * au_km,
                "d_cat_km": d_cat * au_km,
                "sigma_over_dcat": sigma / d_cat if d_cat > 0 else None,
                "n_valid": k,
            }
        )
        logger.info(
            "%d×%d [%s]: d_cat=%.3e AU  σ_dist=%.3e AU (%.0f km)  σ/d=%.1f",
            n1,
            n2,
            p.get("grp"),
            d_cat,
            sigma,
            sigma * au_km,
            (sigma / d_cat) if d_cat > 0 else float("nan"),
        )
    if args.out:
        Path(args.out).write_text(json.dumps(results, indent=1, default=float))
        logger.info("Escrito %s (%d pares)", args.out, len(results))
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--pairs", type=Path, required=True, help="JSON con [{n1,n2,jd_tdb,dist_au,grp}]"
    )
    ap.add_argument("--n-mc", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", type=Path, default=None)
    return run(ap.parse_args())


if __name__ == "__main__":
    raise SystemExit(main())
