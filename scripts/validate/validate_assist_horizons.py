"""Valida el backend ASSIST de orbdet contra JPL Horizons (gate de exactitud T2/T8).

Propaga un asteroide numerado desde sus elementos osculadores de Horizons en una
época dentro de la ventana Gaia y compara la trayectoria baricéntrica eclíptica
contra los vectores de Horizons, para **ambos** backends:

- ``rebound``: planetas integrados libremente (modelo previo, sin GR).
- ``assist``: efeméride DE440 + 16 perturbadores + relatividad EIH (state-of-the-art).

Reporta el residuo máximo en AU y su equivalente angular (mas) a la distancia del
objeto, que es la métrica relevante para astrometría de Gaia (ruido AL ~0.5–2 mas).

Uso:
    docker run --rm -v "$PWD":/app -w /app gaia-asteroid-encounters \
        python3 -m scripts.validate.validate_assist_horizons --target 8 --span-days 900
"""

from __future__ import annotations

import argparse
import math

import numpy as np

from src.orbdet.dynamics import propagate as propagate_rebound
from src.orbdet.dynamics_assist import big_asteroid_perturbers, propagate_assist
from src.orbdet.kepler import KeplerElements

_FULL_PLANETS = (
    "sun", "mercury", "venus", "earth", "mars", "jupiter", "saturn", "uranus", "neptune",
)
_MAS_PER_RAD = math.degrees(1.0) * 3_600_000.0


def _horizons_elements(target: str, epoch_jd_tdb: float) -> KeplerElements:
    from astroquery.jplhorizons import Horizons

    tab = Horizons(
        id=target, id_type="smallbody", location="@sun", epochs=epoch_jd_tdb
    ).elements(refplane="ecliptic")
    return KeplerElements(
        a=float(tab["a"][0]),
        e=float(tab["e"][0]),
        i=math.radians(float(tab["incl"][0])),
        Omega=math.radians(float(tab["Omega"][0])),
        omega=math.radians(float(tab["w"][0])),
        M=math.radians(float(tab["M"][0])),
    )


def _horizons_vectors(target: str, epochs: np.ndarray) -> np.ndarray:
    from astroquery.jplhorizons import Horizons

    vec = Horizons(
        id=target, id_type="smallbody", location="@0", epochs=list(epochs)
    ).vectors(refplane="ecliptic")
    return np.column_stack([vec["x"], vec["y"], vec["z"]]).astype(float)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", default="8", help="número/designación MPC (no perturber)")
    ap.add_argument("--epoch", type=float, default=2_457_000.5, help="JD TDB de la época")
    ap.add_argument("--span-days", type=float, default=900.0)
    ap.add_argument("--n-samples", type=int, default=10)
    args = ap.parse_args()

    el = _horizons_elements(args.target, args.epoch)
    out_epochs = args.epoch + np.linspace(30.0, args.span_days, args.n_samples)
    ref = _horizons_vectors(args.target, out_epochs)
    helio_au = float(np.mean(np.linalg.norm(ref, axis=1)))

    # rebound (planetas libres, sin GR)
    got_reb = propagate_rebound(
        el, args.epoch, out_epochs, perturbers=_FULL_PLANETS, integrator="ias15"
    )
    res_reb = np.linalg.norm(got_reb - ref, axis=1)

    # assist (DE440 + 16 perturbadores + GR)
    bg = big_asteroid_perturbers(args.epoch, exclude=())
    got_ass = propagate_assist(el, args.epoch, out_epochs, asteroid_perturbers=bg, gr=True)
    res_ass = np.linalg.norm(got_ass - ref, axis=1)

    def mas(au: float) -> float:
        return au / helio_au * _MAS_PER_RAD

    print(f"target=({args.target}) helio≈{helio_au:.3f} AU  span={args.span_days:.0f} d")
    print(f"  rebound (libres, sin GR): max {res_reb.max():.3e} AU = {mas(res_reb.max()):8.2f} mas")
    print(f"  assist  (DE440+GR+16):    max {res_ass.max():.3e} AU = {mas(res_ass.max()):8.2f} mas")
    print(f"  mejora: {res_reb.max() / max(res_ass.max(), 1e-30):.1f}×")
    print("  residuos assist por época (mas):", np.array2string(
        np.array([mas(r) for r in res_ass]), precision=3, suppress_small=True))


if __name__ == "__main__":
    main()
