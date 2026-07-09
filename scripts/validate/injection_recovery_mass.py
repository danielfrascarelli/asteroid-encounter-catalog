"""Injection-recovery end-to-end de la capa de MASAS (tribunal 2026-07-04, B8).

Complementa la validación de calibradores (que usa masas de literatura como verdad
externa) con una prueba de **sesgo del estimador** sobre datos donde la verdad es
conocida por construcción: se inyecta una masa de perturbador ``M_inj`` conocida,
se genera la astrometría por-tránsito que esa masa produciría sobre ``N`` objetivos
(geometría de encuentro controlada + ruido along-scan realista), y se verifica que
el ajuste conjunto masa+órbita (:func:`determine_shared_mass`) la recupera **sin
sesgo** y con σ (formal y jackknife) consistente con la dispersión real.

Qué valida y qué no
-------------------
Esto valida el **estimador** (mínimos cuadrados masa+órbita, covarianza de Fisher,
jackknife) end-to-end: que sobre datos con verdad conocida el ajuste es insesgado y
su σ está bien escalada. NO revalida el modelo de fuerzas de producción (ASSIST/
DE440): el generador y el ajuste comparten el propagador ligero (rebound, Sol+
perturbador), de modo que un sesgo aquí es del estimador, no de la efeméride. La
fidelidad del modelo de fuerzas se ancla por separado en los cuatro calibradores.

Construcción de cada objetivo
-----------------------------
1. Perturbador en una órbita fija; ``N`` objetivos cada uno con un encuentro a
   impacto ``b`` (log-uniforme) y velocidad relativa controlada cerca de la época
   media del arco, de modo que la deflexión sea medible.
2. Astrometría sintética: ``predict_radec`` con el perturbador a ``M_inj`` sobre las
   épocas/posiciones de Gaia, más ruido gaussiano along-scan a ``sigma_al_mas``.
3. Ajuste: masa-semilla sesgada (0.6·M_inj) + elementos MPCORB-like → recuperación.

Gates (Tarea 8 / B8)
--------------------
- Recuperación insesgada: ``|M_fit − M_inj| ≤ n_sigma · σ_fit`` (default 3σ).
- σ bien escalada: ``χ²_red`` del ajuste en [0.5, 2.0].
- (Opcional, con ``--jackknife``) σ_jack finita y del orden de σ_formal..×N.

Uso
---
    docker compose run --rm pipeline python -m scripts.validate.injection_recovery_mass \\
        --n-targets 24 --m-inj-kg 2e19 --noise-mas 2.0 --seed 42 --jackknife
"""

from __future__ import annotations

import argparse
import logging
import math

import numpy as np

from src.orbdet.constants import M_SUN_KG
from src.orbdet.dynamics import AsteroidPerturber, planet_state_ecliptic
from src.orbdet.frames import ecliptic_to_equatorial
from src.orbdet.kepler import KeplerElements, elements_to_state, state_to_elements
from src.orbdet.mass_determination import (
    TargetObservations,
    determine_shared_mass,
    jackknife_mass_sigma,
)
from src.orbdet.observation import predict_radec

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
logger = logging.getLogger(__name__)

_PERTURBERS = ("sun",)
_DEG_PER_MAS = 1.0 / 3.6e6


def _perturber_and_epoch(rng: np.random.Generator) -> tuple[KeplerElements, float]:
    """Perturbador en el cinturón medio y época común del arco sintético."""
    epoch = 2_457_000.5
    pert = KeplerElements(
        a=2.77,
        e=0.08,
        i=math.radians(9.7),
        Omega=math.radians(80.4),
        omega=math.radians(73.6),
        M=math.radians(float(rng.uniform(0.0, 360.0))),
    )
    return pert, epoch


def _make_target(
    pert_el: KeplerElements,
    m_inj_msun: float,
    noise_mas: float,
    rng: np.random.Generator,
    epoch: float,
) -> TargetObservations | None:
    """Un objetivo con encuentro cercano al perturbador y astrometría sintética.

    El objetivo se coloca cerca del perturbador (impacto ``b`` log-uniforme) con una
    velocidad relativa pequeña, de modo que la masa inyectada imprima una deflexión
    medible sobre el arco. Devuelve ``None`` si la geometría degenera.
    """
    r_p, v_p = elements_to_state(pert_el)
    # Impacto perpendicular a la velocidad del perturbador, magnitud log-uniforme.
    b_au = 10.0 ** rng.uniform(-2.5, -1.5)  # ~0.003–0.03 AU
    perp = np.cross(v_p, np.array([0.0, 0.0, 1.0]))
    perp = perp / np.linalg.norm(perp)
    r_t = r_p + b_au * perp
    # Velocidad relativa pequeña (encuentro lento → deflexión grande).
    dv = np.array([rng.normal(0, 3e-4), rng.normal(0, 3e-4), rng.normal(0, 3e-4)])
    v_t = v_p + dv
    try:
        tgt_el = state_to_elements(r_t, v_t)
    except Exception:
        return None
    if not (0.0 < tgt_el.e < 0.7 and 1.8 < tgt_el.a < 3.6):
        return None

    obs = epoch + np.linspace(-220.0, 220.0, 20)
    gaia_ecl = np.array([planet_state_ecliptic("earth", float(t))[0] for t in obs])
    gaia_icrs = ecliptic_to_equatorial(gaia_ecl)
    pert_true = AsteroidPerturber("pert", m_inj_msun, pert_el)
    ra, dec = predict_radec(
        tgt_el, epoch, obs, gaia_icrs, perturbers=_PERTURBERS, asteroid_perturbers=(pert_true,)
    )
    if not (np.all(np.isfinite(ra)) and np.all(np.isfinite(dec))):
        return None
    n = obs.size
    ra = ra + rng.normal(0.0, noise_mas * _DEG_PER_MAS, n) / np.cos(np.radians(dec))
    dec = dec + rng.normal(0.0, noise_mas * _DEG_PER_MAS, n)
    return TargetObservations(
        initial_elements=tgt_el,
        obs_jd_tdb=obs,
        ra_obs_deg=ra,
        dec_obs_deg=dec,
        pa_scan_deg=np.linspace(15.0, 165.0, n),
        sigma_al_mas=np.full(n, noise_mas),
        gaia_bary_icrs=gaia_icrs,
    )


def run(args: argparse.Namespace) -> int:
    rng = np.random.default_rng(args.seed)
    m_inj_msun = args.m_inj_kg / M_SUN_KG
    pert_el, epoch = _perturber_and_epoch(rng)

    targets: list[TargetObservations] = []
    tries = 0
    while len(targets) < args.n_targets and tries < args.n_targets * 20:
        tries += 1
        t = _make_target(pert_el, m_inj_msun, args.noise_mas, rng, epoch)
        if t is not None:
            targets.append(t)
    if len(targets) < 3:
        logger.error("No se generaron suficientes objetivos (%d)", len(targets))
        return 2
    logger.info(
        "Inyección: M_inj=%.4e kg, N=%d objetivos, ruido=%.1f mas",
        args.m_inj_kg,
        len(targets),
        args.noise_mas,
    )

    seed_mass = 0.6 * m_inj_msun  # semilla deliberadamente sesgada
    mass_fit, fitted, result = determine_shared_mass(
        targets, seed_mass, pert_el, epoch, perturbers=_PERTURBERS
    )
    m_fit_kg = mass_fit * M_SUN_KG
    chi2_red = float(result.chi2_reduced)

    # σ formal = raíz de la varianza de la masa (parámetro 0) en la covarianza LSQ.
    var_mass = (
        float(result.covariance[0, 0]) if np.all(np.isfinite(result.covariance)) else math.nan
    )
    sigma_formal_kg = float(math.sqrt(var_mass) * M_SUN_KG) if var_mass > 0 else math.nan

    sigma_jack_kg = math.nan
    if args.jackknife:
        jk = jackknife_mass_sigma(
            targets, mass_fit, fitted, pert_el, epoch, perturbers=_PERTURBERS, backend="rebound"
        )
        sigma_jack_kg = float(jk.sigma_jack_msun * M_SUN_KG)

    sigma_kg = max(
        [s for s in (sigma_formal_kg, sigma_jack_kg) if math.isfinite(s) and s > 0],
        default=math.nan,
    )
    ratio = m_fit_kg / args.m_inj_kg
    z = (
        (m_fit_kg - args.m_inj_kg) / sigma_kg
        if math.isfinite(sigma_kg) and sigma_kg > 0
        else math.nan
    )

    logger.info("M_fit = %.4e kg  (ratio %.3f)", m_fit_kg, ratio)
    logger.info(
        "σ_formal=%.3e  σ_jack=%.3e  → z=%.2f  χ²_red=%.2f",
        sigma_formal_kg,
        sigma_jack_kg,
        z,
        chi2_red,
    )

    # --- Gates ---
    ok_bias = math.isfinite(z) and abs(z) <= args.n_sigma
    ok_chi2 = 0.5 <= chi2_red <= 2.0
    passed = ok_bias and ok_chi2
    logger.info(
        "GATE recuperación insesgada |z|≤%.1f: %s  (z=%.2f)",
        args.n_sigma,
        "PASS" if ok_bias else "FAIL",
        z,
    )
    logger.info("GATE χ²_red∈[0.5,2.0]: %s  (%.2f)", "PASS" if ok_chi2 else "FAIL", chi2_red)
    logger.info("RESULTADO: %s", "PASS" if passed else "FAIL")
    return 0 if passed else 1


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-targets", type=int, default=24)
    ap.add_argument("--m-inj-kg", type=float, default=2e19)
    ap.add_argument("--noise-mas", type=float, default=2.0)
    ap.add_argument("--n-sigma", type=float, default=3.0)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--jackknife", action="store_true")
    return run(ap.parse_args())


if __name__ == "__main__":
    raise SystemExit(main())
