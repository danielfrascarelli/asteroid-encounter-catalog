"""Diagnóstico: magnitud del error along-scan sistemático vs aleatorio.

El error sistemático (actitud/calibración) está correlacionado entre los CCDs de
un mismo tránsito FOV; el aleatorio (fotónico) es independiente por CCD. La fuerza
de la correlación intra-FOV la fija el cociente σ_sys/σ_rand proyectados sobre la
dirección de barrido. Si σ_sys domina, los ~7 CCDs de un FOV son casi un único
grado de libertad → σ de la masa subestimada ~√7.
"""

from __future__ import annotations

import sys

import numpy as np

from scripts.mass.fit_mass_gaia_loo import fetch_gaia_full
from src.utils.config import load_config


def _sigma_al(s_ra, s_dec, rho, pa_deg):
    pa = np.radians(pa_deg)
    e_ra, e_dec = np.sin(pa), np.cos(pa)
    var = e_ra**2 * s_ra**2 + 2 * e_ra * e_dec * rho * s_ra * s_dec + e_dec**2 * s_dec**2
    return np.sqrt(np.maximum(var, 0.0))


def main() -> None:
    targets = [int(x) for x in (sys.argv[1:] or ["18937"])]
    release = "fpr"
    cfg = load_config("config.yaml")
    gaia = cfg.sources.gaia_sso
    gaia.release = release
    for target in targets:
        df = fetch_gaia_full(gaia.archive_url, target, gaia.active())
        pa = df["position_angle_scan"].to_numpy().astype(float)
        sig_sys = _sigma_al(
            df["ra_error_systematic"].to_numpy().astype(float),
            df["dec_error_systematic"].to_numpy().astype(float),
            df["ra_dec_correlation_systematic"].to_numpy().astype(float),
            pa,
        )
        sig_rand = _sigma_al(
            df["ra_error_random"].to_numpy().astype(float),
            df["dec_error_random"].to_numpy().astype(float),
            df["ra_dec_correlation_random"].to_numpy().astype(float),
            pa,
        )
        sig_tot = np.sqrt(sig_sys**2 + sig_rand**2)
        frac_sys = sig_sys**2 / sig_tot**2  # fracción de varianza correlacionada
        print(
            f"target {target}: n={df.height}  "
            f"σ_sys={np.median(sig_sys):.3f} mas  σ_rand={np.median(sig_rand):.3f} mas  "
            f"σ_tot={np.median(sig_tot):.3f} mas  "
            f"var_sys/var_tot(med)={np.median(frac_sys):.2f}"
        )


if __name__ == "__main__":
    main()
