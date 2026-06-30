"""Tests de la covarianza en bloques por FOV (correlación intra-tránsito).

Verifican, sin red ni efemérides, que:
  * ``fov_groups_from_epochs`` agrupa los CCDs de un mismo cruce y respeta el orden.
  * ``_block_whiten`` reproduce el blanqueo diagonal clásico cuando no hay piso/grupos,
    y en presencia de piso coincide con la forma cuadrática explícita ``rᵀ C⁻¹ r`` y
    ``Jᵀ C⁻¹ J`` para ``C = diag(σ²) + s² 11ᵀ`` por bloque (y baja el χ² de residuos
    correlacionados, que es el efecto buscado).
"""

from __future__ import annotations

import numpy as np

from src.orbdet.gaia_adapter import fov_groups_from_epochs
from src.orbdet.mass_determination import _block_whiten


def test_fov_groups_cluster_and_preserve_order() -> None:
    # Dos cruces FOV: CCDs a ~5 s (~6e-5 d), separados por >0.05 d entre cruces.
    base = 100.0
    ep = np.array(
        [base, base + 6e-5, base + 1.2e-4, base + 0.08, base + 0.08006, base + 0.5],
        dtype=float,
    )
    g = fov_groups_from_epochs(ep, gap_days=0.01)
    assert g.tolist() == [0, 0, 0, 1, 1, 2]
    # Orden de entrada no monótono → etiquetas en el orden original.
    perm = np.array([5, 0, 3, 1, 4, 2])
    g2 = fov_groups_from_epochs(ep[perm], gap_days=0.01)
    # mismas etiquetas tras deshacer la permutación
    assert g2[np.argsort(perm)].tolist() == [0, 0, 0, 1, 1, 2]


def test_block_whiten_diagonal_without_floor() -> None:
    rng = np.random.default_rng(0)
    n = 10
    r = rng.normal(size=n)
    jac = rng.normal(size=(n, 7))
    sig = rng.uniform(0.5, 2.0, size=n)
    fov = np.array([0, 0, 0, 1, 1, 2, 2, 2, 3, 3])
    # piso 0 → diagonal exacto
    rw, jw = _block_whiten(r, jac, sig, fov, 0.0)
    np.testing.assert_allclose(rw, r / sig)
    np.testing.assert_allclose(jw, jac / sig[:, None])
    # sin grupos → diagonal exacto
    rw2, jw2 = _block_whiten(r, jac, sig, None, 3.0)
    np.testing.assert_allclose(rw2, r / sig)


def test_block_whiten_matches_explicit_inverse() -> None:
    rng = np.random.default_rng(1)
    n = 12
    r = rng.normal(size=n)
    jac = rng.normal(size=(n, 7))
    sig = rng.uniform(0.5, 2.0, size=n)
    fov = np.array([0, 0, 0, 0, 1, 1, 2, 2, 2, 3, 3, 3])
    s_c = 1.7

    rw, jw = _block_whiten(r, jac, sig, fov, s_c)
    chi2_block = float(rw @ rw)
    jtj_block = jw.T @ jw
    jtr_block = jw.T @ rw

    # Forma cuadrática explícita con C = diag(σ²) + s² 11ᵀ por bloque.
    chi2_ref = 0.0
    jtj_ref = np.zeros((7, 7))
    jtr_ref = np.zeros(7)
    for g in np.unique(fov):
        idx = np.where(fov == g)[0]
        cov = np.diag(sig[idx] ** 2) + s_c**2
        cinv = np.linalg.inv(cov)
        rb = r[idx]
        jb = jac[idx]
        chi2_ref += float(rb @ cinv @ rb)
        jtj_ref += jb.T @ cinv @ jb
        jtr_ref += jb.T @ cinv @ rb

    np.testing.assert_allclose(chi2_block, chi2_ref, rtol=1e-10)
    np.testing.assert_allclose(jtj_block, jtj_ref, rtol=1e-9, atol=1e-12)
    np.testing.assert_allclose(jtr_block, jtr_ref, rtol=1e-9, atol=1e-12)


def test_block_whiten_downweights_common_mode() -> None:
    # Residuos idénticos dentro de un bloque (totalmente correlacionados): el piso
    # debe bajar la χ² respecto al conteo independiente.
    sig = np.full(8, 1.0)
    r = np.concatenate([np.full(4, 2.0), np.full(4, -1.5)])  # 2 FOV, residuo común
    fov = np.array([0, 0, 0, 0, 1, 1, 1, 1])
    jac = np.zeros((8, 1))
    rw0, _ = _block_whiten(r, jac, sig, fov, 0.0)
    rw1, _ = _block_whiten(r, jac, sig, fov, 3.0)
    assert float(rw1 @ rw1) < float(rw0 @ rw0)
