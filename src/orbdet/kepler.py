"""Mecánica de dos cuerpos (kepleriana): elementos ↔ estado, propagación.

Elementos orbitales clásicos (marco eclíptico, radianes, AU):
    a      semieje mayor (AU)
    e      excentricidad
    i      inclinación
    Omega  longitud del nodo ascendente (Ω)
    omega  argumento del perihelio (ω)
    M      anomalía media en la época

GM en AU³/día² (por defecto ``GM_SUN``). Solo órbitas elípticas (0 ≤ e < 1).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .constants import GM_SUN
from .frames import drotation_x, drotation_z, rotation_x, rotation_z


@dataclass(frozen=True)
class KeplerElements:
    """Elementos keplerianos clásicos (radianes, AU)."""

    a: float
    e: float
    i: float
    Omega: float
    omega: float
    M: float

    def as_array(self) -> np.ndarray:
        return np.array([self.a, self.e, self.i, self.Omega, self.omega, self.M], dtype=float)


def solve_kepler(M: float, e: float, tol: float = 1e-13, max_iter: int = 100) -> float:
    """Resuelve la ecuación de Kepler ``M = E - e·sin(E)`` para la anomalía
    excéntrica E (rad), por Newton-Raphson.

    Converge para 0 ≤ e < 1. El valor inicial ``E0 = M`` (o ``π`` si e alto)
    da convergencia cuadrática en pocas iteraciones.
    """
    if not 0.0 <= e < 1.0:
        raise ValueError(f"solve_kepler solo soporta órbitas elípticas (0<=e<1); e={e}")
    m = math.fmod(M, 2.0 * math.pi)
    e_anom = m if e < 0.8 else math.pi
    for _ in range(max_iter):
        f = e_anom - e * math.sin(e_anom) - m
        fp = 1.0 - e * math.cos(e_anom)
        delta = f / fp
        e_anom -= delta
        if abs(delta) < tol:
            break
    return e_anom


def true_anomaly_from_eccentric(E: float, e: float) -> float:
    """Anomalía verdadera ν (rad) a partir de la excéntrica E."""
    return 2.0 * math.atan2(
        math.sqrt(1.0 + e) * math.sin(E / 2.0), math.sqrt(1.0 - e) * math.cos(E / 2.0)
    )


def eccentric_from_true_anomaly(nu: float, e: float) -> float:
    """Anomalía excéntrica E (rad) a partir de la verdadera ν."""
    return 2.0 * math.atan2(
        math.sqrt(1.0 - e) * math.sin(nu / 2.0), math.sqrt(1.0 + e) * math.cos(nu / 2.0)
    )


def mean_motion(a: float, mu: float = GM_SUN) -> float:
    """Movimiento medio n = sqrt(mu / a³) (rad/día)."""
    return math.sqrt(mu / (a**3))


def period(a: float, mu: float = GM_SUN) -> float:
    """Período orbital 2π / n (días)."""
    return 2.0 * math.pi / mean_motion(a, mu)


def elements_to_state(el: KeplerElements, mu: float = GM_SUN) -> tuple[np.ndarray, np.ndarray]:
    """Elementos → estado cartesiano (r, v) en el marco de referencia (eclíptico).

    Devuelve ``(r_vec, v_vec)`` en AU y AU/día.
    """
    a, e = el.a, el.e
    E = solve_kepler(el.M, e)
    cosE, sinE = math.cos(E), math.sin(E)
    r = a * (1.0 - e * cosE)

    # Posición y velocidad en el plano perifocal (x hacia el perihelio).
    r_pf = np.array([a * (cosE - e), a * math.sqrt(1.0 - e * e) * sinE, 0.0])
    # v = sqrt(mu·a)/r · [-sinE, sqrt(1-e²)·cosE, 0]
    fac = math.sqrt(mu * a) / r
    v_pf = np.array([-fac * sinE, fac * math.sqrt(1.0 - e * e) * cosE, 0.0])

    # Perifocal → referencia: Rz(Ω) · Rx(i) · Rz(ω)
    rot = rotation_z(el.Omega) @ rotation_x(el.i) @ rotation_z(el.omega)
    return rot @ r_pf, rot @ v_pf


def dstate_delements(el: KeplerElements, mu: float = GM_SUN) -> np.ndarray:
    """Jacobiano analítico ``∂[r, v] / ∂[a, e, i, Ω, ω, M]`` en la época.

    Devuelve la matriz ``(6, 6)`` cuyas filas son ``(rx, ry, rz, vx, vy, vz)`` y
    cuyas columnas son las parciales respecto a ``(a, e, i, Ω, ω, M)`` (mismo
    orden que :meth:`KeplerElements.as_array`). Es el mapa estático de
    reparametrización elementos→estado que, compuesto con la matriz de transición
    de estado de las ecuaciones variacionales (``variational.py``), da
    ``∂[r, v](t) / ∂elementos``.

    Derivadas en marco perifocal propagadas al marco de referencia con
    ``rot = Rz(Ω)·Rx(i)·Rz(ω)``; las parciales respecto a ``i, Ω, ω`` actúan solo
    sobre ``rot`` (el estado perifocal no depende de los ángulos), y las parciales
    respecto a ``a, e, M`` solo sobre el estado perifocal. Las dependencias vía la
    anomalía excéntrica ``E(M, e)`` se incluyen por regla de la cadena
    (``∂E/∂M = 1/(1-e·cosE)``, ``∂E/∂e = sinE/(1-e·cosE)``).
    """
    a, e = el.a, el.e
    E = solve_kepler(el.M, e)
    cosE, sinE = math.cos(E), math.sin(E)
    beta = math.sqrt(1.0 - e * e)  # √(1-e²)
    s = 1.0 - e * cosE  # r/a
    g = math.sqrt(mu / a)
    fac = g / s  # √(μa)/r

    r_pf = np.array([a * (cosE - e), a * beta * sinE, 0.0])
    v_pf = np.array([-fac * sinE, fac * beta * cosE, 0.0])

    dE_de = sinE / s
    dE_dM = 1.0 / s

    # --- ∂r_pf / ∂(a, e, M) ---
    drpf_da = r_pf / a
    drpf_de = np.array(
        [
            a * (-sinE * dE_de - 1.0),
            a * (-(e / beta) * sinE + beta * cosE * dE_de),
            0.0,
        ]
    )
    drpf_dM = np.array([a * (-sinE * dE_dM), a * beta * cosE * dE_dM, 0.0])

    # --- ∂v_pf / ∂(a, e, M) ---
    dvpf_da = -v_pf / (2.0 * a)
    ds_de = -cosE + e * sinE * dE_de  # ∂(1-e·cosE)/∂e total
    ds_dM = e * sinE * dE_dM
    dfac_de = -fac / s * ds_de
    dfac_dM = -fac / s * ds_dM
    dbeta_de = -e / beta
    dvpf_de = np.array(
        [
            -(dfac_de * sinE + fac * cosE * dE_de),
            dfac_de * beta * cosE + fac * dbeta_de * cosE - fac * beta * sinE * dE_de,
            0.0,
        ]
    )
    dvpf_dM = np.array(
        [
            -(dfac_dM * sinE + fac * cosE * dE_dM),
            dfac_dM * beta * cosE - fac * beta * sinE * dE_dM,
            0.0,
        ]
    )

    # Rotación perifocal → referencia y sus derivadas respecto a los ángulos.
    rz_O, rx_i, rz_w = rotation_z(el.Omega), rotation_x(el.i), rotation_z(el.omega)
    rot = rz_O @ rx_i @ rz_w
    drot_dOmega = drotation_z(el.Omega) @ rx_i @ rz_w
    drot_di = rz_O @ drotation_x(el.i) @ rz_w
    drot_domega = rz_O @ rx_i @ drotation_z(el.omega)

    jac = np.zeros((6, 6))
    jac[0:3, 0], jac[3:6, 0] = rot @ drpf_da, rot @ dvpf_da  # a
    jac[0:3, 1], jac[3:6, 1] = rot @ drpf_de, rot @ dvpf_de  # e
    jac[0:3, 2], jac[3:6, 2] = drot_di @ r_pf, drot_di @ v_pf  # i
    jac[0:3, 3], jac[3:6, 3] = drot_dOmega @ r_pf, drot_dOmega @ v_pf  # Ω
    jac[0:3, 4], jac[3:6, 4] = drot_domega @ r_pf, drot_domega @ v_pf  # ω
    jac[0:3, 5], jac[3:6, 5] = rot @ drpf_dM, rot @ dvpf_dM  # M
    return jac


def state_to_elements(r_vec: np.ndarray, v_vec: np.ndarray, mu: float = GM_SUN) -> KeplerElements:
    """Estado cartesiano (r, v) → elementos keplerianos.

    Asume órbita elíptica acotada (energía < 0). Para casos degenerados (e≈0 o
    i≈0) los ángulos no definidos se fijan a 0 de forma consistente.
    """
    r_vec = np.asarray(r_vec, dtype=float)
    v_vec = np.asarray(v_vec, dtype=float)
    r = float(np.linalg.norm(r_vec))
    v = float(np.linalg.norm(v_vec))

    h_vec = np.cross(r_vec, v_vec)
    h = float(np.linalg.norm(h_vec))
    # Vector nodo = ẑ × h
    n_vec = np.cross(np.array([0.0, 0.0, 1.0]), h_vec)
    n = float(np.linalg.norm(n_vec))

    # Vector excentricidad.
    e_vec = ((v * v - mu / r) * r_vec - float(np.dot(r_vec, v_vec)) * v_vec) / mu
    e = float(np.linalg.norm(e_vec))

    energy = v * v / 2.0 - mu / r
    a = -mu / (2.0 * energy)

    i = math.acos(max(-1.0, min(1.0, h_vec[2] / h)))

    if n > 1e-12:
        Omega = math.acos(max(-1.0, min(1.0, n_vec[0] / n)))
        if n_vec[1] < 0.0:
            Omega = 2.0 * math.pi - Omega
    else:  # órbita ecuatorial: nodo indefinido
        Omega = 0.0

    if n > 1e-12 and e > 1e-12:
        omega = math.acos(max(-1.0, min(1.0, float(np.dot(n_vec, e_vec)) / (n * e))))
        if e_vec[2] < 0.0:
            omega = 2.0 * math.pi - omega
    elif e > 1e-12:  # ecuatorial con e>0: usar longitud del perihelio desde X
        omega = math.atan2(e_vec[1], e_vec[0])
        if h_vec[2] < 0.0:
            omega = 2.0 * math.pi - omega
    else:
        omega = 0.0

    if e > 1e-12:
        nu = math.acos(max(-1.0, min(1.0, float(np.dot(e_vec, r_vec)) / (e * r))))
        if float(np.dot(r_vec, v_vec)) < 0.0:
            nu = 2.0 * math.pi - nu
    else:  # circular: anomalía verdadera desde el nodo (o X si ecuatorial)
        ref = n_vec if n > 1e-12 else np.array([1.0, 0.0, 0.0])
        denom = (n if n > 1e-12 else 1.0) * r
        nu = math.acos(max(-1.0, min(1.0, float(np.dot(ref, r_vec)) / denom)))
        if r_vec[2] < 0.0:
            nu = 2.0 * math.pi - nu

    E = eccentric_from_true_anomaly(nu, e)
    M = E - e * math.sin(E)
    return KeplerElements(
        a=a, e=e, i=i, Omega=Omega, omega=omega, M=math.fmod(M + 2.0 * math.pi, 2.0 * math.pi)
    )


def propagate(el: KeplerElements, dt_days: float, mu: float = GM_SUN) -> KeplerElements:
    """Propaga los elementos *dt_days* días (dos cuerpos: solo avanza M)."""
    n = mean_motion(el.a, mu)
    M_new = math.fmod(el.M + n * dt_days, 2.0 * math.pi)
    if M_new < 0.0:
        M_new += 2.0 * math.pi
    return KeplerElements(a=el.a, e=el.e, i=el.i, Omega=el.Omega, omega=el.omega, M=M_new)
