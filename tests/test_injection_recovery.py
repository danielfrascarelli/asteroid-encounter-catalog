"""Injection-recovery de detección como regresión permanente (tribunal M10/B1).

Versión compacta (N=25) del harness de
``scripts.validate.injection_recovery_detection``: pares sintéticos con mínimos
en fase uniforme respecto de la grilla gruesa de 12 h. Falla si un mínimo
inter-sample se pierde o llega sesgado (> 1 μAU o 0.1 %) — la firma exacta de B1.
"""

from __future__ import annotations

from scripts.validate.injection_recovery_detection import run_injection_recovery


def test_injection_recovery_detection_gates() -> None:
    s = run_injection_recovery(n_pairs=25, seed=42, window_days=20.0)
    assert s["n_evaluable"] >= 15  # sanity del generador
    assert s["gate_recovery_ge_99pct"], f"recovery {s['recovery_frac']:.3f}, misses={s['misses']}"
    assert s["gate_distance_ratio_1"], f"max |Δd| = {s['max_abs_d_error_au']:.3e} AU"
    assert s["gate_epoch_within_fine_step"], f"max |Δt| = {s['max_abs_t_error_s']:.1f} s"
