"""orbdet — motor de determinación de órbitas y masas (autocontenido).

Librería de astrodinámica para determinación de órbitas por mínimos cuadrados y,
sobre ella, determinación de masas de asteroides perturbadores (estrategia tipo
Fuentes-Muñoz). Ver ``planning/MASS_DETERMINATION_PLAN.md``.

Contrato de aislamiento
-----------------------
``orbdet`` NO importa ningún otro módulo del proyecto (``src.detect``,
``src.mass``, ``src.propagate``, etc.). Depende solo de la librería estándar +
numpy/scipy/astropy/rebound. Internamente usa imports relativos. Esto está
verificado por ``tests/orbdet/test_isolation.py`` y permite testear el motor de
forma totalmente independiente del pipeline. Ver ``src/orbdet/README.md``.

Convenciones de unidades (internas, hot-path)
---------------------------------------------
- Distancias: AU.
- Tiempos: días (JD en escala TDB salvo que se indique).
- Ángulos: radianes.
- Velocidades: AU/día.
- GM: AU³/día².
"""

from __future__ import annotations

__all__ = [
    "constants",
    "dynamics",
    "frames",
    "kepler",
    "observation",
    "time_scales",
    "variational",
]
