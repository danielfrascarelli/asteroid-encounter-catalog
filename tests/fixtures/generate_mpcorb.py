"""Helper to generate minimal MPCORB.DAT fixture lines for tests.

Not called during tests — run once to inspect output if needed.
"""

from __future__ import annotations


def mpcorb_line(  # noqa: N803
    number: int | str,
    a: float,
    designation: str,
    H: float = 5.0,  # noqa: N803
    G: float = 0.15,  # noqa: N803
    epoch: str = "K205A",
    M: float = 45.0,  # noqa: N803
    omega: float = 60.0,
    Omega: float = 80.0,  # noqa: N803
    i: float = 10.0,
    e: float = 0.08,
) -> str:
    """Return a 202-char MPCORB-format line with values at correct byte offsets."""
    n = 0.9856076 / (a**1.5)  # Kepler's third law (deg/day)

    chars = [" "] * 202

    def place(start: int, s: str, width: int, ljust: bool = False) -> None:
        s = (s.ljust(width) if ljust else s.rjust(width))[:width]
        for j, c in enumerate(s):
            chars[start + j] = c

    if isinstance(number, int):
        place(0, str(number), 7)
    else:
        place(0, str(number), 7, ljust=True)  # provisional: left-justify

    place(8, f"{H:.2f}", 5)
    place(14, f"{G:.3f}", 5)
    place(20, epoch, 5, ljust=True)
    place(26, f"{M:.5f}", 9)
    place(36, f"{omega:.5f}", 10)
    place(47, f"{Omega:.5f}", 10)
    place(58, f"{i:.5f}", 10)
    place(69, f"{e:.7f}", 10)
    place(80, f"{n:.8f}", 11)
    place(92, f"{a:.7f}", 11)
    place(104, "0", 1)
    place(106, "E2020-J83 ", 10, ljust=True)
    place(117, "12345", 5)
    place(123, " 10", 3)
    place(127, "2014-2020", 9, ljust=True)
    place(137, "0.6", 3)
    place(141, "M-vv", 4)
    place(146, "3M4", 3)
    place(166, designation, 28, ljust=True)

    return "".join(chars)


# Known reference values used by tests (a in AU)
CERES_A = 2.7691652
PALLAS_A = 2.7726856
VESTA_A = 2.3615491

HEADER = (
    "Des'n     H     G   Epoch     M        Peri.      Node       Incl."
    "     e            n           a        Reference #Obs #Opp    Arc"
    "    rms  Perts   Computer\n"
    + "-" * 160
    + "\n"
)

LINES = {
    "ceres": mpcorb_line(1, CERES_A, "(1) Ceres", H=3.34, G=0.12,
                         omega=73.597, Omega=80.306, i=10.587, e=0.0785),
    "pallas": mpcorb_line(2, PALLAS_A, "(2) Pallas", H=4.13, G=0.11,
                          omega=310.052, Omega=173.085, i=34.838, e=0.2299),
    "vesta": mpcorb_line(4, VESTA_A, "(4) Vesta", H=3.20, G=0.32,
                         omega=151.199, Omega=103.851, i=7.140, e=0.0887),
    # Provisional designation (should be filtered by only_numbered=True)
    "provisional": mpcorb_line("J98S00A", 2.50, "2014 AA1",
                               H=15.0, G=0.15, e=0.15, i=5.0),
    # Out-of-range semimajor axis (beyond main belt)
    "trojan": mpcorb_line(588, 5.2038, "(588) Achilles",
                          H=8.67, G=0.15, e=0.148, i=10.33),
}
