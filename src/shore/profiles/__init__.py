"""2D airfoil profile generators (NACA 4-digit, 5-digit).

Pure NumPy.  No trimesh dependency — :mod:`shore.profiles` produces
``(N, 2)`` Selig-ordered coordinates only.  For 3D STL extrusion
into wings, see :mod:`shore.shapes`.

Selig ordering convention
-------------------------
Points are ordered around the closed section starting at the upper
trailing edge, walking forward over the suction (upper) surface to
the leading edge, then back over the pressure (lower) surface to
the lower trailing edge::

    out[0]    = upper TE        (x ≈ 1, y ≥ 0)
    out[N//2] = leading edge    (x ≈ 0, y ≈ 0)
    out[-1]   = lower TE        (x ≈ 1, y ≤ 0)

For a **sharp** TE (``closed_te=True``, default) ``out[0] == out[-1]``
in the limit; the polynomial residual is forced to zero so the curve
is exactly closed.  For a **blunt** TE (``closed_te=False``), the two
endpoints are separated by the analytic NACA polynomial residual
(``≈ 0.1% chord`` for the standard 4-digit family).
"""

from __future__ import annotations

from shore.profiles._naca import naca, naca4, naca5
from shore.profiles.io import read_section, write_section

__all__ = ["naca", "naca4", "naca5", "read_section", "write_section"]
