"""3D body STL builders backed by trimesh.

This module produces ``trimesh.Trimesh`` objects (and writes ``.stl``
files via the standard ``mesh.export(path)`` call).  Two families:

- **Analytic primitives** for testing and quick CFD studies:
  :func:`sphere`, :func:`cylinder`, :func:`torus`, :func:`box`,
  :func:`capsule`.  Thin wrappers over ``trimesh.creation``.

- **Section extrusion**: :func:`extrude_section` takes a 2D
  Selig-ordered profile (e.g. from :mod:`shore.profiles`) plus
  span / taper / sweep / twist parameters and produces a wing STL.

For the Selig 2D format that :func:`extrude_section` consumes, see
:mod:`shore.profiles`.
"""

from __future__ import annotations

from shore.shapes._extrude import extrude_section
from shore.shapes._primitives import box, capsule, cylinder, sphere, torus

__all__ = [
    "box",
    "capsule",
    "cylinder",
    "extrude_section",
    "sphere",
    "torus",
]
