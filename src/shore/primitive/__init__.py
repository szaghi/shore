"""Parametric structured-grid primitives for Chimera/overset assembly.

A *primitive* is a mesh built from simple geometric inputs (corners,
radii, axis) — no STL, no projection, no marching.  Each primitive
returns a :class:`~shore.primitive.base.PrimitiveMesh` whose blocks
follow the same :class:`~shore.mesh.block.HexBlock` contract used by
body-fitted meshes, so downstream tooling (``.geo`` writer, VTK
export, Jacobian inspection) is shared.

Available primitives
--------------------

- :class:`~shore.primitive.box.BoxMesh` — single-block axis-aligned
  Cartesian box.  Use case: outer background, simple buffer regions.
- :class:`~shore.primitive.annulus.AnnularCylinderMesh` — single-block
  annular cylinder, optionally a partial sector.  Use case: buffer
  around a near-body sphere/airfoil mesh; cylindrical pipe sections.
- :class:`~shore.primitive.flat_caps.FlatCapsCylinder` — 5-block
  flat-caps cylinder (central square + 4 trapezoidal sub-blocks).
  Use case: cylindrical background without Chimera fringes (high-order
  schemes, conservation, GPU memory locality).
"""

from __future__ import annotations

from shore.primitive.annulus import AnnularCylinderMesh, AnnulusExtent
from shore.primitive.base import PrimitiveMesh, axis_frame, coords_1d, resolve_axis
from shore.primitive.box import BoxExtent, BoxMesh
from shore.primitive.flat_caps import FlatCapsCylinder, FlatCapsExtent

__all__ = [
    "AnnularCylinderMesh",
    "AnnulusExtent",
    "BoxExtent",
    "BoxMesh",
    "FlatCapsCylinder",
    "FlatCapsExtent",
    "PrimitiveMesh",
    "axis_frame",
    "coords_1d",
    "resolve_axis",
]
