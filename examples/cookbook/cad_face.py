"""Rebuild a block face on a parametric CAD surface (ADR 0007).

Writes ``cylinder.step`` and ``sector.geo`` in the current directory.
"""

# requires: cad

from __future__ import annotations

import build123d as bd
import numpy as np

from shore.geometry import BODY, EntityKind, EntityRef
from shore.project import Project

# #region find_surface
bd.export_step(bd.Cylinder(radius=0.5, height=4.0), "cylinder.step")

p = Project()
p.load_geometry("cylinder.step")
# an annulus sector whose inner face (i_lo) lies on the cylinder, clear of its seam
p.build_primitive(
    "annulus", r_in=0.5, r_out=1.5, z0=-1.5, z1=1.5, ni=9, nj=25, nk=21,
    theta_start=0.5, theta_end=2.0, j_law="tanh2", j_beta=2.0, k_law="tanh", k_beta=2.5,
)  # fmt: skip

geometry = p.body.geometry
surfaces = geometry.entities(EntityKind.SURFACE)
probe = 0.5 * np.array([np.cos(1.25), np.sin(1.25), 0.0])  # a point on the wanted face
k = int(geometry.project(probe, BODY).entity_index)
print(f"{len(surfaces)} surfaces; the face under the probe is surface {k}")
# #endregion find_surface

# #region fill_face
block = p.mesh.blocks[0]
before = block.nodes.copy()
p.surface_grid(block.label, "i_lo", surface=k)  # TFI + elliptic smoothing in (u, v)

face = block.nodes[0]
dist = geometry.project(face, EntityRef(EntityKind.SURFACE, k), exact=True).distances
print(f"largest distance from the CAD face: {dist.max():.1e}")
print(f"largest node move: {np.abs(face - before[0]).max():.1e}")
print("other nodes unchanged:", np.array_equal(block.nodes[1:], before[1:]))
print(p.journal().splitlines()[-1])
p.write_geo("sector")
# #endregion fill_face
