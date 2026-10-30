"""Associate block faces, edges and vertices with CAD entities, then snap (ADR 0039).

Writes ``part.step`` and ``snapped.geo`` in the current directory.
"""

# requires: cad

from __future__ import annotations

import build123d as bd
import numpy as np

from shore.geometry import EntityKind, EntityRef
from shore.project import Project

# #region find_entities
# a 2 x 1 x 1 box with a hole; its top face is z = 0.5
bd.export_step(bd.Box(2.0, 1.0, 1.0) - bd.Cylinder(0.3, 2.0), "part.step")

p = Project()
p.load_geometry("part.step", name="part")
p.build_primitive("box", vmin=(-0.8, -0.4, 0.55), vmax=(0.8, 0.4, 0.9), ni=5, nj=4, nk=3)
g = p.geometries["part"].geometry


def nearest(kind: EntityKind, point: tuple[float, float, float]) -> int:
    """Index of the entity of *kind* closest to *point*."""
    xyz = np.array([point])
    d = [g.project(xyz, ref).distances[0] for ref in g.entities(kind)]
    return int(np.argmin(d))


top = nearest(EntityKind.SURFACE, (0.7, 0.3, 0.5))  # the top face
edge = nearest(EntityKind.CURVE, (0.2, -0.5, 0.5))  # the top edge along x at y = -0.5
corner = nearest(EntityKind.POINT, (-1.0, -0.5, 0.5))  # the corner it starts from
print(f"surface {top}, curve {edge}, point {corner}")
# #endregion find_entities

# #region associate_snap
p.associate("box", "k_lo", "part", "surface", [top])  # face -> surface
p.associate("box", "j_lo_k_lo", "part", "curve", [edge])  # edge -> curve
p.associate("box", "v0", "part", "point", [corner])  # vertex -> point
print(p.journal().splitlines()[-1])

p.snap()  # nodes move only now: vertex first, then edges, then faces
n = p.mesh.blocks[0].nodes
print("v0         :", n[0, 0, 0])
print("edge y, z  :", np.unique(n[:, 0, 0, 1:].round(12), axis=0))
print("face z     :", np.unique(n[:, :, 0, 2].round(12)))
dist = g.project(n[:, :, 0].reshape(-1, 3), EntityRef(EntityKind.SURFACE, top)).distances
print(f"k_lo off the surface by at most {dist.max():.1e}")
p.write_geo("snapped")
# #endregion associate_snap

# #region clear
p.associate("box", "v0", "part", "none")  # clear the vertex association
print(p.mesh.blocks[0].vertices()[0].association)
p.undo()  # and bring it back
print(p.mesh.blocks[0].vertices()[0].association is not None)
# #endregion clear
