"""Find the feature curves and points of an STL and snap a grid edge onto one.

Writes wing.stl to the current directory.
"""

from __future__ import annotations

import numpy as np

from shore import profiles, shapes
from shore.geometry import EntityKind, load_geometry
from shore.geometry.stl import StlGeometry
from shore.project import Project

# #region list
# A blunt-TE NACA 0012 wing, span 2 along z, with flat end caps.
wing = shapes.extrude_section(
    profiles.naca("0012", 60, closed_te=False), span=2.0, n_stations=21, end_caps=True
)
wing.export("wing.stl")

geo = load_geometry("wing.stl")  # StlGeometry, crease_deg = 30
curves, points = geo.features()
print(f"{len(geo.entities(EntityKind.CURVE))} curves, {len(points)} points")
for k, c in enumerate(curves):
    lo, hi = c.min(axis=0), c.max(axis=0)
    print(f"  curve {k}: {len(c):3d} nodes, x {lo[0]:.3f}..{hi[0]:.3f}, z {lo[2]:.1f}..{hi[2]:.1f}")
# #endregion list

# #region crease
# A 32-sided cylinder: neighbouring side facets meet at 360/32 = 11.25 degrees.
cylinder = shapes.cylinder(radius=0.5, height=1.0, sections=32)
for crease in (30.0, 10.0):
    curves_c, points_c = StlGeometry(cylinder, crease_deg=crease).features()
    print(f"crease {crease:4.1f} deg: {len(curves_c)} curves, {len(points_c)} points")
# #endregion crease

# #region snap
# The upper trailing-edge line: x = 1 along the whole span, y > 0.
upper = next(
    k
    for k, c in enumerate(curves)
    if np.allclose(c[:, 0], 1.0) and np.all(c[:, 1] > 0) and np.ptp(c[:, 2]) > 1.0
)
p = Project()
p.load_geometry("wing.stl", name="wing")
p.build_primitive("box", vmin=(1.0, 0.0, 0.1), vmax=(1.5, 0.5, 1.9), ni=3, nj=3, nk=9)
p.associate("box", "i_lo_j_lo", "wing", "curve", [upper])  # the k-edge at i = 0, j = 0
p.snap()
edge = p.mesh.blocks[0].nodes[0, 0, :]
dist = geo.project(edge, geo.entities(EntityKind.CURVE)[upper]).distances
print(f"edge on curve {upper}: max distance {dist.max():.1e}")
# #endregion snap
