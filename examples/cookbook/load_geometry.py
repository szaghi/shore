"""Load STL and CAD geometry, list its entities and project points onto it.

Writes wing_body.step, wing_body.brep, wing_body.igs, sphere.step and sphere.stl to the
current directory.
"""

# requires: cad

from __future__ import annotations

import build123d as bd
import numpy as np
from OCP.IGESControl import IGESControl_Writer

from shore import shapes
from shore.geometry import EntityKind, EntityRef, load_geometry
from shore.mesh.body import Body
from shore.profiles import naca
from shore.project import Project

# #region make_cad
# A half wing-body (fuselage + NACA 0012 wing, cut at the symmetry plane z = 0),
# written as STEP, BREP and IGES. Any CAD tool's export works the same way.
radius, length, span = 0.5, 5.0, 3.0
fuselage = bd.Cylinder(radius, length, rotation=(0, 90, 0))
fuselage += bd.Pos(-length / 2, 0, 0) * bd.Sphere(radius) + bd.Pos(length / 2, 0, 0) * bd.Sphere(
    radius
)
sec = naca("0012", 40, closed_te=False)
le = int(np.argmin(sec[:, 0]))
upper = [(float(x), float(y), 0.0) for x, y in sec[: le + 1]]
lower = [(float(x), float(y), 0.0) for x, y in sec[le:]]
wire = bd.Wire([bd.Spline(*upper), bd.Spline(*lower), bd.Line(lower[-1], upper[0])])
wing = bd.extrude(bd.make_face(wire), amount=span)
model = (fuselage + wing) & (bd.Pos(0, 0, 10) * bd.Box(40, 40, 20))
bd.export_step(model, "wing_body.step")
bd.export_brep(model, "wing_body.brep")
writer = IGESControl_Writer()
writer.AddShape(model.wrapped)
writer.ComputeModel()
writer.Write("wing_body.igs")
# #endregion make_cad

# #region entities
for path in ("wing_body.step", "wing_body.brep", "wing_body.igs"):
    geo = load_geometry(path)  # backend chosen by suffix
    counts = {kind.value: len(geo.entities(kind)) for kind in EntityKind}
    print(f"{path:15s} {counts}")
geo = load_geometry("wing_body.step")
print("bounds", np.round(geo.bounds, 3).tolist())
display = geo.display()  # what the GUI draws: triangles per face, edge polylines
print(f"tessellation: {len(display.triangles)} triangles over {len(geo.faces)} faces")
# #endregion entities

# #region exact
# A STEP sphere of radius 1 and the same sphere as a 1280-triangle STL.
bd.export_step(bd.Sphere(1.0), "sphere.step")
shapes.sphere(radius=1.0, subdivisions=3).export("sphere.stl")
points = np.random.default_rng(0).normal(size=(200, 3))
points *= 1.3 / np.linalg.norm(points, axis=1, keepdims=True)  # on a sphere of radius 1.3

cad = load_geometry("sphere.step")
stl = load_geometry("sphere.stl")
for label, geom, exact in (
    ("STL", stl, False),
    ("CAD tier 1", cad, False),
    ("CAD exact", cad, True),
):
    proj = geom.project(points, exact=exact)  # closest point on the whole body
    err = np.abs(np.linalg.norm(proj.points, axis=1) - 1.0).max()
    print(f"{label:10s} max |r - 1| = {err:.1e}")
face = EntityRef(EntityKind.SURFACE, 0)
on_face = cad.project(points, face, exact=True)
print("(u, v) of the first point on face 0:", np.round(on_face.params[0], 4))
# #endregion exact

# #region named
p = Project()
p.load_geometry("wing_body.step", name="wing_body")
p.load_geometry("sphere.stl", name="ball")
for name, body in p.geometries.items():
    print(f"{name}: {type(body.geometry).__name__}, CAD={body.is_cad}")
body = Body.from_file("sphere.step")  # what the marching topologies take
print("anchor", np.round(body.anchor, 6), "CAD", body.is_cad)
# #endregion named
