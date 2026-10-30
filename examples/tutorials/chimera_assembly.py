"""Tutorial 3: a sphere near-body grid inside a cylindrical background grid.

Writes sphere.stl, the project directory assembly/ (components "sphere" and
"background"), the blocks wall_<label>.geo and background_<label>.geo with
their sidecars wall.adjacency.json and background.adjacency.json, to the
current directory.
"""

from __future__ import annotations

import math

import numpy as np

# #region body
import shore.shapes

shore.shapes.sphere(radius=0.5, subdivisions=4).export("sphere.stl")
# #endregion body

# #region nearbody
from shore.project import Project  # noqa: E402

p = Project()
p.use_component("sphere")  # a named component: one grid of the assembly
p.load_geometry("sphere.stl", name="sphere")
p.build_topology(
    "cubed_sphere",
    nk=13,
    ds=0.01,
    growth=1.1,
    geometry="sphere",
    ni=30,
    nj=40,
    theta_cap=math.radians(30.0),
)
p.march()
p.write_geo("wall")  # wall_sub0..3.geo, wall_cap_north.geo, wall_cap_south.geo
p.write_adjacency("wall.adjacency.json")
# #endregion nearbody

# #region background
p.use_component("background")
p.build_primitive(
    "flat_caps",
    r_out=10.0,  # cylinder radius
    z0=-10.0,  # axial extent along z
    z1=50.0,
    inner_frac=0.4,  # corner radius of the central square / r_out
    ni=16,  # radial nodes of the four outer blocks
    nc=24,  # side nodes of the central square
    nk=60,  # axial nodes
)
p.write_geo("background")  # background_center.geo, background_sub_{e,n,w,s}.geo
p.write_adjacency("background.adjacency.json")
# #endregion background

# #region overlap
sphere = p.components["sphere"].blocks
background = p.components["background"].blocks
outer = np.concatenate([b.nodes[:, :, -1].reshape(-1, 3) for b in sphere])
r_outer = np.linalg.norm(outer, axis=1)
print(f"near-body outer layer: r = {r_outer.min():.3f} .. {r_outer.max():.3f}")
last = max(np.linalg.norm(b.nodes[:, :, -1] - b.nodes[:, :, -2], axis=-1).max() for b in sphere)
along = max(
    np.linalg.norm(np.diff(b.nodes[:, :, -1], axis=a), axis=-1).max()
    for b in sphere
    for a in (0, 1)
)
print(f"near-body outer cells: normal {last:.3f}, along the layer up to {along:.3f}")
centre = next(b for b in background if b.label == "center").nodes
dx = np.linalg.norm(np.diff(centre, axis=0), axis=-1).max()
dz = np.linalg.norm(np.diff(centre, axis=2), axis=-1).max()
print(f"background centre block: {centre.shape[:3]} nodes, cells {dx:.3f} across, {dz:.3f} along z")
# #endregion overlap

# #region save
p.save("assembly")  # project.toml, journal.py, inputs/, grids/<component>/
print(p.journal())
# #endregion save
