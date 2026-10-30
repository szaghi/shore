"""Collar grids at the junction of two bodies: full, half (symmetry plane) and split.

Writes ``box.stl``, ``post.stl``, ``hbox.stl``, ``hpost.stl``, ``ellipse.stl``,
``plate.stl``, ``wing.stl`` and the collars ``collar.geo``, ``half.geo``
(+ ``half.adjacency.json``), ``ellipse.geo`` and ``wing_<label>.geo`` (one per
split block) in the current directory. Takes about 20 s.
"""

from __future__ import annotations

import warnings

import numpy as np
import trimesh

from shore.profiles import naca
from shore.project import Project
from shore.quality import scaled_jacobian
from shore.shapes import extrude_section


def min_jacobian(blocks, periodic=False):
    """Smallest scaled Jacobian over *blocks* (with the j-seam cells of a periodic block)."""
    out = []
    for b in blocks:
        nodes = np.concatenate([b.nodes, b.nodes[:, :1]], axis=1) if periodic else b.nodes
        out.append(float(scaled_jacobian(nodes).min()))
    return min(out)


# #region bodies
# body A: a 4 x 4 x 1 box, top face at z = 0; body B: a post of radius 0.5 through it
box = trimesh.creation.box(extents=(4.0, 4.0, 1.0))
box.apply_translation((0.0, 0.0, -0.5))
post = trimesh.creation.cylinder(radius=0.5, height=1.6, sections=96)
post.apply_translation((0.0, 0.0, 0.2))
box.export("box.stl")
post.export("post.stl")
# #endregion bodies

# #region full
p = Project()
p.load_geometry("box.stl")  # geometry "body": body A, where m_a layers march
p.build_collar(
    "post.stl",  # body B: a file path (or the name of a loaded geometry)
    seed=(0.55, 0.0, 0.02),  # a point near the junction: selects the curve
    n_loop=48,  # nodes round the junction curve
    ds_surface=0.01, growth_surface=1.12, m_a=12, m_b=12,  # surface layers on A and B
    ds_wall=1e-4, growth_wall=1.2, nk=25,  # wall-normal layers
)  # fmt: skip
block = p.mesh.blocks[0]
print("full collar:", block.label, block.nodes.shape,
      f"min scaled Jacobian {min_jacobian([block], periodic=True):.3f}")  # fmt: skip
print("faces:", {f.name: f.bc.name for f in block.faces()})
p.write_geo("collar")
# #endregion full

# #region half
# the half model: the box and the post cut at y = 0, the plane y = 0 a symmetry plane
hbox = trimesh.creation.box(extents=(4.0, 2.0, 1.0))
hbox.apply_translation((0.0, 1.0, -0.5))
hbox.update_faces(~np.all(np.isclose(hbox.vertices[hbox.faces][..., 1], 0.0), axis=1))
hpost = post.copy()
hpost.update_faces(hpost.triangles_center[:, 1] > 0.0)
for mesh in (hbox, hpost):
    mesh.remove_unreferenced_vertices()
hbox.export("hbox.stl")
hpost.export("hpost.stl")

h = Project()
with warnings.catch_warnings():
    warnings.simplefilter("ignore")  # open half-model STL files are not watertight
    h.load_geometry("hbox.stl", skip_hygiene=True)
    h.build_collar(
        "hpost.stl",
        seed=(0.45, 0.2, 0.02),  # on the y > 0 side
        n_loop=25,
        ds_surface=0.01, growth_surface=1.12, m_a=12, m_b=12,
        ds_wall=1e-4, growth_wall=1.2, nk=25,
        plane=((0.0, 0.0, 0.0), (0.0, 1.0, 0.0)),  # (point, normal)
    )  # fmt: skip
half = h.mesh.blocks[0]
print("half collar:", half.nodes.shape, f"min scaled Jacobian {min_jacobian([half]):.3f}",
      "j-faces:", half.face("j_lo").bc.name, half.face("j_hi").bc.name,
      f"off the plane {np.abs(half.nodes[:, [0, -1], :, 1]).max():.1e}")  # fmt: skip
h.write_geo("half")
h.write_adjacency("half.adjacency.json")
# #endregion half

# #region spacing
# an elliptic post (1 x 0.2): uniform nodes leave too few at its narrow ends
ellipse = trimesh.creation.cylinder(radius=1.0, height=1.6, sections=256)
ellipse.apply_scale((0.5, 0.1, 1.0))
ellipse.apply_translation((0.0, 0.0, 0.2))
ellipse.export("ellipse.stl")

e = Project()
e.load_geometry("box.stl")
steps = {
    "seed": (0.0, 0.11, 0.02), "n_loop": 96, "ds_surface": 0.005, "growth_surface": 1.12,
    "m_a": 15, "m_b": 15, "ds_wall": 5e-5, "growth_wall": 1.2, "nk": 30,
}  # fmt: skip
e.build_collar("ellipse.stl", **steps)
print(f"ellipse, uniform loop: min scaled Jacobian {min_jacobian(e.mesh.blocks, True):.3f}")
e.undo()
e.build_collar("ellipse.stl", **steps, loop_law="tanh2", loop_beta=2.0, loop_curvature=1.0)
print(f"ellipse, curvature 1 : min scaled Jacobian {min_jacobian(e.mesh.blocks, True):.3f}")
e.write_geo("ellipse")
# #endregion spacing

# #region split
# a NACA 0012 wing with a sharp trailing edge through a plate: the junction turns 164 deg there
section = naca("0012", n_pts=80, closed_te=True)
wing = extrude_section(section, span=1.0, n_stations=2, span_axis="z", end_caps=True)
wing.apply_translation((0.0, 0.0, -0.3))  # chord along x from 0 to 1, through z = 0
plate = trimesh.creation.box(extents=(4.0, 4.0, 1.0))
plate.apply_translation((0.5, 0.0, -0.5))
plate.export("plate.stl")
wing.export("wing.stl")

options = {
    "seed": (0.5, 0.06, 0.0), "n_loop": 96,
    "ds_surface": 0.005, "growth_surface": 1.12, "m_a": 15, "m_b": 15,
    "ds_wall": 5e-5, "growth_wall": 1.2, "nk": 20,
}  # fmt: skip
w = Project()
w.load_geometry("plate.stl")
w.build_collar("wing.stl", **options)  # one block: folds at the sharp trailing edge
one = w.mesh.blocks[0]
folded = int((scaled_jacobian(np.concatenate([one.nodes, one.nodes[:, :1]], axis=1)) <= 0).sum())
print(f"single block: {folded} folded cells")
w.undo()
w.build_collar("wing.stl", **options, split=True)  # cut at the trailing edge
print("split collar:", len(w.mesh.blocks), "blocks,",
      f"min scaled Jacobian {min_jacobian(w.mesh.blocks):.3f}")  # fmt: skip
print("labels:", [b.label for b in w.mesh.blocks])
w.write_geo("wing")
# #endregion split
