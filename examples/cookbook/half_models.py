"""Half models on a symmetry plane: half bodies, a half-collar, SYMMETRY faces.

Writes ``half_box.stl``, ``half_post.stl``, ``half_collar.geo``,
``half_collar.adjacency.json`` and ``cc.par`` in the current directory.
"""

from __future__ import annotations

import warnings

import numpy as np

from shore import shapes
from shore.io.adjacency import read_adjacency_json
from shore.io.cc_par import write_cc_par
from shore.project import Project

# #region half_bodies
# keep the y >= 0 half of each body and drop the faces lying on the plane y = 0
box = shapes.box(extents=(4.0, 2.0, 1.0))
box.apply_translation((0.0, 1.0, -0.5))
box.update_faces(~np.all(np.isclose(box.vertices[box.faces][..., 1], 0.0), axis=1))
post = shapes.cylinder(radius=0.5, height=1.6, sections=96)
post.apply_translation((0.0, 0.0, 0.2))
post.update_faces(post.triangles_center[:, 1] > 0.0)
for mesh in (box, post):
    mesh.remove_unreferenced_vertices()
box.export("half_box.stl")
post.export("half_post.stl")
# #endregion half_bodies

# #region half_collar
p = Project()
with warnings.catch_warnings():
    warnings.simplefilter("ignore")  # an open half body is not watertight: expected
    p.load_geometry("half_box.stl", skip_hygiene=True)
    p.build_collar(
        "half_post.stl",
        seed=(0.45, 0.2, 0.02),  # a point near the junction
        n_loop=25,
        ds_surface=0.01,
        growth_surface=1.12,
        m_a=12,
        m_b=12,
        ds_wall=1e-4,
        growth_wall=1.2,
        nk=25,
        plane=((0.0, 0.0, 0.0), (0.0, 1.0, 0.0)),  # (point, normal) of the symmetry plane
    )
block = p.mesh.blocks[0]
print("nodes:", block.nodes.shape)
print("j-faces off the plane by", np.abs(block.nodes[:, [0, -1], :, 1]).max())
print({f.name: f.bc.name for f in block.faces()})
p.write_geo("half_collar")
# #endregion half_collar

# #region symmetry_bc
p.write_adjacency("half_collar.adjacency.json")
write_cc_par(
    read_adjacency_json("half_collar.adjacency.json"), "cc.par", grd_filename="half_collar.grd"
)
# #endregion symmetry_bc
