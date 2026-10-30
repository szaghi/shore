"""Tutorial 1: a cubed-sphere grid around a sphere, end to end.

Writes sphere.stl, the six blocks sphere_<label>.geo and the journal
sphere_journal.py to the current directory.
"""

from __future__ import annotations

import math

import numpy as np

# #region body
import shore.shapes

body = shore.shapes.sphere(radius=0.5, subdivisions=4)  # a trimesh.Trimesh
body.export("sphere.stl")
# #endregion body

# #region mesh
from shore.project import Project  # noqa: E402

p = Project()
p.load_geometry("sphere.stl")
p.build_topology(
    "cubed_sphere",
    nk=25,  # wall-normal nodes, the wall included
    ds=0.002,  # first cell height
    growth=1.15,  # geometric ratio between successive cells
    ni=33,  # equator-block i nodes, pole to pole
    nj=64,  # nodes round the equator ring: 64 // 4 = 16 cells per block
    theta_cap=math.radians(35.0),  # polar cap half-angle
)
p.march()
p.write_geo("sphere")  # sphere_sub0..3.geo, sphere_cap_north.geo, sphere_cap_south.geo
# #endregion mesh


# #region quality
from shore import quality  # noqa: E402


def report(blocks) -> None:
    """Print the quality summary and the range of the first cell height."""
    for key, value in quality.summary(blocks).items():
        print(f"{key:28s} {value:.4g}")
    first = np.concatenate([quality.k_spacing(b.nodes)[:, :, 0].ravel() for b in blocks])
    print(f"{'first_cell_min':28s} {first.min():.4g}")
    print(f"{'first_cell_max':28s} {first.max():.4g}")


for block in p.mesh.blocks:
    print(f"{block.label:10s} nodes {block.nodes.shape[:3]}")
report(p.mesh.blocks)
# #endregion quality

# #region hyperbolic
p.undo()  # reverts the march (and the write_geo after it)
p.march(march_scheme="hyperbolic")
p.write_geo("sphere")
report(p.mesh.blocks)
# #endregion hyperbolic

# #region journal
p.save_journal("sphere_journal.py")
print(p.journal())
# #endregion journal
