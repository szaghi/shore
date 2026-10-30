"""Override the spacing of single block edges: i/j on a box, k on a cubed sphere and an O-grid.

Writes sphere.stl to the current directory.
"""

from __future__ import annotations

import numpy as np

from shore import shapes
from shore.mesh import CubedSphereTopology, Mesh, OGridTopology
from shore.mesh.spacing import Spacing1D
from shore.primitive import BoxMesh

K_EDGES = ("i_lo_j_lo", "i_hi_j_lo", "i_lo_j_hi", "i_hi_j_hi")


def steps(line: np.ndarray) -> np.ndarray:
    return np.linalg.norm(np.diff(line, axis=0), axis=-1)


# #region box
corners = {"vmin": (0.0, 0.0, 0.0), "vmax": (10.0, 5.0, 3.0), "ni": 41, "nj": 31, "nk": 21}
baseline = BoxMesh.build(**corners)
pinned = BoxMesh.build(**corners, edge_pins={"j_lo_k_lo": Spacing1D(law="tanh", beta=4.0)})
nodes = pinned.blocks[0].nodes
for label, j, k in (("pinned edge (j=0, k=0)", 0, 0), ("opposite edge", -1, -1), ("mid", 15, 10)):
    print(f"{label:22s} first i-step {steps(nodes[:, j, k])[0]:.3g}")
print(f"baseline               first i-step {steps(baseline.blocks[0].nodes[:, 0, 0])[0]:.3g}")
# #endregion box

# #region cubed_sphere
shapes.sphere(radius=1.0, subdivisions=4).export("sphere.stl")
mesh = Mesh.from_stl("sphere.stl", CubedSphereTopology(), ni=33, nj=64, nk=33, theta_cap_deg=35.0)
sp_k = Spacing1D(law="tanh", beta=4.0, total_thickness=2.0)  # one instance for all
for block in mesh.blocks:
    for name in K_EDGES:
        block.edge(name).spacing = sp_k
mesh.march()  # the tanh schedule replaces the (ds, growth) defaults
s = steps(mesh.blocks[0].nodes[16, 8, :])
print(f"cubed sphere: {len(mesh.blocks)} blocks, steps {s[0]:.2e}..{s[-1]:.2e}, sum {s.sum():.3f}")
# #endregion cubed_sphere

# #region ogrid
mesh = Mesh.from_stl("sphere.stl", OGridTopology(), ni=33, nj=49, nk=33, theta_cap_deg=10.0)
sp_k = Spacing1D(law="tanh2", beta=3.0, total_thickness=2.0)
for name in K_EDGES:
    mesh.blocks[0].edge(name).spacing = sp_k
mesh.march()
s = steps(mesh.blocks[0].nodes[16, 0, :])
print(f"O-grid: steps {s[0]:.2e}..{s[-1]:.2e}, sum {s.sum():.3f}")
# #endregion ogrid
