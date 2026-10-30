"""Redistribute the nodes along a block edge, through the seams (ADR 0041).

Writes ``sphere.stl`` and ``sphere_redistributed.geo`` (one ``.geo`` per block)
in the current directory.
"""

from __future__ import annotations

import numpy as np

from shore import shapes
from shore.exceptions import ParameterError
from shore.project import Project

# #region one_block
p = Project()
p.build_primitive("box", vmin=(0.0, 0.0, 0.0), vmax=(2.0, 1.0, 1.0), ni=9, nj=5, nk=5)

# cluster the i-lines towards i = 0: the law is taken over the edge j_lo_k_lo
p.redistribute("box", "j_lo_k_lo", "tanh", beta=3.0)
print("tanh  :", np.round(p.mesh.blocks[0].nodes[:, 0, 0, 0], 4))

# or hold the first and last cells (lengths along the edge) around a uniform interior
p.undo()
p.redistribute("box", "j_lo_k_lo", "pinned", first_step=0.05, last_step=0.4)
x = p.mesh.blocks[0].nodes[:, 0, 0, 0]
print("pinned:", np.round(x, 4), "first", x[1] - x[0], "last", x[-1] - x[-2])
# #endregion one_block

# #region through_seams
shapes.sphere(radius=1.0, subdivisions=3).export("sphere.stl")

s = Project()
s.load_geometry("sphere.stl")
s.build_topology(
    "cubed_sphere", nk=5, ds=0.05, growth=1.1, ni=33, nj=64, theta_cap=np.deg2rad(35.0)
)
s.march()
before = {b.label: b.nodes.copy() for b in s.mesh.blocks}

# along i of sub0: the four equatorial blocks follow through their SHARED j-seams
s.redistribute("sub0", "j_lo_k_lo", "tanh2", beta=2.0)
for b in s.mesh.blocks:
    print(f"{b.label:10s} changed: {not np.array_equal(b.nodes, before[b.label])}")

# the wall nodes moved along their old polylines: put them back on the body
s.snap()
s.write_geo("sphere_redistributed")
# #endregion through_seams

# #region refusals
journal = p.journal()
for args, kwargs in [
    (("box", "k_lo", "tanh"), {}),  # a face, not an edge
    (("box", "j_lo_k_lo", "spline"), {}),  # an unknown law
    (("box", "j_lo_k_lo", "pinned"), {"first_step": 1.5, "last_step": 1.5}),  # longer than the edge
]:
    try:
        p.redistribute(*args, **kwargs)
    except ParameterError as err:
        print("refused:", err)
assert p.journal() == journal  # nothing recorded, nothing moved

# #endregion refusals
