"""Export a two-component project as CGNS and Plot3D.

Writes sphere.stl, the project directory case.shore/, case.cgns, sphere.cgns and
case.xyz to the current directory.
"""

# requires: cgns

from __future__ import annotations

import h5py
import numpy as np

from shore import shapes
from shore.project import Project

# #region project
shapes.sphere(radius=1.0, subdivisions=4).export("sphere.stl")
p = Project()
p.load_geometry("sphere.stl")
p.build_topology("cubed_sphere", nk=17, ds=1e-2, growth=1.15, ni=33, nj=64,
                 theta_cap=np.deg2rad(35.0))  # fmt: skip
p.march()
p.use_component("background")
p.build_primitive("box", vmin=(-4.0, -4.0, -4.0), vmax=(4.0, 4.0, 4.0), ni=33, nj=33, nk=33)
p.save("case.shore")
# #endregion project

# #region cgns
p.write_cgns("case.cgns")  # every component with a mesh: a base each
p.write_cgns("sphere.cgns", components=["main"])
with h5py.File("case.cgns", "r") as f:
    for base in ("main", "background"):
        zones = [f[base][k] for k in f[base] if f[base][k].attrs.get("label") == b"Zone_t"]
        conn = sum(len(z.get("ZoneGridConnectivity", ())) for z in zones)
        bcs = sum(len(z.get("ZoneBC", ())) for z in zones)
        print(f"base {base}: {len(zones)} zones, {conn} 1-to-1 connections, {bcs} BCs")
# #endregion cgns

# #region plot3d
p.write_plot3d("case.xyz")
with open("case.xyz", "rb") as f:
    raw = f.read()
nblocks = int(np.frombuffer(raw, "<i4", 1, 4)[0])
dims = np.frombuffer(raw, "<i4", 3 * nblocks, 16).reshape(nblocks, 3)
print(f"case.xyz: {nblocks} blocks, dims {dims.tolist()}")
# #endregion plot3d
