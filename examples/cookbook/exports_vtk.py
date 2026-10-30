"""Export a marched grid to VTK for ParaView, with the cell Jacobian.

Writes sphere.stl and sphere.vtm (with its block files) to the current directory.
"""

# requires: vis

from __future__ import annotations

import numpy as np
import pyvista as pv

from shore import shapes
from shore.project import Project

# #region vtk
shapes.sphere(radius=1.0, subdivisions=4).export("sphere.stl")
p = Project()
p.load_geometry("sphere.stl")
p.build_topology("cubed_sphere", nk=17, ds=1e-2, growth=1.15, ni=33, nj=64,
                 theta_cap=np.deg2rad(35.0))  # fmt: skip
p.march()
p.mesh.write_vtk("sphere")  # sphere.vtm: one block per HexBlock, cell data "Jacobian"
grid = pv.read("sphere.vtm")
print([(name, block.dimensions) for name, block in zip(grid.keys(), grid, strict=True)])
# #endregion vtk
