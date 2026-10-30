"""Write one Xall run from two components: a cubed sphere and a cylindrical background.

Writes sphere.stl, boxes.json, xall.toml, the project case.shore/ and the run
directory run/ to the current directory.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from shore import shapes
from shore.boxes import boxes_from_stl
from shore.io.boxes import box_to_json
from shore.project import Project

# #region components
shapes.sphere(radius=1.0, subdivisions=4).export("sphere.stl")
p = Project()
p.load_geometry("sphere.stl")
p.build_topology("cubed_sphere", nk=17, ds=1e-2, growth=1.15, ni=33, nj=64,
                 theta_cap=np.deg2rad(35.0))  # fmt: skip
p.march()  # outer shell near r = 1.55
p.use_component("background")
p.build_primitive("flat_caps", r_out=4.0, z0=-4.0, z1=4.0, ni=9, nc=17, nk=33)
p.save("case.shore")
print({name: len(mesh.blocks) for name, mesh in p.components.items()})
# #endregion components

# #region options
# Marker boxes that cut the hole in the background (hole cutting is Xall's job).
boxes = boxes_from_stl("sphere.stl", "sub0", max_outward_gap=0.2)
Path("boxes.json").write_text(json.dumps([box_to_json(b) for b in boxes]))
Path("xall.toml").write_text("""\
[run]
grd = "run.grd"
np = 8

[components.main]
body = 1
boxes = "boxes.json"

[components.background]
free_face_bc = 9
""")
print(f"{len(boxes)} marker boxes")
# #endregion options

# #region write
p.write_xall(
    "run",
    grd="run.grd",
    np=8,
    components={
        "main": {"body": 1, "boxes": [box_to_json(b) for b in boxes]},
        "background": {"free_face_bc": 9},
    },
)
print(sorted(str(f) for f in Path("run").iterdir()))
# #endregion write
