"""Tutorial 4: the Xall files of the sphere-in-cylinder assembly, by hand and in one call.

Writes sphere.stl, the blocks wall_<label>.geo / background_<label>.geo with
their sidecars, sphere_boxes.json, the per-component wall.{grd,cc.par,proc.input}
and background.{grd,cc.par,proc.input}, the merged assembly.{grd,cc.par,proc.input},
and the same run written by Project.write_xall under run/, in the current directory.
"""

from __future__ import annotations

import filecmp
import json
import math

import shore.shapes
from shore.project import Project

# #region grids
shore.shapes.sphere(radius=0.5, subdivisions=4).export("sphere.stl")
p = Project()
p.use_component("sphere")
p.load_geometry("sphere.stl", name="sphere")
p.build_topology(
    "cubed_sphere", nk=13, ds=0.01, growth=1.1, geometry="sphere",
    ni=30, nj=40, theta_cap=math.radians(30.0),
)  # fmt: skip
p.march()
p.write_geo("wall")
p.write_adjacency("wall.adjacency.json")
p.use_component("background")
p.build_primitive("flat_caps", r_out=10.0, z0=-10.0, z1=50.0, inner_frac=0.4, ni=16, nc=24, nk=60)
p.write_geo("background")
p.write_adjacency("background.adjacency.json")
wall_geo = [f"wall_{b.label}.geo" for b in p.components["sphere"].blocks]
background_geo = [f"background_{b.label}.geo" for b in p.components["background"].blocks]
# #endregion grids

# #region boxes
from shore.boxes import boxes_from_stl  # noqa: E402
from shore.io.boxes import box_to_json  # noqa: E402

boxes = boxes_from_stl(
    "sphere.stl",
    "sub0",  # the block the boxes belong to (any block of the near-body grid)
    max_outward_gap=0.10,  # < the near-body thickness, ~0.21
)
with open("sphere_boxes.json", "w") as f:
    json.dump([box_to_json(b) for b in boxes], f)
print(f"{len(boxes)} marker boxes")
# #endregion boxes

# #region grd
from shore.io.grd import write_grd  # noqa: E402

write_grd("wall.grd", wall_geo)  # the block order of the .grd is the order of the files
write_grd("background.grd", background_geo)
# #endregion grd

# #region farfield
# The flat-caps primitive leaves k_lo (its z0 end) as WALL, a viscous wall in
# cc.par. In a far-field background it is an open boundary: make it FREE.
with open("background.adjacency.json") as f:
    sidecar = json.load(f)
for block in sidecar["blocks"]:
    if block["faces"]["k_lo"]["bc"] == "WALL":
        block["faces"]["k_lo"]["bc"] = "FREE"
with open("background.adjacency.json", "w") as f:
    json.dump(sidecar, f, indent=2)
# #endregion farfield

# #region ccpar
from shore.io.cc_par import BlockMeta, write_cc_par  # noqa: E402

near = {b: BlockMeta(free_face_bc=40) for b in ("sub0", "sub1", "sub2", "sub3")}
near |= {b: BlockMeta(level=2, free_face_bc=40) for b in ("cap_north", "cap_south")}
write_cc_par(
    "wall.adjacency.json", "wall.cc.par", grd_filename="wall.grd",
    block_metadata=near, boxes=boxes,
)  # fmt: skip
far = {b.label: BlockMeta(free_face_bc=9) for b in p.components["background"].blocks}
write_cc_par(
    "background.adjacency.json", "background.cc.par", grd_filename="background.grd",
    block_metadata=far,
)  # fmt: skip
# #endregion ccpar

# #region procinput
from shore.io.proc_input import write_proc_input  # noqa: E402

_, summary = write_proc_input(
    "wall.grd", "wall.proc.input", np=4, body_assignment=dict.fromkeys(range(1, 7), 1)
)
write_proc_input("background.grd", "background.proc.input", np=4)
print(summary)
# #endregion procinput

# #region merge
from shore.io.cc_par import merge_cc_par  # noqa: E402
from shore.io.grd import merge_grd  # noqa: E402
from shore.io.proc_input import merge_proc_input  # noqa: E402

order = ["wall", "background"]  # one order for all three merges
merge_grd([f"{c}.grd" for c in order], "assembly.grd")
merge_cc_par([f"{c}.cc.par" for c in order], "assembly.cc.par", grd_filename="assembly.grd")
_, merged = merge_proc_input(
    [f"{c}.proc.input" for c in order], "assembly.proc.input", grd_path="assembly.grd"
)
print(merged)
# #endregion merge

# #region xall
p.write_xall(
    "run",
    grd="assembly.grd",
    np=4,
    components={
        "sphere": {
            "body": 1,
            "boxes": [box_to_json(b) for b in boxes],
            "meta": {"cap_north": {"level": 2}, "cap_south": {"level": 2}},
        },
        "background": {"free_face_bc": 9},
    },
)
print(p.last_result[1])  # the balance of the one proc.input over all 11 blocks
print("same .grd as the manual merge:", filecmp.cmp("assembly.grd", "run/assembly.grd", False))
# #endregion xall
