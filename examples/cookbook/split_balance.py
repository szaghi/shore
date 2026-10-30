"""Plan a balanced split of a cubed sphere for 10 MPI ranks and apply it.

Writes sphere.stl, sphere_*.geo, sphere.adjacency.json, sphere.grd, splits.toml,
split_*.geo, split.grd and split.proc.input to the current directory.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from shore import shapes
from shore.balance import plan_balance
from shore.project import Project

# #region mesh
shapes.sphere(radius=1.0, subdivisions=4).export("sphere.stl")
p = Project()
p.load_geometry("sphere.stl")
p.build_topology("cubed_sphere", nk=33, ds=1e-3, growth=1.15, ni=33, nj=64,
                 theta_cap=np.deg2rad(35.0))  # fmt: skip
p.march()
p.write_geo("sphere")
p.write_adjacency("sphere.adjacency.json")
labels = [b.label for b in p.mesh.blocks]
p.write_grd("sphere.grd", [f"sphere_{label}.geo" for label in labels])
# #endregion mesh

# #region plan
plan = plan_balance(
    "sphere.grd", np=10, ngr=3, min_coarse_cells=4, adjacency_path="sphere.adjacency.json"
)
Path("splits.toml").write_text(plan.to_toml())
before, after = plan.initial_summary, plan.final_summary
print(f"{before.n_blocks} blocks, imbalance {before.imbalance:.0%} -> "
      f"{after.n_blocks} blocks, imbalance {after.imbalance:.0%}")  # fmt: skip
print("cuts:", plan.splits[:3], "...")
# #endregion plan

# #region apply
p.split([(label, axis, [at]) for label, axis, at in plan.splits])  # cuts follow the seams
p.write_geo("split", per_block=True)
p.write_adjacency("split.adjacency.json")
chunks = [b.label for b in p.mesh.blocks]
p.write_grd("split.grd", [f"split_{label}.geo" for label in chunks])
p.write_proc_input("split.grd", "split.proc.input", np=10)
print(f"{len(chunks)} blocks:", ", ".join(chunks[:4]), "...")
# #endregion apply
