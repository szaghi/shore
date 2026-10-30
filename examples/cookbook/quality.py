"""Measure grid quality: summary over blocks, per-block metrics, the worst cell.

Writes sphere.stl and the six sphere_*.geo blocks of a marched cubed sphere to the
current directory.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from shore import quality, shapes
from shore.io.geo import read_geo
from shore.project import Project

# #region build
shapes.sphere(radius=1.0, subdivisions=4).export("sphere.stl")
p = Project()
p.load_geometry("sphere.stl")
p.build_topology("cubed_sphere", nk=33, ds=1e-3, growth=1.15, ni=33, nj=64,
                 theta_cap=np.deg2rad(35.0))  # fmt: skip
p.march()
p.write_geo("sphere")
blocks = p.mesh.blocks
# #endregion build

# #region summary
for key, value in quality.summary(blocks).items():
    print(f"{key:28s} {value:.4g}")
# #endregion summary

# #region per_block
print(f"{'block':10s} {'min SJ':>7s} {'max wall deg':>12s} {'max AR':>8s} {'max growth':>10s}")
for b in blocks:
    sj = quality.scaled_jacobian(b.nodes)
    wall = quality.k_orthogonality_deg(b.nodes, 0)
    ar = quality.aspect_ratio(b.nodes)
    g = quality.k_growth(b.nodes)
    print(f"{b.label:10s} {sj.min():7.3f} {wall.max():12.2f} {ar.max():8.0f} {g.max():10.3f}")
# #endregion per_block

# #region worst
worst = min(blocks, key=lambda b: quality.scaled_jacobian(b.nodes).min())
sj = quality.scaled_jacobian(worst.nodes)
i, j, k = np.unravel_index(np.argmin(sj), sj.shape)
centre = worst.nodes[i : i + 2, j : j + 2, k : k + 2].reshape(-1, 3).mean(axis=0)
print(f"worst cell: block {worst.label}, cell ({i}, {j}, {k}), SJ {sj[i, j, k]:.3f}")
print(f"  centre {np.round(centre, 4)}, {np.count_nonzero(sj < 0.5)} cells below 0.5")
# #endregion worst

# #region compare
# The same grid marched three ways: compare the wall-normal lines.
for options in ({}, {"smooth_iters": 0}, {"march_scheme": "hyperbolic"}):
    q = Project()
    q.load_geometry("sphere.stl")
    q.build_topology("cubed_sphere", nk=33, ds=1e-3, growth=1.15, ni=33, nj=64,
                     theta_cap=np.deg2rad(35.0))  # fmt: skip
    q.march(**options)
    s = quality.summary(q.mesh.blocks)
    print(
        f"{options!s:30s} wall mean {s['wall_orthogonality_mean_deg']:5.2f} deg, "
        f"max {s['wall_orthogonality_max_deg']:5.2f} deg, min SJ {s['scaled_jacobian_min']:.3f}"
    )
# #endregion compare

# #region geo_files
# The same metrics on .geo files written by SHORE or any other tool.
nodes = [read_geo(path) for path in sorted(Path().glob("sphere_*.geo"))]
print(f"{len(nodes)} .geo files, min SJ {quality.summary(nodes)['scaled_jacobian_min']:.3f}")
# #endregion geo_files
