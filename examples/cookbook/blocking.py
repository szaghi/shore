"""Blocking recipes: a sphere in a box, from a lattice to a wall layer.

Writes ``sphere.stl`` and ``sphere_blocking.journal.py`` (the replayable
journal of the sphere case) to the current directory.
"""

from __future__ import annotations

import numpy as np

from shore.project import Project
from shore.quality import scaled_jacobian
from shore.shapes import sphere

FACES = ("i_lo", "i_hi", "j_lo", "j_hi", "k_lo", "k_hi")


def min_jacobian(blocks):
    """Smallest scaled Jacobian over every cell of *blocks*."""
    return min(float(scaled_jacobian(b.nodes).min()) for b in blocks)


def wall_angle(blocks):
    """Largest angle (deg) between a wall node's first grid segment and the radius."""
    worst = 0.0
    for b in blocks:
        for f in b.faces():
            if f.bc.name != "WALL":
                continue
            axis, lo = "ijk".index(f.name[0]), f.name.endswith("lo")
            wall = np.take(b.nodes, 0 if lo else -1, axis=axis)
            seg = np.take(b.nodes, 1 if lo else -2, axis=axis) - wall
            cos = np.einsum("...i,...i", seg, wall) / (
                np.linalg.norm(seg, axis=-1) * np.linalg.norm(wall, axis=-1)
            )
            worst = max(worst, float(np.degrees(np.arccos(np.clip(cos, -1.0, 1.0))).max()))
    return worst


def block_at(blocking, point):
    """Label of the blocking block whose skeleton box holds *point*."""
    for label in blocking.labels():
        pts = blocking.points[list(blocking.block(label))]
        if np.all(pts.min(axis=0) < point) and np.all(point < pts.max(axis=0)):
            return label
    raise LookupError(f"no block holds {point}")


sphere(radius=1.0, subdivisions=3).export("sphere.stl")

# #region split
p = Project()
p.load_geometry("sphere.stl")
p.create_blocking(13, 13, 13, vmin=(-3, -3, -3), vmax=(3, 3, 3))
print(p.blocking.shape("b0"))  # (13, 13, 13): nodes along i, j, k

# Cut at x = -1 and x = +1, then the same along j and k: a 3 x 3 x 3 lattice.
corner = np.full(3, -2.99)
for d, axis in enumerate("ijk"):
    for c in (-1.0, 1.0):
        point = corner.copy()
        point[d] = c
        label = block_at(p.blocking, point)
        pts = p.blocking.points[list(p.blocking.block(label))]
        lo, hi = pts[:, d].min(), pts[:, d].max()
        p.split_blocking(label, axis, (c - lo) / (hi - lo))
print(len(p.blocking.labels()), "blocks")
# #endregion split

# #region hole
centre = block_at(p.blocking, np.zeros(3))
for face in FACES:
    p.associate_blocking(centre, face)  # geometry "body", kind "body"
hole_keys = {p.blocking.entity(centre, face)[0] for face in FACES}
p.delete_blocks([centre])
p.premesh()
blocks = p.mesh.blocks
walls = sum(f.bc.name == "WALL" for b in blocks for f in b.faces())
print(f"{len(blocks)} blocks, {walls} WALL faces")
print(f"min scaled Jacobian {min_jacobian(blocks):.3f}, wall angle {wall_angle(blocks):.1f} deg")
# #endregion hole

# #region ogrid
hole = [
    (label, face)
    for label in p.blocking.labels()
    for face in FACES
    if p.blocking.entity(label, face)[0] in hole_keys
]
p.ogrid(faces=hole, factor=0.5, nodes=17)
p.premesh()
blocks = p.mesh.blocks
print(f"{len(blocks)} blocks")
print(f"min scaled Jacobian {min_jacobian(blocks):.3f}, wall angle {wall_angle(blocks):.1f} deg")
# #endregion ogrid

# #region orthogonalize
p.orthogonalize()  # strength=1.0
blocks = p.mesh.blocks
print(f"min scaled Jacobian {min_jacobian(blocks):.3f}, wall angle {wall_angle(blocks):.1f} deg")
# #endregion orthogonalize

# #region wall-layer
p.wall_layer(1e-3)
b = next(b for b in p.mesh.blocks if any(f.bc.name == "WALL" for f in b.faces()))
wall = next(f.name for f in b.faces() if f.bc.name == "WALL")
axis = "ijk".index(wall[0])
x = np.moveaxis(b.nodes, axis, 0)
line = x[:, x.shape[1] // 2, x.shape[2] // 2]  # one line across the wall
line = line if wall.endswith("lo") else line[::-1]
steps = np.linalg.norm(np.diff(line, axis=0), axis=1)
print(f"{b.label}: first cell {steps[0]:.2e}, ratio {steps[1] / steps[0]:.3f}")
print(f"min scaled Jacobian {min_jacobian(p.mesh.blocks):.3f}")
p.save_journal("sphere_blocking.journal.py")
# #endregion wall-layer

# #region ogrid-box
q = Project()
q.create_blocking(9, 9, 9, vmin=(-1, -1, -1), vmax=(1, 1, 1))
q.ogrid(["b0"], factor=0.5, nodes=4)  # a core b0 and six layer blocks b1..b6
print(q.blocking.labels(), q.blocking.shape("b1"))
q.undo()
q.ogrid(["b0"], exclude=[("b0", "k_lo")], factor=0.5, nodes=4)  # half-O, open on k_lo
print(q.blocking.labels())
# #endregion ogrid-box

# #region move
r = Project()
r.create_blocking(9, 5, 5, vmin=(0, 0, 0), vmax=(2, 1, 1))
r.split_blocking("b0", "i", 0.5)
r.move_vertex("b0", "v7", (1.2, 1.3, 1.1))  # shared with b1: both follow
print(r.blocking.points[r.blocking.block("b1")[3]])

s = Project()
s.load_geometry("sphere.stl")
s.create_blocking(5, 5, 5, vmin=(-0.5, -0.5, -0.5), vmax=(0.5, 0.5, 0.5))
s.associate_blocking("b0", "k_hi")  # v7 inherits the face's association
s.move_vertex("b0", "v7", (0.3, 0.9, 0.9))
moved = s.blocking.points[s.blocking.block("b0")[7]]
print(moved, f"|x| = {np.linalg.norm(moved):.4f}")
# #endregion move
