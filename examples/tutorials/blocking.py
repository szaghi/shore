"""Block a half wing-body from scratch, then replay the journal on an edited CAD model.

The phase-7 acceptance case (ADR 0051): a capsule fuselage and a blunt-TE
NACA 0012 wing fused into one STEP model (half model, symmetry plane z = 0),
blocked with a box, splits, an O-grid round the hole, an exact premesh,
orthogonal wall lines and a 2e-3 wall layer. Writes ``wb.step``,
``wing_body.cgns``, ``xall/`` (the Xall run), ``journal.py``, ``edited/wb.step``
and ``wing_body_edited.cgns`` in the current directory. Takes about 10 s.
"""

# requires: cad
# requires: cgns

from __future__ import annotations

import time
from pathlib import Path

import numpy as np

from shore.geometry.base import BODY, EntityKind
from shore.project import Project
from shore.project.project import run_journal
from shore.quality import scaled_jacobian

FACES = ("i_lo", "i_hi", "j_lo", "j_hi", "k_lo", "k_hi")
start = time.perf_counter()


# #region model
def wing_body(path, span=3.0, radius=0.5, length=5.0):
    """Write the half wing-body as STEP: a capsule along x and a NACA 0012 wing along +z."""
    from build123d import (
        Box,
        Cylinder,
        Line,
        Pos,
        Sphere,
        Spline,
        Wire,
        export_step,
        extrude,
        make_face,
    )

    from shore.profiles import naca

    body = Cylinder(radius, length, rotation=(0, 90, 0))
    body = body + Pos(-length / 2, 0, 0) * Sphere(radius) + Pos(length / 2, 0, 0) * Sphere(radius)
    sec = naca("0012", 40, closed_te=False)  # blunt trailing edge, chord 1 from x = 0
    le = int(np.argmin(sec[:, 0]))
    upper = [(float(x), float(y), 0.0) for x, y in sec[: le + 1]]
    lower = [(float(x), float(y), 0.0) for x, y in sec[le:]]
    wire = Wire([Spline(*upper), Spline(*lower), Line(lower[-1], upper[0])])
    wing = extrude(make_face(wire), amount=span)
    # fuse, then keep the half z >= 0: the plane z = 0 is the symmetry plane
    export_step((body + wing) & (Pos(0, 0, 10) * Box(40, 40, 20)), str(path))
    return path


wing_body("wb.step")
# #endregion model

# #region box
p = Project()
p.load_geometry("wb.step", name="wb")
g = p.geometries["wb"].geometry
print("CAD model:", len(g.entities(EntityKind.SURFACE)), "faces,",
      len(g.entities(EntityKind.CURVE)), "edges")  # fmt: skip

# one block of 53 x 33 x 25 nodes over the far field, the symmetry plane at z = 0
p.create_blocking(53, 33, 25, vmin=(-6, -4, 0), vmax=(7, 4, 6))
# #endregion box


# #region splits
def block_at(point):
    """The label of the blocking block that holds *point*."""
    b = p.blocking
    for label in b.labels():
        pts = b.points[list(b.block(label))]
        if np.all(pts.min(axis=0) < point) and np.all(point < pts.max(axis=0)):
            return label
    raise LookupError(f"no block holds {point}")


def cut(axis, c):
    """Split the blocking across *axis* at coordinate *c* (through the whole sheet)."""
    d = "ijk".index(axis)
    point = np.array([-5.99, -3.99, 0.01])  # a corner block: every sheet passes through it
    point[d] = c
    label = block_at(point)
    pts = p.blocking.points[list(p.blocking.block(label))]
    lo, hi = pts.min(axis=0)[d], pts.max(axis=0)[d]
    p.split_blocking(label, axis, (c - lo) / (hi - lo))


for c in (-3.0, -0.3, 1.3, 3.0):  # x: nose, wing LE, wing TE, tail
    cut("i", c)
for c in (-0.8, -0.3, 0.3, 0.8):  # y: fuselage width, wing thickness
    cut("j", c)
for c in (0.8, 3.3):  # z: fuselage top, beyond the wing tip
    cut("k", c)
print("blocks after the splits:", len(p.blocking.labels()))
# #endregion splits

# #region hole
# the 9 blocks round the fuselage and the one round the wing: the hole
centres = [(x, y, 0.4) for x in (-1.5, 0.5, 2.0) for y in (-0.55, 0.0, 0.55)]
hole = [block_at(np.array(c)) for c in [*centres, (0.5, 0.0, 2.0)]]
wing_block = hole[-1]

# faces of the hole's blocks that bound the hole: owned by one hole block only,
# and not on the symmetry plane
b = p.blocking
owners: dict = {}
for label in hole:
    for f in FACES:
        owners.setdefault(b.entity(label, f)[0], []).append((label, f))
walls = {
    key: faces[0]
    for key, faces in owners.items()
    if len(faces) == 1 and not np.allclose(b.points[list(key[1:])][:, 2], 0.0)
}
print("hole:", len(hole), "blocks,", len(walls), "wall faces")
# #endregion hole

# #region associate
for label, f in walls.values():
    p.associate_blocking(label, f, "wb", "body")


def junction_edges(geometry, radius):
    """The CAD edges lying on the fuselage (radius *radius*) next to the wing."""
    out = []
    for k, pts in enumerate(geometry.display().curves):
        on_fuselage = np.allclose(np.hypot(pts[:, 1], pts[:, 2]), radius, atol=1e-3)
        near_root = (np.abs(pts[:, 1]) < 0.07) & (pts[:, 0] > -0.01) & (pts[:, 0] < 1.01)
        at_wing = bool(np.all(near_root))
        if on_fuselage and at_wing:
            out.append(k)
    return out


junction = junction_edges(g, 0.5)
print("junction curves:", junction)
# the 4 edges at the base of the wing's block run round the wing root
for edge in ("j_lo_k_lo", "j_hi_k_lo", "i_lo_k_lo", "i_hi_k_lo"):
    p.associate_blocking(wing_block, edge, "wb", "curve", junction)
# #endregion associate

# #region ogrid
p.delete_blocks(hole)
# the hole's faces, now owned by the blocks round it
layered = [
    (label, f)
    for label in p.blocking.labels()
    for f in FACES
    if p.blocking.entity(label, f)[0] in walls
]
p.ogrid(faces=layered, factor=0.7, nodes=17)
print("blocks after the O-grid:", len(p.blocking.labels()))
# #endregion ogrid

# #region premesh
p.premesh(exact=True)  # tier-2 projection onto the CAD faces and edges
p.orthogonalize(0.5)  # wall-normal lines; 1.0 is refused at the concave junction
p.wall_layer(2e-3)  # first cell 2e-3 at every wall
print(f"premeshed in {time.perf_counter() - start:.0f} s")
# #endregion premesh


# #region quality
def report(project, label):
    """Minimum scaled Jacobian, wall distance to the body, first cells and wall angles."""
    blocks = project.mesh.blocks
    geometry = project.geometries["wb"].geometry
    sj = min(float(scaled_jacobian(blk.nodes).min()) for blk in blocks)
    pts, seg = [], []
    for blk in blocks:
        for face in blk.faces():
            if face.bc.name == "WALL":
                axis, side = "ijk".index(face.name[0]), face.name.endswith("hi")
                wall = np.take(blk.nodes, -1 if side else 0, axis=axis).reshape(-1, 3)
                nxt = np.take(blk.nodes, -2 if side else 1, axis=axis).reshape(-1, 3)
                pts.append(wall)
                seg.append(nxt - wall)
    pts, seg = np.concatenate(pts), np.concatenate(seg)
    proj = geometry.project(pts, BODY, exact=True)
    first = np.linalg.norm(seg, axis=1)
    cos = np.abs(np.einsum("ij,ij->i", seg, proj.normals)) / first
    angle = np.degrees(np.arccos(np.clip(cos, 0.0, 1.0)))
    off = pts[:, 2] > 0.05  # off the symmetry plane
    print(f"{label}: {len(blocks)} blocks, min scaled Jacobian {sj:.3f}, "
          f"wall off the body {proj.distances.max():.1e}, first cell {first.min():.4g}-"
          f"{first.max():.4g}, wall angle median {np.median(angle[off]):.1f} deg")  # fmt: skip


report(p, "original")
# #endregion quality

# #region write
p.save_journal("journal.py")  # the commands so far: the replay below rebuilds the grid
p.write_cgns("wing_body.cgns")  # CGNS/HDF5: a zone per block, 1-to-1 seams, BCs
p.write_xall("xall")  # .geo, .grd, cc.par, proc.input
# #endregion write

# #region replay
Path("edited").mkdir(exist_ok=True)
wing_body("edited/wb.step", span=2.7, radius=0.55)  # shorter wing, fatter fuselage
# the journal reads "wb.step"; input_dir resolves it in edited/
q = run_journal("journal.py", input_dir="edited")
report(q, "edited")
q.write_cgns("wing_body_edited.cgns")
print(f"total {time.perf_counter() - start:.0f} s")
# #endregion replay
