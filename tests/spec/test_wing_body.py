"""FROZEN acceptance spec: a wing-body blocking replayed after a CAD edit (ADR 0051).

Roadmap phase 7's acceptance. Frozen 2026-10-30. No implementation is
pending: the spec checks that the blocking commands of ADR 0042-0050 meet
the acceptance together. It may not be edited (``tests/spec/SHA256SUMS``,
``scripts/check_spec.py``); changing a threshold is a re-freeze.

Verdict for acceptance (CAD job, needs ``shore[cad]``):
``pytest tests/spec/test_wing_body.py`` reports 6 passed, 0 skipped, 0
xfailed. At freeze: 6 passed.

The case
--------
A half model (z >= 0, the plane z = 0 a symmetry plane), a STEP file made
with build123d: a capsule fuselage along x (radius 0.5, cylinder length 5)
and a blunt-TE NACA 0012 wing (chord 1, from x = 0) along +z to the tip at
z = 3, fused. The edit: span 2.7 and fuselage radius 0.55.

The blocking (``build``), every step a journaled command:
- a box [-6, 7] x [-4, 4] x [0, 6] (53 x 33 x 25 nodes), cut at x = -3,
  -0.3, 1.3, 3, y = +-0.8, +-0.3, z = 0.8, 3.3;
- the hole: the 9 blocks round the fuselage (x in [-3, 3], |y| < 0.8,
  z < 0.8) and the block round the wing (x in [-0.3, 1.3], |y| < 0.3,
  z in [0.8, 3.3]); their faces bounding the hole, except on the symmetry
  plane, associated with the body; the 4 edges at the base of the wing's
  block (z = 0.8) associated with the junction (the CAD edges on both the
  fuselage and the wing); the hole's blocks deleted;
- an O-grid layer on the hole's faces (factor 0.7, 17 nodes across), a
  premesh with exact projection, ``orthogonalize(0.5)``, ``wall_layer(2e-3)``.

Wall angle: between a wall node's first grid segment and the CAD normal
(exact projection), at nodes off the symmetry plane, the TE, the tip and the
junction (edges and a concave corner have no single normal).

Cases
-----
- A1 the model: 8 faces; the junction found as CAD edges.
- A2 the blocking: 90 blocks, valid, seams closed, walls on the body, the
  junction ring on the junction, first cells, wall angles.
- A3 the journal replayed: bit for bit.
- A4 the journal replayed on the edited model: valid, on the new body and
  junction, wall angles.
- A5 without the junction associations the ring leaves the junction.
- A6 the GUI journal (``tests_gui/spec/test_gui_wing_body.py``) is the
  same command sequence: its lines appear in this journal.

Out of scope of this spec: full-span models, hyperbolic marching, edits that
leave the hole (a longer span than the wing block).
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("build123d")

FACES = ("i_lo", "i_hi", "j_lo", "j_hi", "k_lo", "k_hi")
ALONG = {"i": ("j", "k"), "j": ("i", "k"), "k": ("i", "j")}
FUSELAGE_BLOCKS = [(x, y, 0.4) for x in (-1.5, 0.5, 2.0) for y in (-0.55, 0.0, 0.55)]
WING_BLOCK = (0.5, 0.0, 2.0)
BASE_EDGES = ("j_lo_k_lo", "j_hi_k_lo", "i_lo_k_lo", "i_hi_k_lo")


def wing_body(path, span=3.0, radius=0.5, length=5.0):
    """Write the half model as STEP to *path*."""
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
    sec = naca("0012", 40, closed_te=False)
    le = int(np.argmin(sec[:, 0]))
    upper = [(float(x), float(y), 0.0) for x, y in sec[: le + 1]]
    lower = [(float(x), float(y), 0.0) for x, y in sec[le:]]
    wire = Wire([Spline(*upper), Spline(*lower), Line(lower[-1], upper[0])])
    wing = extrude(make_face(wire), amount=span)
    export_step((body + wing) & (Pos(0, 0, 10) * Box(40, 40, 20)), str(path))
    return path


def junction_edges(geometry, radius):
    """Indices of the CAD edges lying on the fuselage and on the wing."""
    from shore.geometry.base import EntityKind

    out = []
    for k in range(len(geometry.entities(EntityKind.CURVE))):
        pts = geometry.display().curves[k]
        on_fuselage = np.allclose(np.hypot(pts[:, 1], pts[:, 2]), radius, atol=1e-3)
        at_wing = np.all(np.abs(pts[:, 1]) < 0.07) and np.all(
            (pts[:, 0] > -0.01) & (pts[:, 0] < 1.01)
        )
        if on_fuselage and at_wing:
            out.append(k)
    return out


def block_at(p, point):
    b = p.blocking
    for label in b.labels():
        pts = b.points[list(b.block(label))]
        if np.all(pts.min(axis=0) < point) and np.all(point < pts.max(axis=0)):
            return label
    raise AssertionError(f"no block holds {point}")


def cut(p, axis, c, other):
    d = "ijk".index(axis)
    point = np.array(other, float)
    point[d] = c
    label = block_at(p, point)
    pts = p.blocking.points[list(p.blocking.block(label))]
    a, z = pts.min(axis=0)[d], pts.max(axis=0)[d]
    p.split_blocking(label, axis, (c - a) / (z - a))


def build(step, junction=True):
    """The wing-body blocking of the module docstring, as journaled commands."""
    from shore.project import Project

    p = Project()
    p.load_geometry(step, name="wb")
    p.create_blocking(53, 33, 25, vmin=(-6, -4, 0), vmax=(7, 4, 6))
    corner = (-5.99, -3.99, 0.01)
    for axis, cuts in (("i", (-3.0, -0.3, 1.3, 3.0)), ("j", (-0.8, -0.3, 0.3, 0.8)),
                       ("k", (0.8, 3.3))):  # fmt: skip
        for c in cuts:
            cut(p, axis, c, corner)
    hole = [block_at(p, np.array(c)) for c in [*FUSELAGE_BLOCKS, WING_BLOCK]]
    b = p.blocking
    owners: dict = {}
    for label in hole:
        for f in FACES:
            owners.setdefault(b.entity(label, f)[0], []).append((label, f))
    walls = {k: v[0] for k, v in owners.items()
             if len(v) == 1 and not np.allclose(b.points[list(k[1:])][:, 2], 0.0)}  # fmt: skip
    for label, f in walls.values():
        p.associate_blocking(label, f, "wb", "body")
    if junction:
        je = junction_edges(p.geometries["wb"].geometry, _radius(p))
        for e in BASE_EDGES:
            p.associate_blocking(hole[-1], e, "wb", "curve", je)
    p.delete_blocks(hole)
    layered = [(label, f) for label in p.blocking.labels() for f in FACES
               if p.blocking.entity(label, f)[0] in walls]  # fmt: skip
    p.ogrid(faces=layered, factor=0.7, nodes=17)
    p.premesh(exact=True)
    p.orthogonalize(0.5)
    p.wall_layer(2e-3)
    return p


def _radius(p):
    bounds = p.geometries["wb"].geometry.bounds
    return float(bounds[1][1])


def face_nodes(nodes, face):
    return np.take(nodes, 0 if face.endswith("lo") else -1, axis="ijk".index(face[0]))


def next_nodes(nodes, face):
    return np.take(nodes, 1 if face.endswith("lo") else -2, axis="ijk".index(face[0]))


def seam_gap(blocks):
    owner = {id(f): b for b in blocks for f in b.faces()}
    worst = 0.0
    for b in blocks:
        for f in b.faces():
            if f.bc.name != "SHARED":
                continue
            theirs = face_nodes(owner[id(f.partner)].nodes, f.partner.name)
            (a0, flip0), (_, flip1) = f.axis_map.along
            t = theirs if ALONG[f.partner.name[0]].index(a0) == 0 else theirs.transpose(1, 0, 2)
            t = t[::-1] if flip0 else t
            t = t[:, ::-1] if flip1 else t
            worst = max(worst, float(np.abs(face_nodes(b.nodes, f.name) - t).max()))
    return worst


def min_jacobian(blocks):
    from shore.quality import scaled_jacobian

    return min(float(scaled_jacobian(b.nodes).min()) for b in blocks)


def wall_nodes(blocks):
    """(wall node, first segment) of every WALL face."""
    pts, segs = [], []
    for b in blocks:
        for f in b.faces():
            if f.bc.name == "WALL":
                w = face_nodes(b.nodes, f.name).reshape(-1, 3)
                pts.append(w)
                segs.append(next_nodes(b.nodes, f.name).reshape(-1, 3) - w)
    return np.concatenate(pts), np.concatenate(segs)


def check(p, radius, tip, *, median, p95):
    from shore.geometry.base import BODY

    g = p.geometries["wb"].geometry
    blocks = p.mesh.blocks
    assert min_jacobian(blocks) > 0.1
    assert seam_gap(blocks) == 0.0
    pts, seg = wall_nodes(blocks)
    proj = g.project(pts, BODY, exact=True)
    assert float(proj.distances.max()) < 1e-8
    np.testing.assert_allclose(np.linalg.norm(seg, axis=1), 2e-3, rtol=1e-9)
    cos = np.abs(np.einsum("ij,ij->i", seg, proj.normals)) / np.linalg.norm(seg, axis=1)
    angle = np.degrees(np.arccos(np.clip(cos, 0.0, 1.0)))
    junction = ((np.hypot(pts[:, 1], pts[:, 2]) < radius + 0.06) & (np.abs(pts[:, 1]) < 0.12)
                & (pts[:, 2] > radius - 0.1))  # fmt: skip
    sharp = (np.abs(pts[:, 0] - 1.0) < 0.05) | (pts[:, 2] > tip - 0.05)
    keep = (pts[:, 2] > 0.05) & ~sharp & ~junction
    assert np.median(angle[keep]) < median and np.percentile(angle[keep], 95) < p95


def ring(p):
    """The nodes of the 4 junction-associated edges of the premeshed blocks."""
    out = []
    for b in p.mesh.blocks:
        for e in b.edges():
            if e.association is not None and e.association.kind.value == "curve":
                fixed = {"i": slice(None), "j": slice(None), "k": slice(None)}
                for token in (e.name[:4], e.name[5:]):
                    fixed[token[0]] = 0 if token.endswith("lo") else -1
                out.append(b.nodes[fixed["i"], fixed["j"], fixed["k"]])
    return np.concatenate(out)


def junction_distance(p, points):
    from shore.geometry.base import EntityKind, EntityRef

    g = p.geometries["wb"].geometry
    je = junction_edges(g, _radius(p))
    d = np.stack([g.project(points, EntityRef(EntityKind.CURVE, k), exact=True).distances
                  for k in je])  # fmt: skip
    return d.min(axis=0)


@pytest.fixture(scope="module")
def model(tmp_path_factory):
    d = tmp_path_factory.mktemp("wing_body")
    return wing_body(d / "wb.step")


@pytest.fixture(scope="module")
def built(model):
    return build(model)


class TestWingBody:
    def test_a1_model(self, model):
        from shore.geometry import load_geometry
        from shore.geometry.base import EntityKind

        g = load_geometry(model)
        assert len(g.entities(EntityKind.SURFACE)) == 8
        assert len(junction_edges(g, 0.5)) == 3  # the upper, lower and TE-base junction

    def test_a2_blocking(self, built):
        assert len(built.mesh.blocks) == 90
        check(built, 0.5, 3.0, median=12.0, p95=35.0)
        assert float(junction_distance(built, ring(built)).max()) < 1e-8

    def test_a3_replay(self, built, model, tmp_path):
        from shore.project.project import run_journal

        built.save_journal(tmp_path / "journal.py")
        q = run_journal(tmp_path / "journal.py")
        for a, b in zip(built.mesh.blocks, q.mesh.blocks, strict=True):
            assert a.label == b.label and np.array_equal(a.nodes, b.nodes)

    def test_a4_edited_model(self, built, model, tmp_path):
        from shore.project.project import run_journal

        (tmp_path / "in").mkdir()
        wing_body(tmp_path / "in" / "wb.step", span=2.7, radius=0.55)
        journal = built.journal().replace(repr(str(model)), repr("wb.step"))
        (tmp_path / "journal.py").write_text(journal)
        q = run_journal(tmp_path / "journal.py", input_dir=tmp_path / "in")
        assert len(q.mesh.blocks) == 90
        check(q, 0.55, 2.7, median=12.0, p95=35.0)
        assert float(junction_distance(q, ring(q)).max()) < 1e-8

    def test_a5_without_junction(self, model, built):
        bare = build(model, junction=False)
        pts, _ = wall_nodes(bare.mesh.blocks)
        near = junction_distance(bare, pts)
        corner = (np.abs(pts[:, 1]) < 0.12) & (pts[:, 2] > 0.4) & (pts[:, 2] < 0.9)
        assert float(near[corner].min()) > 1e-3  # nothing lands on the junction
        assert float(junction_distance(built, ring(built)).max()) < 1e-8

    def test_a6_gui_sequence(self, built):
        lines = built.journal().splitlines()
        methods = [ln.split("(")[0] for ln in lines if ln.startswith("p.")]
        assert methods[:2] == ["p.load_geometry", "p.create_blocking"]
        assert methods[-4:] == ["p.ogrid", "p.premesh", "p.orthogonalize", "p.wall_layer"]
        assert methods.count("p.split_blocking") == 10
        assert methods.count("p.associate_blocking") == 25 + 4
        assert methods.count("p.delete_blocks") == 1
