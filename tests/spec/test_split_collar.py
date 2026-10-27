"""FROZEN acceptance spec: the split collar at sharp corners (ADR 0024).

Frozen 2026-10-27, before any implementation. Implementing work may not edit
this file, ``tests/spec/data/collar_default.npz``, ``tests/spec/__init__.py``,
``tests/conftest.py``, the pytest configuration in ``pyproject.toml``, or
``src/shore/quality.py``. The ONE permitted edit is deleting an ``@PENDING``
line once that case passes (strict xfail turns a passing pending case into a
failure). ``tests/spec/SHA256SUMS`` holds the hashes, taken with ``@PENDING``
lines removed; ``scripts/check_spec.py`` checks them.

Verdict for acceptance (main CI job): ``pytest tests/spec/test_split_collar.py``
reports 33 passed, 0 skipped, 0 xfailed. At freeze: 5 passed (the golden
``build_collar`` blocks and the uniform marches at sharp edges, which already
hold: regression guards), 28 xfailed.

Interface under test
--------------------
``shore.collar.split_collar(body_a, body_b, seed, *, n_loop, steps_a,
steps_b, steps_k, label="collar", params=CollarParams(), loop_spacing=None,
loop_curvature=0.0, wake_length=None) -> list[HexBlock]``:

- the loop's sharp corners are its corners (ADR 0021) that turn by more than
  ``CollarParams.sharp_deg`` (105 deg); without one the result is
  ``[block]``, ``block`` the one ``build_collar`` gives (label *label*);
- otherwise the loop is cut at every sharp corner into pieces; each piece is
  extended at both ends by a wake line on A, from the corner along the
  bisector of the corner's two loop directions, ``wake_length`` long
  (default ``sum(steps_a)``), with cells from the corner's mean loop
  spacing growing by at most ``CollarParams.wake_ratio`` (1.2); over the
  wake the wall is a sheet: the corner's B-line translated along the wake;
- six blocks per piece *p* (pieces in loop order from node 0),
  ``f"{label}_{p}_{part}_{side}"``, *part* in ``wake_lo``, ``body``,
  ``wake_hi`` (along j), *side* ``b`` (i up to the junction row, on B) or
  ``a`` (from it, on A); ``k_lo`` is WALL except on ``wake_lo_b`` /
  ``wake_hi_b``: the sheet, SHARED with the sheet block on the other side of
  the same corner, AxisMap ``(("i", False), ("j", True))``; the seams inside
  a piece are SHARED; all other faces FREE; every SHARED pair matches node
  for node;
- ``march_surface`` (ADR 0008, 0023): at a sharp crease, where the front
  folds back on itself, the node's tangent is lengthened to its arc spacing
  (``MarchParams.sharp_crease``); other marches unchanged (ADR 0023's golden
  file);
- ``ParameterError`` for a bad ``wake_length`` and for a sharp corner that
  is not an edge of B.

Cases
-----
The plane through z-extruded bodies: the NACA 0012 STL wing with a sharp TE
(ADR 0021's ``stl_wing``, turn 164 deg), the analytic biconvex lens (two
tips, turn 157 deg), wedges with a rounded back (frozen below; tips turning
120 to 160 deg). ADR 0010's surface steps ``0.005 * 1.12**j`` (15 layers),
n_loop 96, wall-normal ``5e-5 * 1.2**k``.

- S1 valid (20 layers, 1.9 x the first surface step): every block, every
  cell, scaled Jacobian >= 0.1: the wing and the lens, uniform and
  curvature 1 (prototype 0.155, 0.228, 0.462, 0.447); the lens with 25
  layers (prototype 0.460); the wedges (prototype >= 0.43). The single
  ``build_collar`` block folds on all of them (2-5500 cells, walls
  included).
- S2 layout and seams: block count, labels, face BCs and maps; every SHARED
  pair equal node for node; ``.adjacency.json`` and ``cc.par`` written.
- S3 walls: WALL faces on their body, B-side blocks on B and A-side on A,
  within 1e-9 x the diagonal (STL 1e-6); the body blocks' junction row on
  both bodies.
- S4 wakes: the A-side wake blocks' junction rows on A, cells from the
  corner spacing, adjacent ratio in [1, 1.2], total ``wake_length``
  (default and 0.1); the sheet rows a translation of the corner's B-line.
- S5 marches at sharp edges (``march_surface`` on the wing / lens, periodic,
  from the reversed loop): no flipped cells and in-section drift <= 0.02,
  uniform and curvature 1 (prototype 0.0045, 0.011, 0.0, 0.011; before, with
  curvature: 66 and 2 flips, drift 0.98 and 0.023).
- S6 unchanged: ``data/collar_default.npz``, written by ``build_collar``
  before this ADR (the rectangle prism with curvature, the blunt STL wing,
  the half ellipse on its plane: corners below 105 deg and a half-collar),
  within 1e-12; ``split_collar`` gives the same single blocks.
- S7 errors.

Out of scope of this spec: collars thicker than these at corners, concave
corners and the blunt TE's volume (ADR 0025); half-collars with sharp
corners; sharp corners on A's edges; the project and CLI layers.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from shore.quality import scaled_jacobian
from tests.spec.test_intersection import Plane
from tests.spec.test_loop_corners import RECT, LensPrism, PolygonPrism, stl_wing
from tests.spec.test_loop_spacing import HALF_ELLIPSE_SEED, SYM, EllipticCylinder

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-27; implementation pending")

GOLDEN = Path(__file__).parent / "data" / "collar_default.npz"
STEPS_S = 0.005 * 1.12 ** np.arange(15)
THIN = 5e-5 * 1.2 ** np.arange(20)
MID = 5e-5 * 1.2 ** np.arange(25)
PARTS = ("wake_lo", "body", "wake_hi")


def kernel():
    from shore.collar import split_collar

    return split_collar


def wedge(tip_deg: float, radius: float = 0.1, n_arc: int = 48) -> PolygonPrism:
    """A z-extruded wedge: a tip of *tip_deg* and a back rounded by a 48-gon arc."""
    a = np.radians(tip_deg) / 2.0
    tip = radius / np.sin(a)
    phi = np.linspace(np.pi / 2.0 - a, 3.0 * np.pi / 2.0 + a, n_arc)
    v = np.vstack([[tip, 0.0], np.c_[radius * np.cos(phi), radius * np.sin(phi)]])
    return PolygonPrism(v - [tip / 2.0, 0.0])


def wedge_seed(body: PolygonPrism) -> np.ndarray:
    """A point just off the side from the tip."""
    v0, v1 = body.v[0], body.v[1]
    side = v1 - v0
    normal = np.array([-side[1], side[0]]) / np.linalg.norm(side)
    return np.r_[0.5 * (v0 + v1) + 0.01 * normal, 0.0]


def sharp_case(name: str):
    if name == "wing":
        plate, wing, _ = stl_wing(True)
        return plate, wing, np.array([0.5, 0.06, 0.0]), 1
    return Plane(), LensPrism(), np.array([0.0, 0.06, 0.0]), 2


def collar(name: str, curvature: float = 0.0, steps_k: np.ndarray = THIN, **kw):
    a, b, seed, pieces = sharp_case(name)
    blocks = kernel()(
        a, b, seed, n_loop=96, steps_a=STEPS_S, steps_b=STEPS_S, steps_k=steps_k,
        loop_curvature=curvature, **kw,
    )  # fmt: skip
    return a, b, pieces, blocks


def min_jacobian(blocks) -> float:
    return min(float(scaled_jacobian(blk.nodes).min()) for blk in blocks)


def diagonal(body) -> float:
    return float(np.linalg.norm(np.diff(np.asarray(body.bounds, dtype=float), axis=0)))


def face_nodes(block, name: str) -> np.ndarray:
    axis, end = "ijk".index(name[0]), (0 if name.endswith("_lo") else -1)
    return np.take(block.nodes, end, axis=axis)  # along-axes in (i, j, k) order


def owners(blocks) -> dict[int, object]:
    """The block of every face, by the face's id (faces do not know their block)."""
    return {id(face): blk for blk in blocks for face in blk.faces()}


def partner_view(face, owner: dict[int, object]) -> np.ndarray:
    """The partner face's nodes arranged along this face's along-axes."""
    theirs = [ax for ax in "ijk" if ax != face.partner.name[0]]
    other = face_nodes(owner[id(face.partner)], face.partner.name)
    order = [theirs.index(pair[0]) for pair in face.axis_map.along]
    other = np.transpose(other, (*order, 2))
    for k, (_, flipped) in enumerate(face.axis_map.along):
        if flipped:
            other = np.flip(other, axis=k)
    return other


# ---------------------------------------------------------------------------
# S1: valid collars at sharp corners
# ---------------------------------------------------------------------------


class TestValid:
    @PENDING
    @pytest.mark.parametrize("curvature", [0.0, 1.0], ids=["uniform", "curvature"])
    @pytest.mark.parametrize("name", ["wing", "lens"])
    def test_s1_valid(self, name, curvature):
        assert min_jacobian(collar(name, curvature)[3]) >= 0.1

    @PENDING
    def test_s1_lens_thicker(self):
        assert min_jacobian(collar("lens", 0.0, MID)[3]) >= 0.1

    @PENDING
    @pytest.mark.parametrize("tip", [20, 30, 45, 60])
    def test_s1_wedges(self, tip):
        b = wedge(tip)
        blocks = kernel()(
            Plane(), b, wedge_seed(b), n_loop=96, steps_a=STEPS_S, steps_b=STEPS_S, steps_k=THIN
        )
        assert len(blocks) == 6
        assert min_jacobian(blocks) >= 0.1


# ---------------------------------------------------------------------------
# S2: blocks, boundary conditions and seams
# ---------------------------------------------------------------------------


class TestLayout:
    @PENDING
    @pytest.mark.parametrize("name", ["wing", "lens"])
    def test_s2_blocks_and_faces(self, name):
        from shore.mesh.face import FaceBC

        _, _, pieces, blocks = collar(name)
        labels = [
            f"collar_{p}_{part}_{side}" for p in range(pieces) for part in PARTS for side in "ba"
        ]
        assert sorted(blk.label for blk in blocks) == sorted(labels)
        owner = owners(blocks)
        for blk in blocks:
            part, side = blk.label.split("_", 2)[2].rsplit("_", 1)
            k_lo = blk.face("k_lo")
            if side == "b" and part != "body":
                assert k_lo.bc == FaceBC.SHARED
                other = owner[id(k_lo.partner)].label
                assert "_wake_" in other and other.endswith("_b") and other != blk.label
                assert k_lo.partner.name == "k_lo"
                assert k_lo.axis_map.along == (("i", False), ("j", True))
            else:
                assert k_lo.bc == FaceBC.WALL
            assert blk.face("k_hi").bc == FaceBC.FREE
            assert blk.face("i_lo").bc == (FaceBC.FREE if side == "b" else FaceBC.SHARED)
            assert blk.face("i_hi").bc == (FaceBC.SHARED if side == "b" else FaceBC.FREE)
            assert blk.face("j_lo").bc == (FaceBC.FREE if part == "wake_lo" else FaceBC.SHARED)
            assert blk.face("j_hi").bc == (FaceBC.FREE if part == "wake_hi" else FaceBC.SHARED)

    @PENDING
    @pytest.mark.parametrize("name", ["wing", "lens"])
    def test_s2_seams_match(self, name):
        from shore.mesh.face import FaceBC

        _, _, _, blocks = collar(name)
        owner = owners(blocks)
        seams = 0
        for blk in blocks:
            for face in blk.faces():
                if face.bc == FaceBC.SHARED:
                    seams += 1
                    assert np.array_equal(face_nodes(blk, face.name), partner_view(face, owner))
        assert seams > 0

    @PENDING
    def test_s2_xall_files(self, tmp_path):
        from shore.io.adjacency import write_adjacency_json
        from shore.io.cc_par import write_cc_par

        _, _, _, blocks = collar("lens")
        path = tmp_path / "collar.adjacency.json"
        write_adjacency_json(blocks, path)
        write_cc_par(path, tmp_path / "cc.par", grd_filename="collar.grd")
        assert (tmp_path / "cc.par").stat().st_size > 0


# ---------------------------------------------------------------------------
# S3: walls on the bodies
# ---------------------------------------------------------------------------


class TestWalls:
    @PENDING
    @pytest.mark.parametrize("name", ["wing", "lens"])
    def test_s3_walls_on_bodies(self, name):
        from shore.mesh.face import FaceBC

        a, b, _, blocks = collar(name)
        tol = (1e-6 if name == "wing" else 1e-9) * diagonal(a)
        for blk in blocks:
            if blk.face("k_lo").bc != FaceBC.WALL:
                continue
            body = b if blk.label.endswith("_b") else a
            wall = blk.nodes[:, :, 0].reshape(-1, 3)
            assert float(np.max(body.project(wall).distances)) <= tol
            if "_body_" in blk.label:
                row = blk.nodes[-1 if blk.label.endswith("_b") else 0, :, 0]
                assert float(np.max(a.project(row).distances)) <= tol
                assert float(np.max(b.project(row).distances)) <= tol


# ---------------------------------------------------------------------------
# S4: wakes and sheets
# ---------------------------------------------------------------------------


def wake_rows(blocks):
    """(wake line from the corner outward, sheet rows) of every wake on the A / B side."""
    by_label = {blk.label: blk for blk in blocks}
    out = []
    for label, blk in by_label.items():
        if "_wake_" not in label or not label.endswith("_a"):
            continue
        line = blk.nodes[0, :, 0]
        sheet = by_label[label[:-1] + "b"].nodes[:, :, 0]
        if "_wake_lo_" in label:  # j runs toward the corner
            line, sheet = line[::-1], sheet[:, ::-1]
        out.append((line, sheet))
    return out


class TestWakes:
    @PENDING
    @pytest.mark.parametrize("length", [None, 0.1], ids=["default", "given"])
    def test_s4_wake_lines(self, length):
        a, _, _, blocks = collar("lens", wake_length=length)
        expected = float(STEPS_S.sum()) if length is None else length
        rows = wake_rows(blocks)
        assert len(rows) == 4
        for line, _ in rows:
            assert float(np.max(a.project(line).distances)) <= 1e-9 * diagonal(a)
            cells = np.linalg.norm(np.diff(line, axis=0), axis=1)
            assert abs(cells.sum() - expected) <= 1e-6 * expected
            ratio = cells[1:] / cells[:-1]
            assert ratio.min() >= 1.0 - 1e-9
            assert ratio.max() <= 1.2 + 1e-9

    @PENDING
    @pytest.mark.parametrize("name", ["wing", "lens"])
    def test_s4_sheets(self, name):
        _, _, _, blocks = collar(name)
        for line, sheet in wake_rows(blocks):
            assert np.array_equal(sheet[-1], line)  # the junction row is the wake line
            offset = sheet - sheet[:, :1]
            same = np.broadcast_to(offset[-1], offset.shape)
            np.testing.assert_allclose(offset, same, atol=1e-12)


# ---------------------------------------------------------------------------
# S5: marching at sharp edges
# ---------------------------------------------------------------------------


def sharp_march(name: str, curvature: float):
    from shore.surface.intersect import intersection_loop
    from shore.surface.march import march_surface
    from tests.spec.test_march_creases import flipped

    a, b, seed, _ = sharp_case(name)
    loop = intersection_loop(a, b, seed, 96, curvature=curvature)
    grid = march_surface(b, loop[::-1], STEPS_S, periodic=True)
    drift = float(np.abs(grid[:, :, :2] - loop[::-1][:, None, :2]).max())
    return flipped(b, grid), drift


class TestSharpMarch:
    @pytest.mark.parametrize("name", ["wing", "lens"])
    def test_s5_uniform(self, name):  # already holds: regression guards
        flips, drift = sharp_march(name, 0.0)
        assert flips == 0
        assert drift <= 0.02

    @PENDING
    @pytest.mark.parametrize("name", ["wing", "lens"])
    def test_s5_curvature(self, name):
        flips, drift = sharp_march(name, 1.0)
        assert flips == 0
        assert drift <= 0.02


# ---------------------------------------------------------------------------
# S6: collars without sharp corners unchanged
# ---------------------------------------------------------------------------

G_STEPS_S = 0.01 * 1.2 ** np.arange(8)
G_STEPS_K = 1e-4 * 1.3 ** np.arange(12)


def golden_case(name: str):
    """(bodies, seed, keywords) of the golden collars."""
    if name == "rectangle":
        kw = {"n_loop": 48, "loop_curvature": 1.0}
        return Plane(), PolygonPrism(RECT), np.array([0.0, 0.06, 0.0]), kw
    if name == "blunt_wing":
        plate, wing, _ = stl_wing(False)
        return plate, wing, np.array([0.5, 0.06, 0.0]), {"n_loop": 48}
    kw = {"n_loop": 33, "plane": SYM, "loop_curvature": 1.0}
    return Plane(), EllipticCylinder(), HALF_ELLIPSE_SEED, kw


class TestUnchanged:
    @pytest.mark.parametrize("name", ["rectangle", "blunt_wing", "half_ellipse"])
    def test_s6_build_collar_unchanged(self, name):
        from shore.collar import build_collar

        a, b, seed, kw = golden_case(name)
        steps = {"steps_a": G_STEPS_S, "steps_b": G_STEPS_S, "steps_k": G_STEPS_K}
        block = build_collar(a, b, seed, **steps, **kw)
        np.testing.assert_allclose(block.nodes, np.load(GOLDEN)[name], rtol=0.0, atol=1e-12)

    @PENDING
    @pytest.mark.parametrize("name", ["rectangle", "blunt_wing"])
    def test_s6_split_collar_without_sharp_corners(self, name):
        a, b, seed, kw = golden_case(name)
        blocks = kernel()(
            a, b, seed, steps_a=G_STEPS_S, steps_b=G_STEPS_S, steps_k=G_STEPS_K, label="c", **kw
        )
        assert len(blocks) == 1
        assert blocks[0].label == "c"
        assert blocks[0].periodic_j
        np.testing.assert_allclose(blocks[0].nodes, np.load(GOLDEN)[name], rtol=0.0, atol=1e-12)


# ---------------------------------------------------------------------------
# S7: errors
# ---------------------------------------------------------------------------


class TestErrors:
    @PENDING
    @pytest.mark.parametrize("length", [0.0, -1.0, float("nan")])
    def test_s7_bad_wake_length(self, length):
        from shore.exceptions import ParameterError

        with pytest.raises(ParameterError, match="wake_length"):
            collar("lens", wake_length=length)

    @PENDING
    def test_s7_sharp_corner_on_a(self):
        from shore.exceptions import ParameterError

        with pytest.raises(ParameterError, match="edge of body B"):
            kernel()(
                LensPrism(), Plane(), np.array([0.0, 0.06, 0.0]), n_loop=96,
                steps_a=STEPS_S, steps_b=STEPS_S, steps_k=THIN,
            )  # fmt: skip
