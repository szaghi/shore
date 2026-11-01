"""FROZEN acceptance spec: an O-H cross-section for the H-grid tube (ADR 0060).

Frozen 2026-11-01, before any implementation. Implementing work may not edit
this file, ``tests/spec/__init__.py``, ``tests/conftest.py``, the pytest
configuration in ``pyproject.toml`` or ``src/shore/quality.py``. The ONE
permitted edit is deleting an ``@PENDING`` line once that case passes (strict
xfail turns a passing pending case into a failure). ``tests/spec/SHA256SUMS``
holds the hashes, taken with ``@PENDING`` lines removed;
``scripts/check_spec.py`` checks them.

Verdict for acceptance (main CI job): ``pytest tests/spec/test_h_grid_tube_oh.py``
reports 12 passed, 0 skipped, 0 xfailed. At freeze: 2 passed, 10 xfailed.

Interface under test
--------------------
``HGridTubeMesh.build`` (``shore primitive h-grid-tube``, ``build_primitive
("h_grid_tube", ...)``) fills a straight circular tube with five blocks,
``ring_0`` .. ``ring_3`` and ``core``: a square core and four ring blocks
joining it to the wall. No singular axis: every cell is a proper, right-handed
hexahedron whatever the streamwise axis. The only walls are the rings' outer
faces, on the tube's radius; the inflow and outflow faces are FREE; every
other face is SHARED. ``n_circ`` counts the cells round the wall, a multiple
of 4; ``n_radial`` the rings' nodes from the core to the wall; ``k_law="tanh"``
clusters them at the wall.

Cases
-----
- T1 (x 3, streamwise x, y, z) five blocks; every cell's minimum corner
  scaled Jacobian > 0.5.
- T2 (x 3) the WALL faces are exactly the four rings' k_hi, on the radius;
  every i_lo / i_hi FREE; 16 SHARED faces matching through their axis maps.
- T3 the cross-section's area at every station is the regular
  ``n_circ``-gon's.
- T4 the wall has ``n_circ`` nodes per station, evenly spaced in angle (a
  guard: true of the polar tube too).
- T5 the stream runs from 0 to ``length``; ``i_law="tanh"`` clusters at 0
  (a guard).
- T6 ``k_law="tanh"`` clusters the rings' nodes at the wall.
- T7 an ``n_circ`` that is not a multiple of 4 is refused.
- T8 ``shore primitive h-grid-tube`` writes five valid ``.geo`` files and a
  face graph whose walls are the four ring faces.

Out of scope of this spec: non-circular or varying bores, the polar layout
(removed, not kept as an option).
"""

from __future__ import annotations

import json

import numpy as np
import pytest
from typer.testing import CliRunner

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-11-01; implementation pending")

LABELS = ["ring_0", "ring_1", "ring_2", "ring_3", "core"]
ALONG = {"i": ("j", "k"), "j": ("i", "k"), "k": ("i", "j")}
AXIS = {"x": 0, "y": 1, "z": 2}
R, L = 0.5, 2.0


def tube(axis="x", n_circ=16, **kw):
    from shore.primitive import HGridTubeMesh

    params = {"n_axial": 6, "n_radial": 5, **kw}
    return HGridTubeMesh.build(length=L, radius=R, n_circ=n_circ, streamwise_axis=axis, **params)


def face_nodes(nodes, face):
    return np.take(nodes, 0 if face.endswith("lo") else -1, axis="ijk".index(face[0]))


def cross(points, axis):
    """The two cross-section coordinates of *points*."""
    keep = [a for a in range(3) if a != AXIS[axis]]
    return points[..., keep]


@pytest.mark.parametrize("axis", ["x", "y", "z"])
class TestBlocks:
    def test_t1_valid_right_handed(self, axis):
        from shore.quality import scaled_jacobian

        blocks = tube(axis).blocks
        assert [b.label for b in blocks] == LABELS
        for b in blocks:
            assert float(scaled_jacobian(b.nodes).min()) > 0.5, b.label

    def test_t2_faces(self, axis):
        blocks = tube(axis).blocks
        walls = [(b.label, f.name) for b in blocks for f in b.faces() if f.bc.name == "WALL"]
        assert sorted(walls) == [(f"ring_{q}", "k_hi") for q in range(4)]
        for label, name in walls:
            b = next(x for x in blocks if x.label == label)
            r = np.linalg.norm(cross(face_nodes(b.nodes, name), axis), axis=-1)
            np.testing.assert_allclose(r, R, rtol=1e-12)
        owner = {id(f): b for b in blocks for f in b.faces()}
        shared = 0
        for b in blocks:
            assert b.face("i_lo").bc.name == "FREE" and b.face("i_hi").bc.name == "FREE"
            for f in b.faces():
                if f.bc.name != "SHARED":
                    continue
                theirs = face_nodes(owner[id(f.partner)].nodes, f.partner.name)
                (a0, flip0), (_, flip1) = f.axis_map.along
                t = theirs if ALONG[f.partner.name[0]].index(a0) == 0 else theirs.transpose(1, 0, 2)
                t = t[::-1] if flip0 else t
                t = t[:, ::-1] if flip1 else t
                np.testing.assert_array_equal(face_nodes(b.nodes, f.name), t)
                shared += 1
        assert shared == 16


class TestCrossSection:
    def test_t3_area(self):
        n = 16
        blocks = tube("z", n_circ=n).blocks
        polygon = 0.5 * n * R**2 * np.sin(2.0 * np.pi / n)
        for i in range(blocks[0].nodes.shape[0]):
            area = 0.0
            for b in blocks:
                f = b.nodes[i]
                a, bb, c, d = f[:-1, :-1], f[1:, :-1], f[1:, 1:], f[:-1, 1:]
                area += 0.5 * np.linalg.norm(np.cross(c - a, d - bb), axis=-1).sum()
            assert abs(area - polygon) < 1e-12 * polygon

    def test_t4_wall_nodes(self):
        n = 24
        blocks = tube("y", n_circ=n).blocks
        for i in (0, 3):
            pts = np.concatenate([b.nodes[i, :, -1] for b in blocks[:4]])
            ang = np.unique(np.round(np.arctan2(*cross(pts, "y").T[::-1]), 12))
            assert ang.size == n
            np.testing.assert_allclose(np.diff(ang), 2.0 * np.pi / n, rtol=1e-9)

    def test_t5_stream(self):
        b = tube("x", i_law="tanh", i_beta=3.0).blocks[0]
        s = b.nodes[:, 0, 0, 0]
        assert abs(s[0]) < 1e-15 and abs(s[-1] - L) < 1e-12
        steps = np.diff(s)
        assert steps[0] < steps[-1]

    def test_t6_clustered_at_the_wall(self):
        for b in tube("x", k_law="tanh", k_beta=3.0).blocks[:4]:
            line = b.nodes[0, b.nodes.shape[1] // 2]
            steps = np.linalg.norm(np.diff(line, axis=0), axis=-1)
            assert steps[-1] < steps[0], b.label

    def test_t7_n_circ_multiple_of_four(self):
        from shore.exceptions import ParameterError

        with pytest.raises(ParameterError, match="n_circ"):
            tube("x", n_circ=10)


class TestCli:
    def test_t8_cli(self, tmp_path):
        from shore.cli import app
        from shore.io.geo import read_geo
        from shore.quality import scaled_jacobian

        args = [
            "primitive", "h-grid-tube", "--length", "2", "--radius", "0.5", "--n-axial", "6",
            "--n-circ", "16", "--n-radial", "5", "-o", str(tmp_path / "t"),
            "--adjacency-json", str(tmp_path / "t.json"),
        ]  # fmt: skip
        r = CliRunner().invoke(app, args)
        assert r.exit_code == 0, r.output
        files = sorted(tmp_path.glob("*.geo"))
        assert len(files) == 5
        for f in files:
            assert float(scaled_jacobian(read_geo(f)).min()) > 0.5, f.name
        graph = json.loads((tmp_path / "t.json").read_text())
        walls = [
            (blk["label"], name)
            for blk in graph["blocks"]
            for name, face in blk["faces"].items()
            if face["bc"] == "WALL"
        ]
        assert sorted(walls) == [(f"ring_{q}", "k_hi") for q in range(4)]
