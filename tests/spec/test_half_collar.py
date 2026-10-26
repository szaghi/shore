"""FROZEN acceptance spec: the half-collar block on a symmetry plane (ADR 0013).

Frozen 2026-10-26, before any implementation. Implementing work may not edit
this file, ``tests/spec/__init__.py``, ``tests/conftest.py``, the pytest
configuration in ``pyproject.toml``, or ``src/shore/quality.py``. The ONE
permitted edit is deleting an ``@PENDING`` line once that case passes (strict
xfail turns a passing pending case into a failure). ``tests/spec/SHA256SUMS``
holds the hashes, taken with ``@PENDING`` lines removed;
``scripts/check_spec.py`` checks them.

Verdict for acceptance (main CI job): ``pytest tests/spec/test_half_collar.py``
reports 23 passed, 0 skipped, 0 xfailed. At freeze: 0 passed, 23 xfailed (the
``plane`` argument, ``extrude_collar``, ``FaceBC.SYMMETRY`` and the command
option do not exist yet).

Interface under test
--------------------
``shore.collar.build_collar(body_a, body_b, seed, *, n_loop, steps_a, steps_b,
steps_k, label="collar", plane=None)``: with ``plane = (point, normal)`` the
half-collar of a half model:

- the wall (``k = 0``) is the two halves marched from the open curve: with
  ``curve = intersection_curve(body_a, body_b, seed, n_loop, plane=plane)``,
  ``on_a = march_surface(body_a, curve, steps_a, plane=plane)`` and
  ``on_b = march_surface(body_b, curve[::-1], steps_b, plane=plane)[::-1]``,
  the wall is ``concatenate([on_b[:, ::-1], on_a[:, 1:]], axis=1)``,
  transposed to ``(i, j)`` and reversed in j (as ADR 0010);
- nodes ``(len(steps_b) + len(steps_a) + 1, n_loop, len(steps_k) + 1, 3)``,
  ``periodic_j`` false; faces ``j_lo`` and ``j_hi`` lie on the plane and are
  ``FaceBC.SYMMETRY``, ``k_lo`` is ``WALL``, the others ``FREE``;
- on mirror-symmetric bodies the block is the y >= 0 half of the collar
  volume extruded from the mirrored closed wall:
  ``shore.collar.extrude_collar(wall, steps_k, junction)`` is the ADR 0010
  volume stage on a j-periodic wall (``junction`` the junction row);
- a ``UserWarning`` when the curve meets the plane obliquely (its end
  tangent more than 1 deg off the plane normal; ADR 0011 review decision 2):
  the block's faces are then still on the plane, but its quality is not
  guaranteed;
- ``ParameterError`` for a *plane* that is not a point and a non-zero
  normal.

``FaceBC.SYMMETRY`` is written to ``.adjacency.json`` and read back, and
``write_cc_par`` emits it as the natural-BC code ``2`` (symmetry) with
family 0. ``Project.build_collar(other, ..., plane=None)`` passes the plane
on, journals it, and replays bit-identically.

Cases (n_loop = 33 on the half models of ADR 0011; 15 surface layers
``0.01 * 1.12**j`` on each side, 30 wall-normal layers ``1e-4 * 1.2**k``,
as ADR 0010):

- H1 wall: the composition above (<= 1e-12).
- H2 mirror: equal to ``extrude_collar`` of the mirrored wall, first
  n_loop stations, within 1e-9 x the collar thickness (prototype: the same
  construction).
- H3 validity at the plane: j-faces on the plane within 1e-9 x diagonal;
  every node at y >= 0; the j-lines meet the plane at a right angle, their
  second-order one-sided tangent within 0.5 deg of the plane normal
  (prototype <= 0.12 deg); 0 folded cells, min scaled Jacobian >= 0.05,
  first cell 0.95-1.05 of ``steps_k[0]``, outer layer within 1 % of the
  thickness (the ADR 0010 bounds; prototype 0.081-0.660, 0.984-1.001,
  0.997-1.002).
- H4 topology: shape, periodicity, face BCs.
- H5 oblique plane (x = 0.5, plane x cylinder): warns; j-faces on the plane.
- H6 interface: ``extrude_collar`` reproduces the accepted whole collar
  (ADR 0010) bit for bit; a bad plane raises.
- H7 output: SYMMETRY in the adjacency sidecar and as cc.par code 2.
- C1/C2 command on the STL half model (half box + half post, ADR 0011):
  one valid block with SYMMETRY j-faces; journal, replay.

Out of scope of this spec: two planes, mirror-asymmetric half models'
quality, the STL vertex-normal symmetry (ADR 0012, separate decision),
run time.
"""

from __future__ import annotations

import functools
import warnings
from pathlib import Path

import numpy as np
import pytest
import trimesh

from shore.exceptions import ParameterError
from shore.quality import k_spacing, scaled_jacobian
from tests.spec.test_intersection import CylinderZ, Plane
from tests.spec.test_open_curve import HALF, SYMMETRY

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-26; implementation pending")

N_LOOP = 33
STEPS_S = 0.01 * 1.12 ** np.arange(15)
STEPS_K = 1e-4 * 1.2 ** np.arange(30)
MB = len(STEPS_S)  # the junction row
MIRROR = np.diag([1.0, -1.0, 1.0])  # reflection across y = 0


def kernel():
    """The function under test (imported late)."""
    from shore.collar import build_collar

    return build_collar


@functools.cache
def half_collar(name: str):
    """(body_a, body_b, seed, block) for one half model, built once per session."""
    make_a, make_b, seed, _, _ = HALF[name]
    a, b = make_a(), make_b()
    seed = np.asarray(seed, dtype=np.float64)
    block = kernel()(
        a, b, seed, n_loop=N_LOOP, steps_a=STEPS_S, steps_b=STEPS_S, steps_k=STEPS_K,
        plane=SYMMETRY,
    )  # fmt: skip
    return a, b, seed, block


def expected_wall(a, b, seed) -> np.ndarray:
    from shore.surface.intersect import intersection_curve
    from shore.surface.march import march_surface

    curve = intersection_curve(a, b, seed, N_LOOP, plane=SYMMETRY)
    on_a = march_surface(a, curve, STEPS_S, plane=SYMMETRY)
    on_b = march_surface(b, curve[::-1], STEPS_S, plane=SYMMETRY)[::-1]
    wall = np.concatenate([on_b[:, ::-1], on_a[:, 1:]], axis=1)
    return np.transpose(wall, (1, 0, 2))[:, ::-1]


def j_angle_deg(nodes: np.ndarray) -> float:
    """Largest angle between the plane normal (y) and the second-order one-sided
    j-tangent at either j-face, over all (i, k)."""
    worst = 0.0
    for r0, r1, r2 in (
        (nodes[:, 0], nodes[:, 1], nodes[:, 2]),
        (nodes[:, -1], nodes[:, -2], nodes[:, -3]),
    ):
        t = -3.0 * r0 + 4.0 * r1 - r2
        t /= np.linalg.norm(t, axis=-1, keepdims=True)
        worst = max(worst, float(np.degrees(np.arccos(np.clip(np.abs(t[..., 1]), 0, 1))).max()))
    return worst


# ---------------------------------------------------------------------------
# H1-H4: the half-collar, per half model
# ---------------------------------------------------------------------------


class TestHalfCollar:
    @PENDING
    @pytest.mark.parametrize("name", list(HALF))
    def test_h1_wall_is_the_composed_half_surface(self, name):
        a, b, seed, block = half_collar(name)
        np.testing.assert_allclose(block.nodes[:, :, 0], expected_wall(a, b, seed), atol=1e-12)

    @PENDING
    @pytest.mark.parametrize("name", list(HALF))
    def test_h2_half_of_the_mirrored_collar(self, name):
        from shore.collar import extrude_collar

        *_, block = half_collar(name)
        wall = block.nodes[:, :, 0]
        closed = np.concatenate([wall, wall[:, -2:0:-1] @ MIRROR], axis=1)
        whole = extrude_collar(closed, STEPS_K, MB)
        assert whole.shape == (2 * MB + 1, 2 * (N_LOOP - 1), len(STEPS_K) + 1, 3)
        assert np.abs(block.nodes - whole[:, :N_LOOP]).max() <= 1e-9 * STEPS_K.sum()

    @PENDING
    @pytest.mark.parametrize("name", list(HALF))
    def test_h3_valid_at_the_plane(self, name):
        a, b, _, block = half_collar(name)
        nodes = block.nodes
        diag = float(np.linalg.norm(np.ptp(np.asarray(a.bounds), axis=0)))
        assert np.abs(nodes[:, [0, -1], :, 1]).max() <= 1e-9 * diag
        assert nodes[..., 1].min() >= -1e-9 * diag
        assert j_angle_deg(nodes) <= 0.5
        sj = scaled_jacobian(nodes)
        assert (sj <= 0.0).sum() == 0 and sj.min() >= 0.05
        first = k_spacing(nodes)[..., 0] / STEPS_K[0]
        assert first.min() >= 0.95 and first.max() <= 1.05
        outer = nodes[:, :, -1]
        dist = np.minimum(a.body.distance(outer), b.body.distance(outer)) / STEPS_K.sum()
        assert dist.min() >= 0.99 and dist.max() <= 1.01

    @PENDING
    @pytest.mark.parametrize("name", list(HALF))
    def test_h4_topology(self, name):
        from shore.mesh.face import FaceBC

        *_, block = half_collar(name)
        assert block.nodes.shape == (2 * MB + 1, N_LOOP, len(STEPS_K) + 1, 3)
        assert not block.periodic_j and not block.periodic_i
        bcs = {f.name: f.bc for f in block.faces()}
        assert bcs["j_lo"] is FaceBC.SYMMETRY and bcs["j_hi"] is FaceBC.SYMMETRY
        assert bcs["k_lo"] is FaceBC.WALL
        assert all(bcs[f] is FaceBC.FREE for f in ("i_lo", "i_hi", "k_hi"))


# ---------------------------------------------------------------------------
# H5: oblique plane
# ---------------------------------------------------------------------------


class TestObliquePlane:
    @PENDING
    def test_h5_oblique_plane_warns(self):
        plane = ((0.5, 0.0, 0.0), (1.0, 0.0, 0.0))
        with pytest.warns(UserWarning):
            block = kernel()(
                Plane(), CylinderZ(), np.array([1.1, 0.1, 0.05]), n_loop=N_LOOP,
                steps_a=STEPS_S, steps_b=STEPS_S, steps_k=STEPS_K, plane=plane,
            )  # fmt: skip
        assert np.abs(block.nodes[:, [0, -1], :, 0] - 0.5).max() <= 1e-9


# ---------------------------------------------------------------------------
# H6: interface
# ---------------------------------------------------------------------------


class TestInterface:
    @PENDING
    def test_h6_extrude_collar_is_the_accepted_volume_stage(self):
        from shore.collar import extrude_collar
        from tests.spec.test_collar import STEPS_K as K
        from tests.spec.test_collar import collar

        *_, block = collar("sphere_cylinder")
        np.testing.assert_array_equal(extrude_collar(block.nodes[:, :, 0], K, MB), block.nodes)

    @PENDING
    def test_h6_bad_plane(self):
        make_a, make_b, seed, _, _ = HALF["plane_cylinder"]
        with pytest.raises(ParameterError):
            kernel()(
                make_a(), make_b(), np.asarray(seed), n_loop=N_LOOP, steps_a=STEPS_S,
                steps_b=STEPS_S, steps_k=STEPS_K, plane=((0.0, 0.0, 0.0), (0.0, 0.0, 0.0)),
            )  # fmt: skip


# ---------------------------------------------------------------------------
# H7: output
# ---------------------------------------------------------------------------


class TestOutput:
    @PENDING
    def test_h7_adjacency_round_trip(self, tmp_path):
        from shore.io.adjacency import read_adjacency_json, write_adjacency_json

        *_, block = half_collar("plane_cylinder")
        write_adjacency_json([block], tmp_path / "a.json")
        faces = read_adjacency_json(tmp_path / "a.json")["blocks"][0]["faces"]
        assert faces["j_lo"]["bc"] == "SYMMETRY" and faces["j_hi"]["bc"] == "SYMMETRY"

    @PENDING
    def test_h7_cc_par_symmetry_code(self, tmp_path):
        from shore.io.adjacency import read_adjacency_json, write_adjacency_json
        from shore.io.cc_par import write_cc_par

        *_, block = half_collar("plane_cylinder")
        write_adjacency_json([block], tmp_path / "a.json")
        write_cc_par(read_adjacency_json(tmp_path / "a.json"), tmp_path / "cc.par",
                     grd_filename="c.grd")  # fmt: skip
        patches = {}
        for line in (tmp_path / "cc.par").read_text().splitlines():
            tokens = line.split("!")[0].split()
            if len(tokens) == 10:  # block face bc family imin imax jmin jmax kmin kmax
                patches[int(tokens[1])] = (int(tokens[2]), int(tokens[3]))
        assert patches[3] == (2, 0) and patches[4] == (2, 0)  # j_lo, j_hi: symmetry
        assert patches[5] == (1, 0)  # k_lo: wall


# ---------------------------------------------------------------------------
# C1-C2: the Project command on the STL half model
# ---------------------------------------------------------------------------

COMMAND = {
    "seed": (0.45, 0.2, 0.02), "n_loop": 25, "ds_surface": 0.01, "growth_surface": 1.12,
    "m_a": 12, "m_b": 12, "ds_wall": 1e-4, "growth_wall": 1.2, "nk": 25,
    "plane": ((0.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
}  # fmt: skip


@pytest.fixture(scope="module")
def stl_half(tmp_path_factory) -> dict[str, Path]:
    """The ADR 0011 half box and half post, cut at y = 0."""
    d = tmp_path_factory.mktemp("half_collar")
    box = trimesh.creation.box(extents=(4.0, 2.0, 1.0))
    box.apply_translation((0.0, 1.0, -0.5))
    box.update_faces(~np.all(np.isclose(box.vertices[box.faces][..., 1], 0.0), axis=1))
    post = trimesh.creation.cylinder(radius=0.5, height=1.6, sections=96)
    post.apply_translation((0.0, 0.0, 0.2))
    post.update_faces(post.triangles_center[:, 1] > 0.0)
    for mesh in (box, post):
        mesh.remove_unreferenced_vertices()
    box.export(str(d / "box.stl"))
    post.export(str(d / "post.stl"))
    return {"box": d / "box.stl", "post": d / "post.stl"}


class TestCommand:
    @PENDING
    def test_c1_builds_a_valid_half_collar(self, stl_half):
        from shore.mesh.face import FaceBC
        from shore.project import Project

        p = Project()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # open half-model STL files are not watertight
            p.load_geometry(stl_half["box"], skip_hygiene=True)
            p.build_collar(str(stl_half["post"]), **COMMAND)
        block = p.mesh.blocks[0]
        assert block.nodes.shape == (25, 25, 26, 3)
        assert (scaled_jacobian(block.nodes) <= 0.0).sum() == 0
        assert np.abs(block.nodes[:, [0, -1], :, 1]).max() <= 1e-9
        bcs = {f.name: f.bc for f in block.faces()}
        assert bcs["j_lo"] is FaceBC.SYMMETRY and bcs["j_hi"] is FaceBC.SYMMETRY

    @PENDING
    def test_c2_journal_and_replay(self, stl_half, tmp_path):
        from shore.project import Project

        p = Project()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            p.load_geometry(stl_half["box"], skip_hygiene=True)
            p.build_collar(str(stl_half["post"]), **COMMAND)
            assert "plane=" in p.journal()
            built = p.mesh.blocks[0].nodes.copy()
            p.save(tmp_path / "proj")
            q = Project.load(tmp_path / "proj")
        np.testing.assert_array_equal(q.mesh.blocks[0].nodes, built)
