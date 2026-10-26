"""FROZEN acceptance spec: the collar block of two intersecting bodies (ADR 0010).

Frozen 2026-10-26, before any implementation. Implementing work may not edit
this file, ``tests/spec/__init__.py``, ``tests/conftest.py``, the pytest
configuration in ``pyproject.toml``, or ``src/shore/quality.py``. The ONE
permitted edit is deleting an ``@PENDING`` line once that case passes (strict
xfail turns a passing pending case into a failure). ``tests/spec/SHA256SUMS``
holds the hashes, taken with ``@PENDING`` lines removed;
``scripts/check_spec.py`` checks them.

Verdict for acceptance (main CI job): ``pytest tests/spec/test_collar.py``
reports 26 passed, 0 skipped, 0 xfailed. At freeze: 0 passed, 26 xfailed (the
kernel and the command do not exist yet).

Interface under test
--------------------
``shore.collar.build_collar(body_a, body_b, seed, *, n_loop, steps_a, steps_b,
steps_k, label="collar")`` returns a ``HexBlock``:

- the wall (``k = 0``) is the two collar halves: with
  ``loop = intersection_loop(body_a, body_b, seed, n_loop)``,
  ``on_a = march_surface(body_a, loop, steps_a, periodic=True)`` and
  ``on_b = march_surface(body_b, loop[::-1], steps_b, periodic=True)[::-1]``,
  the wall is ``concatenate([on_b[:, ::-1], on_a[:, 1:]], axis=1)``,
  transposed to ``(i, j)`` and reversed in j (the accepted kernels of
  ADR 0008 and ADR 0009, composed);
- nodes ``(len(steps_b) + len(steps_a) + 1, n_loop, len(steps_k) + 1, 3)``:
  i runs from B's far edge through the junction row ``i = len(steps_b)`` to
  A's far edge; j runs around the loop and is periodic (``periodic_j=True``,
  the seam node not repeated, the O-grid convention); k leaves the walls;
- face ``k_lo`` is ``WALL``, the others ``FREE``;
- the volume is valid at a collar thickness many times the junction's surface
  spacing (the concave junction is the hard part, ADR 0010), the first layer
  follows ``steps_k[0]``, the outer layer lies at ``sum(steps_k)`` from the
  walls;
- ``ParameterError`` for malformed ``steps_*`` and ``n_loop < 3``.

``Project.build_collar(other, *, seed, n_loop, ds_surface, growth_surface,
m_a, m_b, ds_wall, growth_wall, nk, label="collar")``: body A is the project's
loaded geometry, body B is loaded from the file *other*; the surface steps are
``ds_surface * growth_surface**j`` (``m_a`` / ``m_b`` of them), the wall-normal
steps ``ds_wall * growth_wall**k`` (``nk`` of them). The collar replaces the
project's mesh (one block, like ``build_primitive``), is journaled, undoable
and replays bit-identically.

Cases
-----
Every kernel case uses the exact analytic bodies of the intersection spec
(ADR 0009): plane x cylinder, tilted plane x cylinder (junction angle 63 to
117 deg), cylinder x cylinder and sphere x offset cylinder; ``n_loop = 64``,
15 surface layers ``0.01 * 1.12**j`` on each side, 30 wall-normal layers
``1e-4 * 1.2**k`` (collar thickness 0.118 = 12 x the junction spacing 0.01).

- K1 wall: the composition above (<= 1e-12), rows on their bodies.
- K2 validity: 0 folded cells, min scaled Jacobian >= 0.05 (prototype 0.081
  on the tilted junction, 0.48-0.62 on the others; plain hyperbolic
  extrusion folds 384-1408 cells on every pair).
- K3 first cell: within 0.95-1.05 of ``steps_k[0]`` (prototype 0.984-1.001).
- K4 envelope: outer layer within 1 % of ``sum(steps_k)`` from the walls
  (prototype 0.997-1.002).
- K5 topology: shape, periodicity, label, face BCs.
- K6 interface: determinism, errors.
- C1/C2 command on STL files (a box with a cylindrical post through its top):
  one valid collar block, journal line, replay, undo.

Out of scope of this spec: collars much thicker than 12 x the junction spacing
at junctions sharper than 90 deg (the prototype folds the tilted junction at
32 x), open intersection curves and symmetry planes, spacing laws along the
loop, run time.
"""

from __future__ import annotations

import functools
from pathlib import Path

import numpy as np
import pytest
import trimesh

from shore.exceptions import ParameterError
from shore.mesh.face import FaceBC
from shore.quality import k_spacing, scaled_jacobian
from tests.spec.test_intersection import PAIRS

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-26; implementation pending")

N_LOOP = 64
STEPS_S = 0.01 * 1.12 ** np.arange(15)
STEPS_K = 1e-4 * 1.2 ** np.arange(30)
MB = len(STEPS_S)  # the junction row


def kernel():
    """The function under test (imported late: absent until implemented)."""
    from shore.collar import build_collar

    return build_collar


@functools.cache
def collar(name: str):
    """(body_a, body_b, seed, block) for one pair, built once per session."""
    make_a, make_b, seed, _, _ = PAIRS[name]
    a, b = make_a(), make_b()
    seed = np.asarray(seed, dtype=np.float64)
    block = kernel()(a, b, seed, n_loop=N_LOOP, steps_a=STEPS_S, steps_b=STEPS_S, steps_k=STEPS_K)
    return a, b, seed, block


def expected_wall(a, b, seed) -> np.ndarray:
    """The wall as the frozen composition of the accepted kernels (ADR 0008, 0009)."""
    from shore.surface.intersect import intersection_loop
    from shore.surface.march import march_surface

    loop = intersection_loop(a, b, seed, N_LOOP)
    on_a = march_surface(a, loop, STEPS_S, periodic=True)
    on_b = march_surface(b, loop[::-1], STEPS_S, periodic=True)[::-1]
    wall = np.concatenate([on_b[:, ::-1], on_a[:, 1:]], axis=1)
    return np.transpose(wall, (1, 0, 2))[:, ::-1]


def wrapped_jacobian(nodes: np.ndarray) -> np.ndarray:
    """Scaled Jacobian of every cell, including the periodic j-seam cells."""
    return scaled_jacobian(np.concatenate([nodes, nodes[:, :1]], axis=1))


# ---------------------------------------------------------------------------
# K1-K5: the collar block, per junction
# ---------------------------------------------------------------------------


class TestCollar:
    @pytest.mark.parametrize("name", list(PAIRS))
    def test_k1_wall_is_the_composed_collar_surface(self, name):
        a, b, seed, block = collar(name)
        wall = block.nodes[:, :, 0]
        np.testing.assert_allclose(wall, expected_wall(a, b, seed), atol=1e-12)
        diag = float(np.linalg.norm(np.ptp(a.bounds, axis=0)))
        assert b.distance(wall[: MB + 1]).max() <= 1e-9 * diag
        assert a.distance(wall[MB:]).max() <= 1e-9 * diag

    @pytest.mark.parametrize("name", list(PAIRS))
    def test_k2_valid_volume(self, name):
        *_, block = collar(name)
        sj = wrapped_jacobian(block.nodes)
        assert (sj <= 0.0).sum() == 0
        assert sj.min() >= 0.05

    @pytest.mark.parametrize("name", list(PAIRS))
    def test_k3_first_cell(self, name):
        *_, block = collar(name)
        first = k_spacing(block.nodes)[..., 0] / STEPS_K[0]
        assert first.min() >= 0.95 and first.max() <= 1.05

    @pytest.mark.parametrize("name", list(PAIRS))
    def test_k4_envelope(self, name):
        a, b, _, block = collar(name)
        outer = block.nodes[:, :, -1]
        dist = np.minimum(a.distance(outer), b.distance(outer)) / STEPS_K.sum()
        assert dist.min() >= 0.99 and dist.max() <= 1.01

    @pytest.mark.parametrize("name", list(PAIRS))
    def test_k5_topology(self, name):
        *_, block = collar(name)
        assert block.nodes.shape == (2 * MB + 1, N_LOOP, len(STEPS_K) + 1, 3)
        assert block.periodic_j and not block.periodic_i
        assert block.label == "collar"
        bcs = {f.name: f.bc for f in block.faces()}
        assert bcs["k_lo"] is FaceBC.WALL
        assert all(bcs[f] is FaceBC.FREE for f in ("i_lo", "i_hi", "j_lo", "j_hi", "k_hi"))


# ---------------------------------------------------------------------------
# K6: interface
# ---------------------------------------------------------------------------


class TestInterface:
    def test_k6_deterministic(self):
        a, b, seed, block = collar("sphere_cylinder")
        again = kernel()(
            a, b, seed, n_loop=N_LOOP, steps_a=STEPS_S, steps_b=STEPS_S, steps_k=STEPS_K
        )
        np.testing.assert_array_equal(again.nodes, block.nodes)

    @pytest.mark.parametrize("bad", ["steps_a", "steps_k"])
    def test_k6_bad_steps(self, bad):
        make_a, make_b, seed, _, _ = PAIRS["plane_cylinder"]
        steps = {"steps_a": STEPS_S, "steps_b": STEPS_S, "steps_k": STEPS_K}
        steps[bad] = np.array([0.01, -0.01])
        with pytest.raises(ParameterError):
            kernel()(make_a(), make_b(), np.asarray(seed), n_loop=N_LOOP, **steps)

    def test_k6_too_few_loop_nodes(self):
        make_a, make_b, seed, _, _ = PAIRS["plane_cylinder"]
        with pytest.raises(ParameterError):
            kernel()(
                make_a(), make_b(), np.asarray(seed), n_loop=2,
                steps_a=STEPS_S, steps_b=STEPS_S, steps_k=STEPS_K,
            )  # fmt: skip


# ---------------------------------------------------------------------------
# C1-C2: the Project command on STL files
# ---------------------------------------------------------------------------

COMMAND = {
    "seed": (0.55, 0.0, 0.02), "n_loop": 48, "ds_surface": 0.01, "growth_surface": 1.12,
    "m_a": 12, "m_b": 12, "ds_wall": 1e-4, "growth_wall": 1.2, "nk": 25,
}  # fmt: skip


@pytest.fixture(scope="module")
def stl_files(tmp_path_factory) -> dict[str, Path]:
    """A 4 x 4 x 1 box (top face at z = 0) pierced by a cylindrical post (r = 0.5)."""
    d = tmp_path_factory.mktemp("collar")
    box = trimesh.creation.box(extents=(4.0, 4.0, 1.0))
    box.apply_translation((0.0, 0.0, -0.5))
    post = trimesh.creation.cylinder(radius=0.5, height=1.6, sections=96)
    post.apply_translation((0.0, 0.0, 0.2))
    box.export(str(d / "box.stl"))
    post.export(str(d / "post.stl"))
    return {"box": d / "box.stl", "post": d / "post.stl"}


class TestCommand:
    def test_c1_builds_a_valid_collar_component(self, stl_files):
        from shore.project import Project

        p = Project()
        p.load_geometry(stl_files["box"])
        p.build_collar(str(stl_files["post"]), **COMMAND)
        assert [b.label for b in p.mesh.blocks] == ["collar"]
        nodes = p.mesh.blocks[0].nodes
        assert nodes.shape == (25, 48, 26, 3)
        assert (wrapped_jacobian(nodes) <= 0.0).sum() == 0

    def test_c2_journal_replay_and_undo(self, stl_files, tmp_path):
        from shore.project import Project

        p = Project()
        p.load_geometry(stl_files["box"])
        p.build_collar(str(stl_files["post"]), **COMMAND)
        assert "build_collar(" in p.journal()
        built = p.mesh.blocks[0].nodes.copy()
        p.save(tmp_path / "proj")
        q = Project.load(tmp_path / "proj")  # replays and checks every grid bit for bit
        np.testing.assert_array_equal(q.mesh.blocks[0].nodes, built)
        p.undo()
        assert p.mesh is None
