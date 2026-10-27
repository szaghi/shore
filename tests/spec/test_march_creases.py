"""FROZEN acceptance spec: surface marching at creases (ADR 0023).

Frozen 2026-10-27, before any implementation. Implementing work may not edit
this file, ``tests/spec/data/march_default.npz``, ``tests/spec/__init__.py``,
``tests/conftest.py``, the pytest configuration in ``pyproject.toml``, or
``src/shore/quality.py``. The ONE permitted edit is deleting an ``@PENDING``
line once that case passes (strict xfail turns a passing pending case into a
failure). ``tests/spec/SHA256SUMS`` holds the hashes, taken with ``@PENDING``
lines removed; ``scripts/check_spec.py`` checks them.

Verdict for acceptance (main CI job): ``pytest tests/spec/test_march_creases.py``
reports 22 passed, 0 skipped, 0 xfailed. At freeze: 15 passed (the 13 golden
marches and the two prism cases with curvature, which already hold:
regression guards), 7 xfailed.

Interface under test
--------------------
``march_surface`` (ADR 0008, 0012), unchanged in signature: where a front
node's surface normal leans along the front tangent (a corner of the front
on an edge of the surface: ``|n . t| > MarchParams.crease_tolerance``, 0.1),
the marching frame uses the normal made orthogonal to the tangent; elsewhere
the march is unchanged bit for bit.

Cases
-----
The junction of a plane (an STL box plate) with a NACA 0012 wing of blunt
trailing edge (``shore.profiles.naca``, 80 points, extruded along z): the
exact analytic prism of the section (ADR 0021's ``PolygonPrism``) and the
STL (ADR 0021's ``stl_wing``); loop n = 96, uniform and curvature 1
(ADR 0020); surface steps ``0.005 * 1.12**j``, 15 layers.

- C1 no flipped surface cells (cell normal against the body normal) on
  either march, A from the loop and B from the reversed loop (before: 10-11
  uniform, 73 on the STL with curvature).
- C2 on the prismatic wing the march is nearly a translation: in-section
  drift of B's nodes <= 0.08 chord uniform (prototype 0.060; before 0.13)
  and <= 0.02 with curvature (prototype 0.011; before 0.22 on the STL).
- C3 the collar's wall layer (k = 0) of the blunt STL wing, n_loop 96,
  20 wall-normal layers ``5e-5 * 1.2**k``: 0 folded cells (the volume above
  it still folds 6-7 cells at the junction near the TE: out of scope).
- C4 unchanged: ``data/march_default.npz``, written by the marcher before
  this ADR: ADR 0008's seven cases, the elliptic junction of ADR 0020 (both
  bodies, loop and half curve on a plane) and a faceted STL sphere x
  cylinder, within 1e-12.

Out of scope of this spec: sharp trailing edges (the plate march folds at
the 164 deg turn: the split collar), the collar volume at corners
(ADR 0022), marching speed on sliver STL triangles.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import trimesh

from tests.spec.test_intersection import Plane
from tests.spec.test_loop_corners import PolygonPrism, stl_wing

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-27; implementation pending")

GOLDEN = Path(__file__).parent / "data" / "march_default.npz"
STEPS = 0.005 * 1.12 ** np.arange(15)
SEED = np.array([0.5, 0.06, 0.0])


def kernels():
    from shore.surface.intersect import intersection_loop
    from shore.surface.march import march_surface

    return intersection_loop, march_surface


def bodies(kind: str):
    from shore.profiles import naca

    if kind == "prism":
        return Plane(), PolygonPrism(naca("0012", n_pts=80, closed_te=False))
    plate, wing, _ = stl_wing(False)
    return plate, wing


def flipped(body, grid: np.ndarray) -> int:
    """Surface cells whose normal opposes the body normal (sign of the majority)."""
    e1 = np.roll(grid, -1, axis=0) - grid
    e2 = grid[:, 1:] - grid[:, :-1]
    n = np.cross(e1[:, :-1], e2)
    centre = 0.25 * (grid[:, :-1] + np.roll(grid, -1, 0)[:, :-1] + grid[:, 1:]
                     + np.roll(grid, -1, 0)[:, 1:])  # fmt: skip
    normals = np.asarray(body.project(centre.reshape(-1, 3)).normals).reshape(centre.shape)
    o = np.einsum("...k,...k->...", n, normals)
    o *= np.sign(np.median(o))
    return int((o <= 0.0).sum())


def marches(kind: str, curvature: float):
    loop_kernel, march = kernels()
    a, b = bodies(kind)
    loop = loop_kernel(a, b, SEED, 96, curvature=curvature)
    on_a = march(a, loop, STEPS, periodic=True)
    on_b = march(b, loop[::-1], STEPS, periodic=True)
    return a, b, loop, on_a, on_b


# ---------------------------------------------------------------------------
# C1-C2: the marches at the blunt trailing edge
# ---------------------------------------------------------------------------


class TestCreases:
    @pytest.mark.parametrize("kind", ["prism", "stl"])
    @PENDING
    def test_c1_uniform_no_flips(self, kind):
        a, b, _, on_a, on_b = marches(kind, 0.0)
        assert flipped(a, on_a) == 0
        assert flipped(b, on_b) == 0

    @PENDING
    def test_c1_stl_curvature_no_flips(self):
        a, b, _, on_a, on_b = marches("stl", 1.0)
        assert flipped(a, on_a) == 0
        assert flipped(b, on_b) == 0

    def test_c1_prism_curvature_no_flips(self):
        a, b, _, on_a, on_b = marches("prism", 1.0)  # already holds: a regression guard
        assert flipped(a, on_a) == 0
        assert flipped(b, on_b) == 0

    @pytest.mark.parametrize("kind", ["prism", "stl"])
    @PENDING
    def test_c2_uniform_drift(self, kind):
        _, _, loop, _, on_b = marches(kind, 0.0)
        assert np.abs(on_b[:, :, :2] - loop[::-1][:, None, :2]).max() <= 0.08

    @PENDING
    def test_c2_stl_curvature_drift(self):
        _, _, loop, _, on_b = marches("stl", 1.0)
        assert np.abs(on_b[:, :, :2] - loop[::-1][:, None, :2]).max() <= 0.02

    def test_c2_prism_curvature_drift(self):
        _, _, loop, _, on_b = marches("prism", 1.0)  # already holds: a regression guard
        assert np.abs(on_b[:, :, :2] - loop[::-1][:, None, :2]).max() <= 0.02


# ---------------------------------------------------------------------------
# C3: the collar wall
# ---------------------------------------------------------------------------


class TestCollarWall:
    @PENDING
    def test_c3_blunt_wing_collar_wall(self):
        from shore.collar import build_collar
        from shore.quality import scaled_jacobian

        a, b = bodies("stl")
        block = build_collar(
            a, b, SEED, n_loop=96, steps_a=STEPS, steps_b=STEPS,
            steps_k=5e-5 * 1.2 ** np.arange(20),
        )  # fmt: skip
        jac = scaled_jacobian(np.concatenate([block.nodes, block.nodes[:, :1]], axis=1))
        assert (jac[:, :, 0] <= 0.0).sum() == 0


# ---------------------------------------------------------------------------
# C4: smooth marches unchanged
# ---------------------------------------------------------------------------


def _golden_case(name: str):
    from shore.geometry import StlGeometry
    from shore.surface.intersect import intersection_curve
    from tests.spec.test_loop_spacing import (
        ELLIPSE_SEED,
        HALF_ELLIPSE_SEED,
        SYM,
        EllipticCylinder,
    )
    from tests.spec.test_surface_march import run

    loop_kernel, march = kernels()
    if name.startswith("m_"):
        return run(name[2:])[2]
    steps = 0.01 * 1.12 ** np.arange(15)
    a, b = Plane(), EllipticCylinder()
    if name.startswith("ellipse"):
        loop = loop_kernel(a, b, ELLIPSE_SEED, 128, curvature=1.0)
        return (
            march(a, loop, steps, periodic=True)
            if name.endswith("a")
            else march(b, loop[::-1], steps, periodic=True)
        )
    if name.startswith("half"):
        half = intersection_curve(a, b, HALF_ELLIPSE_SEED, 65, plane=SYM, curvature=1.0)
        return (
            march(a, half, steps, plane=SYM)
            if name.endswith("a")
            else march(b, half[::-1], steps, plane=SYM)
        )
    sphere = StlGeometry(trimesh.creation.icosphere(subdivisions=4, radius=1.0))
    cylinder = StlGeometry(trimesh.creation.cylinder(radius=0.4, height=4.0, sections=48))
    loop = loop_kernel(sphere, cylinder, np.array([0.4, 0.0, 0.9]), 48)
    steps = 0.02 * 1.1 ** np.arange(10)
    if name == "stl_sphere":
        return march(sphere, loop, steps, periodic=True)
    return march(cylinder, loop[::-1], steps, periodic=True)


GOLDEN_NAMES = [
    "m_line_plane", "m_ring_cylinder", "m_generator_cylinder", "m_equator_sphere",
    "m_helix_cylinder", "m_concave_l_plane", "m_convex_l_plane", "ellipse_a", "ellipse_b",
    "half_ellipse_a", "half_ellipse_b", "stl_sphere", "stl_cylinder",
]  # fmt: skip


class TestUnchanged:
    @pytest.mark.parametrize("name", GOLDEN_NAMES)
    def test_c4_smooth_marches_unchanged(self, name):
        np.testing.assert_allclose(_golden_case(name), np.load(GOLDEN)[name], rtol=0.0, atol=1e-12)
