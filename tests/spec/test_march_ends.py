"""FROZEN acceptance spec: surface marching with ends on a symmetry plane (ADR 0012).

Frozen 2026-10-26, before any implementation. Implementing work may not edit
this file, ``tests/spec/__init__.py``, ``tests/conftest.py``, the pytest
configuration in ``pyproject.toml``, or ``src/shore/quality.py``. The ONE
permitted edit is deleting an ``@PENDING`` line once that case passes (strict
xfail turns a passing pending case into a failure). ``tests/spec/SHA256SUMS``
holds the hashes, taken with ``@PENDING`` lines removed;
``scripts/check_spec.py`` checks them.

Verdict for acceptance (main CI job): ``pytest tests/spec/test_march_ends.py``
reports 17 passed, 0 skipped, 0 xfailed. At freeze: 0 passed, 17 xfailed
(``march_surface`` has no ``plane`` argument yet).

Interface under test
--------------------
``shore.surface.march.march_surface(geometry, curve, steps, *, periodic=False,
entity=BODY, plane=None)``: with ``plane = (point, normal)`` an open *curve*
whose two ends lie on the plane (from ``intersection_curve``, ADR 0011) is
marched with both ends held on the plane by a symmetry condition: each end
node is solved as an interior node whose missing neighbour is the mirror
image, across the plane, of its inner neighbour; then it is put on the body
and on the plane. ``plane=None`` is the accepted free-end march (ADR 0008).
``ParameterError`` for a *plane* that is not a point and a non-zero normal,
*plane* with ``periodic=True``, and curve ends off the plane by more than
1e-6 x diagonal (the marcher's on-surface tolerance).

Cases (curves from ``intersection_curve`` with n = 33; 15 layers
``0.01 * 1.12**j``, the collar schedule of ADR 0010):

- E1 half models (the four pairs of ADR 0011, cut at y = 0), marched on A
  from the curve and on B from the reversed curve: every node on the body
  and the end nodes on the plane within 1e-9 x diagonal; every node at
  y >= 0; the grid line along the curve meets the plane at a right angle,
  its second-order one-sided tangent at the end within 0.5 deg of the plane
  normal (prototype <= 0.12 deg; marching freely and snapping the end onto
  the plane: 2.3-5.0 deg); equal to the y >= 0 half of the accepted
  periodic march of the mirrored closed curve on the whole body within
  1e-9 x the first step (prototype 3e-14; snapping: 0.6-1.6).
- E2 a plane that is not a symmetry plane (x = 0.5, plane x cylinder): ends
  on the plane and nodes on the body within 1e-9 x diagonal, every cell
  valid.
- E3 STL half model (half box + half post, ADR 0011 O3): ends on the plane
  and nodes on the meshes within 1e-9 x diagonal, every node at y >= 0,
  every cell valid. (No right-angle or mirror requirement: interpolated STL
  vertex normals on the half mesh are not the mirror image of the whole
  mesh's.)
- E4 interface: determinism, errors.

A cell is valid when its signed area along the surface normal at its first
corner is positive.

Out of scope of this spec: ends on two different planes, one end free and
one on a plane, curves meeting the plane obliquely on mirror-symmetric
bodies, the half-collar block, run time.
"""

from __future__ import annotations

import functools

import numpy as np
import pytest
import trimesh

from shore.exceptions import ParameterError
from shore.geometry import BODY
from tests.spec.test_intersection import CylinderZ, Plane
from tests.spec.test_open_curve import HALF, SYMMETRY

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-26; implementation pending")

N = 33
STEPS = 0.01 * 1.12 ** np.arange(15)
MIRROR = np.diag([1.0, -1.0, 1.0])  # reflection across y = 0


def kernel():
    """The function under test."""
    from shore.surface.march import march_surface

    return march_surface


@functools.cache
def half_curve(name: str) -> np.ndarray:
    from shore.surface.intersect import intersection_curve

    make_a, make_b, seed, _, _ = HALF[name]
    return intersection_curve(make_a(), make_b(), np.asarray(seed, dtype=np.float64), N,
                              plane=SYMMETRY)  # fmt: skip


def diagonal(body) -> float:
    return float(np.linalg.norm(np.ptp(np.asarray(body.bounds), axis=0)))


def end_angle_deg(grid: np.ndarray, normal: np.ndarray) -> float:
    """Largest angle between the plane normal and the second-order one-sided
    tangent of the along-curve grid line at either end, over all layers."""
    worst = 0.0
    for r0, r1, r2 in ((grid[0], grid[1], grid[2]), (grid[-1], grid[-2], grid[-3])):
        t = -3.0 * r0 + 4.0 * r1 - r2
        t /= np.linalg.norm(t, axis=-1, keepdims=True)
        worst = max(worst, float(np.degrees(np.arccos(np.clip(np.abs(t @ normal), 0, 1))).max()))
    return worst


def min_signed_area(grid: np.ndarray, body) -> float:
    e1 = grid[1:, :-1] - grid[:-1, :-1]
    e2 = grid[:-1, 1:] - grid[:-1, :-1]
    normals = body.project(grid[:-1, :-1].reshape(-1, 3)).normals.reshape(e1.shape)
    return float(np.einsum("...j,...j->...", np.cross(e1, e2), normals).min())


def plane_offset(points: np.ndarray, plane) -> np.ndarray:
    point, normal = (np.asarray(v, dtype=np.float64) for v in plane)
    return (points - point) @ (normal / np.linalg.norm(normal))


# ---------------------------------------------------------------------------
# E1: half models
# ---------------------------------------------------------------------------


class TestHalfModel:
    @PENDING
    @pytest.mark.parametrize("side", ["A", "B"])
    @pytest.mark.parametrize("name", list(HALF))
    def test_e1_symmetry_plane(self, name, side):
        make_a, make_b, *_ = HALF[name]
        curve = half_curve(name)
        body, start = (make_a(), curve) if side == "A" else (make_b(), curve[::-1])
        grid = kernel()(body, start, STEPS, plane=SYMMETRY)
        diag = diagonal(body)
        assert grid.shape == (N, len(STEPS) + 1, 3)
        assert body.distance(grid.reshape(-1, 3)).max() <= 1e-9 * diag
        assert np.abs(grid[[0, -1], :, 1]).max() <= 1e-9 * diag
        assert grid[..., 1].min() >= -1e-9 * diag
        assert end_angle_deg(grid, np.array([0.0, 1.0, 0.0])) <= 0.5
        loop = np.vstack([start, start[-2:0:-1] @ MIRROR])  # the mirrored closed curve
        whole = kernel()(body.body, loop, STEPS, periodic=True)[:N]
        assert np.abs(grid - whole).max() <= 1e-9 * STEPS[0]


# ---------------------------------------------------------------------------
# E2: a plane that is not a symmetry plane
# ---------------------------------------------------------------------------

OBLIQUE = ((0.5, 0.0, 0.0), (1.0, 0.0, 0.0))


class TestObliquePlane:
    @PENDING
    @pytest.mark.parametrize("side", ["A", "B"])
    def test_e2_ends_on_an_oblique_plane(self, side):
        from shore.surface.intersect import intersection_curve

        a, b = Plane(), CylinderZ()
        curve = intersection_curve(a, b, np.array([1.1, 0.1, 0.05]), N, plane=OBLIQUE)
        body, start = (a, curve) if side == "A" else (b, curve[::-1])
        grid = kernel()(body, start, STEPS, plane=OBLIQUE)
        diag = diagonal(body)
        assert body.distance(grid.reshape(-1, 3)).max() <= 1e-9 * diag
        assert np.abs(plane_offset(grid[[0, -1]], OBLIQUE)).max() <= 1e-9 * diag
        assert min_signed_area(grid, body) > 0.0


# ---------------------------------------------------------------------------
# E3: STL half model
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def stl_half():
    from shore.geometry import StlGeometry
    from shore.surface.intersect import intersection_curve

    box = trimesh.creation.box(extents=(4.0, 2.0, 1.0))
    box.apply_translation((0.0, 1.0, -0.5))
    box.update_faces(~np.all(np.isclose(box.vertices[box.faces][..., 1], 0.0), axis=1))
    post = trimesh.creation.cylinder(radius=0.5, height=1.6, sections=96)
    post.apply_translation((0.0, 0.0, 0.2))
    post.update_faces(post.triangles_center[:, 1] > 0.0)
    for mesh in (box, post):
        mesh.remove_unreferenced_vertices()
    a, b = StlGeometry(box), StlGeometry(post)
    return a, b, intersection_curve(a, b, np.array([0.45, 0.2, 0.02]), N, plane=SYMMETRY)


class TestStlHalfModel:
    @PENDING
    @pytest.mark.parametrize("side", ["box", "post"])
    def test_e3_stl_half_model(self, stl_half, side):
        a, b, curve = stl_half
        body, start, steps = (a, curve, STEPS) if side == "box" else (b, curve[::-1], STEPS[:12])
        grid = kernel()(body, start, steps, plane=SYMMETRY)
        diag = diagonal(a)
        assert body.project(grid.reshape(-1, 3)).distances.max() <= 1e-9 * diag
        assert np.abs(grid[[0, -1], :, 1]).max() <= 1e-9 * diag
        assert grid[..., 1].min() >= -1e-9 * diag
        assert min_signed_area(grid, body) > 0.0


# ---------------------------------------------------------------------------
# E4: interface
# ---------------------------------------------------------------------------


class TestInterface:
    @PENDING
    def test_e4_deterministic(self):
        make_a, *_ = HALF["sphere_cylinder"]
        curve = half_curve("sphere_cylinder")
        first = kernel()(make_a(), curve, STEPS, plane=SYMMETRY)
        again = kernel()(make_a(), curve, STEPS, plane=SYMMETRY, entity=BODY)
        np.testing.assert_array_equal(first, again)

    @PENDING
    def test_e4_periodic_with_plane(self):
        make_a, *_ = HALF["plane_cylinder"]
        with pytest.raises(ParameterError):
            kernel()(make_a(), half_curve("plane_cylinder"), STEPS, periodic=True, plane=SYMMETRY)

    @PENDING
    def test_e4_ends_off_the_plane(self):
        make_a, *_ = HALF["plane_cylinder"]
        curve = half_curve("plane_cylinder")[1:]  # starts one node off the plane
        with pytest.raises(ParameterError):
            kernel()(make_a(), curve, STEPS, plane=SYMMETRY)

    @PENDING
    @pytest.mark.parametrize(
        "plane", [((0.0, 0.0, 0.0), (0.0, 0.0, 0.0)), ((0.0, 0.0), (0.0, 1.0))]
    )
    def test_e4_bad_plane(self, plane):
        make_a, *_ = HALF["plane_cylinder"]
        with pytest.raises(ParameterError):
            kernel()(make_a(), half_curve("plane_cylinder"), STEPS, plane=plane)
