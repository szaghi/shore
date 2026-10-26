"""FROZEN acceptance spec: open intersection curves ending on a symmetry plane (ADR 0011).

Frozen 2026-10-26, before any implementation; RE-FROZEN 2026-10-26 (user review,
ADR 0011): O3 also bounds the run time (5 s), so a trace that stalls at the
plane fails fast. Implementing work may not edit
this file, ``tests/spec/__init__.py``, ``tests/conftest.py``, the pytest
configuration in ``pyproject.toml``, or ``src/shore/quality.py``. The ONE
permitted edit is deleting an ``@PENDING`` line once that case passes (strict
xfail turns a passing pending case into a failure). ``tests/spec/SHA256SUMS``
holds the hashes, taken with ``@PENDING`` lines removed;
``scripts/check_spec.py`` checks them.

Verdict for acceptance (main CI job): ``pytest tests/spec/test_open_curve.py``
reports 20 passed, 0 skipped, 0 xfailed. At freeze: 0 passed, 20 xfailed (the
kernel does not exist yet).

Interface under test
--------------------
``shore.surface.intersect.intersection_curve(body_a, body_b, seed, n, *,
plane, entity_a=BODY, entity_b=BODY)`` returns ``(n, 3)``:

- the intersection curve of the two bodies through the junction nearest
  *seed*, on the side of ``plane = (point, normal)`` where *seed* lies, from
  the plane back to the plane; *n* >= 3 nodes equally spaced in arc length;
- ``curve[0]`` and ``curve[-1]`` on both bodies and on the plane (the
  curve's end points), every node on both bodies and on the seed's side;
- ordered along ``t = n_b x n_a`` (the bodies' outward normals), as
  ``intersection_loop`` (ADR 0009): ``march_surface(body_a, curve, ...)``
  marches on A away from B and ``march_surface(body_b, curve[::-1], ...)``
  on B away from A;
- the bodies may end at the plane (a half model): tracing must not need
  either body beyond it;
- computed with ``Geometry.project`` (with normals) only;
- ``ParameterError`` for n < 3, a *seed* that is not a 3-vector, a *plane*
  that is not a point and a non-zero normal, a *seed* on the plane, no
  intersection near *seed*, and a curve that closes without meeting the
  plane (within 5 s).

Cases (n = 33; the bodies of ADR 0009's spec):

- O1 half models: the four pairs of ADR 0009, each cut at the symmetry
  plane y = 0 (``Half``: the y >= 0 part of the body, exact projection
  onto it), seeds at y > 0: every node on both bodies within
  1e-9 x diagonal (prototype 7e-13), the ends also on the plane, every node
  at y >= 0, chord lengths max / min <= 1.01 (prototype 1.0003), polygon
  length within 0.5 % of the exact half curve (prototype -0.044 %), every
  node on the seed's branch.
- O2 whole bodies: the side of the seed is the one kept (plane x cylinder,
  seed at y < 0); a plane that is not a symmetry plane (x = 0.5, the arc
  of length 2 pi / 3).
- O3 STL half model: a half box and a half post (96 sections) cut at y = 0:
  on both meshes within 1e-9 x diagonal, ends at (+-0.5, 0, 0), within 5 s
  (prototype 0.13 s; a strict crossing test stalls on the cut edge).
- O4 orientation, each half model: one ``march_surface`` layer on A from the
  curve lies outside B; on B, from the reversed curve, outside A.
- O5 interface: determinism, errors.

Out of scope of this spec: ends on anything but a plane (a third body, open
body edges), curves tangent to the plane at their ends, curves meeting the
plane more than twice on the seed's side, the surface marcher's end
constraint (the marched layer's ends are free; ADR 0011 measures their drift
off the plane), spacing laws other than uniform, run time.
"""

from __future__ import annotations

import time

import numpy as np
import pytest
import trimesh

from shore.exceptions import ParameterError
from shore.geometry import BODY, Projection
from tests.spec.test_intersection import (
    CylinderX,
    CylinderZ,
    Plane,
    Sphere,
    _length,
    exact_ellipse,
    exact_sphere_cylinder,
    exact_steinmetz,
)

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-26; implementation pending")

N = 33
SYMMETRY = ((0.0, 0.0, 0.0), (0.0, 1.0, 0.0))


def kernel():
    """The function under test (imported late: absent until implemented)."""
    from shore.surface.intersect import intersection_curve

    return intersection_curve


# ---------------------------------------------------------------------------
# Half bodies (frozen with the spec)
# ---------------------------------------------------------------------------


class Half:
    """The y >= 0 half of a body symmetric about y = 0 (a half model).

    Exact projection: the whole body's foot when it lies at y >= 0, else the
    nearest point of the cut curve, which for these bodies is the foot of
    the point dropped onto y = 0.
    """

    def __init__(self, body) -> None:
        self.body = body
        self.bounds = body.bounds

    def entities(self, kind):
        return self.body.entities(kind)

    def foot(self, p: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        f, n = self.body.foot(p)
        cut = f[..., 1] < 0.0
        if cut.any():
            q = p[cut] * [1.0, 0.0, 1.0]
            fq, nq = self.body.foot(q)
            f[cut], n[cut] = fq * [1.0, 0.0, 1.0], nq
        return f, n

    def project(self, points, entity=None, *, exact=False) -> Projection:
        pts = np.asarray(points, dtype=np.float64)
        foot, normal = self.foot(pts)
        return Projection(foot, np.linalg.norm(pts - foot, axis=-1), normal, None, None)

    def distance(self, p: np.ndarray) -> np.ndarray:
        return np.linalg.norm(p - self.foot(p)[0], axis=-1)

    def side(self, p: np.ndarray) -> np.ndarray:
        """Signed distance to the whole body, > 0 outside."""
        return self.body.side(p)


def _all(p):
    return np.ones(len(p), bool)


# name: (body A, body B, seed, exact length, on the seed's branch)
HALF = {
    "plane_cylinder": (lambda: Half(Plane()), lambda: Half(CylinderZ()), [1.1, 0.1, 0.05],
                       lambda: np.pi, _all),
    "tilted_plane_cylinder": (lambda: Half(Plane((-0.5, 0.0, 1.0))), lambda: Half(CylinderZ()),
                              [1.0, 0.1, 0.6], lambda: 0.5 * _length(exact_ellipse()), _all),
    "cylinder_cylinder": (lambda: Half(CylinderZ()), lambda: Half(CylinderX()),
                          [1.05, 0.05, 0.45], lambda: 0.5 * _length(exact_steinmetz()),
                          lambda p: p[:, 0] > 0.0),
    "sphere_cylinder": (lambda: Half(Sphere()), lambda: Half(CylinderZ(0.5, 0.3)),
                        [0.66, 0.36, 0.68], lambda: 0.5 * _length(exact_sphere_cylinder()),
                        lambda p: p[:, 2] > 0.0),
}  # fmt: skip


def check_curve(a, b, curve, plane, length, on_branch, seed) -> None:
    point, normal = (np.asarray(v, dtype=np.float64) for v in plane)
    normal = normal / np.linalg.norm(normal)
    diag = float(np.linalg.norm(np.ptp(np.asarray(a.bounds), axis=0)))
    assert curve.shape == (N, 3) and curve.dtype == np.float64
    assert a.distance(curve).max() <= 1e-9 * diag
    assert b.distance(curve).max() <= 1e-9 * diag
    signed = (curve - point) @ normal * np.sign((np.asarray(seed) - point) @ normal)
    assert np.abs(signed[[0, -1]]).max() <= 1e-9 * diag
    assert signed.min() >= -1e-9 * diag
    chords = np.linalg.norm(np.diff(curve, axis=0), axis=1)
    assert chords.max() / chords.min() <= 1.01
    assert abs(chords.sum() / length - 1.0) <= 0.005
    assert on_branch(curve).all()


# ---------------------------------------------------------------------------
# O1-O3: the curve
# ---------------------------------------------------------------------------


class TestCurve:
    @pytest.mark.parametrize("name", list(HALF))
    def test_o1_half_model(self, name):
        make_a, make_b, seed, length, on_branch = HALF[name]
        a, b = make_a(), make_b()
        curve = kernel()(a, b, np.asarray(seed, dtype=np.float64), N, plane=SYMMETRY)
        check_curve(a, b, curve, SYMMETRY, length(), on_branch, seed)

    def test_o2_seed_side_is_kept(self):
        a, b, seed = Plane(), CylinderZ(), np.array([1.1, -0.1, 0.05])
        curve = kernel()(a, b, seed, N, plane=SYMMETRY)
        check_curve(a, b, curve, SYMMETRY, np.pi, _all, seed)
        assert (curve[1:-1, 1] < 0.0).all()

    def test_o2_plane_not_a_symmetry_plane(self):
        a, b, seed = Plane(), CylinderZ(), np.array([1.1, 0.1, 0.05])
        plane = ((0.5, 0.0, 0.0), (1.0, 0.0, 0.0))
        curve = kernel()(a, b, seed, N, plane=plane)
        check_curve(a, b, curve, plane, 2.0 * np.pi / 3.0, _all, seed)

    def test_o3_stl_half_model(self):
        from shore.geometry import StlGeometry

        box = trimesh.creation.box(extents=(4.0, 2.0, 1.0))
        box.apply_translation((0.0, 1.0, -0.5))
        on_cut = np.all(np.isclose(box.vertices[box.faces][..., 1], 0.0), axis=1)
        box.update_faces(~on_cut)
        post = trimesh.creation.cylinder(radius=0.5, height=1.6, sections=96)
        post.apply_translation((0.0, 0.0, 0.2))
        post.update_faces(post.triangles_center[:, 1] > 0.0)
        for mesh in (box, post):
            mesh.remove_unreferenced_vertices()
        a, b = StlGeometry(box), StlGeometry(post)
        start = time.perf_counter()
        curve = kernel()(a, b, np.array([0.45, 0.2, 0.02]), N, plane=SYMMETRY)
        assert time.perf_counter() - start <= 5.0
        diag = float(np.linalg.norm(np.ptp(np.asarray(a.bounds), axis=0)))
        assert a.project(curve).distances.max() <= 1e-9 * diag
        assert b.project(curve).distances.max() <= 1e-9 * diag
        ends = sorted(curve[[0, -1]].tolist())
        np.testing.assert_allclose(ends, [[-0.5, 0.0, 0.0], [0.5, 0.0, 0.0]], atol=1e-9 * diag)
        assert curve[:, 1].min() >= -1e-9 * diag
        chords = np.linalg.norm(np.diff(curve, axis=0), axis=1)
        assert chords.max() / chords.min() <= 1.01


# ---------------------------------------------------------------------------
# O4: orientation, composed with the surface marcher (ADR 0008)
# ---------------------------------------------------------------------------


class TestCollarSides:
    @pytest.mark.parametrize("name", list(HALF))
    def test_o4_marching_leaves_the_other_body(self, name):
        from shore.surface.march import march_surface

        make_a, make_b, seed, _, _ = HALF[name]
        a, b = make_a(), make_b()
        curve = kernel()(a, b, np.asarray(seed, dtype=np.float64), N, plane=SYMMETRY)
        steps = np.full(2, 0.02)
        on_a = march_surface(a, curve, steps)
        on_b = march_surface(b, curve[::-1], steps)
        assert (b.side(on_a[:, 1]) > 0.0).all()
        assert (a.side(on_b[:, 1]) > 0.0).all()


# ---------------------------------------------------------------------------
# O5: interface
# ---------------------------------------------------------------------------

SEED = np.array([1.1, 0.1, 0.05])


class TestInterface:
    def test_o5_deterministic(self):
        a, b = Half(Sphere()), Half(CylinderZ(0.5, 0.3))
        seed = np.array([0.66, 0.36, 0.68])
        first = kernel()(a, b, seed, N, plane=SYMMETRY)
        again = kernel()(a, b, seed, N, plane=SYMMETRY, entity_a=BODY, entity_b=BODY)
        np.testing.assert_array_equal(first, again)

    def test_o5_too_few_nodes(self):
        with pytest.raises(ParameterError):
            kernel()(Plane(), CylinderZ(), SEED, 2, plane=SYMMETRY)

    @pytest.mark.parametrize("seed", [np.zeros(2), np.zeros((2, 3))])
    def test_o5_bad_seed(self, seed):
        with pytest.raises(ParameterError):
            kernel()(Plane(), CylinderZ(), seed, N, plane=SYMMETRY)

    @pytest.mark.parametrize(
        "plane", [((0.0, 0.0, 0.0), (0.0, 0.0, 0.0)), ((0.0, 0.0), (0.0, 1.0))]
    )
    def test_o5_bad_plane(self, plane):
        with pytest.raises(ParameterError):
            kernel()(Plane(), CylinderZ(), SEED, N, plane=plane)

    def test_o5_seed_on_the_plane(self):
        with pytest.raises(ParameterError):
            kernel()(Plane(), CylinderZ(), np.array([1.0, 0.0, 0.0]), N, plane=SYMMETRY)

    def test_o5_no_intersection(self):
        with pytest.raises(ParameterError):
            kernel()(Plane(), Sphere((0.0, 0.0, 5.0)), np.array([0.0, 0.1, 2.0]), N, plane=SYMMETRY)

    def test_o5_closed_curve_off_the_plane(self):
        start = time.perf_counter()
        with pytest.raises(ParameterError):
            kernel()(
                Sphere(), CylinderZ(0.5, 0.3), np.array([0.8, 0.0, 0.6]), N,
                plane=((0.0, 0.0, -0.5), (0.0, 0.0, 1.0)),
            )  # fmt: skip
        assert time.perf_counter() - start <= 5.0
