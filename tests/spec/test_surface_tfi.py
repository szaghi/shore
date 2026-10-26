"""FROZEN acceptance spec: surface TFI + elliptic smoothing in (u, v) (ADR 0007, roadmap phase 4).

Frozen 2026-10-26, before any implementation. Implementing work may not edit
this file, ``tests/spec/test_surface_tfi_cad.py``, ``tests/spec/__init__.py``,
``tests/conftest.py``, the pytest configuration in ``pyproject.toml``, or
``src/shore/quality.py``. The ONE permitted edit is deleting an ``@PENDING``
line once that case passes (strict xfail turns a passing pending case into a
failure). ``tests/spec/SHA256SUMS`` holds the hashes, taken with ``@PENDING``
lines removed; ``scripts/check_spec.py`` checks them.

Verdict for acceptance (main CI job, no CAD extra):
``pytest tests/spec/test_surface_tfi.py`` reports 23 passed, 0 skipped,
0 xfailed. The CAD half is ``test_surface_tfi_cad.py`` (CAD job). At freeze:
0 passed, 23 xfailed (the kernel module does not exist yet).

Interface under test
--------------------
``shore.surface.tfi.surface_grid(geometry, surface, boundary, smoothing="elliptic")``
returns ``SurfaceGrid(points, params)``:

- *geometry* is any ``Geometry`` with parametric SURFACE entities; *surface*
  the ``EntityRef`` of one of them.
- *boundary* is ``(ni, nj, 3)``: only its four boundary rows are read (the
  interior is NaN in every case below). They lie on the surface; corners are
  shared by construction.
- ``smoothing="none"`` is plain transfinite interpolation of the boundary's
  (u, v) parameters, evaluated on the surface; ``"elliptic"`` (default) then
  smooths in (u, v) with the surface metric (Winslow-type) and control
  functions that keep the boundary's clustering (Thomas-Middlecoff-type).
- ``points`` ``(ni, nj, 3)`` lie on the surface, the boundary rows equal the
  input, ``params`` ``(ni, nj, 2)`` evaluate to ``points``.
- ``ParameterError`` for: an unknown *smoothing*, a malformed *boundary*
  (not ``(ni, nj, 3)`` with ni, nj >= 3), boundary nodes farther than
  1e-6 x the model diagonal from the surface, a geometry or entity without a
  (u, v) parametrisation (an STL body).
- Smoothing happens in (u, v): the geometry is asked to evaluate or project
  at most 10 points per grid node in total, however many sweeps are run
  (ADR 0002 throughput decision).

The surfaces are exact analytic ``Geometry`` implementations defined here
(numpy only), so every reference is exact and the call count is observable:
the identity plane, a plane and a cylinder with a distorted parametrisation
(``u + K u (1 - u)(2 v - 1)``, K = 0.6, mapping the unit square onto itself),
and a latitude-longitude sphere patch.

Thresholds and references (prototype measurements, ADR 0007; "h" = the
reference grid's smallest local spacing at the node; "TFI" = ``"none"``):

- S1 tensor grid, identity plane: TFI reproduces it (<= 1e-12 diag; exact).
- S1 (u, v)-tensor grid, distorted cylinder: TFI reproduces S(f(xi), g(eta))
  (<= 1e-12 diag); physical TFI is 0.046 off (checked in the case).
- S2 four surfaces x both modes: boundary = input, interior on the surface,
  evaluate(params) = points (all <= 1e-12 diag).
- S3 uniform physical boundaries, distorted plane / cylinder 21x21: elliptic
  <= 0.02 h from the uniform physical grid (prototype (u, v) Winslow
  0.0007 / 0.0008 h); TFI >= 0.3 h (0.65 / 1.03 h).
- S4 clustered tensor rectangle 41x31 (tanh, ds = 8.0e-4): first cell
  0.95-1.05 ds along the clustered edge, <= 0.05 h (Thomas-Middlecoff
  1.000, 0.0000 h; plain Winslow 1.66-13.3, 33 h).
- S5 S-shaped overhang 81x61, arc-length-uniform boundary: TFI folds (29
  cells); elliptic 0 folded, min corner angle >= 2 deg (prototype 5.4 deg).
- S6 sphere patch 31x31 (lat 0-1.2 rad): elliptic 0 folded.
- S7 distorted cylinder 41x41: <= 10 points evaluated + projected per node.
- S8 interface: default, determinism, errors, types (as above).

Out of scope of this spec: faces whose (u, v) range crosses a periodic seam
or contains a singular point (pole, apex); block faces spanning several
surfaces; clustered boundaries on strongly concave boundaries (S4 and S5
test clustering and fold repair separately); run time. Elliptic parameters
(iterations, control-function form, relaxation) may be chosen freely but
must not be passed per case: every case runs with the defaults.
"""

from __future__ import annotations

import numpy as np
import pytest
import trimesh

from shore.exceptions import ParameterError
from shore.geometry import BODY, EntityKind, EntityRef, Projection, StlGeometry

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-26; implementation pending")

K = 0.6  # parametrisation distortion
FACE = EntityRef(EntityKind.SURFACE, 0)


def kernel():
    """The function under test (imported late: absent until implemented)."""
    from shore.surface.tfi import surface_grid

    return surface_grid


# ---------------------------------------------------------------------------
# Exact analytic surfaces (frozen with the spec)
# ---------------------------------------------------------------------------


def _unwarp(w: np.ndarray, v: np.ndarray) -> np.ndarray:
    """Solve ``w = u + K u (1 - u)(2 v - 1)`` for u in [0, 1]."""
    a = -K * (2.0 * v - 1.0)
    b = 1.0 + K * (2.0 * v - 1.0)
    small = np.abs(a) < 1e-12
    safe = np.where(small, 1.0, a)
    root = (-b + np.sqrt(b * b + 4.0 * safe * w)) / (2.0 * safe)
    return np.where(small, w / b, root)


def _warp(u: np.ndarray, v: np.ndarray) -> np.ndarray:
    return u + K * u * (1.0 - u) * (2.0 * v - 1.0)


class AnalyticSurface:
    """One exact parametric SURFACE behind the ``Geometry`` protocol.

    Counts the points it is asked to evaluate or project (``work``).
    """

    def __init__(self) -> None:
        self.work = 0

    @property
    def bounds(self) -> np.ndarray:
        return self._bounds

    def entities(self, kind: EntityKind) -> tuple[EntityRef, ...]:
        return (FACE,) if kind is EntityKind.SURFACE else ()

    def ray_cast(self, *args: object, **kwargs: object) -> object:
        raise NotImplementedError("not used by surface gridding")

    def fingerprint(self, entity: EntityRef) -> object:
        raise NotImplementedError("not used by surface gridding")

    def evaluate(self, entity: EntityRef, params: np.ndarray) -> np.ndarray:
        assert entity == FACE, entity
        prm = np.asarray(params, dtype=np.float64)
        self.work += prm[..., 0].size
        return self.S(prm[..., 0], prm[..., 1])

    def project(
        self, points: np.ndarray, entity: EntityRef | None = None, *, exact: bool = False
    ) -> Projection:
        pts = np.asarray(points, dtype=np.float64)
        self.work += pts[..., 0].size
        uv, foot, normal = self.inverse(pts)
        dist = np.linalg.norm(pts - foot, axis=-1)
        return Projection(foot, dist, normal, uv if entity == FACE else None, None)

    # exact helpers for the checks (not counted)
    def distance(self, pts: np.ndarray) -> np.ndarray:
        return np.linalg.norm(pts - self.inverse(pts)[1], axis=-1)

    def normal(self, pts: np.ndarray) -> np.ndarray:
        return self.inverse(pts)[2]


class Plane(AnalyticSurface):
    """z = 0, (u, v) = (x, y), optionally with the distorted u."""

    def __init__(self, distorted: bool = False) -> None:
        super().__init__()
        self.distorted = distorted
        self._bounds = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 0.0]])

    def S(self, u: np.ndarray, v: np.ndarray) -> np.ndarray:
        x = _warp(u, v) if self.distorted else u
        return np.stack([x, v, np.zeros_like(u)], axis=-1)

    def inverse(self, p: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        x, y = p[..., 0], p[..., 1]
        u = _unwarp(x, y) if self.distorted else x
        foot = np.stack([x, y, np.zeros_like(x)], axis=-1)
        normal = np.broadcast_to([0.0, 0.0, 1.0], p.shape).copy()
        return np.stack([u, y], axis=-1), foot, normal


class Cylinder(AnalyticSurface):
    """Unit cylinder about z: phi = (pi / 2) * warp(u, v), z = v (a developable surface)."""

    def __init__(self) -> None:
        super().__init__()
        self._bounds = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]])

    def S(self, u: np.ndarray, v: np.ndarray) -> np.ndarray:
        phi = 0.5 * np.pi * _warp(u, v)
        return np.stack([np.cos(phi), np.sin(phi), v], axis=-1)

    def inverse(self, p: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        phi = np.arctan2(p[..., 1], p[..., 0])
        z = p[..., 2]
        u = _unwarp(phi / (0.5 * np.pi), z)
        normal = np.stack([np.cos(phi), np.sin(phi), np.zeros_like(phi)], axis=-1)
        foot = normal + np.stack([np.zeros_like(z), np.zeros_like(z), z], axis=-1)
        return np.stack([u, z], axis=-1), foot, normal


class SpherePatch(AnalyticSurface):
    """Unit sphere, (u, v) = (longitude, latitude)."""

    def __init__(self) -> None:
        super().__init__()
        self._bounds = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]])

    def S(self, u: np.ndarray, v: np.ndarray) -> np.ndarray:
        return np.stack([np.cos(v) * np.cos(u), np.cos(v) * np.sin(u), np.sin(v)], axis=-1)

    def inverse(self, p: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        r = np.linalg.norm(p, axis=-1)
        normal = p / r[..., None]
        lon = np.arctan2(p[..., 1], p[..., 0])
        lat = np.arcsin(np.clip(normal[..., 2], -1.0, 1.0))
        return np.stack([lon, lat], axis=-1), normal, normal


# ---------------------------------------------------------------------------
# Distributions, grids and checks (frozen with the spec)
# ---------------------------------------------------------------------------


def one_sided(n: int, beta: float) -> np.ndarray:
    """Roberts tanh law on [0, 1], dense at 0."""
    tau = np.linspace(0.0, 1.0, n)
    return 1.0 - np.tanh(beta * (1.0 - tau)) / np.tanh(beta)


def two_sided(n: int, beta: float) -> np.ndarray:
    """Roberts tanh law on [0, 1], dense at both ends."""
    tau = np.linspace(0.0, 1.0, n)
    return 0.5 * (1.0 + np.tanh(beta * (2.0 * tau - 1.0)) / np.tanh(beta))


def arc_resample(curve, n: int, m: int = 20001) -> np.ndarray:
    """*n* points equally spaced in arc length along ``curve(s)``, s in [0, 1]."""
    pts = curve(np.linspace(0.0, 1.0, m))
    length = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))]
    target = np.linspace(0.0, length[-1], n)
    return np.stack([np.interp(target, length, pts[:, k]) for k in range(pts.shape[1])], -1)


def boundary_of(grid: np.ndarray) -> np.ndarray:
    """*grid* with its interior replaced by NaN: only the boundary may be read."""
    b = np.array(grid, dtype=np.float64)
    b[1:-1, 1:-1] = np.nan
    return b


def tfi(b: np.ndarray) -> np.ndarray:
    """Plain linear transfinite interpolation of an array's boundary rows."""
    ni, nj = b.shape[:2]
    s = np.linspace(0.0, 1.0, ni)[:, None, None]
    t = np.linspace(0.0, 1.0, nj)[None, :, None]
    return (
        (1 - s) * b[0][None] + s * b[-1][None] + (1 - t) * b[:, 0][:, None] + t * b[:, -1][:, None]
        - ((1 - s) * (1 - t) * b[0, 0] + s * (1 - t) * b[-1, 0]
           + (1 - s) * t * b[0, -1] + s * t * b[-1, -1])
    )  # fmt: skip


def diag(points: np.ndarray) -> float:
    flat = points.reshape(-1, 3)
    return float(np.linalg.norm(flat.max(axis=0) - flat.min(axis=0)))


def local_spacing(ref: np.ndarray) -> np.ndarray:
    """Smallest of the four edges at every interior node of *ref*, (ni-2, nj-2)."""
    di = np.linalg.norm(np.diff(ref, axis=0), axis=-1)
    dj = np.linalg.norm(np.diff(ref, axis=1), axis=-1)
    return np.minimum.reduce([di[:-1, 1:-1], di[1:, 1:-1], dj[1:-1, :-1], dj[1:-1, 1:]])


def deviation(points: np.ndarray, ref: np.ndarray) -> float:
    """Max interior node distance from *ref*, in units of the local spacing."""
    d = np.linalg.norm(points - ref, axis=-1)[1:-1, 1:-1]
    return float((d / local_spacing(ref)).max())


def corner_sines(points: np.ndarray, normals: np.ndarray) -> np.ndarray:
    """Signed sine of each cell's smallest corner angle, (ni-1, nj-1).

    Corners in (i, j) order P00 -> P10 -> P11 -> P01, oriented by the surface
    normal; <= 0 is a folded cell.
    """
    p = [points[:-1, :-1], points[1:, :-1], points[1:, 1:], points[:-1, 1:]]
    n = [normals[:-1, :-1], normals[1:, :-1], normals[1:, 1:], normals[:-1, 1:]]
    out = []
    for c in range(4):
        a = p[(c + 1) % 4] - p[c]
        b = p[(c - 1) % 4] - p[c]
        s = np.einsum("...j,...j->...", np.cross(a, b), n[c])
        out.append(s / (np.linalg.norm(a, axis=-1) * np.linalg.norm(b, axis=-1)))
    return np.minimum.reduce(out)


# --- the cases' grids -------------------------------------------------------


def uniform_physical(surface: AnalyticSurface, n: int) -> np.ndarray:
    """Uniform grid in physical arc length on the distorted plane / cylinder."""
    w, v = np.meshgrid(np.linspace(0, 1, n), np.linspace(0, 1, n), indexing="ij")
    return surface.S(_unwarp(w, v), v)


def clustered_rectangle() -> np.ndarray:
    x, y = np.meshgrid(one_sided(41, 3.0), 0.5 * two_sided(31, 2.0), indexing="ij")
    return np.stack([x, y, np.zeros_like(x)], axis=-1)


def overhang_boundary(ni: int = 81, nj: int = 61, a: float = 0.16, c: float = 0.3) -> np.ndarray:
    """Unit square whose top edge is an S-shaped overhang, arc-length-uniform nodes."""
    b = np.full((ni, nj, 3), np.nan)
    s, t = np.linspace(0, 1, ni), np.linspace(0, 1, nj)
    b[:, 0] = np.c_[s, np.zeros(ni), np.zeros(ni)]
    b[0] = np.c_[np.zeros(nj), t, np.zeros(nj)]
    b[-1] = np.c_[np.ones(nj), t, np.zeros(nj)]
    top = arc_resample(
        lambda q: np.c_[q + a * np.sin(2 * np.pi * q), 1 - c * np.sin(2 * np.pi * q)], ni
    )
    b[:, -1] = np.c_[top, np.zeros(ni)]
    return b


def sphere_patch_grid(n: int = 31) -> np.ndarray:
    u, v = np.meshgrid(np.linspace(0, 0.5 * np.pi, n), np.linspace(0, 1.2, n), indexing="ij")
    return SpherePatch().S(u, v)


# ---------------------------------------------------------------------------
# S1: plain TFI is exact on tensor-product data, and it is done in (u, v)
# ---------------------------------------------------------------------------


class TestS1Tfi:
    def test_tensor_grid_identity_plane(self):
        grid = clustered_rectangle()
        g = kernel()(Plane(), FACE, boundary_of(grid), smoothing="none")
        assert np.abs(g.points - grid).max() <= 1e-12 * diag(grid)

    def test_interpolates_parameters_not_points(self):
        surf = Cylinder()
        u, v = np.meshgrid(one_sided(25, 2.0), two_sided(21, 1.5), indexing="ij")
        grid = surf.S(u, v)
        # physical TFI of the same boundary cuts chords: the case discriminates
        assert np.abs(tfi(boundary_of(grid)) - grid).max() > 1e-2
        g = kernel()(surf, FACE, boundary_of(grid), smoothing="none")
        assert np.abs(g.points - grid).max() <= 1e-12 * diag(grid)
        assert np.abs(g.params - np.stack([u, v], axis=-1)).max() <= 1e-12


# ---------------------------------------------------------------------------
# S2: invariants on every surface, both modes
# ---------------------------------------------------------------------------

SURFACES = {
    "plane": (Plane, lambda: clustered_rectangle()),
    "distorted_plane": (lambda: Plane(distorted=True), lambda: uniform_physical(Plane(True), 21)),
    "cylinder": (Cylinder, lambda: uniform_physical(Cylinder(), 21)),
    "sphere_patch": (SpherePatch, sphere_patch_grid),
}


class TestS2Invariants:
    @pytest.mark.parametrize("smoothing", ["none", "elliptic"])
    @pytest.mark.parametrize("name", list(SURFACES))
    def test_boundary_surface_params(self, name, smoothing):
        make_surface, make_grid = SURFACES[name]
        surf, grid = make_surface(), make_grid()
        g = kernel()(surf, FACE, boundary_of(grid), smoothing=smoothing)
        tol = 1e-12 * diag(grid)
        assert g.points.shape == grid.shape and g.points.dtype == np.float64
        assert g.params.shape == (*grid.shape[:2], 2) and g.params.dtype == np.float64
        for edge in (np.s_[0], np.s_[-1], np.s_[:, 0], np.s_[:, -1]):
            assert np.abs(g.points[edge] - grid[edge]).max() <= tol
        assert surf.distance(g.points).max() <= tol
        assert np.abs(surf.S(g.params[..., 0], g.params[..., 1]) - g.points).max() <= tol


# ---------------------------------------------------------------------------
# S3: smoothing uses the surface metric
# ---------------------------------------------------------------------------


class TestS3Metric:
    @pytest.mark.parametrize("make_surface", [lambda: Plane(distorted=True), Cylinder])
    def test_recovers_uniform_physical_grid(self, make_surface):
        surf = make_surface()
        ref = uniform_physical(surf, 21)
        g_tfi = kernel()(surf, FACE, boundary_of(ref), smoothing="none")
        assert deviation(g_tfi.points, ref) >= 0.3  # the case discriminates
        g = kernel()(surf, FACE, boundary_of(ref), smoothing="elliptic")
        assert deviation(g.points, ref) <= 0.02


# ---------------------------------------------------------------------------
# S4: elliptic smoothing keeps the boundary's clustering
# ---------------------------------------------------------------------------


class TestS4Clustering:
    def test_tensor_rectangle_kept(self):
        ref = clustered_rectangle()
        g = kernel()(Plane(), FACE, boundary_of(ref))
        ds = ref[1, 0, 0] - ref[0, 0, 0]
        first = np.linalg.norm(g.points[1, 1:-1] - g.points[0, 1:-1], axis=-1) / ds
        assert first.min() >= 0.95 and first.max() <= 1.05
        assert deviation(g.points, ref) <= 0.05


# ---------------------------------------------------------------------------
# S5: elliptic smoothing repairs a folding TFI
# ---------------------------------------------------------------------------


class TestS5FoldRepair:
    def test_overhang(self):
        b = overhang_boundary()
        up = np.broadcast_to([0.0, 0.0, 1.0], b.shape)
        g_tfi = kernel()(Plane(), FACE, b, smoothing="none")
        assert (corner_sines(g_tfi.points, up) <= 0).sum() >= 1  # the case discriminates
        g = kernel()(Plane(), FACE, b)
        sines = corner_sines(g.points, up)
        assert (sines <= 0).sum() == 0
        assert np.degrees(np.arcsin(sines.min())) >= 2.0


# ---------------------------------------------------------------------------
# S6: a non-developable surface
# ---------------------------------------------------------------------------


class TestS6Sphere:
    def test_patch_valid(self):
        surf, grid = SpherePatch(), sphere_patch_grid()
        g = kernel()(surf, FACE, boundary_of(grid))
        assert (corner_sines(g.points, surf.normal(g.points)) <= 0).sum() == 0


# ---------------------------------------------------------------------------
# S7: smoothing in (u, v), not on the surface (ADR 0002)
# ---------------------------------------------------------------------------


class TestS7Work:
    def test_geometry_queries_per_node(self):
        surf = Cylinder()
        ref = uniform_physical(surf, 41)
        surf.work = 0
        kernel()(surf, FACE, boundary_of(ref))
        assert surf.work <= 10 * ref.shape[0] * ref.shape[1]


# ---------------------------------------------------------------------------
# S8: interface
# ---------------------------------------------------------------------------


class TestS8Interface:
    def test_default_is_elliptic_and_deterministic(self):
        surf, grid = Cylinder(), uniform_physical(Cylinder(), 21)
        a = kernel()(surf, FACE, boundary_of(grid))
        b = kernel()(surf, FACE, boundary_of(grid), smoothing="elliptic")
        np.testing.assert_array_equal(a.points, b.points)
        np.testing.assert_array_equal(a.params, b.params)

    def test_unknown_smoothing(self):
        grid = clustered_rectangle()
        with pytest.raises(ParameterError, match="smoothing"):
            kernel()(Plane(), FACE, boundary_of(grid), smoothing="laplace")

    @pytest.mark.parametrize("shape", [(41, 31, 2), (2, 31, 3), (41, 31)])
    def test_malformed_boundary(self, shape):
        with pytest.raises(ParameterError):
            kernel()(Plane(), FACE, np.zeros(shape))

    def test_boundary_off_surface(self):
        b = boundary_of(clustered_rectangle())
        b[0, 5, 2] += 1e-3  # one boundary node lifted off the plane
        with pytest.raises(ParameterError):
            kernel()(Plane(), FACE, b)

    def test_no_parametrisation(self):
        stl = StlGeometry(trimesh.creation.box(extents=(1.0, 1.0, 1.0)))
        grid = clustered_rectangle() - [0.5, 0.25, 0.5]  # on the box's z = -0.5 face
        with pytest.raises(ParameterError):
            kernel()(stl, BODY, boundary_of(grid))
