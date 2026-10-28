"""FROZEN acceptance spec: surface marching at clustered corners (ADR 0027).

Frozen 2026-10-28, before any implementation. Implementing work may not edit
this file, ``tests/spec/__init__.py``, ``tests/conftest.py``, the pytest
configuration in ``pyproject.toml``, or ``src/shore/quality.py``. The ONE
permitted edit is deleting an ``@PENDING`` line once that case passes (strict
xfail turns a passing pending case into a failure). ``tests/spec/SHA256SUMS``
holds the hashes, taken with ``@PENDING`` lines removed;
``scripts/check_spec.py`` checks them.

Verdict for acceptance (main CI job): ``pytest tests/spec/test_march_corners.py``
reports 15 passed, 0 skipped, 0 xfailed. At freeze: 2 passed (collars that
already hold: regression guards), 13 xfailed.

Interface under test
--------------------
``march_surface`` (ADR 0008, 0012, 0023, 0024), unchanged in signature:

- a layer whose step exceeds ``MarchParams.max_step_ratio`` (24) x the
  smallest front spacing, counting only nodes not closer together than on
  the starting curve, is marched in sub-steps of at most
  ``MarchParams.substep_ratio`` (2) x that spacing;
- where a new layer folds against the previous one (a cell reversed, or a
  front segment pointing backwards: converging fronts), its nodes round the
  fold are smoothed along the front and re-projected until none folds;
- marches where neither triggers are unchanged bit for bit (ADR 0023's
  golden marches, ADR 0024's golden collars).

Cases
-----
The junctions of ADR 0021-0024 (plane x L-shape prism; plate x blunt and
sharp NACA 0012 STL wings; plane x lens; plane x a 60 deg wedge with a
rounded back), loops with curvature clustering at their corners, surface
steps ``0.005 * 1.12**j`` (15 layers), periodic marches from the loop (on A)
or the reversed loop (on B).

- M1 plate marches (A) at concave corners and at the blunt base: no flipped
  cell, and every cell's sine (cell normal against the surface normal,
  normalised by its edges) >= 0.1 (prototype >= 0.18; before: 4-34 flipped
  cells): L-shape n_loop 192 curvature 1, n_loop 384 uniform and curvature
  1; blunt wing n_loop 192 curvature 1, n_loop 384 curvature 3.
- M2 marches along sharp creases (B): no flipped cell and in-section drift
  <= 0.02 (prototype <= 0.011; before: 94-1033 flipped cells, drift up to
  1.0): sharp STL wing n_loop 192 curvature 1; lens n_loop 192 curvature 1,
  n_loop 384 curvature 3; wedge n_loop 384 curvature 1.
- M3 collar walls (k = 0, 1.9 x thick, smoothed volume) without folded
  cells: ``build_collar`` on the L-shape and the blunt wing at n_loop 192
  with curvature 1 (before: 20 and 11), ``split_collar`` on the lens at
  n_loop 192 with curvature 1 (before: 259).
- M4 unchanged where nothing triggers: the L-shape and rectangle collars at
  n_loop 192, uniform, 1.9 x: valid, min scaled Jacobian >= 0.1 (0.180,
  0.191; regression guards).

Out of scope of this spec (measured, ADR 0027): the collar *volume* at these
clusterings (thickness 100-300 x the corner spacing: still folds, the
volume's limit); the sharp STL wing at n_loop 384 with curvature (spacing
3e-5 on sliver triangles: 2-119 flipped cells); run time (sub-stepping takes
up to 150-280 s for n_loop 384 with curvature).
"""

from __future__ import annotations

import numpy as np
import pytest

from shore.quality import scaled_jacobian
from tests.spec.test_intersection import Plane
from tests.spec.test_loop_corners import (
    L_SHAPE,
    RECT,
    LensPrism,
    PolygonPrism,
    stl_wing,
)

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-28; implementation pending")

STEPS = 0.005 * 1.12 ** np.arange(15)
THIN = 5e-5 * 1.2 ** np.arange(20)


def wedge(tip_deg: float = 60.0, radius: float = 0.1, n_arc: int = 48) -> PolygonPrism:
    """A z-extruded wedge with a rounded back (as ADR 0024's spec)."""
    a = np.radians(tip_deg) / 2.0
    tip = radius / np.sin(a)
    phi = np.linspace(np.pi / 2.0 - a, 3.0 * np.pi / 2.0 + a, n_arc)
    v = np.vstack([[tip, 0.0], np.c_[radius * np.cos(phi), radius * np.sin(phi)]])
    return PolygonPrism(v - [tip / 2.0, 0.0])


def bodies(name: str):
    if name == "l_shape":
        return Plane(), PolygonPrism(L_SHAPE), np.array([-0.2, -0.55, 0.0])
    if name == "rectangle":
        return Plane(), PolygonPrism(RECT), np.array([0.0, 0.06, 0.0])
    if name in ("blunt", "sharp"):
        plate, wing, _ = stl_wing(name == "sharp")
        return plate, wing, np.array([0.5, 0.06, 0.0])
    if name == "lens":
        return Plane(), LensPrism(), np.array([0.0, 0.06, 0.0])
    b = wedge()
    v0, v1 = b.v[0], b.v[1]
    side = v1 - v0
    normal = np.array([-side[1], side[0]]) / np.linalg.norm(side)
    return Plane(), b, np.r_[0.5 * (v0 + v1) + 0.01 * normal, 0.0]


def march(name: str, n: int, curvature: float, side: str):
    from shore.surface.intersect import intersection_loop
    from shore.surface.march import march_surface

    a, b, seed = bodies(name)
    loop = intersection_loop(a, b, seed, n, curvature=curvature)
    body, front = (a, loop) if side == "A" else (b, loop[::-1])
    return body, front, march_surface(body, front, STEPS, periodic=True)


def cell_sines(body, grid: np.ndarray) -> np.ndarray:
    """Per cell: (e_front x e_march) . n / (|e_front| |e_march|), majority sign positive."""
    e1 = (np.roll(grid, -1, axis=0) - grid)[:, :-1]
    e2 = grid[:, 1:] - grid[:, :-1]
    centre = 0.5 * (grid[:, :-1] + grid[:, 1:])
    normals = np.asarray(body.project(centre.reshape(-1, 3)).normals).reshape(centre.shape)
    o = np.einsum("...k,...k->...", np.cross(e1, e2), normals)
    o = o / (np.linalg.norm(e1, axis=-1) * np.linalg.norm(e2, axis=-1))
    return o * np.sign(np.median(o))


class TestPlate:
    @PENDING
    @pytest.mark.parametrize(
        "name, n, curvature",
        [("l_shape", 192, 1.0), ("l_shape", 384, 0.0), ("l_shape", 384, 1.0),
         ("blunt", 192, 1.0), ("blunt", 384, 3.0)],
    )  # fmt: skip
    def test_m1_plate_marches(self, name, n, curvature):
        from tests.spec.test_march_creases import flipped

        body, _, grid = march(name, n, curvature, "A")
        assert flipped(body, grid) == 0
        assert cell_sines(body, grid).min() >= 0.1


class TestCreases:
    @PENDING
    @pytest.mark.parametrize(
        "name, n, curvature",
        [("sharp", 192, 1.0), ("lens", 192, 1.0), ("lens", 384, 3.0), ("wedge", 384, 1.0)],
    )
    def test_m2_crease_marches(self, name, n, curvature):
        from tests.spec.test_march_creases import flipped

        body, front, grid = march(name, n, curvature, "B")
        assert flipped(body, grid) == 0
        assert np.abs(grid[:, :, :2] - front[:, None, :2]).max() <= 0.02


def collar_blocks(kind: str, name: str, n: int, curvature: float):
    from shore.collar import build_collar, split_collar

    a, b, seed = bodies(name)
    f = build_collar if kind == "build" else split_collar
    out = f(
        a, b, seed, n_loop=n, steps_a=STEPS, steps_b=STEPS, steps_k=THIN,
        loop_curvature=curvature,
    )  # fmt: skip
    return out if isinstance(out, list) else [out]


def jacobians(blocks):
    for blk in blocks:
        x = np.concatenate([blk.nodes, blk.nodes[:, :1]], axis=1) if blk.periodic_j else blk.nodes
        yield scaled_jacobian(x)


class TestCollars:
    @PENDING
    @pytest.mark.parametrize(
        "kind, name", [("build", "l_shape"), ("build", "blunt"), ("split", "lens")]
    )
    def test_m3_collar_walls(self, kind, name):
        blocks = collar_blocks(kind, name, 192, 1.0)
        assert sum(int((jac[:, :, 0] <= 0.0).sum()) for jac in jacobians(blocks)) == 0

    @pytest.mark.parametrize("name", ["l_shape", "rectangle"])
    def test_m4_untriggered_collars(self, name):  # already hold: regression guards
        blocks = collar_blocks("build", name, 192, 0.0)
        assert min(float(jac.min()) for jac in jacobians(blocks)) >= 0.1

    @PENDING
    def test_m4_params(self):
        from shore.surface.march import MarchParams

        p = MarchParams()
        assert (p.max_step_ratio, p.substep_ratio) == (24.0, 2.0)
