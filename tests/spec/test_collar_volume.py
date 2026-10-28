"""FROZEN acceptance spec: the collar volume at loop corners (ADR 0025).

Frozen 2026-10-28, before any implementation. Implementing work may not edit
this file, ``tests/spec/data/collar_default.npz``, ``tests/spec/__init__.py``,
``tests/conftest.py``, the pytest configuration in ``pyproject.toml``, or
``src/shore/quality.py``. The ONE permitted edit is deleting an ``@PENDING``
line once that case passes (strict xfail turns a passing pending case into a
failure). ``tests/spec/SHA256SUMS`` holds the hashes, taken with ``@PENDING``
lines removed; ``scripts/check_spec.py`` checks them.

Verdict for acceptance (main CI job): ``pytest tests/spec/test_collar_volume.py``
reports 25 passed, 0 skipped, 0 xfailed. At freeze: 3 passed (the default
collars, unchanged: regression guards), 22 xfailed.

Interface under test
--------------------
``CollarParams.volume``: ``"stations"`` (the default: ADR 0010's outer
boundary, per loop station, unchanged bit for bit) or ``"smoothed"`` (ADR
0025: each station's outer line on its 3-D offset polyline, the station
offsets smoothed along the loop and wherever the outer boundary steps back,
each line trimmed again), for ``build_collar`` (closed and half collars) and
``split_collar``; ``ParameterError`` for another value.
``Project.build_collar(..., volume=)`` (``BuildCollar.volume``) forwards it.

Cases
-----
The junctions of ADR 0021-0024 (plane x rectangle prism, plane x L-shape,
blunt NACA 0012 STL wing, the sharp wing and the lens split, the box edge x
ball, the elliptic post and its half), surface steps ``0.005 * 1.12**j`` (15
layers; the ellipse ``0.01 * 1.12**j``), wall-normal ``5e-5 * 1.2**k`` with
20 / 25 / 30 layers (1.9, 4.7, 11.8 x the first surface step), uniform loop
or curvature 1.

- V1 valid with ``"smoothed"`` where ``"stations"`` folds: every cell's
  scaled Jacobian at least the bound (prototype in brackets; folded cells
  with ``"stations"`` before them): rectangle n_loop 96 at 11.8 x uniform
  >= 0.05 (0.105; 20), curvature >= 0.1 (0.213; 52); n_loop 192 at 4.7 x
  curvature >= 0.04 (0.083; 54), 11.8 x curvature >= 0.1 (0.205; 150);
  L-shape n_loop 96 at 1.9 x uniform >= 0.1 (0.230; 4), 4.7 x curvature
  >= 0.03 (0.066; 39); n_loop 192 at 4.7 x uniform >= 0.03 (0.064; 271);
  blunt wing n_loop 96 at 11.8 x curvature >= 0.01 (0.023; 71); n_loop 192
  at 11.8 x uniform >= 0.025 (0.059; 29); split collars at 4.7 x curvature:
  sharp wing >= 0.05 (0.134; 187), lens >= 0.2 (0.437; 164).
- V2 no harm where ``"stations"`` is valid: box edge 4.7 x >= 0.1 (0.161);
  ellipse 1.9 x curvature >= 0.1 (0.270); rectangle n_loop 96 at 4.7 x
  uniform >= 0.025 (0.048, ADR 0022's bound); split sharp wing 1.9 x uniform
  >= 0.1 (0.180); half ellipse 4.7 x curvature >= 0.1 (0.309), its ends on
  the plane within 1e-12.
- V3 the default: ``CollarParams().volume == "stations"``; ``build_collar``
  with the default and with ``"stations"`` given gives ADR 0024's golden
  collars (``data/collar_default.npz``) within 1e-12.
- V4 ``ParameterError`` for an unknown volume.
- V5 ``Project.build_collar(volume="smoothed")`` gives the kernel's block
  and journals ``volume='smoothed'``.

Out of scope of this spec (measured, ADR 0025): concave corners at 4.7 x
and thicker (the L-shape keeps 2-19 folded cells), the blunt trailing edge's
base with n_loop 96 (2-15 cells), and loops whose surface grid itself folds
(the L-shape and the blunt wing with curvature at n_loop 192: the plate's
surface march flips 10 and 4 cells, a marcher problem). Making
``"smoothed"`` the default is a later, user-approved re-freeze of ADR 0024's
golden collars.
"""

from __future__ import annotations

import numpy as np
import pytest

from shore.quality import scaled_jacobian
from tests.spec.test_corner_grading import Ball
from tests.spec.test_intersection import Plane
from tests.spec.test_loop_corners import L_SHAPE, RECT, Box, PolygonPrism, stl_wing
from tests.spec.test_loop_spacing import (
    ELLIPSE_SEED,
    HALF_ELLIPSE_SEED,
    SYM,
    EllipticCylinder,
)
from tests.spec.test_split_collar import (
    G_STEPS_K,
    G_STEPS_S,
    GOLDEN,
    golden_case,
    sharp_case,
)

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-28; implementation pending")

STEPS_S = 0.005 * 1.12 ** np.arange(15)
LAYERS = {1.9: 20, 4.7: 25, 11.8: 30}


def steps_k(thick: float) -> np.ndarray:
    return 5e-5 * 1.2 ** np.arange(LAYERS[thick])


def params(volume: str):
    from shore.collar import CollarParams

    return CollarParams(volume=volume)


def bodies(name: str):
    if name == "rectangle":
        return Plane(), PolygonPrism(RECT), np.array([0.0, 0.06, 0.0])
    if name == "l_shape":
        return Plane(), PolygonPrism(L_SHAPE), np.array([-0.2, -0.55, 0.0])
    plate, wing, _ = stl_wing(False)
    return plate, wing, np.array([0.5, 0.06, 0.0])


def min_jacobian(blocks) -> float:
    worst = np.inf
    for blk in blocks:
        x = np.concatenate([blk.nodes, blk.nodes[:, :1]], axis=1) if blk.periodic_j else blk.nodes
        worst = min(worst, float(scaled_jacobian(x).min()))
    return worst


def collar(name: str, n_loop: int, thick: float, curvature: float, volume: str = "smoothed"):
    from shore.collar import build_collar

    a, b, seed = bodies(name)
    return build_collar(
        a, b, seed, n_loop=n_loop, steps_a=STEPS_S, steps_b=STEPS_S, steps_k=steps_k(thick),
        loop_curvature=curvature, params=params(volume),
    )  # fmt: skip


def split(name: str, thick: float, curvature: float, volume: str = "smoothed"):
    from shore.collar import split_collar

    a, b, seed, _ = sharp_case(name)
    return split_collar(
        a, b, seed, n_loop=96, steps_a=STEPS_S, steps_b=STEPS_S, steps_k=steps_k(thick),
        loop_curvature=curvature, params=params(volume),
    )  # fmt: skip


# ---------------------------------------------------------------------------
# V1: valid where the per-station boundary folds
# ---------------------------------------------------------------------------

V1 = [
    # name, n_loop, thickness, curvature, bound
    ("rectangle", 96, 11.8, 0.0, 0.05),
    ("rectangle", 96, 11.8, 1.0, 0.1),
    ("rectangle", 192, 4.7, 1.0, 0.04),
    ("rectangle", 192, 11.8, 1.0, 0.1),
    ("l_shape", 96, 1.9, 0.0, 0.1),
    ("l_shape", 96, 4.7, 1.0, 0.03),
    ("l_shape", 192, 4.7, 0.0, 0.03),
    ("blunt_wing", 96, 11.8, 1.0, 0.01),
    ("blunt_wing", 192, 11.8, 0.0, 0.025),
]


class TestValid:
    @pytest.mark.parametrize("name, n_loop, thick, curvature, bound", V1)
    def test_v1_collars(self, name, n_loop, thick, curvature, bound):
        assert min_jacobian([collar(name, n_loop, thick, curvature)]) >= bound

    @pytest.mark.parametrize("name, bound", [("wing", 0.05), ("lens", 0.2)])
    def test_v1_split_collars(self, name, bound):
        assert min_jacobian(split(name, 4.7, 1.0)) >= bound


# ---------------------------------------------------------------------------
# V2: no harm where the per-station boundary is valid
# ---------------------------------------------------------------------------


class TestNoHarm:
    def test_v2_box_edge(self):
        from shore.collar import build_collar

        block = build_collar(
            Box(), Ball(), np.array([0.7, 0.0, 0.05]), n_loop=96, steps_a=STEPS_S,
            steps_b=STEPS_S, steps_k=steps_k(4.7), params=params("smoothed"),
        )  # fmt: skip
        assert min_jacobian([block]) >= 0.1

    def test_v2_ellipse(self):
        from shore.collar import build_collar

        s = 0.01 * 1.12 ** np.arange(15)
        block = build_collar(
            Plane(), EllipticCylinder(), ELLIPSE_SEED, n_loop=128, steps_a=s, steps_b=s,
            steps_k=steps_k(1.9), loop_curvature=1.0, params=params("smoothed"),
        )  # fmt: skip
        assert min_jacobian([block]) >= 0.1

    def test_v2_rectangle_adr0022(self):
        assert min_jacobian([collar("rectangle", 96, 4.7, 0.0)]) >= 0.025

    def test_v2_split_wing(self):
        assert min_jacobian(split("wing", 1.9, 0.0)) >= 0.1

    def test_v2_half_collar(self):
        from shore.collar import build_collar

        s = 0.01 * 1.12 ** np.arange(15)
        block = build_collar(
            Plane(), EllipticCylinder(), HALF_ELLIPSE_SEED, n_loop=65, steps_a=s, steps_b=s,
            steps_k=steps_k(4.7), plane=SYM, loop_curvature=1.0, params=params("smoothed"),
        )  # fmt: skip
        assert min_jacobian([block]) >= 0.1
        assert float(np.abs(block.nodes[:, [0, -1], :, 1]).max()) <= 1e-12


# ---------------------------------------------------------------------------
# V3-V5: the default, errors, the project layer
# ---------------------------------------------------------------------------


class TestDefault:
    def test_v3_default_is_stations(self):
        from shore.collar import CollarParams

        assert CollarParams().volume == "stations"

    @pytest.mark.parametrize("name", ["rectangle", "blunt_wing", "half_ellipse"])
    def test_v3_default_unchanged(self, name):  # already holds: regression guards
        from shore.collar import build_collar

        a, b, seed, kw = golden_case(name)
        steps = {"steps_a": G_STEPS_S, "steps_b": G_STEPS_S, "steps_k": G_STEPS_K}
        block = build_collar(a, b, seed, **steps, **kw)
        np.testing.assert_allclose(block.nodes, np.load(GOLDEN)[name], rtol=0.0, atol=1e-12)

    @pytest.mark.parametrize("name", ["rectangle", "blunt_wing", "half_ellipse"])
    def test_v3_stations_unchanged(self, name):
        from shore.collar import build_collar

        a, b, seed, kw = golden_case(name)
        steps = {"steps_a": G_STEPS_S, "steps_b": G_STEPS_S, "steps_k": G_STEPS_K}
        block = build_collar(a, b, seed, **steps, params=params("stations"), **kw)
        np.testing.assert_allclose(block.nodes, np.load(GOLDEN)[name], rtol=0.0, atol=1e-12)

    def test_v4_unknown_volume(self):
        from shore.exceptions import ParameterError

        with pytest.raises(ParameterError, match="volume"):
            collar("rectangle", 32, 1.9, 0.0, volume="elliptic")

    def test_v5_project(self, tmp_path):
        import trimesh

        from shore.collar import build_collar
        from shore.geometry import load_geometry
        from shore.project import Project

        plate = trimesh.creation.box(extents=(4.0, 4.0, 1.0))
        plate.apply_translation((0.0, 0.0, -0.5))
        trimesh.creation.box(extents=(1.0, 0.1, 1.0)).export(str(tmp_path / "prism.stl"))
        plate.export(str(tmp_path / "plate.stl"))
        seed = (0.0, 0.06, 0.0)
        kernel = build_collar(
            load_geometry(tmp_path / "plate.stl"), load_geometry(tmp_path / "prism.stl"),
            np.asarray(seed), n_loop=48, steps_a=G_STEPS_S, steps_b=G_STEPS_S,
            steps_k=G_STEPS_K, params=params("smoothed"),
        )  # fmt: skip
        p = Project()
        p.load_geometry(tmp_path / "plate.stl", name="plate")
        p.build_collar(
            tmp_path / "prism.stl", geometry="plate", seed=seed, n_loop=48, ds_surface=0.01,
            growth_surface=1.2, m_a=8, m_b=8, ds_wall=1e-4, growth_wall=1.3, nk=12,
            volume="smoothed",
        )  # fmt: skip
        np.testing.assert_array_equal(p.mesh.blocks[0].nodes, kernel.nodes)
        assert "volume='smoothed'" in p.journal()
