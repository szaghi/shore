"""FROZEN acceptance spec: spacing laws along the junction curve (ADR 0020).

Frozen 2026-10-27, before any implementation. Implementing work may not edit
this file, ``tests/spec/data/loop_default.npz``, ``tests/spec/__init__.py``,
``tests/conftest.py``, the pytest configuration in ``pyproject.toml``, or
``src/shore/quality.py``. The ONE permitted edit is deleting an ``@PENDING``
line once that case passes (strict xfail turns a passing pending case into a
failure). ``tests/spec/SHA256SUMS`` holds the hashes, taken with ``@PENDING``
lines removed; ``scripts/check_spec.py`` checks them.

Verdict for acceptance (main CI job): ``pytest tests/spec/test_loop_spacing.py``
reports 29 passed, 0 skipped, 0 xfailed. At freeze: 5 passed (S2, the
default outputs against the golden file: regression guards), 24 xfailed.

Interface under test
--------------------
``intersection_loop(body_a, body_b, seed, n, *, spacing=None, curvature=0.0,
...)`` and ``intersection_curve(..., plane, spacing=None, curvature=0.0, ...)``
(ADR 0009, 0011 plus):

- ``spacing`` is a ``shore.mesh.spacing.Spacing1D``: ``uniform`` or
  ``tanh2`` (``beta``, optional ``cluster_at`` in [0, 1]) on a loop;
  ``uniform``, ``tanh`` or ``tanh2`` on an open curve. With ``curvature = 0``
  node *k* sits at arc-length fraction ``spacing.build(m)[k]`` of the curve,
  ``m = n + 1`` on a loop (the last coordinate, the seed again, dropped) and
  ``m = n`` on an open curve; node 0 is where the seed lands (loop) or the
  first plane end (curve), as before;
- ``curvature = alpha >= 0`` clusters nodes where the curve bends: the law
  is applied to a coordinate that equidistributes ``1 + alpha * kappa * L /
  (2 pi)`` (``kappa`` the curvature, averaged over half a uniform cell;
  ``L`` the curve's length; on an open curve both measured on its mirror
  image through the plane too), and adjacent cells differ by at most 1.2
  (plus the projection of the nodes onto the bodies);
- ``spacing=None`` and ``curvature=0`` (the defaults) are today's curves,
  unchanged; ``Spacing1D("uniform")`` with ``curvature=0`` is the same;
- ``build_collar(..., loop_spacing=None, loop_curvature=0.0)`` forwards
  them to the curve kernel it calls (ADR 0010, 0013);
- ``Project.build_collar(..., loop_law="uniform", loop_beta=3.0,
  loop_cluster_at=None, loop_curvature=0.0)`` (``BuildCollar``, journaled)
  and ``shore collar --loop-law --loop-beta --loop-cluster-at
  --loop-curvature`` build ``Spacing1D(loop_law, beta=loop_beta,
  cluster_at=loop_cluster_at)``;
- ``ParameterError`` for a law outside those above (``tanh`` on a loop,
  ``geometric``, ``pinned``), ``cluster_at`` outside [0, 1], and a negative
  or non-finite ``curvature``.

Cases
-----
- S1 law matched, ``curvature = 0``: plane x cylinder (the exact circle,
  arc fraction = angle / 2 pi), n = 64, ``tanh2`` beta 3 and beta 4 at 0.3:
  every node within 2e-6 of its law fraction (prototype 2.0e-7), on both
  bodies within 1e-9 x diagonal, node 0 where the default loop's is; open
  half circle (plane y = 0), n = 33, ``tanh`` beta 2 and ``tanh2`` beta 3:
  within 2e-6 of the law (prototype 1.5e-7), ends on the plane.
- S2 default unchanged: the four loops of ADR 0009 (n = 64) and the open
  half circle (n = 33) equal ``data/loop_default.npz`` (written by the code
  before this ADR) within 1e-12; ``Spacing1D("uniform")`` gives the default
  bit for bit.
- S3 curvature: on an elliptic cylinder ``x^2 + (y / 0.2)^2 = 1`` through a
  plane (curvature ratio 125, the leading edge of a thick wing section),
  n = 128, alpha = 1: the largest turning angle at a node at most 0.25 x the
  uniform loop's (prototype 5.5 deg vs 42.6 deg), adjacent cells within
  1.21 (prototype 1.2001); the half curve (plane y = 0), n = 65: end cells
  at most 0.25 x the uniform ones (prototype 0.17), cells within 1.21; on
  the circle (constant curvature), alpha = 3 leaves the loop uniform within
  2e-6 (prototype 8.7e-8).
- S4 collar: the wall composes the spaced curve (as ADR 0010's K1 and
  ADR 0013's H1, <= 1e-12); on the elliptic junction with ADR 0010's
  steps, 0 folded cells and min scaled Jacobian >= 0.1 for alpha = 1
  (prototype 0.23; the uniform loop: 0.008), ``tanh2`` beta 2 with
  alpha = 1 (0.17), and the half-collar (n = 65, alpha = 1: 0.16; uniform
  0.012).
- S5 command and CLI: a box with an elliptic STL post (256 sections),
  n_loop 96, ``tanh2`` beta 2, alpha 1: one valid collar (0 folded, min
  scaled Jacobian >= 0.05, prototype 0.11; the uniform loop folds 4 cells),
  journal line, bit-exact
  replay, equal to the kernel; ``shore collar --loop-*`` equals the API.
- S6 errors (x 7).

Out of scope of this spec: junction curves with corners (a sharp or blunt
trailing edge: the tracer of ADR 0009 stalls there, a separate step), user
clustering at given 3-D points, spacing laws on the surface-marching or
wall-normal directions (they exist), run time.
"""

from __future__ import annotations

import functools
from pathlib import Path

import numpy as np
import pytest
import trimesh

from shore.exceptions import ParameterError
from shore.mesh.spacing import Spacing1D
from shore.quality import scaled_jacobian
from tests.spec.test_intersection import PAIRS, AnalyticBody, CylinderZ, Plane

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-27; implementation pending")

GOLDEN = Path(__file__).parent / "data" / "loop_default.npz"
SYM = ((0.0, 0.0, 0.0), (0.0, 1.0, 0.0))  # the symmetry plane y = 0
CIRCLE_SEED = np.array([1.1, 0.1, 0.05])
HALF_CIRCLE_SEED = np.array([0.2, 1.1, 0.05])
ELLIPSE_SEED = np.array([1.05, 0.0, 0.0])
HALF_ELLIPSE_SEED = np.array([0.0, 0.25, 0.0])
STEPS_S = 0.01 * 1.12 ** np.arange(15)  # ADR 0010's collar steps
STEPS_K = 1e-4 * 1.2 ** np.arange(30)
DIAG = float(np.linalg.norm([4.0, 4.0, 4.0]))  # the analytic bodies' bounds


def loop_kernel():
    from shore.surface.intersect import intersection_loop

    return intersection_loop


def curve_kernel():
    from shore.surface.intersect import intersection_curve

    return intersection_curve


def collar_kernel():
    from shore.collar import build_collar

    return build_collar


class EllipticCylinder(AnalyticBody):
    """``x^2 / a^2 + y^2 / b^2 = 1``, z-parallel; exact foot by bisection."""

    def __init__(self, a: float = 1.0, b: float = 0.2) -> None:
        self.a, self.b = a, b

    def foot(self, p: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        a, b = self.a, self.b
        sx = np.where(p[..., 0] < 0.0, -1.0, 1.0)
        sy = np.where(p[..., 1] < 0.0, -1.0, 1.0)
        y0, y1 = np.abs(p[..., 0]), np.abs(p[..., 1])
        # the root t >= -b^2 of (a y0 / (t + a^2))^2 + (b y1 / (t + b^2))^2 = 1
        lo = np.full_like(y0, -b * b)
        hi = a * np.hypot(y0, y1) + 1.0
        for _ in range(200):
            t = 0.5 * (lo + hi)
            f = (a * y0 / (t + a * a)) ** 2 + (b * y1 / (t + b * b)) ** 2 - 1.0
            lo, hi = np.where(f > 0.0, t, lo), np.where(f > 0.0, hi, t)
        t = 0.5 * (lo + hi)
        x0, x1 = sx * a * a * y0 / (t + a * a), sy * b * b * y1 / (t + b * b)
        n = np.stack([x0 / (a * a), x1 / (b * b), np.zeros_like(x0)], axis=-1)
        n /= np.linalg.norm(n, axis=-1, keepdims=True)
        return np.stack([x0, x1, p[..., 2]], axis=-1), n


# ---------------------------------------------------------------------------
# Frozen measures
# ---------------------------------------------------------------------------


def angle_fraction(nodes: np.ndarray, period: float) -> np.ndarray:
    """|polar angle - angle of node 0| / period, unwrapped (circle about z)."""
    theta = np.unwrap(np.arctan2(nodes[:, 1], nodes[:, 0]))
    return np.abs(theta - theta[0]) / period


def cells(nodes: np.ndarray, closed: bool) -> np.ndarray:
    pts = np.vstack([nodes, nodes[:1]]) if closed else nodes
    return np.linalg.norm(np.diff(pts, axis=0), axis=1)


def cell_ratio(nodes: np.ndarray, closed: bool) -> float:
    h = cells(nodes, closed)
    h2 = np.r_[h, h[:1]] if closed else h
    return float(np.max(np.maximum(h2[1:] / h2[:-1], h2[:-1] / h2[1:])))


def max_turn_deg(nodes: np.ndarray) -> float:
    """Largest angle between consecutive chords of a closed loop."""
    t1 = nodes - np.roll(nodes, 1, axis=0)
    t2 = np.roll(nodes, -1, axis=0) - nodes
    c = np.einsum("ij,ij->i", t1, t2) / np.linalg.norm(t1, axis=1) / np.linalg.norm(t2, axis=1)
    return float(np.degrees(np.arccos(np.clip(c, -1.0, 1.0))).max())


def wrapped_jacobian(nodes: np.ndarray) -> np.ndarray:
    return scaled_jacobian(np.concatenate([nodes, nodes[:, :1]], axis=1))


def composed_wall(a, b, curve, *, plane=None) -> np.ndarray:
    """ADR 0010's K1 / ADR 0013's H1 composition on a given curve."""
    from shore.surface.march import march_surface

    march = {"periodic": True} if plane is None else {"plane": plane}
    on_a = march_surface(a, curve, STEPS_S, **march)
    on_b = march_surface(b, curve[::-1], STEPS_S, **march)[::-1]
    wall = np.concatenate([on_b[:, ::-1], on_a[:, 1:]], axis=1)
    return np.transpose(wall, (1, 0, 2))[:, ::-1]


def ellipse() -> tuple[Plane, EllipticCylinder]:
    return Plane(), EllipticCylinder()


@functools.cache
def ellipse_collar(spacing: tuple | None, curvature: float, half: bool):
    a, b = ellipse()
    law = None if spacing is None else Spacing1D(*spacing[:1], **dict(spacing[1:]))
    if half:
        return collar_kernel()(
            a, b, HALF_ELLIPSE_SEED, n_loop=65, steps_a=STEPS_S, steps_b=STEPS_S,
            steps_k=STEPS_K, plane=SYM, loop_spacing=law, loop_curvature=curvature,
        )  # fmt: skip
    return collar_kernel()(
        a, b, ELLIPSE_SEED, n_loop=128, steps_a=STEPS_S, steps_b=STEPS_S, steps_k=STEPS_K,
        loop_spacing=law, loop_curvature=curvature,
    )  # fmt: skip


# ---------------------------------------------------------------------------
# S1: the law, without curvature
# ---------------------------------------------------------------------------

LOOP_LAWS = {
    "tanh2": Spacing1D("tanh2", beta=3.0),
    "tanh2-at-0.3": Spacing1D("tanh2", beta=4.0, cluster_at=0.3),
}
CURVE_LAWS = {"tanh": Spacing1D("tanh", beta=2.0), "tanh2": Spacing1D("tanh2", beta=3.0)}


class TestLaw:
    @pytest.mark.parametrize("law", list(LOOP_LAWS))
    def test_s1_loop_follows_the_law(self, law):
        a, b = Plane(), CylinderZ()
        spacing = LOOP_LAWS[law]
        loop = loop_kernel()(a, b, CIRCLE_SEED, 64, spacing=spacing)
        frac = np.r_[angle_fraction(loop, 2.0 * np.pi), 1.0]
        assert np.all(np.diff(frac) > 0.0)
        assert np.abs(frac - spacing.build(65, 0.0, 1.0)).max() <= 2e-6
        assert a.distance(loop).max() <= 1e-9 * DIAG
        assert b.distance(loop).max() <= 1e-9 * DIAG
        default = loop_kernel()(a, b, CIRCLE_SEED, 64)
        assert np.array_equal(loop[0], default[0])

    @pytest.mark.parametrize("law", list(CURVE_LAWS))
    def test_s1_open_curve_follows_the_law(self, law):
        a, b = Plane(), CylinderZ()
        spacing = CURVE_LAWS[law]
        curve = curve_kernel()(a, b, HALF_CIRCLE_SEED, 33, plane=SYM, spacing=spacing)
        frac = angle_fraction(curve, np.pi)
        assert np.abs(frac - spacing.build(33, 0.0, 1.0)).max() <= 2e-6
        assert np.abs(curve[[0, -1], 1]).max() <= 1e-9 * DIAG
        assert a.distance(curve).max() <= 1e-9 * DIAG
        assert b.distance(curve).max() <= 1e-9 * DIAG


# ---------------------------------------------------------------------------
# S2: the default is unchanged
# ---------------------------------------------------------------------------


class TestDefault:
    @pytest.mark.parametrize("name", list(PAIRS))
    def test_s2_default_loops_unchanged(self, name):
        make_a, make_b, seed, _, _ = PAIRS[name]
        loop = loop_kernel()(make_a(), make_b(), np.asarray(seed, dtype=np.float64), 64)
        np.testing.assert_allclose(loop, np.load(GOLDEN)[name], rtol=0.0, atol=1e-12)

    def test_s2_default_curve_unchanged(self):
        curve = curve_kernel()(Plane(), CylinderZ(), HALF_CIRCLE_SEED, 33, plane=SYM)
        golden = np.load(GOLDEN)["plane_cylinder_curve"]
        np.testing.assert_allclose(curve, golden, rtol=0.0, atol=1e-12)

    def test_s2_uniform_law_is_the_default(self):
        a, b = Plane(), CylinderZ()
        uniform = Spacing1D("uniform")
        loop = loop_kernel()(a, b, CIRCLE_SEED, 64, spacing=uniform, curvature=0.0)
        assert np.array_equal(loop, loop_kernel()(a, b, CIRCLE_SEED, 64))
        curve = curve_kernel()(a, b, HALF_CIRCLE_SEED, 33, plane=SYM, spacing=uniform)
        assert np.array_equal(curve, curve_kernel()(a, b, HALF_CIRCLE_SEED, 33, plane=SYM))


# ---------------------------------------------------------------------------
# S3: curvature clustering
# ---------------------------------------------------------------------------


class TestCurvature:
    def test_s3_loop_clusters_where_it_bends(self):
        a, b = ellipse()
        uniform = loop_kernel()(a, b, ELLIPSE_SEED, 128)
        loop = loop_kernel()(a, b, ELLIPSE_SEED, 128, curvature=1.0)
        assert max_turn_deg(loop) <= 0.25 * max_turn_deg(uniform)
        assert cell_ratio(loop, closed=True) <= 1.21
        assert a.distance(loop).max() <= 1e-9 * DIAG
        assert b.distance(loop).max() <= 1e-9 * DIAG

    def test_s3_open_curve_clusters_at_the_plane(self):
        a, b = ellipse()
        uniform = curve_kernel()(a, b, HALF_ELLIPSE_SEED, 65, plane=SYM)
        curve = curve_kernel()(a, b, HALF_ELLIPSE_SEED, 65, plane=SYM, curvature=1.0)
        h, h0 = cells(curve, closed=False), cells(uniform, closed=False)
        assert h[0] <= 0.25 * h0[0]
        assert h[-1] <= 0.25 * h0[-1]
        assert cell_ratio(curve, closed=False) <= 1.21
        assert np.abs(curve[[0, -1], 1]).max() <= 1e-9 * DIAG

    def test_s3_constant_curvature_stays_uniform(self):
        loop = loop_kernel()(Plane(), CylinderZ(), CIRCLE_SEED, 64, curvature=3.0)
        frac = np.r_[angle_fraction(loop, 2.0 * np.pi), 1.0]
        assert np.abs(frac - np.linspace(0.0, 1.0, 65)).max() <= 2e-6


# ---------------------------------------------------------------------------
# S4: the collar
# ---------------------------------------------------------------------------


class TestCollar:
    def test_s4_wall_composes_the_spaced_loop(self):
        a, b = ellipse()
        law = Spacing1D("tanh2", beta=2.0)
        block = ellipse_collar(("tanh2", ("beta", 2.0)), 1.0, False)
        loop = loop_kernel()(a, b, ELLIPSE_SEED, 128, spacing=law, curvature=1.0)
        expected = composed_wall(a, b, loop)
        assert np.abs(block.nodes[:, :, 0] - expected).max() <= 1e-12

    def test_s4_half_wall_composes_the_spaced_curve(self):
        a, b = ellipse()
        block = ellipse_collar(None, 1.0, True)
        curve = curve_kernel()(a, b, HALF_ELLIPSE_SEED, 65, plane=SYM, curvature=1.0)
        expected = composed_wall(a, b, curve, plane=SYM)
        assert np.abs(block.nodes[:, :, 0] - expected).max() <= 1e-12

    @pytest.mark.parametrize(
        "spacing", [None, ("tanh2", ("beta", 2.0))], ids=["curvature", "tanh2+curvature"]
    )
    def test_s4_valid_collar(self, spacing):
        jac = wrapped_jacobian(ellipse_collar(spacing, 1.0, False).nodes)
        assert (jac <= 0.0).sum() == 0
        assert jac.min() >= 0.1

    def test_s4_valid_half_collar(self):
        jac = scaled_jacobian(ellipse_collar(None, 1.0, True).nodes)
        assert (jac <= 0.0).sum() == 0
        assert jac.min() >= 0.1


# ---------------------------------------------------------------------------
# S5: command and CLI
# ---------------------------------------------------------------------------

COMMAND = {
    "seed": (0.0, 0.11, 0.02), "n_loop": 96, "ds_surface": 0.005, "growth_surface": 1.12,
    "m_a": 15, "m_b": 15, "ds_wall": 5e-5, "growth_wall": 1.2, "nk": 30,
}  # fmt: skip
LOOP = {"loop_law": "tanh2", "loop_beta": 2.0, "loop_curvature": 1.0}


@pytest.fixture(scope="module")
def stl_files(tmp_path_factory) -> dict[str, Path]:
    """A 4 x 4 x 1 box (top face at z = 0) pierced by an elliptic post (1 x 0.2)."""
    d = tmp_path_factory.mktemp("loop_spacing")
    box = trimesh.creation.box(extents=(4.0, 4.0, 1.0))
    box.apply_translation((0.0, 0.0, -0.5))
    post = trimesh.creation.cylinder(radius=1.0, height=1.6, sections=256)
    post.apply_scale((0.5, 0.1, 1.0))
    post.apply_translation((0.0, 0.0, 0.2))
    round_post = trimesh.creation.cylinder(radius=0.5, height=1.6, sections=48)
    round_post.apply_translation((0.0, 0.0, 0.2))
    for name, mesh in (("box", box), ("post", post), ("round", round_post)):
        mesh.export(str(d / f"{name}.stl"))
    return {k: d / f"{k}.stl" for k in ("box", "post", "round")}


@functools.cache
def command_project(box: Path, post: Path):
    from shore.project import Project

    p = Project()
    p.load_geometry(box)
    p.build_collar(str(post), **COMMAND, **LOOP)
    return p


class TestCommand:
    def test_s5_valid_collar_component(self, stl_files):
        p = command_project(stl_files["box"], stl_files["post"])
        jac = wrapped_jacobian(p.mesh.blocks[0].nodes)
        assert (jac <= 0.0).sum() == 0
        assert jac.min() >= 0.05

    def test_s5_equals_the_kernel(self, stl_files):
        from shore.geometry import load_geometry

        p = command_project(stl_files["box"], stl_files["post"])
        c = COMMAND
        block = collar_kernel()(
            load_geometry(stl_files["box"]), load_geometry(stl_files["post"]),
            np.asarray(c["seed"]), n_loop=c["n_loop"],
            steps_a=c["ds_surface"] * c["growth_surface"] ** np.arange(c["m_a"]),
            steps_b=c["ds_surface"] * c["growth_surface"] ** np.arange(c["m_b"]),
            steps_k=c["ds_wall"] * c["growth_wall"] ** np.arange(c["nk"]),
            loop_spacing=Spacing1D("tanh2", beta=2.0), loop_curvature=1.0,
        )  # fmt: skip
        assert np.array_equal(p.mesh.blocks[0].nodes, block.nodes)

    def test_s5_journal_and_replay(self, stl_files, tmp_path):
        from shore.project import Project

        p = command_project(stl_files["box"], stl_files["post"])
        line = next(ln for ln in p.journal().splitlines() if "build_collar(" in ln)
        assert "loop_law='tanh2'" in line and "loop_curvature=1.0" in line
        p.save(tmp_path / "proj")
        q = Project.load(tmp_path / "proj")  # replays and checks every grid bit for bit
        assert np.array_equal(q.mesh.blocks[0].nodes, p.mesh.blocks[0].nodes)

    def test_s5_cli_equals_the_api(self, stl_files, tmp_path):
        from typer.testing import CliRunner

        from shore.cli import app
        from shore.project import Project

        small = {
            "seed": (0.55, 0.0, 0.02), "n_loop": 16, "ds_surface": 0.02,
            "growth_surface": 1.15, "m_a": 5, "m_b": 5, "ds_wall": 1e-3, "growth_wall": 1.25,
            "nk": 8,
        }  # fmt: skip
        loop = {"loop_law": "tanh2", "loop_beta": 2.0, "loop_cluster_at": 0.5,
                "loop_curvature": 1.0}  # fmt: skip
        args = ["collar", str(stl_files["box"]), str(stl_files["round"]),
                "--seed", *map(str, small["seed"]), "--project", str(tmp_path / "cli")]  # fmt: skip
        for key, value in {**small, **loop}.items():
            if key != "seed":
                args += ["--" + key.replace("_", "-"), str(value)]
        r = CliRunner().invoke(app, args)
        assert r.exit_code == 0, r.output
        api = Project()
        api.load_geometry(stl_files["box"], name="box")
        api.load_geometry(stl_files["round"], name="round")
        api.use_component("collar")
        api.build_collar("round", geometry="box", **small, **loop)
        cli = Project.load(tmp_path / "cli")
        cli.select_component("collar")
        assert np.array_equal(cli.mesh.blocks[0].nodes, api.mesh.blocks[0].nodes)


# ---------------------------------------------------------------------------
# S6: errors
# ---------------------------------------------------------------------------


class TestErrors:
    @pytest.mark.parametrize(
        "kwargs",
        [
            {"spacing": Spacing1D("tanh", beta=2.0)},
            {"spacing": Spacing1D("geometric", ds=0.01)},
            {"spacing": Spacing1D("pinned", first_step=0.01, last_step=0.01)},
            {"spacing": Spacing1D("tanh2", beta=3.0, cluster_at=1.5)},
            {"curvature": -1.0},
            {"curvature": float("nan")},
        ],
        ids=["tanh-on-loop", "geometric", "pinned", "cluster-at", "negative", "nan"],
    )
    def test_s6_kernel_errors(self, kwargs):
        with pytest.raises(ParameterError):
            loop_kernel()(Plane(), CylinderZ(), CIRCLE_SEED, 16, **kwargs)

    def test_s6_project_unknown_law(self, stl_files):
        from shore.project import Project

        p = Project()
        p.load_geometry(stl_files["box"])
        with pytest.raises(ParameterError):
            p.build_collar(str(stl_files["round"]), **{**COMMAND, "n_loop": 16}, loop_law="cubic")
        assert p.mesh is None
