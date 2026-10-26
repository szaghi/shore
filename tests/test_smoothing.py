"""Regression tests for the limited Laplacian smoothing (ADR 0005 finding).

Unlimited smoothing of node positions shrinks a curved layer inward by about
``strength * h**2 * curvature`` regardless of the step, so small first steps
folded (negative Jacobian at k=1). ``smooth_mode='limited'`` caps that
motion per node; ``'legacy'`` keeps the old behaviour.
"""

from __future__ import annotations

import numpy as np
import pytest
import trimesh
from typer.testing import CliRunner

from shore.cli import app
from shore.exceptions import MarchingError, ParameterError
from shore.mesh import CubedSphereTopology, Mesh, OGridTopology
from shore.mesh.spacing import Spacing1D
from shore.topology.sphere_cap import _apply_smooth_mode
from shore.volume.hyperbolic import _limit_smoothing

runner = CliRunner()


def _wall_tilt_deg(mesh: Mesh) -> float:
    """Largest angle between the first k-segment and the radial direction."""
    cos = []
    for b in mesh.blocks:
        e = b.nodes[:, :, 1] - b.nodes[:, :, 0]
        r = b.nodes[:, :, 0] / np.linalg.norm(b.nodes[:, :, 0], axis=-1, keepdims=True)
        cos.append((np.einsum("ijk,ijk->ij", e, r) / np.linalg.norm(e, axis=-1)).ravel())
    return float(np.rad2deg(np.arccos(np.clip(np.concatenate(cos).min(), -1.0, 1.0))))


# ---------------------------------------------------------------------------
# The limiter itself
# ---------------------------------------------------------------------------


class TestLimiter:
    n = np.array([[[0.0, 0.0, 1.0]]])
    raw = np.zeros((1, 1, 3))

    def test_untouched_when_within_limits(self):
        smoothed = np.array([[[0.01, 0.0, -0.01]]])  # step 0.1: within 0.25 * step
        out = _limit_smoothing(self.raw, smoothed, self.n, step=0.1)
        assert out is smoothed

    def test_pullback_capped(self):
        smoothed = np.array([[[0.0, 0.0, -1.0]]])
        out = _limit_smoothing(self.raw, smoothed, self.n, step=0.1)
        np.testing.assert_allclose(out, [[[0.0, 0.0, -0.025]]])

    def test_skew_capped(self):
        smoothed = np.array([[[1.0, 0.0, 0.0]]])
        out = _limit_smoothing(self.raw, smoothed, self.n, step=0.1)
        np.testing.assert_allclose(out, [[[0.025, 0.0, 0.0]]])

    def test_outward_motion_kept(self):
        smoothed = np.array([[[0.0, 0.0, 1.0]]])  # pushes forward: never limited
        out = _limit_smoothing(self.raw, smoothed, self.n, step=0.1)
        assert out is smoothed

    def test_modes(self):
        smoothed = np.array([[[0.0, 0.0, -1.0]]])
        assert _apply_smooth_mode(self.raw, smoothed, self.n, 0.1, "legacy") is smoothed
        with pytest.raises(ParameterError, match="smooth_mode"):
            _apply_smooth_mode(self.raw, smoothed, self.n, 0.1, "magic")


# ---------------------------------------------------------------------------
# The defect: small first steps on a coarse curved grid
# ---------------------------------------------------------------------------

_TANH = Spacing1D(law="tanh", beta=2.0, total_thickness=0.2)
_CASES = {
    "cubed ds=1e-3": (CubedSphereTopology, {"ds": 1e-3}),
    "cubed ds=1e-5": (CubedSphereTopology, {"ds": 1e-5, "growth": 1.3}),
    "cubed tanh": (CubedSphereTopology, {"steps": _TANH.build_steps(8)}),
    "ogrid ds=1e-3": (OGridTopology, {"ds": 1e-3}),
}


class TestSmallFirstStep:
    @pytest.mark.parametrize("case", list(_CASES))
    def test_legacy_folds(self, sphere_stl, case):
        topo, kw = _CASES[case]
        mesh = Mesh.from_file(sphere_stl, topo(), ni=10, nj=12, nk=8)
        with pytest.raises(MarchingError, match="Negative Jacobian"):
            mesh.march(smooth_mode="legacy", **kw)

    @pytest.mark.parametrize("case", list(_CASES))
    def test_limited_marches(self, sphere_stl, case):
        topo, kw = _CASES[case]
        mesh = Mesh.from_file(sphere_stl, topo(), ni=10, nj=12, nk=8)
        mesh.march(**kw)  # default: limited
        # 19-24 deg measured; unsmoothed ~8 deg, the (0.5, 1.0) limits gave 64-68
        assert _wall_tilt_deg(mesh) < 25.0

    def test_limited_equals_legacy_where_safe(self, sphere_stl):
        # the regression grid the limits were chosen on: nothing binds there
        kw = {
            "topology": CubedSphereTopology(),
            "ni": 60,
            "nj": 60,
            "nk": 10,
            "theta_cap": float(np.pi / 6),
            "spacing_k": Spacing1D(law="geometric", ds=0.02, growth=1.1),
        }
        legacy = Mesh.from_file(sphere_stl, **kw)
        legacy.march(ds=0.02, growth=1.1, smooth_mode="legacy")
        limited = Mesh.from_file(sphere_stl, **kw)
        limited.march(ds=0.02, growth=1.1)
        for a, b in zip(legacy.blocks, limited.blocks, strict=True):
            np.testing.assert_array_equal(a.nodes, b.nodes)


class TestCli:
    _ARGS = ("--ni", "10", "--nj", "12", "--nk", "8", "--ds", "1e-3")

    def test_default_marches_small_first_step(self, sphere_stl, tmp_path):
        r = runner.invoke(app, ["mesh", str(sphere_stl), *self._ARGS, "-o", str(tmp_path / "s")])
        assert r.exit_code == 0, r.output

    def test_legacy_flag(self, sphere_stl, tmp_path):
        r = runner.invoke(
            app,
            ["mesh", str(sphere_stl), *self._ARGS, "--smooth-mode", "legacy",
             "-o", str(tmp_path / "s")],
        )  # fmt: skip
        assert r.exit_code == 1
        assert "Negative Jacobian" in r.output

    def test_bad_mode(self, sphere_stl, tmp_path):
        r = runner.invoke(
            app, ["mesh", str(sphere_stl), "--smooth-mode", "x", "-o", str(tmp_path / "s")]
        )
        assert r.exit_code == 1

    def test_config_key(self, sphere_stl, tmp_path):
        cfg = tmp_path / "c.toml"
        cfg.write_text(f'version = 1\ninput = "{sphere_stl}"\n\n[march]\nsmooth_mode = "legacy"\n')
        r = runner.invoke(
            app,
            ["mesh", str(sphere_stl), *self._ARGS, "-c", str(cfg), "-o", str(tmp_path / "s")],
        )
        assert r.exit_code == 1  # legacy from the config folds this case
        assert "Negative Jacobian" in r.output


def test_icosphere_fixture_is_the_documented_one(sphere_stl):
    # the measurements quoted above are on a subdivision-3 unit icosphere
    m = trimesh.load(sphere_stl)
    assert len(m.faces) == 1280
