"""FROZEN acceptance spec: wing-tip caps on the tip face (ADR 0059).

Frozen 2026-10-31, before any implementation. Implementing work may not edit
this file, ``tests/spec/__init__.py``, ``tests/conftest.py``, the pytest
configuration in ``pyproject.toml`` or ``src/shore/quality.py``. The ONE
permitted edit is deleting an ``@PENDING`` line once that case passes (strict
xfail turns a passing pending case into a failure). ``tests/spec/SHA256SUMS``
holds the hashes, taken with ``@PENDING`` lines removed;
``scripts/check_spec.py`` checks them.

Verdict for acceptance (main CI job): ``pytest tests/spec/test_tip_cap.py``
reports 19 passed, 0 skipped, 0 xfailed. At freeze: 0 passed, 19 xfailed.

Interface under test
--------------------
A wing-tip cap closes one span end of a C-grid wing as an overset component:
five blocks per tip, ``cap_<tip>_{core,upper,le,lower,te}`` on the C-grid
(``cap_tips``), ``cap_{core,upper,le,lower,te}`` from ``Project.build_tip_cap``
/ ``shore cap-tip``. Their ``k_lo`` faces are the wall and together cover the
section exactly, on the wing's tip face; k leaves it outwards. The faces
between cap blocks are SHARED (16 per tip), every other face FREE. With
``cap_tips`` the C-grid's tip stations sit on the tip faces and every wall
node of the body's tip station is a cap wall node, bit for bit. Ring k-lines
lean outwards at the section's edge by ``cap_flare_deg`` (default 30).

Cases
-----
- C1 (x 2, sharp and blunt NACA 0012) five blocks per tip, every cell's
  minimum corner scaled Jacobian > 0.2; marching the wrap leaves them alone.
- C2 (x 2) 16 SHARED faces per tip, partners in the same cap, matching
  through their axis maps node for node.
- C3 (x 2) ``k_lo`` is the only WALL face; every other face SHARED or FREE.
- C4 (x 2) each tip's wall is one plane, within 1e-4 of the span of the
  STL's tip, and holds every wall node of the body's tip station.
- C5 (x 2) the wall's area is the section's (the tip station closed across
  the TE base), to 1e-9 relative.
- C6 the core's k-lines run along the wall's outward normal at the wrap's
  layer distances (the wall-normal schedule of ``ds`` and ``growth``).
- C7 (x 2) at the section's edge the first k segments lean outwards by
  ``cap_flare_deg``: 30 deg by default, 0 deg (along the normal) when 0.
- C8 (x 2) a swept, tapered, cambered wing (NACA 4412, sweep 20 deg, taper
  0.6): valid cells, the wall planar and holding the tip-station nodes.
- C9 ``Project.build_tip_cap`` on the body of a C-grid without caps: five
  valid blocks, 16 SHARED faces, the wall on the STL's tip face (its default
  inset: the C-grid's stations sit 2 % of the span inside the tips).
- C10 ``shore mesh --cap-tips both`` writes 14 valid ``.geo`` files;
  ``shore cap-tip`` writes 5.
- C11 a capped project saves and loads back (replay) bit for bit.
- C12 ``cap_te_frac`` is accepted with a ``DeprecationWarning`` and changes
  nothing.

Out of scope of this spec: rounded or non-planar tips, a conformal (face
matched) cap-to-wrap interface, Xall runs on the capped wing.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest
from typer.testing import CliRunner

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-31; implementation pending")

SPAN = 3.0
CGRID = {"nk": 12, "ds": 2e-3, "ni_body": 40, "ni_wake": 10, "n_stations": 4, "wake_length": 5.0}
PARTS = ("core", "upper", "le", "lower", "te")
ALONG = {"i": ("j", "k"), "j": ("i", "k"), "k": ("i", "j")}


@pytest.fixture(scope="module")
def wings(tmp_path_factory):
    from shore.profiles import naca
    from shore.shapes import extrude_section

    d = tmp_path_factory.mktemp("tip_cap")
    for name, code, closed, extra in (
        ("sharp", "0012", True, {}),
        ("blunt", "0012", False, {}),
        ("swept", "4412", False, {"sweep_deg": 20.0, "taper": 0.6}),
    ):
        wing = extrude_section(naca(code, 60, closed_te=closed), span=SPAN, n_stations=5, **extra)
        wing.export(d / f"{name}.stl")
    return d


def build(stl, march=True, **params):
    from shore.project import Project

    p = Project()
    p.load_geometry(stl, skip_hygiene=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        p.build_topology("cgrid", cap_tips="both", **{**CGRID, **params})
        if march:
            p.march()
    return p


@pytest.fixture(scope="module")
def capped(wings):
    return {
        te: build(wings / f"{te}.stl", **({"blocks": 4} if te == "blunt" else {}))
        for te in ("sharp", "blunt")
    }


def caps(blocks, tip):
    return [b for b in blocks if b.label.startswith(f"cap_{tip}_")]


def face_nodes(nodes, face):
    return np.take(nodes, 0 if face.endswith("lo") else -1, axis="ijk".index(face[0]))


def min_sj(blocks):
    from shore.quality import scaled_jacobian

    return min(float(scaled_jacobian(b.nodes).min()) for b in blocks)


def plane_of(points):
    centre = points.mean(axis=0)
    normal = np.linalg.svd(points - centre)[2][-1]
    return centre, normal / np.linalg.norm(normal)


def wall_points(blocks):
    return np.concatenate([b.nodes[:, :, 0].reshape(-1, 3) for b in blocks])


def outward(blocks, normal):
    """The wall normal on the side the cap grows on."""
    b = blocks[0]
    return normal if float((b.nodes[0, 0, 1] - b.nodes[0, 0, 0]) @ normal) > 0 else -normal


def seams_match(blocks):
    owner = {id(f): b for b in blocks for f in b.faces()}
    n = 0
    for b in blocks:
        for f in b.faces():
            if f.bc.name != "SHARED":
                continue
            assert id(f.partner) in owner, f"{b.label}.{f.name}: partner outside the cap"
            theirs = face_nodes(owner[id(f.partner)].nodes, f.partner.name)
            mine = face_nodes(b.nodes, f.name)
            (a0, flip0), (_, flip1) = f.axis_map.along
            t = theirs if ALONG[f.partner.name[0]].index(a0) == 0 else theirs.transpose(1, 0, 2)
            t = t[::-1] if flip0 else t
            t = t[:, ::-1] if flip1 else t
            assert t.shape == mine.shape, f"{b.label}.{f.name}"
            np.testing.assert_allclose(mine, t, rtol=0, atol=1e-12, err_msg=f"{b.label}.{f.name}")
            n += 1
    return n


def inside(ring, point, normal):
    """Whether *point* lies inside the closed polygon *ring* (even-odd rule, in its plane)."""
    e1 = ring[1] - ring[0]
    e1 = e1 - (e1 @ normal) * normal
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(normal, e1)
    x, y = (ring - point) @ e1, (ring - point) @ e2
    x2, y2 = np.roll(x, -1), np.roll(y, -1)
    crosses = (y > 0) != (y2 > 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        at = x + (0 - y) * (x2 - x) / (y2 - y)
    return bool(np.count_nonzero(crosses & (at > 0)) % 2)


def quad_area(face):
    a, b, c, d = face[:-1, :-1], face[1:, :-1], face[1:, 1:], face[:-1, 1:]
    return 0.5 * np.linalg.norm(np.cross(c - a, d - b), axis=-1).sum()


@pytest.mark.parametrize("te", ["sharp", "blunt"])
class TestConformalCaps:
    @PENDING
    def test_c1_blocks_valid(self, wings, capped, te):
        blocks = capped[te].mesh.blocks
        for tip in ("jlo", "jhi"):
            assert sorted(b.label for b in caps(blocks, tip)) == sorted(
                f"cap_{tip}_{n}" for n in PARTS
            )
        cap_blocks = [b for b in blocks if b.label.startswith("cap_")]
        assert len(cap_blocks) == 10
        assert min_sj(cap_blocks) > 0.2
        built = build(wings / f"{te}.stl", march=False, **({"blocks": 4} if te == "blunt" else {}))
        for b in cap_blocks:
            (same,) = [c for c in built.mesh.blocks if c.label == b.label]
            assert np.array_equal(b.nodes, same.nodes), b.label

    @PENDING
    def test_c2_seams(self, capped, te):
        for tip in ("jlo", "jhi"):
            assert seams_match(caps(capped[te].mesh.blocks, tip)) == 16

    @PENDING
    def test_c3_bcs(self, capped, te):
        for b in [b for b in capped[te].mesh.blocks if b.label.startswith("cap_")]:
            bcs = {f.name: f.bc.name for f in b.faces()}
            assert bcs.pop("k_lo") == "WALL", b.label
            assert set(bcs.values()) <= {"SHARED", "FREE"}, (b.label, bcs)

    @PENDING
    def test_c4_wall_on_tip(self, capped, te):
        blocks = capped[te].mesh.blocks
        body = next(b for b in blocks if b.label == "body")
        for tip, j, end in (("jlo", 0, 0.0), ("jhi", -1, SPAN)):
            pts = wall_points(caps(blocks, tip))
            centre, normal = plane_of(pts)
            assert np.abs((pts - centre) @ normal).max() < 1e-9
            assert abs(float(centre[2]) - end) < 1e-4 * SPAN
            keys = {p.tobytes() for p in pts}
            assert all(p.tobytes() in keys for p in body.nodes[:, j, 0])

    @PENDING
    def test_c5_wall_covers_section(self, capped, te):
        blocks = capped[te].mesh.blocks
        body = next(b for b in blocks if b.label == "body")
        for tip, j in (("jlo", 0), ("jhi", -1)):
            cap = caps(blocks, tip)
            ring = body.nodes[:, j, 0]
            _, normal = plane_of(ring)
            section = 0.5 * abs(
                float(np.cross(ring, np.roll(ring, -1, axis=0)).sum(axis=0) @ normal)
            )
            wall = sum(quad_area(b.nodes[:, :, 0]) for b in cap)
            assert abs(wall - section) <= 1e-9 * section

    @PENDING
    def test_c7_flare(self, wings, capped, te):
        extra = {"blocks": 4} if te == "blunt" else {}
        for p, angle in (
            (capped[te], 30.0),
            (build(wings / f"{te}.stl", cap_flare_deg=0.0, **extra), 0.0),
        ):
            blocks = p.mesh.blocks
            body = next(b for b in blocks if b.label == "body")
            for tip, j in (("jlo", 0), ("jhi", -1)):
                cap = caps(blocks, tip)
                _, normal = plane_of(wall_points(cap))
                normal = outward(cap, normal)
                edge = {q.tobytes() for q in body.nodes[:, j, 0]}
                seen = 0
                for b in cap:
                    w, nxt = b.nodes[:, :, 0].reshape(-1, 3), b.nodes[:, :, 1].reshape(-1, 3)
                    for a, c in zip(w, nxt, strict=True):
                        if a.tobytes() not in edge:
                            continue
                        seg = (c - a) / np.linalg.norm(c - a)
                        tilt = np.degrees(np.arccos(np.clip(seg @ normal, -1.0, 1.0)))
                        assert abs(tilt - angle) < 0.5, (b.label, tilt)
                        if angle > 0:  # a short step along the lean leaves the section
                            lean = seg - (seg @ normal) * normal
                            assert not inside(body.nodes[:, j, 0], a + 1e-6 * lean, normal)
                        seen += 1
                assert seen >= body.nodes.shape[0]


class TestLayers:
    @PENDING
    def test_c6_core_layers(self, capped):
        from shore.mesh.spacing import Spacing1D

        blocks = capped["sharp"].mesh.blocks
        steps = Spacing1D(law="geometric", ds=CGRID["ds"], growth=1.15).build_steps(CGRID["nk"])
        expected = np.concatenate([[0.0], np.cumsum(steps)])
        for tip in ("jlo", "jhi"):
            cap = caps(blocks, tip)
            _, normal = plane_of(wall_points(cap))
            normal = outward(cap, normal)
            core = next(b for b in cap if b.label.endswith("_core")).nodes
            rel = core - core[:, :, :1]
            np.testing.assert_allclose(
                np.linalg.norm(rel, axis=-1),
                np.broadcast_to(expected, rel.shape[:3]),
                rtol=1e-12,
                atol=1e-15,
            )
            along = rel[:, :, 1:] / np.linalg.norm(rel[:, :, 1:], axis=-1, keepdims=True)
            np.testing.assert_allclose(along @ normal, 1.0, atol=1e-12)


@pytest.mark.parametrize("tip", ["jlo", "jhi"])
class TestSwept:
    @PENDING
    def test_c8_swept_tapered_cambered(self, wings, tip):
        p = build(wings / "swept.stl", blocks=4)
        blocks = p.mesh.blocks
        cap = caps(blocks, tip)
        assert len(cap) == 5
        assert min_sj(cap) > 0.2
        pts = wall_points(cap)
        centre, normal = plane_of(pts)
        assert np.abs((pts - centre) @ normal).max() < 1e-9
        body = next(b for b in blocks if b.label == "body")
        keys = {q.tobytes() for q in pts}
        assert all(q.tobytes() in keys for q in body.nodes[:, 0 if tip == "jlo" else -1, 0])


@pytest.fixture(scope="module")
def body_geo(wings):
    """The body block of a sharp C-grid without caps, as ``shore mesh`` writes it."""
    from shore.project import Project

    p = Project()
    p.load_geometry(wings / "sharp.stl", skip_hygiene=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        p.build_topology("cgrid", **CGRID)
        p.march()
    p.write_geo(wings / "w", per_block=True)
    (path,) = wings.glob("w*body*.geo")
    return path


class TestStandalone:
    @PENDING
    def test_c9_build_tip_cap(self, body_geo):
        from shore.project import Project

        for tip, end in (("jlo", 0.0), ("jhi", SPAN)):
            p = Project()
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                p.build_tip_cap(body_geo, tip=tip)
            blocks = p.mesh.blocks
            assert sorted(b.label for b in blocks) == sorted(f"cap_{n}" for n in PARTS)
            assert min_sj(blocks) > 0.2
            assert seams_match(blocks) == 16
            pts = wall_points(blocks)
            centre, normal = plane_of(pts)
            assert np.abs((pts - centre) @ normal).max() < 1e-9
            assert abs(float(centre[2]) - end) < 1e-6 * SPAN

    @PENDING
    def test_c10_cli(self, wings, body_geo, tmp_path):
        from shore.cli import app
        from shore.io.geo import read_geo
        from shore.quality import scaled_jacobian

        runner = CliRunner()
        mesh_dir, cap_dir = tmp_path / "mesh", tmp_path / "cap"
        mesh_dir.mkdir()
        cap_dir.mkdir()
        args = [
            "mesh", str(wings / "sharp.stl"), "--topology", "cgrid", "--nk", "12", "--ds", "2e-3",
            "--ni-body", "40", "--ni-wake", "10", "--n-stations", "4", "--wake-length", "5",
            "--cap-tips", "both", "-o", str(mesh_dir / "w"),
        ]  # fmt: skip
        assert runner.invoke(app, args).exit_code == 0
        r = runner.invoke(app, ["cap-tip", str(body_geo), "--tip", "jhi", "-o", str(cap_dir / "c")])
        assert r.exit_code == 0, r.output
        for d, count in ((mesh_dir, 14), (cap_dir, 5)):
            files = sorted(d.glob("*.geo"))
            assert len(files) == count
            for f in files:
                assert float(scaled_jacobian(read_geo(f)).min()) > 0.0, f.name


class TestProjects:
    @PENDING
    def test_c11_save_load(self, capped, tmp_path):
        from shore.project import Project

        p = capped["sharp"]
        p.save(tmp_path / "proj")
        q = Project.load(tmp_path / "proj")
        assert [b.label for b in q.mesh.blocks] == [b.label for b in p.mesh.blocks]
        for a, b in zip(p.mesh.blocks, q.mesh.blocks, strict=True):
            assert np.array_equal(a.nodes, b.nodes), a.label

    @PENDING
    def test_c12_te_frac_deprecated(self, wings, capped):
        from shore.project import Project

        p = Project()
        p.load_geometry(wings / "sharp.stl", skip_hygiene=True)
        with pytest.warns(DeprecationWarning, match="cap_te_frac"):
            p.build_topology("cgrid", cap_tips="both", cap_te_frac=0.9, **CGRID)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            p.march()
        for a, b in zip(p.mesh.blocks, capped["sharp"].mesh.blocks, strict=True):
            assert np.array_equal(a.nodes, b.nodes), a.label
