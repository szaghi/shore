"""FROZEN acceptance spec: background primitives have no wall (ADR 0055).

Frozen 2026-10-31, before any implementation. Implementing work may not edit
this file, ``tests/spec/__init__.py``, ``tests/conftest.py``, the pytest
configuration in ``pyproject.toml`` or ``src/shore/quality.py``. The ONE
permitted edit is deleting an ``@PENDING`` line once that case passes (strict
xfail turns a passing pending case into a failure). ``tests/spec/SHA256SUMS``
holds the hashes, taken with ``@PENDING`` lines removed;
``scripts/check_spec.py`` checks them.

Verdict for acceptance (main CI job): ``pytest tests/spec/test_background_bcs.py``
reports 8 passed, 0 skipped, 0 xfailed. At freeze: 2 passed, 6 xfailed.

Interface under test
--------------------
The background primitives ``box``, ``annulus`` and ``flat_caps``
(``Project.build_primitive``, ``shore primitive``) have no WALL face: every
face that is not a seam is FREE. Before ADR 0055 they inherited
``HexBlock``'s default ``k_lo = WALL`` (meant for marched topologies, whose
body is there), so their far-field ``k_lo`` was written as a viscous wall:
``cc.par`` BC 1, CGNS ``BCWall``; ``free_face_bc`` could not change it.
The duct primitives keep their walls, and marched grids their body walls.

Cases
-----
- B1 (x 3) box, annulus, flat-caps: no WALL face, every non-seam face FREE.
- B2 a flat-caps background through ``write_xall`` with ``free_face_bc`` 9:
  no patch with BC 1, every boundary patch 9.
- B3 a box through CGNS: no ``BCWall``.
- B4 the H-grid duct (rect: j and k faces) and tube (k_hi) keep their walls
  (a guard; holds at freeze).
- B5 a marched cubed sphere keeps its body wall at k_lo (a guard).
- B6 ``shore primitive box --project``: the saved face graph has no WALL
  (a guard against a second path; pending at freeze).

Out of scope of this spec: the H-grid tube's ``k_lo`` (its singular axis,
marked WALL by design there), user-chosen walls on primitives.
"""

from __future__ import annotations

import json

import pytest
from typer.testing import CliRunner

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-31; implementation pending")

BACKGROUNDS = {
    "box": {"vmin": (0.0, 0.0, 0.0), "vmax": (1.0, 2.0, 3.0), "ni": 4, "nj": 4, "nk": 4},
    "annulus": {"r_in": 0.5, "r_out": 2.0, "z0": 0.0, "z1": 1.0, "ni": 4, "nj": 9, "nk": 4},
    "flat_caps": {
        "r_out": 3.0,
        "z0": -1.0,
        "z1": 4.0,
        "ni": 5,
        "nc": 6,
        "nk": 5,
        "inner_frac": 0.4,
    },
}


def primitive(kind, **params):
    from shore.project import Project

    p = Project()
    p.build_primitive(kind, **(params or BACKGROUNDS[kind]))
    return p


def bcs(p):
    return [(b.label, f.name, f.bc.name) for b in p.mesh.blocks for f in b.faces()]


class TestBackgrounds:
    @PENDING
    @pytest.mark.parametrize("kind", sorted(BACKGROUNDS))
    def test_b1_no_wall(self, kind):
        faces = bcs(primitive(kind))
        assert not [f for f in faces if f[2] == "WALL"]
        assert {f[2] for f in faces} <= {"FREE", "SHARED"}

    @PENDING
    def test_b2_cc_par(self, tmp_path):
        from tests.spec.test_cc_par_xall import parse

        p = primitive("flat_caps")
        p.write_xall(tmp_path / "run", components={"main": {"free_face_bc": 9}})
        _, patches = parse((tmp_path / "run" / "cc.par").read_text())
        codes = [patch[2] for patch in patches]  # (block, face, code, ...)
        assert 1 not in codes
        assert 9 in codes

    @PENDING
    def test_b3_cgns(self, tmp_path):
        import CGNS.MAP

        primitive("box").write_cgns(tmp_path / "box.cgns")
        found = []

        def walk(node):
            if node[3] == "BC_t":
                found.append(node[1].tobytes().decode())
            for child in node[2]:
                walk(child)

        walk(CGNS.MAP.load(str(tmp_path / "box.cgns"))[0])
        assert found and "BCWall" not in found

    def test_b4_ducts_keep_walls(self):
        rect = primitive("h_grid_rect", length=4.0, width=1.0, height=1.0, ni=5, nj=4, nk=4)
        assert {rect.mesh.blocks[0].face(f).bc.name for f in ("j_lo", "j_hi", "k_lo", "k_hi")} == {
            "WALL"
        }
        tube = primitive("h_grid_tube", length=4.0, radius=1.0, n_axial=5, n_circ=8, n_radial=4)
        assert tube.mesh.blocks[0].face("k_hi").bc.name == "WALL"

    def test_b5_marched_wall(self, sphere_stl):
        from shore.project import Project

        p = Project()
        p.load_geometry(sphere_stl)
        p.build_topology("cubed_sphere", nk=5, ds=0.05, growth=1.1, ni=9, nj=16)
        assert all(b.face("k_lo").bc.name == "WALL" for b in p.mesh.blocks)

    @PENDING
    def test_b6_cli_project(self, tmp_path):
        from shore.cli import app

        r = CliRunner().invoke(app, [
            "primitive", "box", "--min", "0", "0", "0", "--max", "1", "1", "1",
            "--ni", "4", "--nj", "4", "--nk", "4", "-o", str(tmp_path / "bg"),
            "--project", str(tmp_path / "proj"),
        ])  # fmt: skip
        assert r.exit_code == 0, r.output
        graph = json.loads((tmp_path / "proj" / "grids" / "main" / "adjacency.json").read_text())
        assert "WALL" not in json.dumps(graph)
