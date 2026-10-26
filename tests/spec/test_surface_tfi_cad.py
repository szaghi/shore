"""FROZEN acceptance spec, CAD half: the surface-grid command on a STEP face (ADR 0007).

Frozen 2026-10-26 with ``test_surface_tfi.py`` (same rules: the only
permitted edit is deleting an ``@PENDING`` line once that case passes).
Needs the ``cad`` extra; the CI CAD job runs it.

Verdict for acceptance (CAD job): ``pytest tests/spec/test_surface_tfi_cad.py``
reports 4 passed, 0 skipped, 0 xfailed. At freeze: 0 passed, 4 xfailed.

Interface under test
--------------------
``Project.surface_grid(block, face, *, surface, smoothing="elliptic")``
rebuilds the nodes of one face of one block from that face's boundary, on
SURFACE entity number *surface* of the loaded geometry, with
``shore.surface.tfi.surface_grid``. It is a journaled, undoable command:
``Project.load`` replays it bit-identically. ``ParameterError`` for an
unknown block or face name, and for a geometry without (u, v) (STL).

Case: a STEP cylinder (r = 0.5) and an annulus-sector primitive whose inner
face (``i_lo``) lies on it, clustered in both directions (tanh2 in theta,
tanh in z) and clear of the cylinder's periodic seam. The cylinder is
developable and its (u, v) = (angle, height) is a scaled isometry, so the
primitive's own tensor-product face is the exact clustered elliptic
solution: the command must reproduce it (<= 0.05 local spacing), with every
node on the CAD face (exact projection, <= 1e-9 x diagonal), and must leave
every other node of the block bit-identical.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import trimesh

bd = pytest.importorskip("build123d")

from shore.exceptions import ParameterError  # noqa: E402
from shore.geometry import BODY, EntityKind, EntityRef  # noqa: E402
from shore.project import Project  # noqa: E402

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-26; implementation pending")

ANNULUS = {
    "r_in": 0.5, "r_out": 1.5, "z0": -1.5, "z1": 1.5, "ni": 9, "nj": 25, "nk": 21,
    "theta_start": 0.5, "theta_end": 2.0, "j_law": "tanh2", "j_beta": 2.0,
    "k_law": "tanh", "k_beta": 2.5,
}  # fmt: skip


@pytest.fixture(scope="module")
def files(tmp_path_factory) -> dict[str, Path]:
    d = tmp_path_factory.mktemp("tfi_cad")
    bd.export_step(bd.Cylinder(radius=0.5, height=4.0), str(d / "cylinder.step"))
    trimesh.creation.cylinder(radius=0.5, height=4.0, sections=64).export(str(d / "cylinder.stl"))
    return {"step": d / "cylinder.step", "stl": d / "cylinder.stl"}


def _project(path: Path) -> Project:
    p = Project()
    p.load_geometry(path)
    p.build_primitive("annulus", **ANNULUS)
    return p


def _cylinder_face(p: Project) -> int:
    """Index of the SURFACE entity under the annulus' inner wall."""
    a = 0.5 * (ANNULUS["theta_start"] + ANNULUS["theta_end"])
    probe = np.array([0.5 * np.cos(a), 0.5 * np.sin(a), 0.0])
    return int(p.body.geometry.project(probe, BODY).entity_index)


def _local_spacing(ref: np.ndarray) -> np.ndarray:
    di = np.linalg.norm(np.diff(ref, axis=0), axis=-1)
    dj = np.linalg.norm(np.diff(ref, axis=1), axis=-1)
    return np.minimum.reduce([di[:-1, 1:-1], di[1:, 1:-1], dj[1:-1, :-1], dj[1:-1, 1:]])


class TestCommand:
    def test_rebuilds_face_on_cad_surface(self, files):
        p = _project(files["step"])
        block = p.mesh.blocks[0]
        before = block.nodes.copy()
        k = _cylinder_face(p)
        p.surface_grid(block.label, "i_lo", surface=k)
        face = block.nodes[0]
        ref = before[0]
        dev = np.linalg.norm(face - ref, axis=-1)[1:-1, 1:-1] / _local_spacing(ref)
        assert dev.max() <= 0.05
        proj = p.body.geometry.project(face, EntityRef(EntityKind.SURFACE, k), exact=True)
        diag = float(np.linalg.norm(np.ptp(before.reshape(-1, 3), axis=0)))
        assert proj.distances.max() <= 1e-9 * diag
        np.testing.assert_array_equal(block.nodes[1:], before[1:])

    def test_journal_replay_and_undo(self, files, tmp_path):
        p = _project(files["step"])
        before = p.mesh.blocks[0].nodes.copy()
        p.surface_grid(p.mesh.blocks[0].label, "i_lo", surface=_cylinder_face(p))
        assert "surface_grid(" in p.journal()
        after = p.mesh.blocks[0].nodes.copy()
        p.save(tmp_path / "proj")
        q = Project.load(tmp_path / "proj")  # replays and checks every grid bit for bit
        np.testing.assert_array_equal(q.mesh.blocks[0].nodes, after)
        p.undo()
        np.testing.assert_array_equal(p.mesh.blocks[0].nodes, before)

    def test_rejects_unknown_block_or_face(self, files):
        p = _project(files["step"])
        k = _cylinder_face(p)
        with pytest.raises(ParameterError):
            p.surface_grid("no_such_block", "i_lo", surface=k)
        with pytest.raises(ParameterError):
            p.surface_grid(p.mesh.blocks[0].label, "x_lo", surface=k)

    def test_rejects_stl(self, files):
        p = _project(files["stl"])
        with pytest.raises(ParameterError):
            p.surface_grid(p.mesh.blocks[0].label, "i_lo", surface=0)
