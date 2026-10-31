"""FROZEN acceptance spec: loading a saved project without replay (ADR 0033).

Frozen 2026-10-28, before any implementation. Implementing work may not edit
this file, ``tests/spec/__init__.py``, ``tests/conftest.py``, the pytest
configuration in ``pyproject.toml``, or ``src/shore/quality.py``. The ONE
permitted edit is deleting an ``@PENDING`` line once that case passes (strict
xfail turns a passing pending case into a failure). ``tests/spec/SHA256SUMS``
holds the hashes, taken with ``@PENDING`` lines removed;
``scripts/check_spec.py`` checks them.

Verdict for acceptance (main CI job): ``pytest tests/spec/test_fast_load.py``
reports 10 passed, 0 skipped, 0 xfailed. At freeze: 1 passed (the replayed
load still verifies the grids: a regression guard), 9 xfailed.

Interface under test
--------------------
- ``Project.save`` also writes each non-empty component's face graph,
  ``grids/<component>/adjacency.json`` (project format unchanged, 3).
- ``Project.load(directory, replay=False)``: the input copies are checked
  against their hashes; the journal is recorded, not replayed (geometries
  loaded, components selected, no mesh rebuilt); each component's blocks are
  read from its grids and rewired from its face graph. The result has the
  same components, current component, geometries, blocks (labels, nodes bit
  for bit, face BCs, partners and axis maps) and journal as a replayed
  load; saving it or writing its Xall run gives the same files.
- The recorded history cannot be undone (``ParameterError`` naming the
  replay); a loaded mesh cannot be marched (``ParameterError``); a project
  without the face graphs cannot be loaded without replay
  (``MeshInputError`` naming the replay).
- ``shore xall DIR --no-replay`` (and ``--no-replay`` with ``--project`` on
  the mesh, collar and primitive commands) loads without replay.
- ``Project.load(directory)`` is unchanged: it replays and verifies.

Cases
-----
A project with two components: ``main``, a flat-caps primitive; ``collar``,
the collar of an STL plate (a 4 x 4 x 1 box) and an STL prism (a 1 x 0.1 x 1
box), n_loop 48, surface steps ``0.01 * 1.2**j`` (8), wall-normal
``1e-4 * 1.3**k`` (12). L1-L7 as the interface above; L8 the replayed load
still rejects a changed grid.

Out of scope of this spec: load times (measured in ADR 0033), the copies'
renaming on re-save (``00_00_plate.stl``, an older behaviour).

Re-frozen 2026-10-31 (ADR 0054, approved by the owner): projects are
format 4; L1 checks ``format = 4``.
"""

from __future__ import annotations

import filecmp
import re
from pathlib import Path

import numpy as np
import pytest

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-28; implementation pending")

SEED = (0.0, 0.06, 0.0)


def make_project(d: Path):
    import trimesh

    from shore.project import Project

    plate = trimesh.creation.box(extents=(4.0, 4.0, 1.0))
    plate.apply_translation((0.0, 0.0, -0.5))
    plate.export(str(d / "plate.stl"))
    trimesh.creation.box(extents=(1.0, 0.1, 1.0)).export(str(d / "prism.stl"))
    p = Project()
    p.build_primitive("flat_caps", r_out=3.0, z0=-1.0, z1=4.0, ni=5, nc=6, nk=5, inner_frac=0.4)
    p.use_component("collar")
    p.load_geometry(d / "plate.stl", name="plate")
    p.build_collar(
        d / "prism.stl", geometry="plate", seed=SEED, n_loop=48, ds_surface=0.01,
        growth_surface=1.2, m_a=8, m_b=8, ds_wall=1e-4, growth_wall=1.3, nk=12,
    )  # fmt: skip
    p.save(d / "proj")
    return d / "proj"


@pytest.fixture
def saved(tmp_path) -> Path:
    return make_project(tmp_path)


def load(directory: Path, replay: bool):
    from shore.project import Project

    return Project.load(directory, replay=replay) if not replay else Project.load(directory)


def face_graph(mesh):
    """Per block: (label, [(face, bc, partner label or None, partner face, axis map)])."""
    owner = {id(f): b.label for b in mesh.blocks for f in b.faces()}
    out = []
    for b in mesh.blocks:
        faces = []
        for f in b.faces():
            partner = f.partner
            faces.append(
                (
                    f.name,
                    f.bc.name,
                    None if partner is None else owner.get(id(partner)),
                    None if partner is None else partner.name,
                    None if f.axis_map is None else f.axis_map.along,
                )
            )
        out.append((b.label, faces))
    return out


class TestSave:
    @PENDING
    def test_l1_face_graphs_saved(self, saved):
        for name in ("main", "collar"):
            assert (saved / "grids" / name / "adjacency.json").is_file()
        assert re.search(r"^format = 4$", (saved / "project.toml").read_text(), re.M)


class TestLoad:
    def test_l2_same_state(self, saved):
        full, fast = load(saved, True), load(saved, False)
        assert list(fast.components) == list(full.components)
        assert fast.current_component == full.current_component
        assert sorted(fast.geometries) == sorted(full.geometries)
        for name, mesh in full.components.items():
            other = fast.components[name]
            assert [b.label for b in other.blocks] == [b.label for b in mesh.blocks]
            for a, b in zip(mesh.blocks, other.blocks, strict=True):
                assert np.array_equal(a.nodes, b.nodes)
            assert face_graph(other) == face_graph(mesh)

    def test_l3_no_kernel_runs(self, saved, monkeypatch):
        import shore.collar
        import shore.pipeline

        def boom(*args, **kwargs):
            raise AssertionError("a mesh was rebuilt")

        monkeypatch.setattr(shore.collar, "build_collar", boom)
        monkeypatch.setattr(shore.collar, "split_collar", boom)
        monkeypatch.setattr(shore.pipeline, "build_primitive", boom)
        fast = load(saved, False)
        assert fast.components["collar"] is not None

    def test_l4_journal_and_resave(self, saved, tmp_path):
        full, fast = load(saved, True), load(saved, False)
        assert fast.journal() == full.journal()
        full.save(tmp_path / "full")
        fast.save(tmp_path / "fast")
        cmp = filecmp.dircmp(tmp_path / "full", tmp_path / "fast")
        assert not cmp.diff_files and not cmp.left_only and not cmp.right_only
        for sub in cmp.subdirs.values():
            assert not sub.diff_files and not sub.left_only and not sub.right_only

    def test_l5_xall_run(self, saved, tmp_path):
        load(saved, True).write_xall(tmp_path / "full")
        load(saved, False).write_xall(tmp_path / "fast")
        files = sorted(p.relative_to(tmp_path / "full") for p in (tmp_path / "full").rglob("*"))
        assert files == sorted(
            p.relative_to(tmp_path / "fast") for p in (tmp_path / "fast").rglob("*")
        )
        for rel in files:
            if (tmp_path / "full" / rel).is_file():
                assert (tmp_path / "full" / rel).read_bytes() == (
                    tmp_path / "fast" / rel
                ).read_bytes(), rel


class TestRefusals:
    def test_l6_undo_refused(self, saved):
        from shore.exceptions import ParameterError

        fast = load(saved, False)
        with pytest.raises(ParameterError, match="replay"):
            fast.undo()

    def test_l6_march_refused(self, saved):
        from shore.exceptions import ParameterError

        fast = load(saved, False)
        with pytest.raises(ParameterError):
            fast.march()

    def test_l6_needs_face_graphs(self, saved):
        from shore.exceptions import MeshInputError

        for name in ("main", "collar"):
            (saved / "grids" / name / "adjacency.json").unlink(missing_ok=True)
        with pytest.raises(MeshInputError, match="replay"):
            load(saved, False)


class TestCli:
    def test_l7_xall_no_replay(self, saved, tmp_path):
        from typer.testing import CliRunner

        from shore.cli import app

        runner = CliRunner()
        r = runner.invoke(app, ["xall", str(saved), "--out", str(tmp_path / "a"), "--no-replay"])
        assert r.exit_code == 0, r.output
        r = runner.invoke(app, ["xall", str(saved), "--out", str(tmp_path / "b")])
        assert r.exit_code == 0, r.output
        for f in (tmp_path / "a").rglob("*"):
            if f.is_file():
                assert (
                    f.read_bytes() == (tmp_path / "b" / f.relative_to(tmp_path / "a")).read_bytes()
                )


class TestReplayedLoad:
    def test_l8_replay_still_verifies(self, saved):  # already holds: a regression guard
        from shore.exceptions import MeshInputError

        grid = saved / "grids" / "collar" / "0000.npz"
        with np.load(grid) as g:
            nodes = g["nodes"].copy()
        nodes[1, 1, 1, 0] += 1e-9
        np.savez(grid, nodes=nodes)
        with pytest.raises(MeshInputError, match="differs"):
            load(saved, True)
