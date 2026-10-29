"""FROZEN acceptance spec: seams kept closed by in-place edits (ADR 0040).

Frozen 2026-10-29, before any implementation. Implementing work may not edit
this file, ``tests/spec/__init__.py``, ``tests/conftest.py``, the pytest
configuration in ``pyproject.toml`` or ``src/shore/quality.py``. The ONE
permitted edit is deleting an ``@PENDING`` line once that case passes (strict
xfail turns a passing pending case into a failure). ``tests/spec/SHA256SUMS``
holds the hashes, taken with ``@PENDING`` lines removed;
``scripts/check_spec.py`` checks them.

Verdict for acceptance (main CI job): ``pytest tests/spec/test_seam_sync.py``
reports 9 passed, 0 skipped, 0 xfailed. At freeze: 6 passed (S2 x 5, S4:
guards), 3 xfailed.

Interface under test
--------------------
- SHARED and DIRICHLET seams are one set of points. Built topologies share
  SHARED nodes as numpy views; DIRICHLET seams, and every seam of blocks read
  back from saved grids (ADR 0033), are copies.
- After ``snap`` (and ``fill_face``), each group of aliased nodes takes the
  value of its last member the edit changed, in block order, as views
  would. A group the edit did not change, or changed to the same value on
  every member, is left as it is.
- If a seam is not point-matched, ``snap`` warns (``UserWarning``) and
  snaps without closing seams.

Cases
-----
A cubed sphere (icosphere STL of ``tests/conftest.py``).
- S1 (x 2) one-sided associations across a SHARED seam (``sub0.j_hi``) and
  a DIRICHLET seam (``cap_north.i_lo``): after ``snap`` every seam is
  closed, both on the replayed project and on the one loaded without
  replay, and the two agree within 1e-12.
- S2 (x 5) wall-only snaps of built meshes (cubed sphere, O-grid, C-grid,
  CH, OH) give the same nodes, bit for bit, as projecting each block on its
  own (the snap before ADR 0040).
- S3 a seam that is not point-matched: ``snap`` warns and still snaps.
- S4 undo restores both sides of a seam.

Out of scope of this spec: ``fill_face`` (needs CAD; the same code path),
PERIODIC seams (not coincident), seams across components.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-29; implementation pending")

OG = {"nk": 5, "ds": 0.05, "growth": 1.1, "ni": 9, "nj": 16}
CG = {"nk": 6, "ds": 1e-3, "ni_body": 40, "ni_wake": 10, "n_stations": 4, "wake_length": 5.0}
CH = {"chord_min": -2.0, "chord_max": 15.0, "normal_min": -2.5, "normal_max": 2.5, "nk_channel": 4}
OH = {"nk_channel": 4, "box_vmin": (-4.0, -4.0, -4.0), "box_vmax": (4.0, 4.0, 4.0)}


def built(name, sphere_stl, naca0012_stl):
    from shore.project import Project

    p = Project()
    wing = name in ("cgrid", "ch")
    p.load_geometry(naca0012_stl if wing else sphere_stl, skip_hygiene=wing)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if name == "cgrid":
            p.build_topology("cgrid", **CG)
        elif name == "ch":
            p.build_topology("ch", options=CH, **{**CG, "nk": 8})
        elif name == "oh":
            p.build_topology("oh", options=OH, nk=5, ds=0.05, growth=1.1, ni=16, nj=16)
        else:
            p.build_topology(name, **OG)
        p.march()
    return p


def seam_gap(blocks) -> float:
    """Largest distance between the two sides of any SHARED / DIRICHLET seam."""
    from shore.mesh.face import FaceBC
    from shore.mesh.seams import _map_indices, seams

    owner = {b.label: b for b in blocks}
    worst = 0.0
    for label, sides in seams(blocks, 1.0, bcs=(FaceBC.SHARED, FaceBC.DIRICHLET)).items():
        for s in sides:
            src, dst = owner[label].nodes, owner[s.donor].nodes
            shape = tuple(int(n) for n in (s.rng[1] - s.rng[0] + 1))
            here = src[tuple(slice(lo - 1, hi) for lo, hi in zip(s.rng[0], s.rng[1], strict=True))]
            there = dst[_map_indices(s.donor_rng[0], s.transform, shape)]
            worst = max(worst, float(np.abs(here - there).max()))
    return worst


def one_sided(p) -> None:
    p.associate("sub0", "j_hi", "body", "body")  # SHARED with sub1.j_lo
    p.associate("cap_north", "i_lo", "body", "body")  # DIRICHLET with sub0.i_hi
    p.snap()


class TestOneSided:
    @PENDING
    @pytest.mark.parametrize("replay", [True, False])
    def test_s1_seams_closed(self, sphere_stl, naca0012_stl, tmp_path, replay):
        from shore.project import Project

        built("cubed_sphere", sphere_stl, naca0012_stl).save(tmp_path / "saved")
        q = Project.load(tmp_path / "saved", replay=replay)
        one_sided(q)
        assert seam_gap(q.mesh.blocks) == 0.0
        other = Project.load(tmp_path / "saved", replay=not replay)
        one_sided(other)
        for a, b in zip(q.mesh.blocks, other.mesh.blocks, strict=True):
            np.testing.assert_allclose(a.nodes, b.nodes, rtol=0.0, atol=1e-12)


class TestUnchanged:
    @pytest.mark.parametrize("name", ["cubed_sphere", "ogrid", "cgrid", "ch", "oh"])
    def test_s2_wall_snaps_bit_identical(self, sphere_stl, naca0012_stl, name):  # a guard
        from shore.mesh.association import AssociationResolver, snap_block

        p = built(name, sphere_stl, naca0012_stl)
        reference = built(name, sphere_stl, naca0012_stl)
        resolver = AssociationResolver({"body": reference.mesh.body.geometry})
        for block in reference.mesh.blocks:
            snap_block(block, resolver)
        p.snap()
        for a, b in zip(p.mesh.blocks, reference.mesh.blocks, strict=True):
            assert np.array_equal(a.nodes, b.nodes), a.label


class TestEdges:
    @PENDING
    def test_s3_unmatched_seam_warns(self, sphere_stl, naca0012_stl, tmp_path):
        from shore.project import Project

        built("cubed_sphere", sphere_stl, naca0012_stl).save(tmp_path / "saved")
        q = Project.load(tmp_path / "saved", replay=False)
        q.mesh.blocks[1].nodes[0, 2, 2] += 1e-3  # sub1.j_lo off its partner
        with pytest.warns(UserWarning, match="point-matched"):
            q.snap()
        assert q.journal().splitlines()[-1] == "p.snap()"

    def test_s4_undo_restores_both_sides(self, sphere_stl, naca0012_stl, tmp_path):  # a guard
        from shore.project import Project

        built("cubed_sphere", sphere_stl, naca0012_stl).save(tmp_path / "saved")
        q = Project.load(tmp_path / "saved", replay=False)
        before = [b.nodes.copy() for b in q.mesh.blocks]
        one_sided(q)
        q.undo()
        for b, nodes in zip(q.mesh.blocks, before, strict=True):
            assert np.array_equal(b.nodes, nodes)
