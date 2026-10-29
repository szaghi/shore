"""FROZEN acceptance spec: redistributing the nodes along a block edge (ADR 0041).

Frozen 2026-10-29, before any implementation. Implementing work may not edit
this file, ``tests/spec/__init__.py``, ``tests/conftest.py``, the pytest
configuration in ``pyproject.toml`` or ``src/shore/quality.py``. The ONE
permitted edit is deleting an ``@PENDING`` line once that case passes (strict
xfail turns a passing pending case into a failure). ``tests/spec/SHA256SUMS``
holds the hashes, taken with ``@PENDING`` lines removed;
``scripts/check_spec.py`` checks them.

Verdict for acceptance (main CI job): ``pytest tests/spec/test_redistribute.py``
reports 15 passed, 0 skipped, 0 xfailed. At freeze: 0 passed, 15 xfailed.

Interface under test
--------------------
- ``Project.redistribute(block, edge, law="uniform", *, beta=3.0,
  cluster_at=None, first_step=None, last_step=None, interior="uniform")``,
  the command ``Redistribute``; journal ``p.redistribute('box',
  'j_lo_k_lo', 'tanh')`` (defaults left out).
- The law (``uniform``, ``tanh``, ``tanh2``, ``pinned``) is
  ``Spacing1D(law, ...).build(n, 0, L)`` over the picked edge's arc length
  ``L``, taken as arc-length fractions ``t``.
- Every grid line of the block along the edge's axis is resampled at ``t``
  along its own arc length, its nodes on its old polyline, its two ends
  kept, so the two faces across the axis do not move.
- Through SHARED and DIRICHLET seams, every block sharing a face that
  contains the axis is redistributed along its mapped axis (reversed across
  a flipped seam), transitively; seams stay closed.
- ``ParameterError`` for:
  - an unknown block, edge or law;
  - a law that gives non-increasing nodes;
  - a face shared with several blocks on the way (CH's ``body.k_hi``).

  Project and journal are then unchanged.
- Undo restores the nodes bit for bit.

Cases
-----
- D1 (x 4) each law on a straight box edge: the nodes are the law's; pinned
  holds its first and last cell.
- D2 a curved grid (a marched O-grid along i): every line's new nodes lie on
  its old polyline, at arc-length fractions ``t`` of it.
- D3 the faces across the axis are unchanged, bit for bit.
- D4 a cubed sphere along i: the four equatorial blocks change, the caps do
  not, every seam is closed.
- D5 the C-grid's flipped wake cut: the lower block takes the reversed
  fractions, and the seam opens no further than after the build.
- D6 (x 5) refusals (bad names, law or pinned steps x 4; a shared face).
- D7 undo and the journal line.
- D8 a cubed sphere loaded without replay (seams as copies): closed.

Out of scope of this spec: the GUI (``tests_gui/spec/test_gui_spacing.py``),
snapping associated nodes after a redistribution, periodic seams (refused),
the march's k schedule (``march`` keeps its own).
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

from shore.exceptions import ParameterError

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-29; implementation pending")

BOX = {"vmin": (0.0, 0.0, 0.0), "vmax": (2.0, 1.0, 1.0), "ni": 9, "nj": 5, "nk": 5}
OG = {"nk": 5, "ds": 0.05, "growth": 1.1, "ni": 9, "nj": 16}
CG = {"nk": 6, "ds": 1e-3, "ni_body": 40, "ni_wake": 10, "n_stations": 4, "wake_length": 5.0}
CH = {"chord_min": -2.0, "chord_max": 15.0, "normal_min": -2.5, "normal_max": 2.5, "nk_channel": 4}


def box():
    from shore.project import Project

    p = Project()
    p.build_primitive("box", **BOX)
    return p


def marched(name, sphere_stl, naca0012_stl):
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
        else:
            p.build_topology(name, **OG)
        p.march()
    return p


def arc_params(points: np.ndarray, polyline: np.ndarray) -> np.ndarray:
    """Arc-length fraction along *polyline* (m, 3) of each of *points* (n, 3),
    from the closest point of the polyline."""
    a, b = polyline[:-1], polyline[1:]
    ab = b - a
    seg = np.linalg.norm(ab, axis=1)
    start = np.concatenate([[0.0], np.cumsum(seg)])[:-1]
    out = []
    for p in points:
        t = np.clip(np.einsum("ij,ij->i", p - a, ab) / np.einsum("ij,ij->i", ab, ab), 0, 1)
        i = int(np.argmin(np.linalg.norm(a + t[:, None] * ab - p, axis=1)))
        out.append(start[i] + t[i] * seg[i])
    return np.asarray(out) / seg.sum()


def on_polyline(points: np.ndarray, polyline: np.ndarray) -> float:
    """Largest distance of *points* (n, 3) from the polyline (m, 3)."""
    a, b = polyline[:-1], polyline[1:]
    ab = b - a
    worst = 0.0
    for p in points:
        t = np.clip(np.einsum("ij,ij->i", p - a, ab) / np.einsum("ij,ij->i", ab, ab), 0, 1)
        worst = max(worst, float(np.linalg.norm(a + t[:, None] * ab - p, axis=1).min()))
    return worst


def seam_gap(blocks) -> float:
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


# ---------------------------------------------------------------------------
# D1-D3: one block
# ---------------------------------------------------------------------------

LAWS = [
    ("uniform", {}),
    ("tanh", {"beta": 3.0}),
    ("tanh2", {"beta": 2.0, "cluster_at": 0.3}),
    ("pinned", {"first_step": 0.05, "last_step": 0.4, "interior": "uniform"}),
]


class TestOneBlock:
    @PENDING
    @pytest.mark.parametrize(("law", "params"), LAWS)
    def test_d1_laws(self, law, params):
        from shore.mesh.spacing import Spacing1D

        p = box()
        p.redistribute("box", "j_lo_k_lo", law, **params)
        x = p.mesh.blocks[0].nodes[:, 0, 0, 0]
        expected = Spacing1D(law, **params).build(9, 0.0, 2.0)  # type: ignore[arg-type]
        np.testing.assert_allclose(x, expected, rtol=0.0, atol=1e-12)
        if law == "pinned":
            assert abs((x[1] - x[0]) - 0.05) <= 1e-12 and abs((x[-1] - x[-2]) - 0.4) <= 1e-12

    @PENDING
    def test_d2_curved_lines(self, sphere_stl, naca0012_stl):
        from shore.mesh.spacing import Spacing1D

        p = marched("ogrid", sphere_stl, naca0012_stl)
        old = p.mesh.blocks[0].nodes.copy()
        p.redistribute("ogrid", "j_lo_k_lo", "tanh2", beta=2.0)
        new = p.mesh.blocks[0].nodes
        edge = old[:, 0, 0]
        length = float(np.linalg.norm(np.diff(edge, axis=0), axis=1).sum())
        t = Spacing1D("tanh2", beta=2.0).build(9, 0.0, length) / length
        lines_new = np.moveaxis(new, 0, -2)  # (nj, nk, ni, 3)
        lines_old = np.moveaxis(old, 0, -2)
        for j in range(0, old.shape[1], 3):
            for k in range(old.shape[2]):
                assert on_polyline(lines_new[j, k], lines_old[j, k]) <= 1e-12
                # resampled at t along the OLD line (chords of the new one are shorter)
                np.testing.assert_allclose(
                    arc_params(lines_new[j, k], lines_old[j, k]), t, rtol=0.0, atol=1e-12
                )

    @PENDING
    def test_d3_faces_across_unchanged(self):
        p = box()
        before = p.mesh.blocks[0].nodes.copy()
        p.redistribute("box", "i_lo_k_lo", "tanh", beta=2.5)  # along j
        after = p.mesh.blocks[0].nodes
        assert np.array_equal(after[:, 0], before[:, 0]) and np.array_equal(
            after[:, -1], before[:, -1]
        )
        assert not np.array_equal(after, before)


# ---------------------------------------------------------------------------
# D4-D5, D8: through seams
# ---------------------------------------------------------------------------


class TestSeams:
    @PENDING
    def test_d4_cubed_sphere(self, sphere_stl, naca0012_stl):
        from shore.mesh.spacing import Spacing1D

        p = marched("cubed_sphere", sphere_stl, naca0012_stl)
        before = {b.label: b.nodes.copy() for b in p.mesh.blocks}
        p.redistribute("sub0", "j_lo_k_lo", "tanh2", beta=2.0)
        edge = before["sub0"][:, 0, 0]
        length = float(np.linalg.norm(np.diff(edge, axis=0), axis=1).sum())
        t = Spacing1D("tanh2", beta=2.0).build(len(edge), 0.0, length) / length
        blocks = {b.label: b for b in p.mesh.blocks}
        for label in ("sub0", "sub1", "sub2", "sub3"):
            assert not np.array_equal(blocks[label].nodes, before[label]), label
            new = np.moveaxis(blocks[label].nodes, 0, -2)
            old = np.moveaxis(before[label], 0, -2)
            for j in (0, new.shape[0] // 2):
                for k in (0, new.shape[1] - 1):
                    np.testing.assert_allclose(
                        arc_params(new[j, k], old[j, k]), t, rtol=0.0, atol=1e-12
                    )
        for label in ("cap_north", "cap_south"):
            assert np.array_equal(blocks[label].nodes, before[label])
        assert seam_gap(p.mesh.blocks) == 0.0

    @PENDING
    def test_d5_flipped_wake_cut(self, sphere_stl, naca0012_stl):
        p = marched("cgrid", sphere_stl, naca0012_stl)
        built = seam_gap(p.mesh.blocks)
        before = {b.label: b.nodes.copy() for b in p.mesh.blocks}
        p.redistribute("upper", "j_lo_k_hi", "tanh", beta=2.0)
        blocks = {b.label: b for b in p.mesh.blocks}

        def params(label):  # an i-line off the seam (k = 2): the seam sync cannot set it
            return arc_params(blocks[label].nodes[:, 0, 2], before[label][:, 0, 2])

        np.testing.assert_allclose(
            params("lower"), 1.0 - params("upper")[::-1], rtol=0.0, atol=1e-12
        )
        assert seam_gap(p.mesh.blocks) <= built

    @PENDING
    def test_d8_loaded_without_replay(self, sphere_stl, naca0012_stl, tmp_path):
        from shore.project import Project

        marched("cubed_sphere", sphere_stl, naca0012_stl).save(tmp_path / "saved")
        q = Project.load(tmp_path / "saved", replay=False)
        q.redistribute("sub2", "j_hi_k_hi", "tanh", beta=2.0)
        assert seam_gap(q.mesh.blocks) == 0.0


# ---------------------------------------------------------------------------
# D6-D7: refusals, undo, journal
# ---------------------------------------------------------------------------


class TestRefusals:
    @PENDING
    @pytest.mark.parametrize(
        "args",
        [
            (("nope", "j_lo_k_lo", "tanh"), {}),
            (("box", "k_lo", "tanh"), {}),
            (("box", "j_lo_k_lo", "spline"), {}),
            (("box", "j_lo_k_lo", "pinned"), {"first_step": 1.5, "last_step": 1.5}),
        ],
    )
    def test_d6_bad_arguments(self, args):
        p = box()
        before = p.mesh.blocks[0].nodes.copy()
        journal = p.journal()
        with pytest.raises(ParameterError):
            p.redistribute(*args[0], **args[1])
        assert p.journal() == journal and np.array_equal(p.mesh.blocks[0].nodes, before)

    @PENDING
    def test_d6_shared_face(self, sphere_stl, naca0012_stl):
        p = marched("ch", sphere_stl, naca0012_stl)
        before = [b.nodes.copy() for b in p.mesh.blocks]
        journal = p.journal()
        with pytest.raises(ParameterError, match="several"):
            p.redistribute("body", "j_lo_k_lo", "tanh")
        assert p.journal() == journal
        assert all(np.array_equal(b.nodes, n) for b, n in zip(p.mesh.blocks, before, strict=True))


class TestUndo:
    @PENDING
    def test_d7_undo_and_journal(self, sphere_stl, naca0012_stl):
        p = marched("cubed_sphere", sphere_stl, naca0012_stl)
        before = [b.nodes.copy() for b in p.mesh.blocks]
        p.redistribute("sub0", "j_lo_k_lo", "tanh", beta=2.0)
        assert (
            p.journal().splitlines()[-1] == "p.redistribute('sub0', 'j_lo_k_lo', 'tanh', beta=2.0)"
        )
        p.undo()
        assert all(np.array_equal(b.nodes, n) for b, n in zip(p.mesh.blocks, before, strict=True))
