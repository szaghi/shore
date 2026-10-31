"""FROZEN acceptance spec: right-handed CH outflow blocks (ADR 0057).

Frozen 2026-10-31, before any implementation. Implementing work may not edit
this file, ``tests/spec/__init__.py``, ``tests/conftest.py``, the pytest
configuration in ``pyproject.toml`` or ``src/shore/quality.py``. The ONE
permitted edit is deleting an ``@PENDING`` line once that case passes (strict
xfail turns a passing pending case into a failure). ``tests/spec/SHA256SUMS``
holds the hashes, taken with ``@PENDING`` lines removed;
``scripts/check_spec.py`` checks them.

Verdict for acceptance (main CI job): ``pytest tests/spec/test_ch_handedness.py``
reports 6 passed, 0 skipped, 0 xfailed. At freeze: 2 passed, 4 xfailed.

Interface under test
--------------------
Every block of the CH hybrid is right-handed: positive scaled Jacobian in
every cell (``shore.quality``). Before ADR 0057 the outflow blocks above
the wake cut (``chan_outflow_top_upper``, ``chan_outflow_wake_upper``) and
the blunt TE's ``chan_outflow_wake_base`` were left-handed in every cell:
channel i along the chord, j along the span, k up. Their j now runs
against their partners', with the seams' axis maps saying so; their nodes
are the old ones with j reversed, every other block is unchanged, and the
seams still match node for node.

Cases
-----
- H1 (x 2) sharp and blunt TE: every cell of every block positive.
- H2 (x 2) every SHARED face matches its partner through its axis map,
  node for node, faces with several partners aside (a guard: holds before
  and after).
- H3 (x 2) the grid equals the one before ADR 0057
  (``data/ch_handedness_before.npz``, made at the base commit): the
  flipped blocks with j reversed, the others bit for bit.

Out of scope of this spec: the C-grid's wing-tip cap blocks (several are
mirrored, the LE ones degenerate, the TE ones without cells: a redesign
recorded in ADR 0057, not done), solver acceptance of left-handed blocks.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pytest

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-31; implementation pending")

DATA = Path(__file__).parent / "data"
CG = {"nk": 12, "ds": 1e-3, "ni_body": 40, "ni_wake": 10, "n_stations": 4, "wake_length": 5.0}
CH = {"chord_min": -2.0, "chord_max": 15.0, "normal_min": -2.5, "normal_max": 2.5, "nk_channel": 4}
FLIPPED = {"chan_outflow_top_upper", "chan_outflow_wake_upper", "chan_outflow_wake_base"}
ALONG = {"i": ("j", "k"), "j": ("i", "k"), "k": ("i", "j")}


@pytest.fixture(scope="module")
def grids(tmp_path_factory):
    from shore.profiles import naca
    from shore.project import Project
    from shore.shapes import extrude_section

    d = tmp_path_factory.mktemp("ch_hand")
    out = {}
    for te, closed in (("sharp", True), ("blunt", False)):
        stl = d / f"{te}.stl"
        extrude_section(naca("0012", 60, closed_te=closed), span=3.0, n_stations=5).export(stl)
        p = Project()
        p.load_geometry(stl, skip_hygiene=True)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            p.build_topology("ch", options=CH, **CG)
            p.march()
        out[te] = p.mesh.blocks
    return out


def face_nodes(nodes, face):
    return np.take(nodes, 0 if face.endswith("lo") else -1, axis="ijk".index(face[0]))


@pytest.mark.parametrize("te", ["sharp", "blunt"])
class TestHandedness:
    def test_h1_right_handed(self, grids, te):
        from shore.quality import scaled_jacobian

        bad = {b.label: float(scaled_jacobian(b.nodes).min()) for b in grids[te]}
        assert {k: v for k, v in bad.items() if v <= 0.0} == {}

    def test_h2_seams_match(self, grids, te):
        blocks = grids[te]
        owner = {id(f): b for b in blocks for f in b.faces()}
        compared = 0
        for b in blocks:
            for f in b.faces():
                if f.bc.name != "SHARED":
                    continue
                theirs = face_nodes(owner[id(f.partner)].nodes, f.partner.name)
                mine = face_nodes(b.nodes, f.name)
                (a0, flip0), (_, flip1) = f.axis_map.along
                t = theirs if ALONG[f.partner.name[0]].index(a0) == 0 else theirs.transpose(1, 0, 2)
                t = t[::-1] if flip0 else t
                t = t[:, ::-1] if flip1 else t
                if t.shape != mine.shape:  # a face shared with several partners
                    continue
                np.testing.assert_allclose(
                    mine, t, rtol=0, atol=1e-9, err_msg=f"{b.label}.{f.name}"
                )
                compared += 1
        assert compared >= 20

    def test_h3_same_grid(self, grids, te):
        with np.load(DATA / "ch_handedness_before.npz") as gold:
            for b in grids[te]:
                old = gold[f"{te}/{b.label}"]
                expect = old[:, ::-1] if b.label in FLIPPED else old
                assert np.array_equal(b.nodes, expect), b.label
