"""FROZEN acceptance spec: cc.par connections as Xall reads them (ADR 0036).

Frozen 2026-10-28, before any implementation. Implementing work may not edit
this file, ``tests/spec/__init__.py``, ``tests/conftest.py``, the pytest
configuration in ``pyproject.toml``, ``src/shore/quality.py`` or the data
file below. The ONE permitted edit is deleting an ``@PENDING`` marker once
that case passes (strict xfail turns a passing pending case into a failure).
``tests/spec/SHA256SUMS`` holds the hashes, taken with ``@PENDING`` lines
removed; ``scripts/check_spec.py`` checks them.

Verdict for acceptance (main CI job): ``pytest tests/spec/test_cc_par_xall.py``
reports 45 passed, 0 skipped, 0 xfailed. At freeze: 31 passed (guards: Q0,
Q1 and Q2 guards, Q3, Q5), 14 xfailed.

What Xall reads
---------------
The reader of ``cc.par`` is the overset pre-processor (``Xnavis-projects/
overset/src``; the solvers read only its output). For a connection patch
``p`` (``BC > 99``) with partner ``q`` (``family``):

- ``readpar.f90`` requires ``q``'s family to be ``p`` (reciprocity) and
  extents within the block; nothing compares the two patches' sizes.
- ``metrix.f90`` ``adiaface`` / ``tangind`` / ``normind`` map the nodes. The
  digits ``l_i l_j l_k`` of ``BC`` name, for this block's i, j, k, the
  partner axis ``(l + 1) // 2`` (1 = I, 2 = J, 3 = K).
  - **Along the face:** an odd digit runs up from ``q``'s lower extent, an
    even one down from its upper extent, by the node's offset from ``p``'s
    lower extent.
  - **Across it:** with ``m = 1`` on a lo face and ``0`` on a hi face, an
    odd digit puts the partner face at ``m * (q's upper extent)`` and an
    even one at ``(1 - m) * (q's upper extent)``. So odd pairs opposite
    sides (lo <-> hi) and even pairs the same sides.

This reading checked all 834 connection patches of a solver-run INSEAN case
(``xnavis_c/test/CC``: ``assetto.grd`` and its ``cc.par``) to 0 errors. The
two hand-built seams of Q0 pin it to that file's codes: 135 on both sides
of a j lo <-> hi seam (89 patch pairs there), and of a k lo <-> hi seam (98).

Interface under test
--------------------
- ``write_cc_par`` (``shore cc-par``, ``write_xall``) writes connection
  codes that Xall reads as the face graph means them.
- A face shared with several partners (CH's ``body.k_hi`` with its three
  channel blocks) gets a patch per partner, over that partner's sub-range.
  Its adjacency sidecar entry carries ``patches``: per partner a 0-based
  node ``range``, the ``partner`` and the ``axis_map``. Other faces and
  sidecars are unchanged.
- Nothing else in ``cc.par`` changes: header, blocks, patch order, extents
  and families are those of the writer before ADR 0036, except the
  connection codes of lo faces (``data/xall_before.zip``, made before it).

Cases
-----
- Q0 (x 2) the interpreter on hand-built seams with the real file's codes.
- Q1 (x 11) every connection patch, as Xall reads it, lands on coincident
  nodes (1e-12 x the bounding-box diagonal), reciprocal, of equal length,
  inside the block.
- Q2 (x 11) the patches of each block tile its six faces: every boundary
  cell face once.
- Q3 (x 9) against the writer before ADR 0036: only lo-face connection codes
  differ, and only in their normal digit.
- Q4 (x 2) CH (sharp, blunt): ``write_xall`` succeeds; ``body.k_hi`` has
  three patches, one per channel block, and is the only face with
  ``patches`` in the sidecar.
- Q5 (x 9) the sidecars of the other block sets are byte-identical to
  before.
- Q6 a CH project saves, loads without replay and saves the same files.

Out of scope of this spec: multigrid divisibility of extents (Xall's
``readpar`` stops on it; a mesh-size choice), splitting CH blocks, Xall's
cell-centre mapping (``calcorn.f90``), running Xall.
"""

from __future__ import annotations

import filecmp
import json
import warnings
import zipfile
from pathlib import Path

import numpy as np
import pytest
from typer.testing import CliRunner

from tests.spec.test_exports import CASES, CG, CH, _files, build

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-28; implementation pending")

DATA = Path(__file__).parent / "data"
runner = CliRunner()
SETS = [*CASES, "ch_blunt"]
SYMMETRIC = [n for n in CASES if n != "ch"]
FACES = ("i_lo", "i_hi", "j_lo", "j_hi", "k_lo", "k_hi")


# ---------------------------------------------------------------------------
# cc.par as Xall's overset pre-processor reads it
# ---------------------------------------------------------------------------


def parse(text: str) -> tuple[int, list[list[int]]]:
    """(number of blocks, patch lines ``[block, face, bc, family, I0, I1, J0, J1, K0, K1]``)."""
    lines = [ln.split("!")[0].split() for ln in text.splitlines()]
    lines = [ln for ln in lines if ln]
    i = 0
    while True:  # the blocks section: a count, then that many 'level group priority'
        if len(lines[i]) == 1 and lines[i][0].isdigit():
            n = int(lines[i][0])
            rows = lines[i + 1 : i + 1 + n]
            if len(rows) == n and all(len(r) == 3 for r in rows) and len(lines[i + 1 + n]) == 1:
                break
        i += 1
    nblocks = int(lines[i][0])
    i += 1 + nblocks
    npatch = int(lines[i][0])
    return nblocks, [list(map(int, lines[i + 1 + k][:10])) for k in range(npatch)]


def xall_errors(text: str, nodes: list[np.ndarray], tol: float) -> list[str]:
    """Connection patches checked as ``readpar`` / ``adiaface`` read them."""
    nblocks, patches = parse(text)
    errors = [] if nblocks == len(nodes) else [f"{nblocks} blocks for {len(nodes)}"]
    for p, (blk, nf, bc, fam, *ext) in enumerate(patches, start=1):
        cells = [n - 1 for n in nodes[blk - 1].shape[:3]]
        lo, hi = ext[0::2], ext[1::2]
        if any(a < 0 or b > c or a > b for a, b, c in zip(lo, hi, cells, strict=True)):
            errors.append(f"p{p}: extents {ext} outside block {blk} ({cells} cells)")
            continue
        if bc <= 99:
            continue
        if not 1 <= fam <= len(patches):
            errors.append(f"p{p}: family {fam}")
            continue
        q = patches[fam - 1]
        if q[3] != p:
            errors.append(f"p{p}: partner p{fam} points at p{q[3]}")
        if q[2] <= 99:
            errors.append(f"p{p}: partner p{fam} is not a connection")
        digits = (bc // 100, bc // 10 % 10, bc % 10)
        if sorted((d + 1) // 2 for d in digits) != [1, 2, 3]:
            errors.append(f"p{p}: code {bc} is not a permutation")
            continue
        nla = q[4:10]
        m = nf % 2
        normal = (nf - 1) // 2
        for a in range(3):
            pa = (digits[a] - 1) // 2
            if a != normal and hi[a] - lo[a] != nla[2 * pa + 1] - nla[2 * pa]:
                errors.append(f"p{p}: length {hi[a] - lo[a]} on axis {a}, partner's differs")
        here, there = nodes[blk - 1], nodes[q[0] - 1]
        worst = 0.0
        for i in range(lo[0], hi[0] + 1):
            for j in range(lo[1], hi[1] + 1):
                for k in range(lo[2], hi[2] + 1):
                    idx = [0, 0, 0]
                    for a, r in enumerate((i - lo[0], j - lo[1], k - lo[2])):
                        d = digits[a]
                        pa = (d - 1) // 2
                        top = nla[2 * pa + 1]
                        if a == normal:
                            idx[pa] = m * top if d % 2 else (1 - m) * top
                        else:
                            idx[pa] = nla[2 * pa] + r if d % 2 else top - r
                    if any(v < 0 or v >= n for v, n in zip(idx, there.shape[:3], strict=True)):
                        errors.append(f"p{p}: maps outside block {q[0]}")
                        break
                    worst = max(worst, float(np.abs(here[i, j, k] - there[tuple(idx)]).max()))
        if worst > tol:
            errors.append(f"p{p} (block {blk} face {nf} code {bc}): nodes {worst:.3g} off")
    return errors


def tiling_errors(text: str, nodes: list[np.ndarray]) -> list[str]:
    _, patches = parse(text)
    cover = {}
    for blk, nf, _, _, *ext in patches:
        a = (nf - 1) // 2
        cells = [n - 1 for n in nodes[blk - 1].shape[:3]]
        grid = cover.setdefault((blk, nf), np.zeros([cells[b] for b in range(3) if b != a]))
        lo, hi = ext[0::2], ext[1::2]
        grid[tuple(slice(lo[b], hi[b]) for b in range(3) if b != a)] += 1
    errors = []
    for blk in range(1, len(nodes) + 1):
        for nf in range(1, 7):
            grid = cover.get((blk, nf))
            if grid is None or grid.min() != 1 or grid.max() != 1:
                errors.append(f"block {blk} face {nf} is not tiled once")
    return errors


# ---------------------------------------------------------------------------
# block sets
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def sets(tmp_path_factory) -> dict:
    d = tmp_path_factory.mktemp("ccpar")
    f = _files(d)
    out = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for name in CASES:
            out[name] = build(name, f)
        from shore.project import Project

        p = Project()
        p.load_geometry(f["blunt"], skip_hygiene=True)
        p.build_topology("ch", options=CH, **{**CG, "nk": 8})
        p.march()
        out["ch_blunt"] = (p, p.mesh.blocks)
    return {"dir": d, "files": f, **out}


def write(sets, name: str) -> tuple[str, str, list[np.ndarray]]:
    """(cc.par text, sidecar text, block nodes) of a block set."""
    from shore.io.adjacency import write_adjacency_json
    from shore.io.cc_par import write_cc_par

    blocks = list(sets[name][1])
    d = sets["dir"]
    write_adjacency_json(blocks, d / f"{name}.json")
    payload = json.loads((d / f"{name}.json").read_text())
    write_cc_par(payload, d / f"{name}.par", grd_filename="mesh.grd")
    nodes = [b.nodes for b in blocks]
    return (d / f"{name}.par").read_text(), (d / f"{name}.json").read_text(), nodes


def diagonal(nodes) -> float:
    pts = np.concatenate([n.reshape(-1, 3) for n in nodes])
    return float(np.linalg.norm(pts.max(axis=0) - pts.min(axis=0)))


def before(name: str, ext: str) -> str:
    with zipfile.ZipFile(DATA / "xall_before.zip") as z:
        return z.read(f"{name}.{ext}").decode()


# ---------------------------------------------------------------------------
# Q0: the interpreter on the real file's codes
# ---------------------------------------------------------------------------


def two_blocks(axis: int) -> list[np.ndarray]:
    """Two 3 x 4 x 5 blocks abutting along *axis* (block 2 after block 1)."""
    shape = (3, 4, 5)
    grid = np.stack(np.meshgrid(*(np.arange(n, dtype=float) for n in shape), indexing="ij"), -1)
    other = grid.copy()
    other[..., axis] += shape[axis] - 1
    return [grid, other]


def two_block_cc_par(face_lo: int, face_hi: int, code_lo: int, code_hi: int, axis: int) -> str:
    cells = [2, 3, 4]
    ext = [0, cells[0], 0, cells[1], 0, cells[2]]
    hi_ext = list(ext)
    hi_ext[2 * axis] = hi_ext[2 * axis + 1] = cells[axis]
    lo_ext = list(ext)
    lo_ext[2 * axis] = lo_ext[2 * axis + 1] = 0
    lines = [
        "'x.grd'", "'cc.'", ".true.", ".true.", ".false.", "", "1 1", "", "-0.1", "-1.0", "",
        "2", "1 0 0 ! a", "1 0 0 ! b", "2",
        " ".join(map(str, [1, face_hi, code_hi, 2, *hi_ext])),
        " ".join(map(str, [2, face_lo, code_lo, 1, *lo_ext])),
        "0", "0",
    ]  # fmt: skip
    return "\n".join(lines) + "\n"


class TestInterpreter:
    @pytest.mark.parametrize(("axis", "face_lo", "face_hi"), [(1, 3, 4), (2, 5, 6)])
    def test_q0_real_codes(self, axis, face_lo, face_hi):  # already holds: a guard
        nodes = two_blocks(axis)
        assert xall_errors(two_block_cc_par(face_lo, face_hi, 135, 135, axis), nodes, 0.0) == []
        wrong = 145 if axis == 1 else 136  # the lo side's code before ADR 0036
        assert xall_errors(two_block_cc_par(face_lo, face_hi, wrong, 135, axis), nodes, 0.0)


# ---------------------------------------------------------------------------
# Q1-Q3: every block set
# ---------------------------------------------------------------------------

NO_CONNECTIONS = ["ogrid", "half_collar"]
CONNECTED = [n for n in SETS if n not in NO_CONNECTIONS]


def connection_errors(sets, name: str) -> list[str]:
    text, _, nodes = write(sets, name)
    return xall_errors(text, nodes, 1e-12 * diagonal(nodes))


def tiling(sets, name: str) -> list[str]:
    text, _, nodes = write(sets, name)
    return tiling_errors(text, nodes)


class TestBlockSets:
    @pytest.mark.parametrize("name", NO_CONNECTIONS)
    def test_q1_no_connections(self, sets, name):  # already holds: a guard
        assert connection_errors(sets, name) == []

    @PENDING
    @pytest.mark.parametrize("name", CONNECTED)
    def test_q1_connections(self, sets, name):
        assert connection_errors(sets, name) == []

    @pytest.mark.parametrize("name", SYMMETRIC)
    def test_q2_tiling(self, sets, name):  # already holds: a guard
        assert tiling(sets, name) == []

    @PENDING
    @pytest.mark.parametrize("name", ["ch", "ch_blunt"])
    def test_q2_tiling_ch(self, sets, name):
        assert tiling(sets, name) == []

    @pytest.mark.parametrize("name", SYMMETRIC)
    def test_q3_only_lo_codes_change(self, sets, name):  # holds before and after: a guard
        text, _, _ = write(sets, name)
        old, new = before(name, "par").splitlines(), text.splitlines()
        assert len(old) == len(new)
        _, old_patches = parse(before(name, "par"))
        _, new_patches = parse(text)
        changed = {
            i for i, (a, b) in enumerate(zip(old_patches, new_patches, strict=True)) if a != b
        }
        for i in changed:
            a, b = old_patches[i], new_patches[i]
            assert a[:2] == b[:2] and a[3:] == b[3:], (a, b)  # only the code
            assert a[1] % 2 == 1, (a, b)  # on a lo face
            normal = (a[1] - 1) // 2
            da, db = str(a[2]), str(b[2])
            assert [x for n, x in enumerate(da) if n != normal] == [
                x for n, x in enumerate(db) if n != normal
            ]
            assert (
                int(da[normal]) != int(db[normal])
                and (int(da[normal]) + 1) // 2 == (int(db[normal]) + 1) // 2
            )
        patch_lines = {i for i, ln in enumerate(old) if ln != new[i]}
        assert len(patch_lines) == len(changed)


# ---------------------------------------------------------------------------
# Q4-Q6: CH and the sidecars
# ---------------------------------------------------------------------------


class TestCH:
    @PENDING
    @pytest.mark.parametrize("name", ["ch", "ch_blunt"])
    def test_q4_ch(self, sets, name, tmp_path):
        project, blocks = sets[name]
        _, sidecar, _ = write(sets, name)
        project.write_xall(tmp_path / "run")
        payload = json.loads(sidecar)
        shared = [
            (b["label"], f)
            for b in payload["blocks"]
            for f, face in b["faces"].items()
            if "patches" in face
        ]
        assert shared == [("body", "k_hi")]
        body = next(b for b in payload["blocks"] if b["label"] == "body")
        partners = [p["partner"] for p in body["faces"]["k_hi"]["patches"]]
        channels = {"chan_body_top", "chan_body_front", "chan_body_bot"}
        assert {p["label"] for p in partners} == channels
        assert all(p["face"] == "k_lo" for p in partners)
        _, patches = parse((tmp_path / "run" / "cc.par").read_text())
        index = [b.label for b in blocks].index("body") + 1
        assert len([p for p in patches if p[0] == index and p[1] == 6]) == 3

    @pytest.mark.parametrize("name", SYMMETRIC)
    def test_q5_sidecars_unchanged(self, sets, name):  # holds before and after: a guard
        _, sidecar, _ = write(sets, name)
        assert sidecar == before(name, "json")

    @PENDING
    def test_q6_ch_project_round_trip(self, sets, tmp_path):
        from shore.project import Project

        project = sets["ch"][0]
        project.save(tmp_path / "a")
        assert '"patches"' in (tmp_path / "a" / "grids" / "main" / "adjacency.json").read_text()
        fast = Project.load(tmp_path / "a", replay=False)
        fast.save(tmp_path / "b")
        pending = [filecmp.dircmp(tmp_path / "a", tmp_path / "b")]
        while pending:
            cmp = pending.pop()
            assert not cmp.left_only and not cmp.right_only, (cmp.left, cmp.left_only)
            _, mismatch, errors = filecmp.cmpfiles(
                cmp.left, cmp.right, cmp.common_files, shallow=False
            )
            assert not mismatch and not errors, (cmp.left, mismatch)
            pending.extend(cmp.subdirs.values())
        fast.write_xall(tmp_path / "run")
        assert (tmp_path / "run" / "cc.par").is_file()
