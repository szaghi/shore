"""FROZEN acceptance spec: CGNS and Plot3D exports (ADR 0035, roadmap phase 5).

Frozen 2026-10-28, before any implementation. Implementing work may not edit
this file, ``tests/spec/__init__.py``, ``tests/conftest.py``, the pytest
configuration in ``pyproject.toml`` or ``src/shore/quality.py``. The ONE
permitted edit is deleting an ``@PENDING`` line once that case passes (strict
xfail turns a passing pending case into a failure). ``tests/spec/SHA256SUMS``
holds the hashes, taken with ``@PENDING`` lines removed;
``scripts/check_spec.py`` checks them.

Verdict for acceptance (main CI job): ``pytest tests/spec/test_exports.py``
reports 61 passed, 0 skipped, 0 xfailed. At freeze: 0 passed, 61 xfailed.

Interface under test
--------------------
- ``shore.io.cgns.write_cgns(components, path)``: ``components`` maps a base
  name to its blocks. One CGNS/HDF5 file: a ``CGNSBase_t`` per component
  (cell and physical dimension 3), a structured ``Zone_t`` per block named
  by its label, ``GridCoordinates`` in double precision.
  - Connected faces (SHARED, DIRICHLET, PERIODIC) become
    ``GridConnectivity1to1_t`` in both zones, with ``PointRange``,
    ``PointRangeDonor`` and ``Transform`` as the SIDS define them. A face
    shared with several smaller faces (CH's ``body.k_hi``) gets one
    connection per partner, over the partner's sub-range.
  - WALL, SYMMETRY and FREE faces become ``BC_t`` over the whole face:
    ``BCWall``, ``BCSymmetryPlane`` and ``BCTypeUserDefined`` with a
    ``Descriptor_t`` reading ``shore:FREE``.
  - A connection whose nodes do not coincide, or a zone name over 32
    characters, raises ``ParameterError``.
- ``shore.io.plot3d.write_plot3d(blocks, path)``: one multi-block, whole,
  3-D, double-precision Plot3D grid file.
- ``Project.write_cgns(path, components=())`` and
  ``Project.write_plot3d(path, components=())``: output commands
  (journaled, never undone) over the named components, all those with a mesh
  when empty, in project order. An unknown or empty component raises
  ``ParameterError``.
- ``shore cgns DIR -o FILE [-c C ...] [--no-replay]`` and ``shore plot3d``
  likewise write a saved project.

Independent readers: VTK's ``vtkCGNSReader`` and
``vtkMultiBlockPLOT3DReader`` (bases, zones, coordinates), and pyCGNS's
``CGNS.MAP`` loader (the CGNS/Python tree that the SIDS checks below read),
with its ``CGNS.VAL`` DEFAULT grammar. pyCGNS 6.3.5's SIDS grammar crashes
under numpy 2 even on trees built by pyCGNS itself, so the SIDS rules for
1-to-1 connections and boundary coverage are checked here instead.

Cases
-----
Ten block sets: cubed sphere, O-grid, sharp and blunt C-grids, CH, OH, a
cubed sphere split along k, flat caps, a split collar at a lens's two sharp
corners, a half collar on a symmetry plane.

- X1 (x 10) VTK reads the CGNS file: one base, the zones, coordinates bit
  for bit.
- X2 (x 10) SIDS: zone sizes; every 1-to-1 maps each node onto a donor node
  within 1e-12 x the bounding-box diagonal, with the donor range's end the
  image of its own end and a reciprocal connection in the donor; every cell
  face on a zone's boundary described exactly once.
- X3 (x 10) the file says what the face graph says: BC types per face; the
  donors of each connected face are the blocks it is wired to.
- X4 (x 10) pyCGNS's DEFAULT grammar reports no error.
- X5 (x 10) VTK reads the Plot3D file: blocks, dimensions, coordinates bit
  for bit.
- X6 (x 3) projects: two components give two bases in project order, and
  ``components=`` selects; journal lines; unknown and empty components
  raise.
- X7 (x 4) the CLI writes the same bytes as the API, also with
  ``--no-replay``; ``-c`` selects.
- X8 (x 2) refusals: nodes that do not coincide, a long zone name.
- X9 (x 2) writing twice gives the same bytes (CGNS, Plot3D).

Out of scope of this spec: the O-grid's periodic j seam (not a face of the
block graph, as in ``cc.par``), overset (chimera) hole cutting and
connectivity in CGNS, flow solutions, CGNS/ADF files, ASCII or
single-precision Plot3D, IBLANK, reading CGNS back into SHORE.
"""

from __future__ import annotations

import filecmp
import importlib.util  # noqa: F401  (pyCGNS 6.3.5 uses importlib.util without importing it)
import warnings
from pathlib import Path

import CGNS.MAP
import CGNS.VAL.simplecheck
import numpy as np
import pytest
import trimesh
import vtk
from typer.testing import CliRunner
from vtk.util.numpy_support import vtk_to_numpy

from shore.exceptions import ParameterError
from tests.conftest import _make_naca0012_stl

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-28; implementation pending")

runner = CliRunner()
OG = {"nk": 5, "ds": 0.05, "growth": 1.1, "ni": 9, "nj": 16}
CG = {"nk": 6, "ds": 1e-3, "ni_body": 40, "ni_wake": 10, "n_stations": 4, "wake_length": 5.0}
CH = {"chord_min": -2.0, "chord_max": 15.0, "normal_min": -2.5, "normal_max": 2.5, "nk_channel": 4}
OH = {"nk_channel": 4, "box_vmin": (-4.0, -4.0, -4.0), "box_vmax": (4.0, 4.0, 4.0)}
COLLAR = {
    "seed": (0.45, 0.2, 0.02), "n_loop": 16, "ds_surface": 0.02, "growth_surface": 1.15,
    "m_a": 5, "m_b": 5, "ds_wall": 1e-3, "growth_wall": 1.25, "nk": 8,
    "plane": ((0.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
}  # fmt: skip
CASES = [
    "cubed_sphere", "ogrid", "cgrid", "cgrid_blunt", "ch", "oh", "cubed_split", "flat_caps",
    "split_collar", "half_collar",
]  # fmt: skip
SIX = ("sub0", "sub1", "sub2", "sub3", "cap_north", "cap_south")


# ---------------------------------------------------------------------------
# block sets
# ---------------------------------------------------------------------------


def _files(d: Path) -> dict[str, Path]:
    trimesh.creation.icosphere(subdivisions=3).export(str(d / "sphere.stl"))
    _make_naca0012_stl(d / "wing.stl")
    _make_naca0012_stl(d / "blunt.stl", blunt_te_thickness=0.01)
    hbox = trimesh.creation.box(extents=(4.0, 2.0, 1.0))
    hbox.apply_translation((0.0, 1.0, -0.5))
    hbox.update_faces(~np.all(np.isclose(hbox.vertices[hbox.faces][..., 1], 0.0), axis=1))
    post = trimesh.creation.cylinder(radius=0.5, height=1.6, sections=48)
    post.apply_translation((0.0, 0.0, 0.2))
    post.update_faces(post.triangles_center[:, 1] > 0.0)
    for mesh in (hbox, post):
        mesh.remove_unreferenced_vertices()
    hbox.export(str(d / "hbox.stl"))
    post.export(str(d / "hpost.stl"))
    return {k: d / f"{k}.stl" for k in ("sphere", "wing", "blunt", "hbox", "hpost")}


def build(name: str, f: dict[str, Path]):
    """The project of case *name* (its current component holds the blocks)."""
    from shore.project import Project

    p = Project()
    if name == "split_collar":
        from tests.spec.test_split_collar import collar

        blocks = collar("lens")[3]
        return None, blocks
    if name == "flat_caps":
        p.build_primitive("flat_caps", r_out=3.0, z0=-1.0, z1=4.0, ni=5, nc=6, nk=5, inner_frac=0.4)
        return p, p.mesh.blocks
    if name == "half_collar":
        p.load_geometry(f["hbox"], name="hbox")
        p.load_geometry(f["hpost"], name="hpost")
        p.build_collar("hpost", geometry="hbox", **COLLAR)
        return p, p.mesh.blocks
    stl = {"cgrid": "wing", "ch": "wing", "cgrid_blunt": "blunt"}.get(name, "sphere")
    p.load_geometry(f[stl], skip_hygiene=stl != "sphere")
    if name in ("cubed_sphere", "cubed_split"):
        p.build_topology("cubed_sphere", **OG)
    elif name == "ogrid":
        p.build_topology("ogrid", **OG)
    elif name == "cgrid":
        p.build_topology("cgrid", **CG)
    elif name == "cgrid_blunt":
        p.build_topology("cgrid", **CG, blocks=4)
    elif name == "ch":
        p.build_topology("ch", options=CH, **{**CG, "nk": 8})
    elif name == "oh":
        p.build_topology("oh", options=OH, nk=5, ds=0.05, growth=1.1, ni=16, nj=16)
    p.march()
    if name == "cubed_split":
        p.split([(label, "k", (2,)) for label in SIX])
    return p, p.mesh.blocks


@pytest.fixture(scope="module")
def written(tmp_path_factory) -> dict[str, tuple[list, Path, Path]]:
    """Per case: (blocks, CGNS file, Plot3D file); files None until implemented."""
    d = tmp_path_factory.mktemp("exports")
    f = _files(d)
    out = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for name in CASES:
            blocks = list(build(name, f)[1])
            out[name] = (blocks, d / f"{name}.cgns", d / f"{name}.xyz")
    return out


def files(written, name: str) -> tuple[list, Path, Path]:
    from shore.io.cgns import write_cgns
    from shore.io.plot3d import write_plot3d

    blocks, cgns, xyz = written[name]
    if not cgns.exists():
        write_cgns({name: blocks}, cgns)
        write_plot3d(blocks, xyz)
    return blocks, cgns, xyz


# ---------------------------------------------------------------------------
# independent readers
# ---------------------------------------------------------------------------


def vtk_cgns(path: Path) -> list[tuple[str, list[np.ndarray]]]:
    """Per base: (name, [zone points (n, 3), i fastest])."""
    r = vtk.vtkCGNSReader()
    r.SetFileName(str(path))
    r.UpdateInformation()  # fills the base list that EnableAllBases selects
    r.EnableAllBases()
    r.Update()
    out = r.GetOutput()
    bases = []
    for b in range(out.GetNumberOfBlocks()):
        zones = out.GetBlock(b)
        pts = [
            vtk_to_numpy(zones.GetBlock(z).GetPoints().GetData())
            for z in range(zones.GetNumberOfBlocks())
        ]
        bases.append((out.GetMetaData(b).Get(vtk.vtkCompositeDataSet.NAME()), pts))
    return bases


def kids(node, label=None):
    return [c for c in node[2] if label is None or c[3] == label]


def text(node) -> str:
    return node[1].tobytes().decode()


def tree(path: Path):
    return CGNS.MAP.load(str(path))[0]


def sids(path: Path) -> dict[str, dict[str, dict]]:
    """``{base: {zone: {size, nodes, conns, bcs, type}}}`` from pyCGNS's tree."""
    bases = {}
    for base in kids(tree(path), "CGNSBase_t"):
        zones = {}
        for z in kids(base, "Zone_t"):
            (gc,) = kids(z, "GridCoordinates_t")
            xyz = {c[0]: c[1] for c in kids(gc, "DataArray_t")}
            nodes = np.stack([xyz[f"Coordinate{a}"] for a in "XYZ"], axis=-1)
            conns = []
            for zgc in kids(z, "ZoneGridConnectivity_t"):
                for c in kids(zgc, "GridConnectivity1to1_t"):
                    get = {k[0]: k[1] for k in c[2]}
                    conns.append(
                        (c[0], text(c), get["PointRange"].T, get["PointRangeDonor"].T,
                         tuple(int(v) for v in get["Transform"]))
                    )  # fmt: skip
            bcs = []
            for zbc in kids(z, "ZoneBC_t"):
                for bc in kids(zbc, "BC_t"):
                    rng = {k[0]: k[1] for k in bc[2]}["PointRange"].T
                    bcs.append((bc[0], text(bc), rng, [text(d) for d in kids(bc, "Descriptor_t")]))
            ztype = text(kids(z, "ZoneType_t")[0])
            zones[z[0]] = {"size": z[1], "nodes": nodes, "conns": conns, "bcs": bcs, "type": ztype}
        bases[base[0]] = zones
    return bases


def image(t, begin, begin_donor, index) -> np.ndarray:
    """SIDS: ``Index2 = T.(Index1 - Begin1) + Begin2``."""
    out = np.array(begin_donor, dtype=np.int64).copy()
    for a, v in enumerate(t):
        out[abs(v) - 1] += int(np.sign(v)) * (int(index[a]) - int(begin[a]))
    return out


def face_of(rng: np.ndarray, shape) -> str | None:
    lo, hi = np.minimum(rng[0], rng[1]), np.maximum(rng[0], rng[1])
    axes = [a for a in range(3) if lo[a] == hi[a]]
    if len(axes) != 1 or lo[axes[0]] not in (1, shape[axes[0]]):
        return None
    a = axes[0]
    return "ijk"[a] + ("_lo" if lo[a] == 1 else "_hi")


def sids_errors(path: Path, tol: float) -> list[str]:
    errors = []
    for zones in sids(path).values():
        for zname, z in zones.items():
            shape = z["nodes"].shape[:3]
            if z["type"] != "Structured":
                errors.append(f"{zname}: ZoneType {z['type']}")
            if tuple(z["size"][:, 0]) != shape or tuple(z["size"][:, 1]) != tuple(
                n - 1 for n in shape
            ):
                errors.append(f"{zname}: sizes {z['size'].tolist()} for nodes {shape}")
            cover: dict[str, np.ndarray] = {}

            def mark(rng, what, shape=shape, cover=cover, zname=zname):
                face = face_of(rng, shape)
                if face is None:
                    errors.append(f"{zname}/{what}: range {rng.tolist()} is not on a face")
                    return
                a = "ijk".index(face[0])
                lo, hi = np.minimum(rng[0], rng[1]), np.maximum(rng[0], rng[1])
                cells = [shape[b] - 1 for b in range(3) if b != a]  # cell faces, not nodes
                grid = cover.setdefault(face, np.zeros(cells))
                grid[tuple(slice(lo[b] - 1, hi[b] - 1) for b in range(3) if b != a)] += 1

            for name, donor, rng, rngd, t in z["conns"]:
                mark(rng, name)
                if donor not in zones:
                    errors.append(f"{zname}/{name}: donor {donor} not in the base")
                    continue
                if not np.array_equal(image(t, rng[0], rngd[0], rng[1]), rngd[1]):
                    errors.append(f"{zname}/{name}: donor range end is not the image of its end")
                dz = zones[donor]
                lo, hi = np.minimum(rng[0], rng[1]), np.maximum(rng[0], rng[1])
                axes = [np.arange(lo[a], hi[a] + 1) for a in range(3)]
                idx = np.stack(np.meshgrid(*axes, indexing="ij"), -1).reshape(-1, 3)
                img = np.array([image(t, rng[0], rngd[0], i) for i in idx])
                if (img < 1).any() or (img > np.array(dz["nodes"].shape[:3])).any():
                    errors.append(f"{zname}/{name}: maps outside {donor}")
                    continue
                err = np.abs(z["nodes"][tuple((idx - 1).T)] - dz["nodes"][tuple((img - 1).T)])
                if float(err.max()) > tol:
                    errors.append(f"{zname}/{name}: nodes {float(err.max()):.3g} off {donor}")
                back = [
                    c for c in dz["conns"]
                    if c[1] == zname and np.array_equal(np.sort(c[3], axis=0), np.sort(rng, axis=0))
                ]  # fmt: skip
                if len(back) != 1:
                    errors.append(f"{zname}/{name}: {len(back)} reciprocal connections")
            for name, _, rng, _ in z["bcs"]:
                mark(rng, name)
            for face in ("i_lo", "i_hi", "j_lo", "j_hi", "k_lo", "k_hi"):
                grid = cover.get(face)
                if grid is None or grid.min() < 1:
                    errors.append(f"{zname}: face {face} not fully described")
                elif grid.max() > 1:
                    errors.append(f"{zname}: face {face} partly described twice")
    return errors


class _Array(np.ndarray):
    """pyCGNS 6.3.5 calls ``tostring`` and array truthiness, gone in numpy 2."""

    def tostring(self, order="C"):
        return self.tobytes(order)

    def __bool__(self):
        return self.size > 0


def default_grammar_errors(path: Path) -> list[str]:
    t = tree(path)
    pending = [t]
    while pending:
        node = pending.pop()
        if isinstance(node[1], np.ndarray):
            node[1] = node[1].view(_Array)
        pending.extend(node[2])
    log = CGNS.VAL.simplecheck.run(t, False, ["DEFAULT"])
    return [f"{p}: {m}" for p in log.keys() for m in log[p] if '"E"' in str(m)]


def vtk_plot3d(path: Path) -> list[tuple[list[int], np.ndarray]]:
    r = vtk.vtkMultiBlockPLOT3DReader()
    r.SetXYZFileName(str(path))
    r.AutoDetectFormatOn()
    r.Update()
    out = r.GetOutput()
    blocks = []
    for b in range(out.GetNumberOfBlocks()):
        grid = out.GetBlock(b)
        dims = [0, 0, 0]
        grid.GetDimensions(dims)
        blocks.append((dims, vtk_to_numpy(grid.GetPoints().GetData())))
    return blocks


def diagonal(blocks) -> float:
    pts = np.concatenate([b.nodes.reshape(-1, 3) for b in blocks])
    return float(np.linalg.norm(pts.max(axis=0) - pts.min(axis=0)))


# ---------------------------------------------------------------------------
# X1-X5: per block set
# ---------------------------------------------------------------------------

FREE_DESCRIPTOR = "shore:FREE"
BC_TYPE = {"WALL": "BCWall", "SYMMETRY": "BCSymmetryPlane", "FREE": "BCTypeUserDefined"}


class TestBlockSets:
    @pytest.mark.parametrize("name", CASES)
    def test_x1_vtk_reads_cgns(self, written, name):
        blocks, cgns, _ = files(written, name)
        bases = vtk_cgns(cgns)
        assert [b[0] for b in bases] == [name]
        points = bases[0][1]
        assert len(points) == len(blocks)
        for pts, block in zip(points, blocks, strict=True):
            assert np.array_equal(pts, block.nodes.reshape(-1, 3, order="F"))

    @pytest.mark.parametrize("name", CASES)
    def test_x2_sids(self, written, name):
        blocks, cgns, _ = files(written, name)
        assert list(sids(cgns)[name]) == [b.label for b in blocks]
        assert sids_errors(cgns, 1e-12 * diagonal(blocks)) == []

    @pytest.mark.parametrize("name", CASES)
    def test_x3_face_graph(self, written, name):
        blocks, cgns, _ = files(written, name)
        zones = sids(cgns)[name]
        owner = {id(f): b.label for b in blocks for f in b.faces()}
        wired: dict[tuple[str, str], set[str]] = {}
        for b in blocks:
            for f in b.faces():
                if f.partner is not None and f.bc.name in ("SHARED", "DIRICHLET", "PERIODIC"):
                    wired.setdefault((b.label, f.name), set()).add(owner[id(f.partner)])
                    wired.setdefault((owner[id(f.partner)], f.partner.name), set()).add(b.label)
        for b in blocks:
            z = zones[b.label]
            shape = b.nodes.shape[:3]
            bcs = {face_of(rng, shape): (bctype, desc) for _, bctype, rng, desc in z["bcs"]}
            donors: dict[str, set[str]] = {}
            for _, donor, rng, _, _ in z["conns"]:
                donors.setdefault(face_of(rng, shape), set()).add(donor)
            for f in b.faces():
                if f.bc.name in BC_TYPE:
                    bctype, desc = bcs[f.name]
                    assert bctype == BC_TYPE[f.bc.name], (b.label, f.name)
                    assert (FREE_DESCRIPTOR in desc) == (f.bc.name == "FREE")
                else:
                    assert donors.get(f.name) == wired[(b.label, f.name)], (b.label, f.name)

    @pytest.mark.parametrize("name", CASES)
    def test_x4_pycgns_default_grammar(self, written, name):
        _, cgns, _ = files(written, name)
        assert default_grammar_errors(cgns) == []

    @pytest.mark.parametrize("name", CASES)
    def test_x5_vtk_reads_plot3d(self, written, name):
        blocks, _, xyz = files(written, name)
        read = vtk_plot3d(xyz)
        assert len(read) == len(blocks)
        for (dims, pts), block in zip(read, blocks, strict=True):
            assert dims == list(block.nodes.shape[:3])
            assert np.array_equal(pts, block.nodes.reshape(-1, 3, order="F"))


# ---------------------------------------------------------------------------
# X6-X7: projects and the CLI
# ---------------------------------------------------------------------------


@pytest.fixture
def saved(tmp_path) -> Path:
    from tests.spec.test_fast_load import make_project

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return make_project(tmp_path)


def load(directory: Path, replay: bool = True):
    from shore.project import Project

    return Project.load(directory) if replay else Project.load(directory, replay=False)


class TestProject:
    def test_x6_bases_per_component(self, saved, tmp_path):
        p = load(saved)
        p.write_cgns(tmp_path / "all.cgns")
        p.write_cgns(tmp_path / "one.cgns", components=["collar"])
        names = [n for n, m in p.components.items() if m is not None]
        assert [b[0] for b in vtk_cgns(tmp_path / "all.cgns")] == names == ["main", "collar"]
        assert [b[0] for b in vtk_cgns(tmp_path / "one.cgns")] == ["collar"]
        for name in names:
            labels = [b.label for b in p.components[name].blocks]
            assert list(sids(tmp_path / "all.cgns")[name]) == labels
        p.write_plot3d(tmp_path / "all.xyz")
        blocks = [b for n in names for b in p.components[n].blocks]
        read = vtk_plot3d(tmp_path / "all.xyz")
        assert len(read) == len(blocks)
        for (_, pts), block in zip(read, blocks, strict=True):
            assert np.array_equal(pts, block.nodes.reshape(-1, 3, order="F"))

    def test_x6_journal(self, saved, tmp_path):
        p = load(saved)
        before = p.journal()
        p.write_cgns(tmp_path / "a.cgns", components=["collar"])
        p.write_plot3d(tmp_path / "a.xyz")
        lines = p.journal()[len(before) :].strip().splitlines()
        assert lines[0].startswith("p.write_cgns(") and "components=('collar',)" in lines[0]
        assert lines[1].startswith("p.write_plot3d(")

    def test_x6_errors(self, saved, tmp_path):
        from shore.project import Project

        p = load(saved)
        with pytest.raises(ParameterError):
            p.write_cgns(tmp_path / "x.cgns", components=["nope"])
        with pytest.raises(ParameterError):
            p.write_plot3d(tmp_path / "x.xyz", components=["nope"])
        with pytest.raises(ParameterError):
            Project().write_cgns(tmp_path / "y.cgns")


class TestCli:
    @pytest.mark.parametrize("fmt", ["cgns", "plot3d"])
    def test_x7_cli_equals_api(self, saved, tmp_path, fmt):
        from shore.cli import app

        write = getattr(load(saved), f"write_{fmt}")
        write(tmp_path / "api")
        r = runner.invoke(app, [fmt, str(saved), "-o", str(tmp_path / "cli")])
        assert r.exit_code == 0, r.output
        r = runner.invoke(app, [fmt, str(saved), "-o", str(tmp_path / "fast"), "--no-replay"])
        assert r.exit_code == 0, r.output
        assert filecmp.cmp(tmp_path / "api", tmp_path / "cli", shallow=False)
        assert filecmp.cmp(tmp_path / "api", tmp_path / "fast", shallow=False)

    @pytest.mark.parametrize("fmt", ["cgns", "plot3d"])
    def test_x7_cli_component(self, saved, tmp_path, fmt):
        from shore.cli import app

        getattr(load(saved), f"write_{fmt}")(tmp_path / "api", components=["collar"])
        r = runner.invoke(app, [fmt, str(saved), "-o", str(tmp_path / "cli"), "-c", "collar"])
        assert r.exit_code == 0, r.output
        assert filecmp.cmp(tmp_path / "api", tmp_path / "cli", shallow=False)


# ---------------------------------------------------------------------------
# X8-X9: refusals and determinism
# ---------------------------------------------------------------------------


class TestRefusals:
    def test_x8_unmatched_nodes(self, written, tmp_path):
        from shore.io.cgns import write_cgns

        blocks = list(written["cgrid"][0])
        nodes = blocks[0].nodes.copy()
        nodes[-1, 1, 1] += 1e-3  # upper's i_hi face, SHARED with body's i_lo
        moved = type(blocks[0]).__new__(type(blocks[0]))
        moved.__dict__.update(blocks[0].__dict__)
        moved.nodes = nodes
        with pytest.raises(ParameterError):
            write_cgns({"c": [moved, *blocks[1:]]}, tmp_path / "x.cgns")

    def test_x8_long_zone_name(self, written, tmp_path):
        from shore.io.cgns import write_cgns

        (block,) = written["ogrid"][0]
        long = type(block).__new__(type(block))
        long.__dict__.update(block.__dict__)
        long.label = "z" * 33
        with pytest.raises(ParameterError):
            write_cgns({"o": [long]}, tmp_path / "x.cgns")


class TestDeterminism:
    def test_x9_cgns_twice(self, written, tmp_path):
        from shore.io.cgns import write_cgns

        blocks = written["ch"][0]
        write_cgns({"ch": blocks}, tmp_path / "a.cgns")
        write_cgns({"ch": blocks}, tmp_path / "b.cgns")
        assert filecmp.cmp(tmp_path / "a.cgns", tmp_path / "b.cgns", shallow=False)

    def test_x9_plot3d_twice(self, written, tmp_path):
        from shore.io.plot3d import write_plot3d

        blocks = written["ch"][0]
        write_plot3d(blocks, tmp_path / "a.xyz")
        write_plot3d(blocks, tmp_path / "b.xyz")
        assert filecmp.cmp(tmp_path / "a.xyz", tmp_path / "b.xyz", shallow=False)
