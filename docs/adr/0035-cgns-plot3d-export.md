# 0035 — CGNS and Plot3D exports

- **Status:** Accepted
- **Date:** 2026-10-28

## Context

Roadmap phase 5 ([ADR 0001](./0001-scope.md)) delivers the standard
structured-grid formats besides Xall's:
- CGNS, with structured zones and 1-to-1 connectivity from the face graph;
- Plot3D.

Its acceptance is a round trip through an independent reader, such as
ParaView's CGNS reader.

Phase 4 was closed on 2026-10-28 with its open limits recorded on the
roadmap. The user was away, and the night protocol was in force: go ahead on
my recommendations, which were to close phase 4 and start phase 5. The
choices below are therefore mine, and the **review questions** at the end
are for the user.

## Investigation (2026-10-28, prototype in the session scratchpad `cg/`; T1)

**Writing CGNS needs HDF5.** CGNS/HDF5 is the current CGNS file format (ADF
is legacy).
- `h5py` writes it directly, following the CGNS HDF5 file mapping:
  - a group per node, with `name`, `label`, `type` and `flags`
    attributes;
  - data in a `" data"` dataset, with the SIDS dimensions reversed;
  - children in creation order.
- The CGNS C library has no maintained pip wheel. VTK has a CGNS reader but
  no writer.

**Independent readers available from pip:**
- **VTK 9.6** (already a dev dependency through pyvista):
  - `vtkCGNSReader` reads bases, zones and coordinates;
  - it skips structured `BC_t` patches with "not supported yet";
  - `vtkMultiBlockPLOT3DReader` reads Plot3D.
- **pyCGNS 6.3.5** (ONERA, LGPL-2, wheels for CPython 3.11-3.13): its
  `CGNS.MAP` loads the file into a CGNS/Python tree.
  - Its `CGNS.VAL` DEFAULT grammar validates names, labels, data types and
    the library version; it accepts versions ≤ 3.4 only.
  - Its SIDS grammar crashes under numpy 2 even on trees built by pyCGNS
    itself: `IndexRange_t` and `GridConnectivity1to1_t` raise "too many
    values to unpack". The spec checks those SIDS rules itself, on the tree
    pyCGNS loads.
  - Two numpy-2 shims are needed: `ndarray.tostring`, and array
    truthiness.

**A sweep over every builder** (prototype; T1). Sixteen block sets:
- the five topologies, plus blunt C-grid and blunt CH;
- a cubed sphere split along k;
- the five primitives;
- a collar, a half collar and a split collar.

Results:
- VTK reads every zone's coordinates bit for bit.
- The SIDS checks pass with a tolerance of 1e-12 × the bounding-box
  diagonal: every 1-to-1 node maps onto its donor's, the connections are
  reciprocal, and every boundary cell face is described exactly once.
  - The C-grid wake cuts and the flat-caps centre seams match to 2e-15,
    not bit for bit.
- **CH shares a face with several partners.** `body.k_hi` (40 nodes along
  i) is shared with three channel blocks of 14, 14 and 14 nodes. The face
  graph keeps one partner per face, so `body.k_hi` names only
  `chan_body_front`; CH's builder calls this its "Phase-1 limitation".
  - Whole-face connections would be wrong there. The writer instead takes
    each connected pair of faces once and locates the smaller face inside
    the larger one by its nodes. `body.k_hi` thus gets three sub-range
    connections, one per channel block.
  - Whether `cc.par` describes these three seams correctly is a separate
    question, and not checked here.
- **The O-grid's j seam** is periodic in the marcher but is not a face of
  the block graph (its j faces are FREE, with a gap of one cell). As in
  `cc.par`, it is not written as a connection.
- **Files are deterministic** (same bytes when written twice). The CLI, the
  Python API and a load without replay give identical files.

## Decision

- **`shore.io.cgns.write_cgns(components, path)`**: one CGNS/HDF5 file.
  - **Bases:** one `CGNSBase_t` per component (cell and physical dimension
    3), named after the component.
  - **Zones:** one structured `Zone_t` per block, named by its label, with
    `GridCoordinates` in double precision.
  - **Connections:** SHARED, DIRICHLET and PERIODIC faces become
    `GridConnectivity1to1_t` in both zones: `PointRange`,
    `PointRangeDonor` (begin = the image of `PointRange`'s begin) and
    `Transform`.
    - Along the face, `Transform` comes from the face's `AxisMap`; across
      it, a step into the zone is a step out of the donor.
    - The smaller face of each pair is located inside the larger one, so
      faces with several partners work.
    - Nodes that do not coincide within 1e-9 × the diagonal are refused
      (`ParameterError`).
  - **Boundary conditions:** WALL, SYMMETRY and FREE faces become `BC_t`
    over the whole face:
    - WALL is `BCWall`, SYMMETRY is `BCSymmetryPlane`;
    - FREE is `BCTypeUserDefined` with a `Descriptor_t` reading
      `shore:FREE`. A FREE face is an overset fringe or a natural boundary,
      which only the solver setup decides.
  - **Library version** 3.4: nothing beyond the 3.4 SIDS is written.
  - **Names** longer than 32 characters, or repeated, are refused.
- **`shore.io.plot3d.write_plot3d(blocks, path)`**: one multi-block, whole
  (no IBLANK), 3-D Plot3D grid file. It is Fortran unformatted, double
  precision, little endian, with node counts.
- **Projects.** `Project.write_cgns(path, components=())` and
  `Project.write_plot3d(path, components=())` are output commands
  (`WriteCgns`, `WritePlot3d`): journaled, never undone, dropped by `save`.
  - With no components named, they write every component with a mesh, in
    project order.
  - For Plot3D, the components' blocks are concatenated in that order.
- **CLI.** `shore cgns DIR -o FILE [-c C ...] [--no-replay]` and
  `shore plot3d ...` write a saved project.
- **Dependencies.**
  - `h5py` is a new optional extra, `shore-mesh[cgns]`. Without it,
    `write_cgns` raises a `ParameterError` naming the extra.
  - `h5py` and `pycgns` join the `dev` extra.

## Alternatives considered

| Alternative | Why not |
|---|---|
| The CGNS C library through pyCGNS's writer | pyCGNS's write path and SIDS grammar are unmaintained under numpy 2; `h5py` is a common, maintained wheel. |
| CGNS/ADF | Legacy; HDF5 is the default since CGNS 3.0. |
| One base for all components, zone names prefixed by component | Zone names would exceed 32 characters (collar pieces), and chimera components have no 1-to-1 connections between them anyway. |
| FREE as `BCFarfield` | Wrong for overset fringes, the common case in Xall runs. |
| Close the O-grid's j seam by duplicating it | Changes the zone's node count against the `.geo` / `.grd` files of the same block. |
| Whole-face 1-to-1 connections only | Wrong for CH's `body.k_hi`. |

## Consequences

- **Phase 5's acceptance** is met by VTK's readers (coordinates, bases,
  zones) and by pyCGNS's loader plus the spec's SIDS checks (connections and
  BCs). VTK skips our BC patches, so the BC check rests on pyCGNS.
- **Not exported:** overset hole cutting and donor search (that stays
  Xall's job), flow solutions, CGNS/ADF, ASCII or single-precision Plot3D,
  IBLANK.
- **CI.** The `dev` extra gains `h5py` and `pycgns`. CI minutes are
  exhausted for October, so the export spec has run only locally.

## Verification

### Acceptance spec (frozen 2026-10-28, before implementation)

`tests/spec/test_exports.py` runs in the main CI job.

**Verdicts.**
- Required: 61 passed, 0 skipped, 0 xfailed.
- At freeze: 61 xfailed.
- With the prototype behind the interface: 61 passed. The full suite with
  the prototype (before the spec existed): 1701 passed.

| Case | Required |
|---|---|
| X1 (× 10) | VTK reads the CGNS file: the base, the zones, coordinates bit for bit |
| X2 (× 10) | SIDS: zone sizes; 1-to-1 nodes onto donor nodes (1e-12 × diagonal), donor range end, reciprocity; every boundary cell face described once |
| X3 (× 10) | BC types per face and the donors of each connected face match the face graph |
| X4 (× 10) | pyCGNS's DEFAULT grammar reports no error |
| X5 (× 10) | VTK reads the Plot3D file: blocks, dimensions, coordinates bit for bit |
| X6 (× 3) | a base per component in project order, `components=` selects; journal lines; errors |
| X7 (× 4) | the CLI writes the API's bytes, also with `--no-replay`; `-c` selects |
| X8 (× 2) | unmatched nodes and a 33-character zone name are refused |
| X9 (× 2) | the same bytes when written twice |

The ten block sets: cubed sphere, O-grid, sharp and blunt C-grids, CH, OH,
a cubed sphere split along k, flat caps, a split collar (a lens's two sharp
corners), and a half collar.

### Acceptance (2026-10-28): export spec met

**Implementation.** As decided:
- `shore/io/cgns.py`: `write_cgns`; `_links` pairs the connected faces
  and `_locate` places the smaller face inside the larger one;
  `transform_vector`; the CGNS/HDF5 node writer.
- `shore/io/plot3d.py`: `write_plot3d`.
- `WriteCgns` / `WritePlot3d` and `Project.write_cgns` /
  `Project.write_plot3d`; `shore cgns`, `shore plot3d`.
- The `[cgns]` extra; `h5py` and `pycgns` in `dev`.

**Verdicts** (T2 locally; `scripts/check_spec.py` reports all files ok):
- `tests/spec/test_exports.py`: **61 passed, 0 skipped, 0 xfailed**. All
  61 `@PENDING` lines are removed.
- Unit tests (`tests/test_cgns_plot3d.py`, 12): Transform cases and its
  inverse, the Plot3D record layout, and the errors (missing `h5py`, a
  connection outside its base, repeated names).

**Not verified:** a run in CI (minutes exhausted until November); CGNS
files opened in ParaView's GUI or read by a CFD solver.

### Review questions (for the user)

1. **FREE faces as `BCTypeUserDefined` + `shore:FREE`.** Should they
   instead be written as overset connectivity (`GridConnectivity_t` of type
   Overset), or as a family per component?
2. **Dependencies.** Accept the new `[cgns]` extra (`h5py`), and `pycgns`
   (LGPL-2, test-only) in `dev`?
3. **CH's multi-partner face.** The CGNS writer works around it. Should the
   face graph (and `cc.par`) get real multi-partner faces? That would be a
   later ADR.
