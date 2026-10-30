# File formats

Every file SHORE reads or writes, with the module that owns it. The
excerpts below were produced by running SHORE 0.4.0 on a unit sphere
(`shore body sphere`, then `shore mesh --ni 33 --nj 64 --nk 25 --ds 0.01
--theta-cap 0.6109`) and a 9 × 17 × 5 `shore primitive annulus`
background in a second component.

Node counts are used everywhere except in the `.grd` header and the
`cc.par` patch extents, which hold cell counts (= nodes − 1).

| File | Written by | Read by | Module |
|---|---|---|---|
| [`.geo`](#geo-block-coordinates) | `shore mesh`, `primitive`, `collar`, `split`, `cap-tip`, `Project.write_geo` | `shore grd`, `split`, `info`, `check`, `export`, `view`, `Project.load_blocks` | `shore.io.geo` |
| [`.grd`](#grd-xall-grid-pack) | `shore grd`, `grd-merge`, `xall` | `shore proc-input`, `balance`, Xall | `shore.io.grd` |
| [`cc.par`](#cc-par-chimera-connectivity) | `shore cc-par`, `cc-par-merge`, `xall` | Xall's overset pre-processor | `shore.io.cc_par` |
| [`proc.input`](#proc-input-mpi-distribution) | `shore proc-input`, `proc-input-merge`, `xall` | Xall | `shore.io.proc_input` |
| [`.adjacency.json`](#adjacency-json-face-graph) | `--adjacency-json`, `shore split`, `Project.write_adjacency`, `Project.save` | `shore cc-par`, `split --adjacency`, `balance --adjacency` | `shore.io.adjacency` |
| [surface `.npz`](#surface-npz) | `shore project`, `Body.save_surface` | `shore info`, `export`, `view`, `Body.load_surface` | `shore.mesh.body` |
| [project directory](#project-directory) | `Project.save`, `--project` | `Project.load`, `shore xall`, `cgns`, `plot3d`, `gui` | `shore.project.project` |
| [CGNS](#cgns) | `shore cgns`, `Project.write_cgns` | CGNS tools | `shore.io.cgns` |
| [Plot3D](#plot3d) | `shore plot3d`, `Project.write_plot3d` | Plot3D tools | `shore.io.plot3d` |
| marker boxes JSON | `shore boxes-from-stl` | `shore cc-par --boxes`, `boxes-export` | `shore.io.boxes`, see [Marker boxes](boxes) |
| journal `.py` | `--journal`, `Project.save_journal` | `shore run` | see [Project scripts](run) |
| TOML configs | you | `shore mesh -c`, `primitive -c`, `split -c` | see [Config files](config) |

## `.geo`: block coordinates

ASCII, one block per file, the Xall / `geogrd` input. Line 1 holds the
node counts `ni nj nk`; then every x, then every y, then every z, each in
Fortran order (i fastest, then j, then k), six values per line in
`%.10e`:

```
33 17 25
5.7166017204e-01  6.0943135842e-01  6.5618365045e-01  7.0249805148e-01  7.4284455485e-01  7.8161536013e-01
8.1793167500e-01  8.5344217193e-01  8.8174580937e-01  9.1004127350e-01  9.3203607998e-01  9.5262912749e-01
...
```

A file of `n = ni·nj·nk` nodes has `1 + 3·⌈n/6⌉` lines (7015 for the
33 × 17 × 25 block above). The reader accepts any whitespace layout and
checks the value count. The 11 significant digits are the precision
every downstream file packed from `.geo` files carries (`.grd`,
`proc.input` weights). Naming rules: [`shore mesh` output](mesh#output).

## `.grd`: Xall grid pack

Fortran unformatted sequential, native byte order (little endian on
x86-64): every record is wrapped in 4-byte `int32` length markers.

| Record | Content |
|---|---|
| 1 | `nb`, the block count (`int32`) |
| 2 … nb+1 | `ni−1, nj−1, nk−1` of block n (3 × `int32`, **cell** counts) |
| then, per block | X, Y, Z: three records of `ni·nj·nk` `float64`, i fastest |

The header of the 6-block sphere (`od -t d4`: marker, value, marker…):

```
4   6   4                  # record 1: 6 blocks
12  32 16 24  12           # block 1: 33 x 17 x 25 nodes -> 32 16 24 cells
12  32 16 24  12           # block 2
...
```

`shore grd` packs from `.geo` files in the order given; `cc.par` and
`proc.input` index blocks by that order. `read_grd_metadata` reads only
the header (the weight source of `proc-input` and `balance`).
`grd-merge` concatenates the block records byte for byte. Details:
[`shore grd`](grd).

## `cc.par`: chimera connectivity

ASCII, read by Xall's overset pre-processor with Fortran list-directed
reads: a 10-line header, the blocks section (`level group priority !
comment`), the patches section (one line per block face, or per
partner on a face shared with several), an edges section and the
marker-boxes section. Excerpt for the sphere:

```
's35.grd'
'cc.'
.true.
.false.
.false.

  4   1

-1.0000E-01
-1.0000E+00

6 ! blocks number

        1         0         1 ! 1 : sub0
...
36 ! patches number

        1         1       415        33         0         0         0        16         0        24 !         1
        1         2       135        25        32        32         0        16         0        24 !         2
...
        1         5         1         0         0        32         0        16         0         0 !         5
        1         6        40         0         0        32         0        16        24        24 !         6
...
0 ! edges number

0 ! boxes number
```

Patch columns: block, face (1–6 = Imin, Imax, Jmin, Jmax, Kmin, Kmax),
BC code, family (the partner patch for a connection), then the i, j, k
extents in cells. SHORE maps WALL → `1`, FREE → `40` (chimera
cell-centre recipient) or a positive natural BC in `[1, 19]` (`9`
extrapolation), SYMMETRY → `2`, and SHARED / DIRICHLET → a three-digit
orientation code `li lj lk` derived from the face's `AxisMap`. The
digit of the axis across the face is odd when the two faces lie on
opposite sides of their blocks (lo ↔ hi) and even on the same side
(lo ↔ lo, hi ↔ hi), as the pre-processor reads it
([ADR 0036](/adr/0036-cc-par-as-xall-reads-it)): a plain lo ↔ hi seam is
`135` on both sides. Full field and code reference:
[`shore cc-par`](cc-par).

## `proc.input`: MPI distribution

ASCII, read by Xall at start-up with list-directed reads: a free comment
line, two blank lines, the block count, two blank lines, a column-header
line, then one `block group body proc` record per block (`proc` is the
0-based MPI rank; the trailing `!` comment is ignored by the reader).
`shore proc-input sphere.grd proc.input --np 2`:

```
generated by SHORE


6 number of blocks


 block - group - body - processor [comment]
     1      0      0      0 ! weight=12288
     2      0      0      1 ! weight=12288
     3      0      0      0 ! weight=12288
     4      0      0      1 ! weight=12288
     5      0      0      0 ! weight=6144
     6      0      0      1 ! weight=6144
```

The weight is the block's cell count from the `.grd` header; blocks go
greedily to the least-loaded rank. Details: [`shore proc-input`](proc-input).

## `.adjacency.json`: face graph

SHORE's own sidecar, the bridge from the block model to `cc.par`, to
`shore split` and to a project loaded without replay. JSON, version 3
(readers accept 2 and 3; version 1 is rejected with a regenerate
message). Per block: `label`, node `shape`, and the six faces with `bc`
(`WALL`, `FREE`, `SHARED`, `DIRICHLET`, `PERIODIC`, `SYMMETRY`),
`partner` (`label`, `face`) and `axis_map` (per along-face axis of this
face, the partner axis and a flip flag):

```json
{
  "version": 3,
  "associations": [{"geometry": "body", "kind": "body", "keys": []}],
  "blocks": [
    {
      "label": "sub0",
      "shape": [33, 17, 25],
      "faces": {
        "i_lo": {"bc": "DIRICHLET", "partner": {"label": "cap_south", "face": "j_lo"},
                 "axis_map": [["i", false], ["k", false]]},
        "j_lo": {"bc": "SHARED", "partner": {"label": "sub3", "face": "j_hi"},
                 "axis_map": [["i", false], ["k", false]]},
        "k_lo": {"bc": "WALL", "partner": null, "axis_map": null},
        "k_hi": {"bc": "FREE", "partner": null, "axis_map": null}
      },
      "associated": {
        "faces": {"k_lo": 0},
        "edges": {"j_lo_k_lo": 0, "j_hi_k_lo": 0, "i_lo_k_lo": 0, "i_hi_k_lo": 0},
        "vertices": {"0": 0, "2": 0, "4": 0, "6": 0}
      }
    }
  ]
}
```

(two of the six faces left out). Optional keys:

- `associations` (top level) and `associated` (per block): the
  block-to-geometry associations ([ADR 0003](/adr/0003-association)),
  each distinct one stored once. An entry names the `geometry`, the
  entity `kind` (`body`, `surface`, `curve`, `point`) and the entity
  fingerprints `keys`, so they are found again in a re-imported file.
  A CAD face key:
  `{"kind": "surface", "type": "Plane", "centroid": [0.436, 0.0, 0.109],
  "measure": 0.190, "extent": [0.873, 0.0, 0.218]}` (coordinates
  normalised by the bounding box). Vertices are keyed by corner index
  (`4i + 2j + k`). Absent when nothing is associated.
- `bc_code` on a FREE face: a per-face natural BC that overrides the
  block's `free_face_bc` in `cc.par`.
- `patches` on a face shared with several partners (CH's `body.k_hi`,
  [ADR 0036](/adr/0036-cc-par-as-xall-reads-it)): per partner, the
  0-based node `range` `[[i0, j0, k0], [i1, j1, k1]]`, the `partner` and
  its `axis_map`.

## Surface `.npz`

The k=0 layer of `shore project` (and `Body.save_surface`): a NumPy
`.npz` with one array `k0` of shape `(ni, nj, 3)`, float64:

```python
>>> np.load("k0.npz")["k0"].shape     # shore project sphere.stl --ni 9 --nj 16 -o k0.npz
(9, 16, 3)
```

## Project directory

`Project.save(DIR)` (or `--project DIR` on `mesh`, `collar`,
`primitive *`) writes a portable directory
([ADR 0016](/adr/0016-multi-component-projects),
[ADR 0033](/adr/0033-load-without-replay)):

```
proj/
├── project.toml          # format, version, components and their block labels, input hashes
├── journal.py            # the state-changing commands, paths into inputs/
├── inputs/00_sphere.stl  # every input file, copied and numbered
├── adjacency.json        # face graph of the current component's mesh
└── grids/
    ├── main/             # one folder per component with a mesh
    │   ├── 0000.npz      # block 0: array "nodes", (ni, nj, nk, 3) float64
    │   ├── ...
    │   └── adjacency.json   # the component's face graph
    └── background/
        ├── 0000.npz
        └── adjacency.json
```

```toml
format = 3
shore_version = "0.4.0"
commands = 5

[components]
"main" = ["sub0", "sub1", "sub2", "sub3", "cap_north", "cap_south"]
"background" = ["annulus"]

[inputs]
"inputs/00_sphere.stl" = "b44c7c37f9338682ead3a23ca966591ab3301383adee6a3dc0b95171403f7ac4"
```

- `journal.py` holds the state-changing commands only; output commands
  (`write_*`) are dropped, since replaying them would write where the
  journal says.
- `Project.load(DIR)` checks every input's SHA-256, replays the
  journal (inputs resolved against `DIR`) and requires every replayed
  grid to equal the saved `.npz` bit for bit.
- `Project.load(DIR, replay=False)` (`--no-replay`) records the journal
  instead, loads only geometries and component switches, and reads each
  component's blocks from `grids/<component>/` and its face graph: no
  march, no undo of the loaded history, no blocking.
- Format 3 is current; formats 1 and 2 still load. Projects saved before
  format 3 replay collars with the `stations` volume
  ([ADR 0026](/adr/0026-smoothed-collar-volume-by-default)).

## Xall run directory

`shore xall PROJECT --out run` (or `Project.write_xall`) writes one
folder per component and the merged files:

```
run/
├── main/        main_<label>.geo …, main.adjacency.json, main.grd, cc.par
├── background/  background_annulus.geo, background.adjacency.json, background.grd, cc.par
├── mesh.grd     # merged, components in project order
├── cc.par       # merged: block and patch indices shifted per component
└── proc.input   # one, balanced over all blocks
```

See [`shore xall`](xall) for the options file.

## CGNS

One CGNS/HDF5 file (`h5py`, the `cgns` extra), CGNS library version 3.4
([ADR 0035](/adr/0035-cgns-plot3d-export)): a `CGNSBase_t` per component
(cell and physical dimension 3), a structured `Zone_t` per block named by
its label, double-precision `GridCoordinates`, and the face graph:

```
CGNSLibraryVersion
main                      CGNSBase_t
  sub0                    Zone_t
    GridCoordinates       CoordinateX, CoordinateY, CoordinateZ
    ZoneBC                k_lo (BCWall), k_hi (BCTypeUserDefined, shore:FREE)
    ZoneGridConnectivity  i_hi_cap_north, i_lo_cap_south, j_hi_sub1, j_lo_sub3
    ZoneType
  ...
background                CGNSBase_t
  annulus                 Zone_t
    ...
```

SHARED, DIRICHLET and PERIODIC faces become `GridConnectivity1to1_t`
nodes in both zones (named `<face>_<partner>`, with `PointRange`,
`PointRangeDonor` and the `Transform` from the `AxisMap`); WALL becomes
`BCWall`, SYMMETRY `BCSymmetryPlane`, FREE `BCTypeUserDefined` with a
`shore:FREE` descriptor. Details: [`shore cgns`](cgns).

## Plot3D

Multi-block, whole, 3-D grid file; Fortran unformatted sequential,
little endian (`<i4` markers and counts, `<f8` coordinates)
([ADR 0035](/adr/0035-cgns-plot3d-export)):

| Record | Content |
|---|---|
| 1 | block count |
| 2 | `ni, nj, nk` of every block (node counts) |
| 2 + n | X, Y, Z of block n, i fastest |

The two-component sphere project above (`od -t d4`):

```
4   7   4                # 7 blocks
84  33 17 25  33 17 25 ...  # 7 x 3 node counts (84 bytes)
```

No names, connectivity or boundary conditions: use CGNS for those.
Details: [`shore plot3d`](plot3d).

## See also

- [Block model](/concepts/block-model): faces, BCs and `AxisMap`
- [Projects and journals](/concepts/projects)
- [Chimera grids](/concepts/chimera)
- [Python API: `shore.io`](/reference/api/shore.io)
