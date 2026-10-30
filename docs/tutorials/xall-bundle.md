---
title: The Xall bundle
---

# The Xall bundle

You turn the two grids of [3. A Chimera assembly](/tutorials/chimera-assembly)
into the files the Xall solver chain reads: marker boxes that cut the body
out of the background, one binary `.grd` with every block, the `cc.par`
connectivity and boundary file, and the `proc.input` MPI distribution. You
do it twice: in one `shore xall` call on the project, and by hand, file by
file, which is what the one call does for you. You also meet one defect of
the background grid that you have to fix by hand today.

::: info You need
- `shore-mesh` (no extras); `shore boxes-export` needs the `vis` extra.
- About 15 minutes. The script runs in about 7 s, most of it in the marker
  boxes.
- The script: [`examples/tutorials/xall_bundle.py`](https://github.com/szaghi/shore/blob/main/examples/tutorials/xall_bundle.py).
  It rebuilds the grids of tutorial 3 first, so it runs on its own; the CLI
  commands below run in the directory where you ran tutorial 3.
:::

| File | Read by | Content |
|---|---|---|
| `.grd` | Xall, overset pre-processor | Fortran unformatted binary: block count, cell counts per block, then X, Y, Z of every block |
| `cc.par` | overset pre-processor | header options, per-block level / group / priority, one patch per block face (boundary code or connection), marker boxes |
| `proc.input` | Xall at MPI start-up | per block: group, body, MPI rank |

The three files share one 1-based block numbering. Whatever writes them must
use one block order for all three. The formats are in
[`shore grd`](/reference/grd), [`shore cc-par`](/reference/cc-par) and
[`shore proc-input`](/reference/proc-input).

<<< @/../examples/tutorials/xall_bundle.py#grids{python}

## 1. Generate the marker boxes

The background covers the body. Marker boxes tell the overset
pre-processor which background cells lie inside the body and must be blanked.
`voxel-fill` (the default) voxelises the body's bounding box and keeps the
voxels inside the body, merged into larger boxes; the voxel size is
`max_outward_gap / sqrt(3)`, so no box surface is farther than the gap
outside the body. Keep the gap below the near-body thickness (about 0.21
here), so the boxes end inside the near-body grid.

::: code-group
```bash [CLI]
shore boxes-from-stl sphere.stl --block sub0 --max-outward-gap 0.10 -o sphere_boxes.json
shore boxes-export sphere_boxes.json -o sphere_boxes.vtm
```
<<< @/../examples/tutorials/xall_bundle.py#boxes{python} [Python]
:::

```
✓ Wrote sphere_boxes.json (101 box(es), method='voxel-fill')
```

`--block sub0` is the near-body block the boxes belong to in `cc.par`; any
block of that grid will do. See [Marker boxes](/concepts/marker-boxes) for
the algorithm and the `aabb` / `obb` alternatives.

<Figure name="marker-boxes" caption="The sphere STL and its voxel-fill marker boxes." />

## 2. Write the run in one call

`shore xall` writes the Xall run of every component of a saved project. Its
options file has one `[run]` table for the run-wide options and one
`[components.<name>]` table per component:

```toml
# xall.toml
[run]
grd = "assembly.grd"
np = 4

[components.background]
free_face_bc = 9   # natural BC 9: extrapolation at the far field
body = 0

[components.sphere]
body = 1
boxes = "sphere_boxes.json"

[components.sphere.meta.cap_north]
level = 2

[components.sphere.meta.cap_south]
level = 2
```

::: code-group
```bash [CLI]
shore xall assembly --out run --options xall.toml
```
<<< @/../examples/tutorials/xall_bundle.py#xall{python} [Python]
:::

```
✓ Wrote the Xall run of ['sphere', 'background'] to run
```

```
run/
  assembly.grd           # the 11 blocks, sphere first (component creation order)
  cc.par                 # merged
  proc.input             # balanced once over all 11 blocks
  sphere/                # per component: .geo files, sidecar, sphere.grd, cc.par
  background/
```

What the options mean:

| Option | Effect |
|---|---|
| `free_face_bc` | the code of the component's FREE faces: 40 (the default) is a chimera fringe that receives from other grids; a far-field background takes a natural boundary code in 1–19, e.g. 9 (extrapolation), 3 (inflow), 4 (outflow) |
| `body` | the Xall body number of all the component's blocks in `proc.input` (forces and moments are integrated per body) |
| `boxes` | marker boxes for the component's `cc.par`, a file relative to the options file |
| `meta.<block>` | per-block `level`, `group`, `priority`, `comment`, `free_face_bc` in `cc.par`; here the caps get overset level 2 |
| `proc_meta.<block>` | per-block `proc` / `group` / `body` in `proc.input` |

Unknown tables and keys, components without a grid and block labels a
component does not have are errors, reported before anything is written. The
[`shore xall`](/reference/xall) page lists the `[run]` options.

::: warning Never a negative far-field code
Older SHORE documentation used `-9` for far-field faces. Xall rejects any
face code below 1 ("icc sulle facce deve essere > 0"): the natural codes are
positive, 1 to 19.
:::

::: warning The background's z0 end is written as a wall
The flat-caps primitive leaves the `k_lo` face of its five blocks (the
`z0 = -10` end of the cylinder) as WALL, and `cc.par` gets code 1, a
viscous wall, there:

```
        7         5         1         0         0        23         0        23         0         0 !        41
```

(block 7, face 5 = `k_lo`, code 1). The `z1` end gets the 9 you asked for.
`shore xall` has no option to change it. In the file-by-file route below,
make the face FREE in the sidecar before writing `cc.par` (step 3).
:::

## 3. The same, file by file

`shore xall` runs the steps below for you. Run them yourself when the
components do not live in one project, or to control each file. The shell
script
[`examples/05-chimera-assembly/chimera_assembly.sh`](https://github.com/szaghi/shore/blob/main/examples/05-chimera-assembly/chimera_assembly.sh)
chains the same commands from the STL on (without the sidecar fix below).

**Pack the blocks.** The order of the files is the block order of the
`.grd` and must be the order `cc.par` uses: the order of the blocks in the
sidecar.

::: code-group
```bash [CLI]
shore grd wall_sub0.geo wall_sub1.geo wall_sub2.geo wall_sub3.geo \
    wall_cap_north.geo wall_cap_south.geo -o wall.grd
shore grd background_center.geo background_sub_e.geo background_sub_n.geo \
    background_sub_w.geo background_sub_s.geo -o background.grd
```
<<< @/../examples/tutorials/xall_bundle.py#grd{python} [Python]
:::

**Open the background's z0 end.** Edit the sidecar so every `k_lo` face of
the background is FREE (see the warning above):

<<< @/../examples/tutorials/xall_bundle.py#farfield{python}

**Write one `cc.par` per component.** Near-body FREE faces default to 40;
the background's take 9. The boxes go with the near-body grid.

::: code-group
```bash [CLI]
shore cc-par wall.adjacency.json wall.cc.par --grd wall.grd \
    --boxes sphere_boxes.json --meta '{"cap_north": {"level": 2}, "cap_south": {"level": 2}}'
shore cc-par background.adjacency.json background.cc.par --grd background.grd --free-bc 9
```
<<< @/../examples/tutorials/xall_bundle.py#ccpar{python} [Python]
:::

**Write one `proc.input` per component.** `--meta` keys are 1-based block
indices of that `.grd`; the near-body blocks are body 1.

::: code-group
```bash [CLI]
shore proc-input wall.grd wall.proc.input --np 4 \
    --meta '{"1": {"body": 1}, "2": {"body": 1}, "3": {"body": 1}, "4": {"body": 1}, "5": {"body": 1}, "6": {"body": 1}}'
shore proc-input background.grd background.proc.input --np 4
```
<<< @/../examples/tutorials/xall_bundle.py#procinput{python} [Python]
:::

```
proc.input: 5 blocks across 4 ranks
  rank 0: 1 blocks, 31.2k cells (27.7%)
  rank 1: 2 blocks, 40.7k cells (36.1%)
  rank 2: 1 blocks, 20.4k cells (18.1%)
  rank 3: 1 blocks, 20.4k cells (18.1%)
  imbalance: 72.3%
Warning: load imbalance is 72.3% (> 5%).  Run `shore balance` to plan a topology-safe
set of splits that brings the imbalance within tolerance.
```

**Merge.** The three merges must take the components in the same order;
each shifts the block (and patch) indices of the later inputs.

::: code-group
```bash [CLI]
shore grd-merge wall.grd background.grd -o assembly.grd
shore cc-par-merge wall.cc.par background.cc.par -o assembly.cc.par --grd assembly.grd
shore proc-input-merge wall.proc.input background.proc.input -o assembly.proc.input \
    --grd assembly.grd
```
<<< @/../examples/tutorials/xall_bundle.py#merge{python} [Python]
:::

```
✓ Wrote assembly.proc.input (2 inputs merged)
proc.input: 11 blocks across 4 ranks
  rank 0: 3 blocks, 35.9k cells (27.8%)
  rank 1: 4 blocks, 45.4k cells (35.2%)
  rank 2: 2 blocks, 23.8k cells (18.5%)
  rank 3: 2 blocks, 23.8k cells (18.5%)
  imbalance: 66.9%
```

The CLI and the script write the same bytes for all nine files.

## 4. Compare the two routes

With the near-body grid first in both, the merged `.grd` of the two routes
is the same file. The other two differ:

- `cc.par`: besides the z0 faces of the background, only the block comments
  differ (`wall.sub0` by hand, `cc.sub0` from `shore xall`).
- `proc.input`: the merge does not re-balance; it keeps each component's
  ranks, balanced on its own. `shore xall` balances once over all 11 blocks:
  imbalance 37.8 % instead of 66.9 %.

Neither is good, and no assignment can be: the background has five large
blocks (31 211 and 4 × 20 355 of the 129 k cells) for four ranks, so one rank
carries two of them, 40 710 cells against a mean of 32 238. Only splitting
blocks helps.
[Split and balance](/cookbook/split-balance) plans the splits
([`shore balance`](/reference/balance), [`shore split`](/reference/split)).

## What you learned

- Marker boxes from the body STL blank the background inside the body; the
  outward gap must be smaller than the near-body thickness.
- `cc.par` codes: 1 wall, 40 chimera fringe, 1–19 natural far-field codes,
  never negative.
- `shore xall` writes the per-component files and the merged run of a
  project in one call, from one options file, and balances `proc.input`
  over the whole run.
- The file-by-file route needs one block order for `.grd`, `cc.par` and
  `proc.input`, and its merged `proc.input` is not re-balanced.
- The flat-caps background's `z0` end comes out as a wall: fix the sidecar
  before `shore cc-par`.

## Next

- [5. Collar at a junction](/tutorials/collar).
- [Multi-component Xall run](/cookbook/multi-component),
  [After-overset splitting](/cookbook/after-overset-split).
- [Chimera grids](/concepts/chimera), [Marker boxes](/concepts/marker-boxes),
  [Load balancing](/concepts/load-balancing).
- [`shore xall`](/reference/xall), [`shore cc-par`](/reference/cc-par),
  [`shore cc-par-merge`](/reference/cc-par-merge),
  [`shore grd-merge`](/reference/grd-merge),
  [`shore proc-input`](/reference/proc-input),
  [Marker boxes reference](/reference/boxes).
