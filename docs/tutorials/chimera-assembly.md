---
title: A Chimera assembly
---

# A Chimera assembly

You build the two grids of a minimal overset (Chimera) case: a six-block
cubed-sphere near-body grid around a sphere and a five-block cylindrical
background grid, generated independently in the same coordinates. You keep
them as two components of one SHORE project, measure how well the background
resolves the overlap, and look at the pair in ParaView. The next tutorial
turns this project into the files the Xall solver reads.

::: info You need
- `shore-mesh` (no extras); `shore export` needs the `vis` extra.
- About 10 minutes. Each command runs in about a second.
- The script: [`examples/tutorials/chimera_assembly.py`](https://github.com/szaghi/shore/blob/main/examples/tutorials/chimera_assembly.py).
- [1. A sphere, end to end](/tutorials/sphere) explains the near-body
  parameters.
:::

| Grid | Topology | Blocks | Extent |
|---|---|---|---|
| Near-body | [cubed sphere](/concepts/topologies/cubed-sphere) | 6 | sphere of radius 0.5, marched to r ≈ 0.7 |
| Background | [flat-caps cylinder](/concepts/topologies/primitives) | 5 | r ≤ 10, −10 ≤ z ≤ 50, no hole |

SHORE generates the grids and the solver files. Hole cutting, fringe
identification and donor search are done by the overset pre-processor of the
solver, guided by the marker boxes and boundary codes of the next tutorial
([Chimera grids](/concepts/chimera)).

## 1. Make the body

::: code-group
```bash [CLI]
shore body sphere --radius 0.5 --subdivisions 4 -o sphere.stl
```
<<< @/../examples/tutorials/chimera_assembly.py#body{python} [Python]
:::

## 2. Build the near-body grid

A cubed sphere marched 12 cells off the wall with a geometric law
(`ds = 0.01`, `growth = 1.1`): total thickness about 0.21, outer layer at
r ≈ 0.7. The 30° caps keep the cap and equator spacings within the 2x
window for C1 pinning (cap step $0.524/10 = 0.052$ rad, equator step
$(\pi - 1.047)/29 = 0.072$ rad).

::: code-group
```bash [CLI]
shore mesh sphere.stl --ni 30 --nj 40 --nk 13 --ds 0.01 --growth 1.1 \
    --theta-cap-deg 30 -o wall --adjacency-json wall.adjacency.json -v
```
<<< @/../examples/tutorials/chimera_assembly.py#nearbody{python} [Python]
:::

```
✓ Wrote 6 .geo files  (sub: 30x11x13, caps: 11x11x13)
✓ Wrote wall.adjacency.json
```

The `.adjacency.json` sidecar records every face's boundary condition and
its partner: WALL at the sphere (`k_lo`), FREE at the outer layer (`k_hi`),
SHARED or DIRICHLET at the block seams. `shore cc-par` reads it in the next
tutorial.

## 3. Build the background grid

The flat-caps cylinder is a central square block and four blocks that fill
the ring between the square and the outer circle, extruded along z. There
is no hole: the background covers the body too, and the solver blanks the
cells inside it.

::: code-group
```bash [CLI]
shore primitive flat-caps --r-out 10 --z0 -10 --z1 50 --inner-frac 0.4 \
    --ni 16 --nc 24 --nk 60 -o background --adjacency-json background.adjacency.json
```
<<< @/../examples/tutorials/chimera_assembly.py#background{python} [Python]
:::

```
wrote background_center.geo
wrote background_sub_e.geo
wrote background_sub_n.geo
wrote background_sub_w.geo
wrote background_sub_s.geo
wrote background.adjacency.json
```

`--inner-frac 0.4` puts the corners of the central square at r = 4: the
whole near-body grid sits inside the `center` block. The CLI and the script
write the same bytes for every `.geo` and sidecar.

## 4. Keep both grids in one project

A project holds several named components, each with its own grid, and saves
them with the journal that built them. On the command line, `--project DIR
--component NAME` adds the result of `shore mesh` or `shore primitive` to the
project in `DIR` (created if missing):

::: code-group
```bash [CLI]
shore mesh sphere.stl --ni 30 --nj 40 --nk 13 --ds 0.01 --growth 1.1 \
    --theta-cap-deg 30 --project assembly --component sphere
shore primitive flat-caps --r-out 10 --z0 -10 --z1 50 --inner-frac 0.4 \
    --ni 16 --nc 24 --nk 60 --project assembly --component background
```
<<< @/../examples/tutorials/chimera_assembly.py#save{python} [Python]
:::

```
✓ Saved project assembly (component 'sphere')
✓ Saved project assembly (component 'background')
```

The project directory holds `project.toml`, `journal.py`, a copy of the STL
under `inputs/` (with its SHA-256) and the grids under
`grids/<component>/`. The CLI also writes the blocks as `.geo` files in the
working directory, under the default stems (`sphere_*.geo`,
`flat_caps_*.geo`). The journal of the CLI run:

```python
p = Project()
p.use_component('sphere')
p.load_geometry('inputs/00_sphere.stl', name='sphere')
p.build_topology('cubed_sphere', nk=13, ds=0.01, growth=1.1, geometry='sphere', ni=30, nj=40,
                 theta_cap=0.5235987755982988, i_spacing='uniform', ...)
p.march()
p.use_component('background')
p.build_primitive('flat_caps', r_out=10.0, z0=-10.0, z1=50.0, ni=16, nc=24, nk=60,
                  inner_frac=0.4, axis='z', ...)
```

See [6. Projects and replay](/tutorials/projects) for loading, undo and
replay.

## 5. Measure the overlap

The solver interpolates between the grids where they overlap. The donor
cells of the background should be about the size of the near-body cells
they serve. The script measures both:

<<< @/../examples/tutorials/chimera_assembly.py#overlap{python}

```
near-body outer layer: r = 0.686 .. 0.714
near-body outer cells: normal 0.029, along the layer up to 0.122
background centre block: (24, 24, 60) nodes, cells 0.246 across, 1.017 along z
```

A background cell here is 2x the largest near-body outer cell across and 8x
along z; one background cell along z is longer than the near-body grid's
diameter (1.43). This layout shows the workflow, not a grid to run on: for a
real case refine the background around the body (more `--nk` and `--nc`, a
shorter domain, or an intermediate refined box as a third component) until
the cell sizes match in the overlap.

## 6. Look at it

```bash
shore export wall_*.geo -o wall.vtm
shore export background_*.geo -o background.vtm
```

Open both `.vtm` files in ParaView. A slice at z = 0 shows the near-body
grid inside the `center` block of the background; the cell Jacobian is
attached as cell data.

<Figure name="chimera-assembly" caption="The sphere's near-body grid inside the cylindrical background, cut away." />

## What you learned

- An overset case is a set of grids built independently in one coordinate
  system; SHORE builds them, the solver's pre-processor connects them.
- `--adjacency-json` writes the face graph the solver files are made from.
- `--project DIR --component NAME` (or `Project.use_component`) keeps
  several grids in one project with one journal.
- The overlap needs comparable cell sizes on both sides; measure them before
  writing solver files.

## Next

- [4. The Xall bundle](/tutorials/xall-bundle): marker boxes, `cc.par`,
  `.grd`, `proc.input`, merged by hand or in one `shore xall` call.
- [Chimera grids](/concepts/chimera),
  [Primitives](/concepts/topologies/primitives),
  [Projects and journals](/concepts/projects).
- [`shore primitive`](/reference/primitive), [`shore mesh`](/reference/mesh),
  [`shore export`](/reference/export).
