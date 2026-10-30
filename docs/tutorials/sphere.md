---
title: A sphere, end to end
---

# A sphere, end to end

You build a six-block cubed-sphere grid around a sphere, starting from
nothing: an analytic STL, the near-body grid, a quality check that finds a
problem, a fix, a VTK file for ParaView and a journal that replays the whole
run bit for bit. Every step is shown from the command line and from Python;
both produce the same bytes.

::: info You need
- `shore-mesh` (no extras); the VTK export in step 6 needs the `vis` extra
  (`pip install "shore-mesh[vis]"`).
- About 5 minutes. Each command runs in a few seconds.
- The script: [`examples/tutorials/sphere.py`](https://github.com/szaghi/shore/blob/main/examples/tutorials/sphere.py).
:::

## 1. Make the body

SHORE needs a closed surface to march from. A sphere of radius 0.5 (an
icosphere with 5120 triangles) is enough to see every part of the
cubed-sphere topology.

::: code-group
```bash [CLI]
shore body sphere --radius 0.5 --subdivisions 4 -o sphere.stl
```
<<< @/../examples/tutorials/sphere.py#body{python} [Python]
:::

```
✓ Wrote sphere.stl  (2562 vertices, 5120 faces, sphere)
```

## 2. Build and march the grid

The cubed sphere is four equatorial blocks (`sub0`..`sub3`) and two polar caps.
The wall layer (k = 0) is ray-cast onto the STL from the body's centroid, then
marched outward layer by layer.

::: code-group
```bash [CLI]
shore mesh sphere.stl --ni 33 --nj 64 --nk 25 --ds 0.002 --growth 1.15 \
    --theta-cap-deg 35 -o sphere -v
```
<<< @/../examples/tutorials/sphere.py#mesh{python} [Python]
:::

```
✓ Loaded sphere.stl  (5120 triangles)
  Marching       ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 24/24 layer 24/24, min J 8.13e-05 0:00:00
✓ Wrote 6 .geo files  (sub: 33x17x25, caps: 17x17x25)
```

The counts are **nodes**. `--nk 25` is 24 cells off the wall. `--nj` is the
node count of the whole equator ring: each of the four sub-blocks gets
`64 // 4 = 16` cells (17 nodes), and so does each side of a cap.

`--theta-cap-deg` is the half-angle of the polar caps (`--theta-cap` takes
radians). With $\theta_{cap} = 35°$ the cap's radial step
$\theta_{cap}/16 = 0.038$ rad and the equator's
$(\pi - 2\theta_{cap})/32 = 0.060$ rad are within a factor of 2, so the
equator meridians are pinned to the cap edges with a continuous first
interval (C1).

::: warning Do not keep the default cap angle
Left at its default (0.157 rad, 9°), the caps are small and their cells
much finer than the equator's. With these counts SHORE warns and falls back
to an unpinned distribution:

```
UserWarning: Cap radial step (~0.0138 rad) is incompatible with equator 'uniform' law
(natural first step ~0.0884 rad); falling back to unpinned spacing. ...
```

Worse, the march then folds 40 cells at the cap corners and `shore mesh`
still exits with code 0. `shore check` on `sphere_cap_north.geo` and
`sphere_cap_south.geo` reports `✗ Inverted cells detected.` Always check
(step 4).
:::

<Figure name="sphere-blocks" caption="The k = 0 wall of the cubed sphere, one colour per block: four equatorial blocks and the north cap." />

## 3. Read the wall spacing

The wall-normal spacing is the default geometric law: the first cell is
`ds = 0.002`, each next cell is `growth = 1.15` times the previous. The total
thickness is

$$
\Delta s \, \frac{g^{\,n_k - 1} - 1}{g - 1} = 0.002 \cdot \frac{1.15^{24} - 1}{0.15} \approx 0.368,
$$

so the outer layer sits at $r \approx 0.87$ (`shore info sphere_sub0.geo`
reports `x bounds [-1.3194e-05, 0.868336]`). To fix the thickness instead of the
growth rate, use a tanh law:

```bash
shore mesh sphere.stl --ni 33 --nj 64 --nk 25 --theta-cap-deg 35 \
    --volume-k-spacing tanh --volume-k-thickness 0.37 -o tanh/sphere
```

In Python the same is `build_topology(..., k_law="tanh", k_thickness=0.37)`.
The first cell then follows from `nk`, the thickness and `--volume-k-beta`
(3.0 by default): the schedule gives $5.2 \times 10^{-4}$ here (step 4
explains why the default marching scheme can make it thinner). See
[Wall-normal spacing](/cookbook/wall-spacing) and
[Spacing laws](/concepts/spacing-laws).

## 4. Check the quality

`shore check` reports the smallest cell Jacobian of every k-layer of a block
and exits with code 1 if any is negative. It is a volume, not a normalised
measure, so the values scale with the body; `warn` marks a layer below
`--warn-threshold` (default `1e-06`), which on this small body is the thin
first layer, not a defect.

```bash
shore check sphere_sub0.geo
```

```
       Jacobian check -
        sphere_sub0.geo
┏━━━━━━━┳━━━━━━━━━━━━┳━━━━━━━━┓
┃ Layer ┃       jmin ┃ Status ┃
┡━━━━━━━╇━━━━━━━━━━━━╇━━━━━━━━┩
│     1 │ 8.2299e-07 │ warn   │
│     2 │ 1.0009e-06 │ ok     │
...
│    24 │ 6.3937e-05 │ ok     │
└───────┴────────────┴────────┘
✓ All Jacobians positive.
```

No folded cell, but positivity is a low bar. `shore.quality.summary` gives
the normalised metrics the acceptance specs are written in (scaled Jacobian,
angle of the first k-segment to the wall normal, aspect ratio, k growth); the
script adds the range of the first cell height:

<<< @/../examples/tutorials/sphere.py#quality{python}

```
scaled_jacobian_min          0.7732
folded_cells                 0
wall_orthogonality_max_deg   18.43
wall_orthogonality_mean_deg  8.676
aspect_ratio_max             41.42
k_growth_max                 1.253
first_cell_min               0.0015
first_cell_max               0.002
```

Two numbers are wrong for a wall-resolved grid. The first k-line leaves the
wall up to 18° off the normal, and the first cell is as thin as
`0.0015 = 0.75 ds`. Both come from the default marching scheme
(`algebraic`): after each normal step, Laplacian smoothing moves the nodes,
and the limiter lets it pull a node back by up to 0.25 of the step and shift
it sideways by up to 0.25 of the step. On the thin first layers of this grid
the limit is reached.

## 5. March again with the hyperbolic scheme

The `hyperbolic` scheme (Steger–Chaussee) solves for orthogonal k-lines and
prescribed cell volumes instead of smoothing. In a project, `undo` reverts
the march (and the output written after it); march again and rewrite:

<<< @/../examples/tutorials/sphere.py#hyperbolic{python}

```
scaled_jacobian_min          0.7395
folded_cells                 0
wall_orthogonality_max_deg   2.74
wall_orthogonality_mean_deg  0.3863
aspect_ratio_max             31.2
k_growth_max                 1.155
first_cell_min               0.001998
first_cell_max               0.002036
```

The first cell is now `ds` within 1.8 %, the worst wall angle is 2.7° and
the mean below half a degree. On the command line, add `--march-scheme hyperbolic`:

```bash
shore mesh sphere.stl --ni 33 --nj 64 --nk 25 --ds 0.002 --growth 1.15 \
    --theta-cap-deg 35 --march-scheme hyperbolic -o sphere --journal sphere_journal.py
```

<Figure name="quality-histogram" ext="svg" caption="Scaled Jacobian and k-orthogonality of the tutorial grid." />

## 6. Export for ParaView

`shore export` bundles several `.geo` files into one VTK MultiBlock, with the
cell Jacobian attached as cell data (`--no-jacobian` to skip it).

```bash
shore export sphere_*.geo -o sphere.vtm
```

```
✓ Wrote sphere.vtm (6 blocks: sphere_cap_north, sphere_cap_south, sphere_sub0,
sphere_sub1, sphere_sub2, sphere_sub3)
```

The `.vtm` is a manifest; the blocks are in the directory `sphere/` next to
it. From Python, `p.mesh.write_vtk("sphere")` writes the same bundle.

<Figure name="marching-layers" caption="The marched volume with a quarter removed, k-layers coloured by depth." />

## 7. Replay the journal

Every project command is one line of a journal: a plain Python script of
`Project` calls. `--journal` on the CLI and `p.save_journal(...)` in Python
write it:

<<< @/../examples/tutorials/sphere.py#journal{python}

```python
# SHORE journal: replay with `shore run <this file>`.
from shore.project import Project

p = Project()
p.load_geometry('sphere.stl')
p.build_topology('cubed_sphere', nk=25, ds=0.002, ni=33, nj=64, theta_cap=0.6108652381980153)
p.march(march_scheme='hyperbolic')
p.write_geo('sphere')
```

The `undo`ne default march is not in it: the journal is the history that
produced the grid. Replaying it rewrites the six `.geo` files bit for bit:

```bash
shore run sphere_journal.py
```

```
✓ Replayed 4 commands from sphere_journal.py
```

The journal written by `shore mesh` lists every surface-spacing keyword
explicitly (`i_spacing='uniform', i_beta=3.0, ...`); the one written from
Python leaves out the defaults. Both replay to the same bytes. The cap angle is recorded in radians. A journal is
executed as Python: replay only journals you trust.

## What you learned

- `shore body` and `shore.shapes` make analytic test bodies.
- `shore mesh` / `Project.build_topology` + `march` build a cubed sphere;
  `--nj` is the equator-ring node count, split four ways.
- `--theta-cap` decides whether the cap and equator spacings join with C1
  continuity; the 9° default folds cap-corner cells without failing the run.
- `ds`, `growth` and `nk` fix the wall spacing and the thickness; a tanh law
  fixes the thickness instead.
- `shore check` finds folded cells; `shore.quality.summary` finds skewed and
  thinned first cells. The `hyperbolic` scheme keeps the first cell at `ds`
  and normal to the wall.
- `undo` reverts a command; a journal replays the history bit for bit.

## Next

- [2. A wing with a C-grid](/tutorials/wing): sharp trailing edges, wakes and
  far-field channels.
- [Cubed-sphere topology](/concepts/topologies/cubed-sphere) and
  [k=0 construction](/concepts/k0-construction): how the caps and the
  equator are built.
- [Marching schemes](/concepts/hyperbolic-marching) and
  [Grid quality](/concepts/quality).
- [`shore mesh`](/reference/mesh), [`shore check`](/reference/check),
  [`shore export`](/reference/export), [`shore run`](/reference/run).
