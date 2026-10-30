# Check grid quality

Check a grid before you hand it to the solver: no folded cells, a bounded
scaled Jacobian, wall-normal lines that leave the wall straight, and a
smooth growth away from it. `shore check` gives the folded-cell test on a
`.geo` file; `shore.quality` gives the metrics on any block list or node
array.

The code is
[`examples/cookbook/quality.py`](https://github.com/szaghi/shore/blob/main/examples/cookbook/quality.py)
(about 5 s). It marches a cubed sphere first:

<<< @/../examples/cookbook/quality.py#build{python}

## Check a `.geo` file for folded cells

```bash
shore check sphere_cap_north.geo
```

```text
       Jacobian check -
     sphere_cap_north.geo
┏━━━━━━━┳━━━━━━━━━━━━┳━━━━━━━━┓
┃ Layer ┃       jmin ┃ Status ┃
┡━━━━━━━╇━━━━━━━━━━━━╇━━━━━━━━┩
│     1 │ 1.6730e-06 │ ok     │
│     2 │ 1.9557e-06 │ ok     │
...
✓ All Jacobians positive.
```

- `jmin` is the smallest triple product `(e_i × e_j) · e_k` of the cells
  between layers `k-1` and `k`: a volume, in geometry units cubed, not a
  normalised measure.
- A layer with `jmin <= 0` holds a folded cell; the command then exits
  with code 1, so it can gate a script.
- `warn` marks layers with `0 < jmin < --warn-threshold` (default
  `1e-6`). Thin wall cells trip it on sound grids: a first cell of `5e-5`
  on a unit sphere gives `jmin` near `1e-7`. Scale the threshold to your
  cell volumes (`--warn-threshold 1e-9`) or read the scaled Jacobian below.
- `Mesh.check()` returns the same layer minima as a dict per block.

## Summarise a whole grid

<<< @/../examples/cookbook/quality.py#summary{python}

```text
scaled_jacobian_min          0.7695
folded_cells                 0
wall_orthogonality_max_deg   18.43
wall_orthogonality_mean_deg  14.64
aspect_ratio_max             165.7
k_growth_max                 1.271
```

| Key | Metric | Good |
|---|---|---|
| `scaled_jacobian_min` | smallest corner scaled Jacobian of any cell (1 = orthogonal corner, <= 0 = folded) | > 0.3 for a marched grid |
| `folded_cells` | cells with scaled Jacobian <= 0 | 0 |
| `wall_orthogonality_*_deg` | angle between the first k-segment and the wall-layer normal | a few degrees |
| `aspect_ratio_max` | longest over shortest edge of a cell | large at the wall by design |
| `k_growth_max` | largest ratio of consecutive k-segments | <= 1.2–1.3 |

The "good" column gives rules of thumb, not SHORE limits; your solver
decides. `summary` takes `HexBlock`s or bare `(ni, nj, nk, 3)` arrays.

<Figure name="quality-jacobian" caption="A grid coloured by the scaled Jacobian of its cells." />

## Compare blocks and find the worst cell

The functions behind `summary` return one value per cell or node, so you
can slice them per block, per layer or per region:

<<< @/../examples/cookbook/quality.py#per_block{python}

```text
block       min SJ max wall deg   max AR max growth
sub0         0.771        18.43      166      1.242
sub1         0.770        18.43      166      1.239
sub2         0.771        18.43      166      1.242
sub3         0.770        18.43      166      1.239
cap_north    0.917        18.24       82      1.271
cap_south    0.917        18.24       82      1.271
```

<<< @/../examples/cookbook/quality.py#worst{python}

```text
worst cell: block sub3, cell (31, 15, 0), SJ 0.770
  centre [ 0.5801 -0.0189 -0.8134], 0 cells below 0.5
```

| Function | Shape of the result |
|---|---|
| `scaled_jacobian(nodes)` | `(ni-1, nj-1, nk-1)`, per cell |
| `aspect_ratio(nodes)` | `(ni-1, nj-1, nk-1)`, per cell |
| `k_orthogonality_deg(nodes, k=0)` | `(ni, nj)`, per node of layer `k` |
| `k_spacing(nodes)` | `(ni, nj, nk-1)`, per k-segment |
| `k_growth(nodes)` | `(ni, nj, nk-2)` |
| `k_angle_deg(nodes, directions, k=0)` | `(ni, nj)`, against given directions |

`k_orthogonality_deg` uses a one-sided layer normal on block edges; on a
multi-block wall `k_angle_deg` against the body normals
(`Body.normals(points)`) avoids that bias.

## Compare march settings

The table above shows wall lines tilted by 14.6° on average. The default
algebraic march with `smooth_mode="limited"` lets the smoothing move a node
sideways by up to a quarter of its step, `atan(0.25) = 14°`, and on a
sphere almost every first-layer node uses that allowance. Measure before
you choose:

<<< @/../examples/cookbook/quality.py#compare{python}

```text
{}                             wall mean 14.64 deg, max 18.43 deg, min SJ 0.770
{'smooth_iters': 0}            wall mean  0.32 deg, max  2.76 deg, min SJ 0.773
{'march_scheme': 'hyperbolic'} wall mean  0.34 deg, max  2.74 deg, min SJ 0.747
```

- On a smooth convex body, `march_scheme="hyperbolic"` (`--march-scheme
  hyperbolic`) or no smoothing gives wall lines within 3° of the normal.
- Smoothing earns its place on concave regions and at sharp features,
  where unsmoothed layers fold; check the folded-cell count after each
  change.

## Check `.geo` files from any source

<<< @/../examples/cookbook/quality.py#geo_files{python}

```text
6 .geo files, min SJ 0.770
```

`read_geo` returns the `(ni, nj, nk, 3)` array of one block, so the metrics
work on grids from other generators too.

## See also

- [Grid quality](/concepts/quality): the definitions.
- [Marching schemes](/concepts/hyperbolic-marching): smoothing, the limiter
  and the hyperbolic scheme.
- [Set the wall-normal spacing](/cookbook/wall-spacing).
- CLI: [`shore check`](/reference/check), [`shore info`](/reference/info),
  [`shore export`](/reference/export) (VTK with the Jacobian, for
  ParaView).
- API: [`shore.quality`](/reference/api/shore.quality).
