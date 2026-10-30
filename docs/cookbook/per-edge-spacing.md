# Override the spacing of one block edge

A `HexBlock` has 12 edges, 4 along each of i, j and k. By default the 4
parallel edges of an axis share one `Spacing1D` law and the block is the
tensor product of the three distributions. Give one edge its own law and the
primitive builders blend the four parallel edges by transfinite
interpolation (TFI); give the k-edges a law and the marcher uses it as the
wall-normal schedule. Use this when one side of a block needs clustering and
the opposite side does not.

The code is
[`examples/cookbook/per_edge_spacing.py`](https://github.com/szaghi/shore/blob/main/examples/cookbook/per_edge_spacing.py)
(about 2 s).

## Name the edges

An edge's name gives the two axes that are constant along it:

| Edges along | Names |
|---|---|
| i | `j_lo_k_lo`, `j_hi_k_lo`, `j_lo_k_hi`, `j_hi_k_hi` |
| j | `i_lo_k_lo`, `i_hi_k_lo`, `i_lo_k_hi`, `i_hi_k_hi` |
| k | `i_lo_j_lo`, `i_hi_j_lo`, `i_lo_j_hi`, `i_hi_j_hi` |

`j_lo_k_lo` runs along i at `j = 0`, `k = 0`. `block.edge(name)` returns
the `Edge`; its `spacing` field holds the law.

## Cluster one edge of a box

Pass `edge_pins={name: Spacing1D}` to `BoxMesh.build`. Once any edge of an
axis is pinned, that axis is built by TFI between its 4 edges; the unpinned
ones keep the axis default (`i_law`, `j_law`, `k_law`).

<<< @/../examples/cookbook/per_edge_spacing.py#box{python}

```text
pinned edge (j=0, k=0) first i-step 0.00149
opposite edge          first i-step 0.25
mid                    first i-step 0.188
baseline               first i-step 0.25
```

The pinned edge clusters at `vmin`, the opposite edge stays uniform
(`10 / 40 = 0.25`) and the interior blends the two.

Which builders take `edge_pins`:

| Builder | i-edges | j-edges | k-edges |
|---|---|---|---|
| `BoxMesh`, `HGridRectMesh` | yes | yes | yes |
| `AnnularCylinderMesh`, sector | yes | yes | yes |
| `AnnularCylinderMesh`, full annulus | yes | refused (j periodic) | yes |
| `FlatCapsCylinder` | refused | refused | refused |

`FlatCapsCylinder` refuses `edge_pins` altogether (its 5 blocks must stay
coincident at every k-layer): use its `i_law` / `c_law` / `k_law`. A
refusal raises `ParameterError`.

## Set the k-schedule of a cubed sphere

The marching topologies read only the k-edges. Assign the **same**
`Spacing1D` instance to the 4 k-edges of every block:

<<< @/../examples/cookbook/per_edge_spacing.py#cubed_sphere{python}

```text
cubed sphere: 6 blocks, steps 2.86e-04..2.44e-01, sum 1.938
```

The layers sum to slightly less than `total_thickness = 2.0` because the
limited smoothing pulls each node back by up to a quarter of its step.

## Set the k-schedule of an O-grid

The same on the single O-grid block, here with the two-sided law:

<<< @/../examples/cookbook/per_edge_spacing.py#ogrid{python}

```text
O-grid: steps 1.78e-02..1.53e-02, sum 1.851
```

## Precedence rules

`Mesh.march()` picks one k-schedule for all blocks of the mesh, in this
order:

1. the `steps=` argument;
2. block 0's k-edges, when all 4 carry the same user-assigned `Spacing1D`;
3. the `ds` / `growth` arguments of `march()` (geometric, defaults
   `1e-3` / `1.15`).

- Overrides are detected by identity: `HexBlock` gives its 4 k-edges one
  default instance, and an edge whose `spacing is not` that instance is an
  override. The 4 overrides of a block must be equal in value (one shared
  instance is simplest).
- A block with some k-edges overridden and others not, or with 4 unequal
  overrides, warns and falls back to its default.
- A block whose override differs from block 0's warns; block 0 wins. All
  blocks share one schedule because the k-coupled seams must stay
  coincident layer by layer.
- `Project.march` (and `shore mesh`) passes the schedule built from the
  `build_topology` spacing arguments as `steps=`, so per-edge k overrides
  apply to the `Mesh` API, not to a project.
- The i and j edges of the marching topologies are not read: their wall
  grid comes from the topology's own arguments (`i_spacing`, `j_spacing`,
  `theta_cap`, ...).

## Pin the end intervals of an edge

`Spacing1D(law="pinned", first_step=..., last_step=..., interior=...)`
fixes the first and last interval lengths and shapes the interior with
`interior` (`uniform`, `tanh`, `tanh2`): the C1 match the cubed sphere uses
at its cap–equator seams. On an existing grid, apply it with
[`Project.redistribute`](/cookbook/redistribute), which carries it through
the seams.

## See also

- [Set the wall-normal spacing](/cookbook/wall-spacing): the k-schedule
  through `shore mesh` and `Project.build_topology`.
- [Spacing laws](/concepts/spacing-laws) and the
  [block model](/concepts/block-model).
- API: [`Spacing1D`](/reference/api/shore.mesh.spacing#Spacing1D),
  [`BoxMesh.build`](/reference/api/shore.primitive.box#BoxMesh-build),
  [`Mesh.march`](/reference/api/shore.mesh.mesh#Mesh-march).
- The longer scripts in `examples/03-per-edge-spacing/` compare the
  overridden grids with the defaults and open a viewer.
