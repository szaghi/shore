# C-grid topology

A C-grid wraps a body with a trailing edge (TE), such as an airfoil,
hydrofoil or wing, and continues downstream in two wake branches. In a
span section the grid lines form a C: the outer envelope opens behind the
TE, and the two wake branches meet along the **wake cut**. In SHORE the
k=0 layer comes from slicing the body at span stations
([C-grid construction](/concepts/c-grid-construction)), and the blocks are
marched like the other topologies. Built by
`shore mesh --topology cgrid`.

<Figure name="cgrid-wing" caption="A 3-block C-grid around a NACA 0012 wing: a span section of the volume over the wall." />

## Block layouts

`--cgrid-blocks` selects the layout:

| Layout | Blocks | TE | Notes |
|---|---|---|---|
| 3 (default) | `upper`, `body`, `lower` | sharp | wake cut as a `SHARED` face between two blocks |
| 4 | + `wake_base` | blunt | a block filling the gap between the wake branches; its `i_hi` is the TE base wall |
| 1 (legacy) | one block | sharp | the wake cut is an overlap of the block with itself, on which Xall cannot apply BCs cleanly |

Per block, j runs along the span (`--n-stations` nodes) and k away from
the wall:

| Block | Shape | i runs |
|---|---|---|
| `upper` | `(ni_wake, n_stations, nk)` | upper wake tip → upper TE |
| `body` | `(ni_body, n_stations, nk)` | upper TE → suction side → LE → pressure side → lower TE |
| `lower` | `(ni_wake, n_stations, nk)` | lower TE → lower wake tip |
| `wake_base` | `(ni_wake, n_stations, nk_w)` | wake tip → TE; k spans the TE thickness |

## Seams and boundary conditions

| Seam | Faces | Axis map |
|---|---|---|
| upper TE column | `upper.i_hi ↔ body.i_lo` | identity (shared by view) |
| lower TE column | `body.i_hi ↔ lower.i_lo` | identity |
| wake cut (3-block) | `upper.k_lo ↔ lower.k_lo` | i flipped: upper runs toward the TE, lower away from it |
| wake base (4-block) | `wake_base.k_hi ↔ upper.k_lo`, `wake_base.k_lo ↔ lower.k_lo` | identity; i flipped |

| Block | i_lo | i_hi | j_lo, j_hi | k_lo | k_hi |
|---|---|---|---|---|---|
| `upper` | `FREE` | `SHARED` | `FREE` | `SHARED` (wake cut) | `FREE` |
| `body` | `SHARED` | `SHARED` | `FREE` | `WALL` | `FREE` |
| `lower` | `SHARED` | `FREE` | `FREE` | `SHARED` (wake cut) | `FREE` |
| `wake_base` | `FREE` | `WALL` (TE base) | `FREE` | `SHARED` | `SHARED` |

In `cc.par` the `FREE` faces get the block's `free_face_bc`: 40 (chimera
fringe) by default, or a positive natural BC such as 9 (extrapolation)
for an outer boundary. The `j` faces are `FREE` because the wrap alone
does not close the wing tips: see [tip caps](#tip-caps).

The `wake_base` block is not marched: it is interpolated in k between the
two wake k=0 strips at build time, and the marcher skips it.

## Marching

`Mesh._march_cgrid_three` marches `upper`, `body` and `lower` together
with ghost rows across the two TE columns, so the TE seams keep one
normal ([Seam-aware normals](/concepts/seam-aware-normals)). Both
`march_scheme` values apply ([Marching](/concepts/hyperbolic-marching));
the hyperbolic scheme's acceptance covers the 3-block sharp-TE C-grid, not
the 4-block blunt one.

## Parameters

| Flag | Default | Meaning |
|---|---|---|
| `--ni-body` | required | body-wrap nodes per station |
| `--ni-wake` | required | nodes per wake branch, TE included |
| `--n-stations` | required | span stations (the blocks' j count) |
| `--wake-length` | required | wake branch length from the TE, body units |
| `--wake-growth` | 1.15 | geometric growth of the wake cells, starting from the TE cell |
| `--span-axis` | `z` | span axis; the chord axis is the longer other extent |
| `--te-dihedral-deg` | 30 | dihedral above which an edge is a TE candidate |
| `--cgrid-blocks` | 3 | layout |

The freestream, which sets the wake direction, defaults to the chord
axis from the LE toward the TE. `Mesh.from_stl` does not accept the
C-grid (its parameters do not fit the `(ni, nj, nk)` signature): in Python
build the blocks with `CGridTopology().build(...)` and wrap them in a
`Mesh`. The [wing tutorial](/tutorials/wing) runs the whole chain from a
NACA section.

## Tip caps

::: warning Known problems of the tip-cap blocks
Measured on a NACA 0012 wing ([ADR 0057](/adr/0057-right-handed-ch-outflow)):
- three cap blocks are mirrored, with every cell's scaled Jacobian
  negative: one tip's `lower` and `wake_base`, the other's `upper`;
- the LE blocks are degenerate (`jlo`) or twisted (`jhi`);
- the four `te_*` blocks have a single node along i, so no cells.

A redesign of the cap is proposed, not done. Check a cap with
`shore check` / `quality.summary` before using it.
:::

The wrap leaves the span ends open. A finite wing is closed at a tip by
an H-style set of cap blocks that fills the airfoil section in the tip
plane:

| Block | Fills |
|---|---|
| `cap_le` | around the LE, between the two LE cuts (`le_frac`, default 0.05 chord) |
| `cap_upper`, `cap_lower` | the chord band between the LE cut and the TE cut (`te_frac`, default 0.95), above and below the chord |
| `cap_te_upper`, `cap_te_lower` | the TE region, between the TE cut and the TE base |
| `cap_wake_base` | (conformal path, or blunt input) the cap of the `wake_base` block |

Each cap block has `k_lo` `WALL` (the tip wall) and `j_hi` `WALL` (the
airfoil surface carried into the cap); the cap blocks are connected to
each other by `SHARED` faces. Caps are volumetric at build time: marching
skips them.

<Figure name="cgrid-tip-cap" caption="Tip cap blocks closing the wing section at the tip, over the C-grid wall." />

**Why the TE is rounded.** At a sharp TE the first and last nodes of the
section coincide, and an H-fill block with a corner there degenerates to
a singular axis. Both cap paths therefore split the TE node into two,
`te_thickness` apart across the airfoil (default 0.2 % of chord,
`te_thickness_frac` 0.002).

There are two paths:

| | Conformal: `shore mesh --cap-tips jlo\|jhi\|both` | Overset: `shore cap-tip BODY_GEO` |
|---|---|---|
| Body TE | rounded over the whole span; the C-grid is promoted to 4 blocks | left sharp; only the cap is rounded |
| Blocks | 4 C-grid + 6 cap blocks per tip (16 for `both`) | 5 cap blocks per tip (sharp body) |
| Cap–body interface | the cap's wall nodes equal the body's tip-station wall nodes | not coincident: the solver interpolates across the fringe |
| When | building the wing from scratch | adding caps to an existing C-grid without remeshing it |

On the conformal path the cap and body nodes are equal at the interface
but live in separate arrays; the face graph declares no cap–body seam,
because the two volumes meet along a curve (the section at the tip), not
a face. The overset path couples them through chimera interpolation. The
rounding is an aerodynamic approximation at the TE (0.2 % of chord by
default). Blunt input sections are used as they are. The CLI of the
overset path builds sharp-TE caps only; blunt caps need the Python API
(`CGridTipCap(is_blunt=True, ...)` with the wake-base section).

## Limits

- Each span station must cut one closed section; the TE must run within
  about 30° of the span axis ([C-grid construction](/concepts/c-grid-construction#limits)).
- The span-wise stations are uniform between the body's span limits,
  2 % inside each end.
- A wing meeting a fuselage needs a [collar](/concepts/collars) or a
  [blocking](/concepts/blocking), not a C-grid.

## See also

- [Wing tutorial](/tutorials/wing)
- [C-grid construction](/concepts/c-grid-construction): TE detection, slicing, wake curves
- [CH hybrid](ch-hybrid): the C-grid wrap extended to a far-field frame
- [`shore mesh`](/reference/mesh), [`shore cap-tip`](/reference/cap-tip)
- [`CGridTopology`](/reference/api/shore.mesh.topology#CGridTopology), [`CGridTipCap`](/reference/api/shore.mesh.topology#CGridTipCap), [`shore.topology.cgrid_cap`](/reference/api/shore.topology.cgrid_cap)
