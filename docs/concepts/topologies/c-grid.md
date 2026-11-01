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

The wrap leaves the span ends open. A cap closes a finite wing's tip
([ADR 0059](/adr/0059-tip-cap-on-the-tip-face)): an overset grid on the flat
tip face, made of five blocks.

| Block | Covers |
|---|---|
| `core` | the inside of the section, inset by `core_frac` (0.5) of the local half-thickness from both sides |
| `upper`, `lower` | from the core to the wall on each side (`upper` is the C-grid's i = 0 side), between the LE cuts and the TE |
| `le` | from the core round the leading edge, between the two LE cuts (`le_frac`, default 0.05 chord) |
| `te` | from the core across the TE base |

Footprint and volume:
- **The footprint is an O-H layout.** The ring's lines run straight from
  the core to the wall nodes, which are the C-grid's own tip-station nodes.
  The core's corners sit on 45° diagonals from the LE cuts and the TE
  corners, so where three blocks meet the angles are about 90 / 135 / 135°,
  and no corner is a straight angle.
- **The volume grows outwards.** Every block's `k_lo` face is the wall;
  together the five cover the section exactly. k leaves the tip face along
  its normal at the wrap's layer distances. The ring's lines lean out by up
  to `flare_deg` (30°) at the section's edge, so the cap overhangs the side
  wall.
- **Faces:** the blocks are right-handed. Faces between cap blocks are
  `SHARED` (16 per tip); the other outer faces are `FREE` fringes, so the
  cap is an overset component of the wrap.

Caps are built, not marched, and a cap with a folded cell or an empty block
is refused ([ADR 0058](/adr/0058-no-silent-tip-cap-folds)).

<Figure name="cgrid-tip-cap" caption="The cap's wall closing the wing section at the tip, inside the C-grid's end face." />

**Why the TE is rounded.** The `te` block spans the TE base, so a sharp TE
gets one: its node is split in two, `te_thickness` apart across the
airfoil (default 0.2 % of the chord, `te_thickness_frac` 0.002).

There are two ways to build a cap:

| | `shore mesh --cap-tips jlo\|jhi\|both` | `shore cap-tip BODY_GEO` |
|---|---|---|
| Body TE | rounded over the whole span; the C-grid is promoted to 4 blocks | left sharp; only the cap is rounded |
| Stations | the first and last sit on the tips | as built: 2 % of the span inside the tips |
| Blocks | 4 C-grid + 5 cap blocks per tip | 5 cap blocks |
| The cap's edge | the body's tip-station wall nodes, bit for bit | those nodes moved out onto the tip face (`--tip-inset`) |
| Layers | the wrap's schedule (`ds`, `growth`) | the body's k-lines at the tip station |

The cap and the wrap meet along the tip edge, not across a face, so the
face graph has no cap–body seam: the solver interpolates between them
through the fringes. With `shore cap-tip` the strip of side wall between
the tip station and the tip face has no grid of its own.

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
