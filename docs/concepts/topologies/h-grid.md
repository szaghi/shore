# H-grid topology

An H-grid fills the inside of a duct: i runs streamwise from the inlet
face to the outlet face, j and k span the cross-section, and the faces on
the duct's surface are its walls. SHORE's H-grids are analytic
primitives without an internal body. A body in a channel is handled the
overset way: an H-grid for the channel plus a body-fitted grid around the
body ([Chimera mesh generation](/concepts/chimera)).

<Figure name="h-grid" caption="H-grid primitives: the rectangular duct, one block, and the circular tube, five blocks in an O-H section." />

## Variants

`shore.primitive` exposes two analytic H-grid primitives:

| Variant | CLI | Cross-section | Block layout | Use case |
|---|---|---|---|---|
| **Rectangular** | `shore primitive h-grid-rect` | rectangle | 1 hex block, no singular axis | prismatic rectangular ducts; channels |
| **Circular** | `shore primitive h-grid-tube` | circle | 5 hex blocks, O-H: a square core and four ring blocks, no singular axis | straight cylindrical ducts; simple SRM chambers |

Ducts of varying or non-circular cross-section are not covered by a
primitive; build them as a [blocking](/concepts/blocking) on the duct
geometry.

## Conventions (both variants)

- Block axes: ``i`` is **streamwise**, ``j`` and ``k`` are the two
  cross-section directions.
- World alignment: ``streamwise_axis`` (``"x"``, ``"y"``, or ``"z"``)
  selects which world axis the duct length runs along.  Default ``x``.
- Default face BCs: streamwise faces (``i_lo``, ``i_hi``) are
  ``FREE`` (inflow / outflow); the faces on the duct's wall are
  ``WALL``.
- On the rectangular duct the ``j`` and ``k`` faces are the four walls; on
  the circular tube the faces between its five blocks are ``SHARED``.

## Rectangular duct

Single hex block of shape ``(ni, nj, nk, 3)``.  ``length`` / ``width``
/ ``height`` set the three Cartesian extents (``length`` along the
streamwise axis, ``width`` along block ``j``, ``height`` along block
``k``).  Tanh / tanh2 spacing on the cross-section axes clusters
nodes near the duct walls — the natural choice for boundary-layer
resolution.

### Face boundary conditions

| Face | BC | Meaning |
|---|---|---|
| `i_lo` | `FREE` | Streamwise inflow |
| `i_hi` | `FREE` | Streamwise outflow |
| `j_lo` / `j_hi` | `WALL` | Two of the four duct walls |
| `k_lo` / `k_hi` | `WALL` | The other two duct walls |

In `cc.par` the streamwise `FREE` faces take the block's
`free_face_bc`; set a positive natural code for them, per block or per face
(3 / 4 for inflow / outflow).

## Circular tube

Five hex blocks in an O-H ("butterfly") cross-section
([ADR 0060](/adr/0060-h-grid-tube-oh)):

- ``core``: a square in the middle, its corners at ``core_frac`` (0.5) of
  the radius, a Cartesian grid of ``n_circ / 4`` cells a side;
- ``ring_0`` .. ``ring_3``: each joins one side of the core to a quarter of
  the wall (``ring_q`` from −45 + 90 q to 45 + 90 q degrees). ``j`` runs
  along the wall, ``k`` from the core (``k = 0``) to the wall
  (``k = n_radial - 1``); the lines are straight, and at the core's
  corners they follow the 45° diagonals.
- ``i`` = streamwise in every block.

Every block is right-handed, whatever ``streamwise_axis``, and there is no
singular axis: the worst cells, at the core's corners, score about 0.71
(sin 135°). Before ADR 0060 the tube was one block on a polar grid whose
centreline row collapsed to a point and was marked ``WALL``.

### Face boundary conditions

| Face | BC | Meaning |
|---|---|---|
| `i_lo` (every block) | `FREE` | Streamwise inflow |
| `i_hi` (every block) | `FREE` | Streamwise outflow |
| the rings' `k_hi` | `WALL` | The tube's wall |
| the rest (16 faces) | `SHARED` | Between the core and the rings, and between neighbouring rings |

## Parameters

| Command | Sizes | Nodes | Spacing |
|---|---|---|---|
| `h-grid-rect` | `--length` (1.0), `--width` (0.1), `--height` (0.1) | `--ni` 40, `--nj` 20, `--nk` 20 | `--i/j/k-spacing` `uniform` / `tanh` / `tanh2`, `--i/j/k-beta` 3.0 |
| `h-grid-tube` | `--length` (1.0), `--radius` (0.05), `--core-frac` (0.5) | `--n-axial` 40, `--n-circ` 32 (a multiple of 4), `--n-radial` 8 | `--i-spacing`, `--k-spacing` (core to wall; `tanh` clusters at the wall), betas 3.0 |

Both take `--streamwise-axis` (`x` default). Two-sided `tanh2` on the
cross-section axes clusters toward the walls. In Python they are
`HGridRectMesh.build(...)` and `HGridTubeMesh.build(...)`, with the same
`PrimitiveMesh` contract as the other [primitives](primitives). Commands:
[`shore primitive`](/reference/primitive).

## See also

- [Topology catalogue](./) — pick the right topology for your case.
- [`shore primitive`](/reference/primitive) — full CLI reference for all primitives.
- [Chimera assembly tutorial](/tutorials/chimera-assembly) — a body-fitted grid inside a background.
- [Per-face natural BC codes](/concepts/load-balancing#per-face-natural-bc-codes) — different codes on the inlet and outlet faces.
