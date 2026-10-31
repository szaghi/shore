---
title: shore.mesh.topology
outline: [2, 3]
---

# `shore.mesh.topology`

Topology ABC and concrete topologies (cubed-sphere, O-grid, C-grid, OH-/CH-hybrid).

## `Topology` <Badge type="tip" text="class" /> {#Topology}

```python
class Topology(ABC)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/topology.py#L28">source</a>

Abstract topology builder.

Receives an already-projected k=0 surface layer and returns a list of
initialised blocks with shared-face views wired.  Marching is not
performed here.

### `Topology.build` <Badge type="info" text="method" /> {#Topology-build}

```python
build(surface: np.ndarray, nk: int, **kwargs) -> list[HexBlock]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/topology.py#L37">source</a>

Build blocks from the k=0 surface layer.

**Parameters**

- **`surface`** — Projected k=0 surface, shape `(ni, nj, 3)`.
- **`nk`** — Total number of wall-normal layers (including k=0).

**Returns**

- `list[HexBlock]` — Initialised blocks with k=0 filled and shared views wired. Volume layers k=1..nk-1 are allocated but not filled.

## `CubedSphereTopology` <Badge type="tip" text="class" /> {#CubedSphereTopology}

```python
class CubedSphereTopology(Topology)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/topology.py#L248">source</a>

6-block cubed-sphere topology (Sadourny 1972 / Ronchi-Iacono-Paolucci 1996).

A six-block decomposition of the body surface: 4 equatorial sub-blocks
plus 2 polar caps.  Each cap k=0 face is built by gnomonic flat-square
projection (a tangent-plane square radially projected onto the body),
matching the construction of cubed-sphere atmospheric models.

Layout: `[sub0, sub1, sub2, sub3, cap_north, cap_south]`

Construction order (gnomonic algorithm):

1. Build the south cap k=0 by gnomonic flat-square projection (4 corners
   at colatitude `theta_cap`, axis-aligned with sub-block junctions).
2. Build the north cap k=0 the same way.
3. Build the equator k=0 by slerp-meridian interpolation between the two
   cap seams, ray-cast onto the STL.
4. Wire the blocks: sub_q ↔ sub_(q+1) periodic SHARED face seams, and
   sub_q ↔ cap DIRICHLET face links at the seam.

Shared faces use numpy views into a contiguous equator-ring allocation.

### `CubedSphereTopology.build` <Badge type="info" text="method" /> {#CubedSphereTopology-build}

```python
build(
    surface: np.ndarray | None,
    nk: int,
    ni: int | None = None,
    nj: int | None = None,
    spacing_k: Spacing1D | None = None,
    mesh: object | None = None,
    projector: object | None = None,
    theta_cap: float = CUBED_SPHERE_THETA_CAP,
    i_spacing: str = 'uniform',
    i_beta: float = 3.0,
    j_spacing: str = 'uniform',
    j_beta: float = 3.0,
    phi_offset: float = 0.0,
    **kwargs,
) -> list[HexBlock]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/topology.py#L271">source</a>

Build 6 blocks from the body's STL via gnomonic cap projection.

**Parameters**

- **`surface`** — Ignored — preserved as the first positional for backwards-compat with the old CubedSphereTopology.build(surface, nk, ...) signature.  The new construction projects directly through the cached Projector, so a pre-projected surface array is not used.  When `ni` or `nj` are missing, they are derived from `surface.shape`.
- **`nk`** — Total number of wall-normal layers (including k=0).
- **`ni, nj`** — Equator resolution.  `nj` must be a multiple of 4.  If omitted, both are derived from `surface.shape`.
- **`spacing_k`** — Wall-normal spacing law shared by all blocks.  Defaults to `Spacing1D(law='geometric', ds=1e-3)`.
- **`mesh`** — `trimesh.Trimesh` instance.  Legacy alternative to *projector*; creates a temporary Projector using `mesh.centroid`.
- **`projector`** — Pre-built `Projector` instance (shares cached BVH with the projection step).  Preferred over *mesh* — avoids a second BVH build.  Required for the gnomonic algorithm; if `None` and no `mesh` is supplied, a ValueError is raised.
- **`theta_cap`** — Polar exclusion angle in radians.  The 4 cap corners sit at colatitude `theta_cap` on a sphere body; cap edges sweep poleward at midpoints (the natural gnomonic-square geometry).
- **`i_spacing`** — Latitude spacing law of the equator (`'uniform'` | `'tanh'` | `'tanh2'`).  Affects only the equator's meridional distribution.
- **`i_beta`** — Clustering strength for tanh/tanh2 i-spacing.
- **`j_spacing`** — Tangent-plane spacing law of the cap edges (`'uniform'` or `'tanh2'`).  Inherited by the equator seam since the equator's seam rows are the cap edges bit-exactly.
- **`j_beta`** — Clustering strength for tanh2 j-spacing.

**Returns**

- `list[HexBlock]` — `[sub0, sub1, sub2, sub3, cap_north, cap_south]`

## `OGridTopology` <Badge type="tip" text="class" /> {#OGridTopology}

```python
class OGridTopology(Topology)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/topology.py#L377">source</a>

Single-block lat-lon O-grid topology.

A single hexahedral block with shape `(ni, nj, nk, 3)`, j-periodic
(the longitude direction wraps), and polar singularities at `i = 0`
(south pole at colatitude `theta_cap`) and `i = ni - 1` (north
pole at colatitude `pi - theta_cap`).  This is the simpler and
older sphere covering — useful as a Chimera background or in
downstream tools that expect a single periodic-j block — and it
coexists with :class:`CubedSphereTopology` as a peer, not a
component.

The k=0 surface is built by ray-cast from the body's anchor through
a structured lat-lon sphere template.  The same `i_spacing` /
`j_spacing` flags as the cubed sphere apply, with the same
semantics: `i` controls the meridional distribution, `j` the
longitude.

Per-edge spacing
----------------
Per-edge **k-spacing** works on O-grid via the same identity-check
mechanism as :class:`CubedSphereTopology`: replace the four k-edges'
`spacing` field with a shared `Spacing1D` instance, and
:meth:`shore.mesh.Mesh.march` reads it instead of its `ds`/`growth`
arguments.  See :func:`shore.mesh.mesh._resolve_k_steps` for the rules.

Per-edge **i-spacing** has no degree of freedom on O-grid: the j
direction is periodic, so the four parallel i-edges
(`j_lo_k_lo` / `j_hi_k_lo` / `j_lo_k_hi` / `j_hi_k_hi`) lie
at *adjacent* azimuths under the j-wrap rather than at distinct
transverse corners.  Use the per-axis `i_spacing` / `i_beta`
arguments instead.

Per-edge **j-spacing** is meaningless on a periodic axis (the
"first" and "last" cells of any j-edge are the same cell).

### `OGridTopology.build` <Badge type="info" text="method" /> {#OGridTopology-build}

```python
build(
    surface: np.ndarray | None,
    nk: int,
    ni: int | None = None,
    nj: int | None = None,
    spacing_k: Spacing1D | None = None,
    mesh: object | None = None,
    projector: object | None = None,
    theta_cap: float = 0.05 * np.pi,
    i_spacing: str = 'uniform',
    i_beta: float = 3.0,
    j_spacing: str = 'uniform',
    j_beta: float = 3.0,
    **kwargs,
) -> list[HexBlock]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/topology.py#L414">source</a>

Build a single-block O-grid topology from the body's STL.

**Parameters**

- **`surface`** — Optional pre-projected k=0 surface, shape `(ni, nj, 3)`. When supplied, *ni* and *nj* are derived from the array and no ray-cast is performed.  When `None`, the surface is built by ray-cast through *projector* (or *mesh*).
- **`nk`** — Total number of wall-normal layers (including k=0).
- **`ni, nj`** — Latitude (meridional) and longitude (azimuthal) resolution. Ignored when *surface* is supplied.
- **`spacing_k`** — Wall-normal spacing law.  Defaults to `Spacing1D(law='geometric', ds=1e-3)`.
- **`mesh`** — `trimesh.Trimesh` instance.  Legacy alternative to *projector*.
- **`projector`** — Pre-built `Projector` instance.  Required when *surface* is None and *mesh* is not supplied.
- **`theta_cap`** — Polar exclusion angle in radians.  `i = 0` sits at colatitude *theta_cap*, `i = ni - 1` at `pi - theta_cap`.
- **`i_spacing, j_spacing`** — Surface spacing laws.  `i_spacing` controls the meridional distribution; `j_spacing` controls the longitude distribution (only `"uniform"` and `"tanh2"` are valid for j when the j direction is periodic).
- **`i_beta, j_beta`** — Clustering strengths for tanh / tanh2.

**Returns**

- `list[HexBlock]` — Single-element list containing the O-grid block.

## `CGridTopology` <Badge type="tip" text="class" /> {#CGridTopology}

```python
class CGridTopology(Topology)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/topology.py#L520">source</a>

C-grid topology for sharp-TE bodies (airfoils, hydrofoils).

Two block layouts are supported, selected via the `blocks`
argument to :meth:`build`:

**3-block (default, Xall-correct).** Three hex blocks `upper`,
`body`, `lower`:

- `upper` `(ni_wake, n_stations, nk)` — sweeps i from the
  upper-wake far-field tip toward the upper-TE.
- `body` `(ni_body, n_stations, nk)` — sweeps i from the
  upper-TE around the LE to the lower-TE.
- `lower` `(ni_wake, n_stations, nk)` — sweeps i from the
  lower-TE toward the lower-wake far-field tip.

Three SHARED seams:

- `upper.i_hi  ↔ body.i_lo`  — upper-TE column (same axis,
  identity AxisMap on along-seam `(j, k)`).
- `body.i_hi   ↔ lower.i_lo` — lower-TE column (same).
- `upper.k_lo  ↔ lower.k_lo` — wake cut.  Same-axis (k↔k) but
  `i` reverses on the partner side: upper's `i = 0` is the
  far-field tip while lower's `i = 0` is the TE.  AxisMap on
  `upper.k_lo`: `from_pairs(("i", True), ("j", False))`.

For sharp TE the wake-cut nodes coincide bit-exactly in space and
SHARED is a true memory-aliased view; for blunt TE the upper and
lower wake k=0 strips are parallel-but-separated and the 3-block
topology is rejected (use the 4-block variant instead — see the
`blocks=4` arm of :meth:`build`).

**1-block (legacy, retained for opt-in compatibility).** A single
hex block of shape `(ni_total, n_stations, nk)` whose i-axis
sweeps ``upper-wake-far-tip → upper-TE → suction → LE → pressure
→ lower-TE → lower-wake-far-tip`` as one continuous parametric
line.  The wake cut is encoded as an i-direction overlap: the
upper wake column at `i = 0..ni_wake - 1` and the lower wake
column at `i = ni_wake + ni_body - 1..ni_total - 1` are
geometrically coincident at k=0 but topologically distinct in i.
Xall cannot apply BCs cleanly on the wake-cut faces in this
layout (the j-counter "flips direction" along the overlapping
edge); use `blocks=3` instead for new work.

Conventions, both layouts:

- `span_axis` defaults to `"z"`.
- Wake far-field tips (`upper.i_lo`, `lower.i_hi` for the
  3-block layout; `i_lo` and `i_hi` for the 1-block layout)
  are `FREE` by default — typically Chimera fringes (BC code
  `40`) or extrapolation (`-9`) in a downstream `cc.par`.
- `j_lo` and `j_hi` on every block are wing-tip `FREE`.
- `k_lo` is `WALL` on the body block; on upper / lower it
  is `FREE` (the wake-cut faces become SHARED in 3-block
  layout, so the FREE default only ever surfaces in the
  legacy 1-block path).
- `k_hi` is `FREE` (outer / Chimera fringe).

**See Also**

:func:`shore.surface.cgrid.build_cgrid_k0` :
    Builds the k=0 layer fed to this topology's :meth:`build`.

### `CGridTopology.build` <Badge type="info" text="method" /> {#CGridTopology-build}

```python
build(
    surface: np.ndarray | None,
    nk: int,
    ni_body: int | None = None,
    ni_wake: int | None = None,
    n_stations: int | None = None,
    wake_length: float | None = None,
    span_axis: str = 'z',
    freestream_world: np.ndarray | None = None,
    te_dihedral_deg: float = 30.0,
    wake_growth: float = 1.15,
    spacing_k: Spacing1D | None = None,
    mesh: object | None = None,
    blocks: int = 3,
    nk_w: int = 3,
    cap_tips: str = 'none',
    cap_te_thickness: float | None = None,
    cap_te_thickness_frac: float = 0.002,
    cap_le_frac: float = 0.05,
    cap_te_frac: float = 0.95,
    cap_tip_inset: float | None = None,
    cap_n_cap_j: int = 12,
    cap_j_spacing: str = 'tanh',
    cap_j_beta: float = 4.0,
    cap_nk: int | None = None,
    **kwargs,
) -> list[HexBlock]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/topology.py#L584">source</a>

Build a C-grid block list from a body STL.

**Parameters**

- **`surface`** — Optional pre-built `(ni_total, n_stations, 3)` k=0 layer from :func:`shore.surface.cgrid.build_cgrid_k0`.  When supplied, *ni_body* / *ni_wake* / *n_stations* / *wake_length* are ignored.  When `None`, the k=0 layer is built from *mesh* using the supplied parameters. For `blocks=3` the surface is split internally into per- block (upper / body / lower) k=0 layers via :meth:`shore.surface.cgrid.CGridK0Result.per_block_layers`; pre-built surfaces work equivalently as long as the `ni_total = 2 * (ni_wake - 1) + ni_body` invariant holds.
- **`nk`** — Total number of wall-normal layers (including k=0).
- **`ni_body`** — Per-station body wrap node count.  Required when *surface* is None.
- **`ni_wake`** — Per-branch wake node count (shared TE node included). Required when *surface* is None.
- **`n_stations`** — Number of span-axis stations.  Required when *surface* is None.  Becomes each block's `nj`.
- **`wake_length`** — Length of each wake branch from the TE.  Required when *surface* is None.
- **`span_axis`** — `"x"` / `"y"` / `"z"`.  Default `"z"`.
- **`freestream_world`** — Length-3 vector giving the freestream direction in world coordinates.  `None` = default to chord-axis +1.
- **`te_dihedral_deg, wake_growth`** — Forwarded to :func:`build_cgrid_k0`.
- **`spacing_k`** — Wall-normal spacing law.  Defaults to `Spacing1D(law='geometric', ds=1e-3)`.
- **`mesh`** — `trimesh.Trimesh` instance.  Required when *surface* is None.
- **`blocks`** — `3` (default, Xall-correct sharp-TE), `4` (Xall-correct blunt-TE: adds a wake-base block W bridging the TE thickness gap), or `1` (legacy single-block).  See class docstring for the trade-offs.
- **`nk_w`** — Wake-base block W's resolution across the TE thickness gap.  Used only for `blocks=4`.  Default `3` is adequate for the thin gaps (≤ 1% chord) typical of blunt-TE airfoils.
- **`cap_tips`** — `"none"` (default), `"jlo"`, `"jhi"`, or `"both"`. When set to anything other than `"none"`, the C-grid mesh is augmented with **conformal** wing-tip cap blocks at the requested spanwise extreme(s).  The cap blocks close the wing-tip airfoil cross-section as a wall, with their airfoil-wall nodes coincident bit-exactly with the body's k=0 wall at the spanwise tip. For sharp-TE bodies the cap-tips path **auto-rounds the body's TE over its full span** by `cap_te_thickness` (or `cap_te_thickness_frac * chord` if not given), so the body wrap and the cap meshes match without any chimera-fringe handover.  `blocks=3` is automatically promoted to `blocks=4` in that case (a stderr note is emitted).  `blocks=1` (legacy) + `cap_tips != "none"` is rejected. For already-blunt bodies the existing TE thickness is used and no rounding is applied.
- **`cap_te_thickness`** — Absolute TE thickness for body-rounding when `cap_tips != "none"` and the body is sharp.  Mutually exclusive with `cap_te_thickness_frac`.
- **`cap_te_thickness_frac`** — Same, as a fraction of chord.  Default `0.002` (0.2%).
- **`cap_le_frac, cap_te_frac`** — Chord fractions for the cap's LE-cut and TE-cut.  Default `0.05` and `0.95`.
- **`cap_tip_inset`** — Cap spanwise wall thickness.  Default = `1e-3 * chord`.
- **`cap_n_cap_j`** — Cap j-axis (airfoil-thickness) node count.  Default `12`.
- **`cap_j_spacing, cap_j_beta`** — Cap j-axis spacing law (`"uniform"` / `"tanh"` / `"tanh2"`) and clustering strength.  Default `"tanh"`, beta `4.0`.
- **`cap_nk`** — Cap k-axis (spanwise) node count.  Default = `nk` (the body's wall-normal layer count).

**Returns**

- `list[HexBlock]` — For `blocks=3`: `[upper, body, lower]`. For `blocks=4`: `[upper, body, lower, wake_base]`. For `blocks=1`: `[cgrid]`. When `cap_tips != "none"`: the C-grid blocks are followed by 5 cap blocks per requested tip (`cap_jlo_le`, `cap_jlo_upper`, `cap_jlo_lower`, `cap_jlo_te_upper`, `cap_jlo_te_lower` for jlo; mirror for jhi), plus 1 extra `cap_jX_wake_base` block per tip when the body wrap is blunt (`blocks=4`).

## `CGridTipCap` <Badge type="tip" text="class" /> {#CGridTipCap}

```python
CGridTipCap(
    le_cut_upper_idx: int,
    le_cut_lower_idx: int,
    span_dir: tuple[float, float, float],
    tip_inset: float,
    n_cap_j: int = 12,
    cap_j_spacing: str = 'tanh',
    cap_j_beta: float = 4.0,
    is_blunt: bool = False,
    wake_base_cross_section: np.ndarray | None = None,
    te_cut_upper_idx: int | None = None,
    te_cut_lower_idx: int | None = None,
    auto_round_te: bool = True,
    te_thickness: float | None = None,
    te_thickness_frac: float | None = 0.002,
)
```
Bases: `Topology`<br>
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/topology.py#L1418">source</a>

C-grid wing-tip cap topology — 3 (sharp) or 4 (blunt) overset blocks.

Builds the H-style tip closure for one tip of a C-grid wing as an
independent chimera component.  The cap is composed of 3 (sharp TE)
or 4 (blunt TE) hexahedral sub-blocks that fill the airfoil-interior
cross-section at the wing-tip plane in an "H3 T-fill":

- `cap_le`    — wraps the leading edge from the upper-LE-cut to the
                   lower-LE-cut.
- `cap_upper` — chord-band on the suction side, between LE-cut and
                   upper-TE-cut, above the chord cut.
- `cap_lower` — mirror of `cap_upper` on the pressure side.
- `cap_wake_base` (blunt only) — closes the wake_base block's
                   spanwise face downstream of the TE base segment.

All sub-blocks share the same axis convention: i = chordwise along
the airfoil profile (or LE wrap), j = airfoil-thickness direction
(chord cut → wall), k = spanwise.  k=0 is the outboard wing-tip
wall (offset from the cross-section by `tip_inset * span_dir`);
k=nk-1 is the inboard chimera-overlap plane (coincident with the
body wrap's spanwise-extreme cross-section).

Default face BCs:

- **k_lo**: WALL (the actual wing-tip wall closure).
- **k_hi**: FREE (chimera fringe receiving from the body wrap).
- **j_hi**: WALL (the airfoil suction / pressure / LE / TE-base
  surface, mirrored into the cap).
- **j_lo**: SHARED between cap_upper and cap_lower (chord cut),
  otherwise FREE.
- **i_lo / i_hi**: SHARED at the LE-cuts between cap_le and
  cap_upper / cap_lower; FREE otherwise (TE-cut chord-perp
  sliver).

The cap is **fully volumetric at build time** — calling
:meth:`shore.mesh.Mesh.march` on a CGridTipCap mesh is a no-op.
The cap is meant to be assembled against a C-grid body wrap as a
chimera component via `shore cc-par`.

**Parameters**

- **`le_cut_upper_idx, le_cut_lower_idx`** — i-indices in the `surface` cross-section polyline at the LE-cut on the suction and pressure sides.  Must satisfy `0 < le_cut_upper_idx < le_cut_lower_idx < ni - 1`.
- **`span_dir`** — Length-3 unit vector pointing **outward** along the span axis at this tip.
- **`tip_inset`** — Spanwise thickness of the cap slab.  &gt; 0.
- **`n_cap_j`** — Node count along the cap's j-axis (airfoil-thickness).  &gt;= 2.
- **`cap_j_spacing`** — `"uniform"` / `"tanh"` / `"tanh2"`.  `"tanh"` (default) clusters near the airfoil wall.
- **`cap_j_beta`** — Clustering strength.
- **`is_blunt`** — Toggles the 4th sub-block `cap_wake_base`.
- **`wake_base_cross_section`** — Required when `is_blunt=True`.  `(ni_wake, nk_w, 3)` — the wake_base block's spanwise-extreme face from the C-grid.
- `Future input modes (not implemented yet)`
- `----------------------------------------`
- `The current ``build()`` accepts a 1D cross-section polyline as`
- ```surface``.  In future, ``build()`` will additionally accept a 2D`
- `sliced STL/STEP cross-section ``(ni, nj_thin, 3)`` so the cap can`
- `be driven directly by a CAD/STL surface input — the same way the`
- `body wrap is.  The signature remains the same.`

### `CGridTipCap.build` <Badge type="info" text="method" /> {#CGridTipCap-build}

```python
build(
    surface: np.ndarray,
    nk: int,
    spacing_k: Spacing1D | None = None,
    **kwargs,
) -> list[HexBlock]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/topology.py#L1561">source</a>

Build the tip cap blocks from a cross-section polyline.

**Parameters**

- **`surface`** — `(ni_body, 3)` cross-section polyline at the wing tip. Sweeps upper-TE-cut → suction → LE → pressure → lower-TE-cut.
- **`nk`** — Spanwise node count for the cap (`nk_cap`).  Must be &gt;= 2.
- **`spacing_k`** — Reserved for future per-block k-spacing overrides. Currently the cap's k-axis is uniformly spaced from the outboard wall to the cross-section plane.

**Returns**

- `list[HexBlock]` — `[cap_le, cap_upper, cap_lower]` for sharp TE, or `[cap_le, cap_upper, cap_lower, cap_wake_base]` for blunt.

## `OHTopology` <Badge type="tip" text="class" /> {#OHTopology}

```python
OHTopology(
    box_vmin: tuple[float, float, float],
    box_vmax: tuple[float, float, float],
    nk_channel: int,
    k_law_channel: str = 'uniform',
    k_beta_channel: float = 3.0,
    c1_seam: bool = True,
)
```
Bases: `Topology`<br>
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/topology.py#L1819">source</a>

Style-2 conforming hybrid: cubed-sphere wrap + 6-block H-channel.

A 12-block decomposition: 6 wrap blocks (cubed-sphere, marched outward
from the STL surface to `nk_wrap - 1`) and 6 channel blocks (frustums
bridging the wrap's outer envelope to the user-supplied bounding box).
Wrap and channel are SHARED-by-view at the seam — they share one
contiguous numpy buffer per pair, so `wrap.k_hi` and `chan.k_lo` are
bit-identical.

The wrap is built with `phi_offset = pi/4` so the four equatorial
sub-block *centres* (rather than corners) face the world cardinal axes:
sub0→+y, sub1→-x, sub2→-y, sub3→+x.  Combined with the two cap blocks
(cap_south→+z, cap_north→-z), each of the 6 wrap blocks faces exactly
one of the box's 6 faces, allowing a 1:1 wrap-block-to-channel-block
mapping.

Block layout returned by :meth:`build`:

    [sub0, sub1, sub2, sub3, cap_north, cap_south,
     chan_sub0, chan_sub1, chan_sub2, chan_sub3, chan_north, chan_south]

Channel construction is deferred until after marching: the wrap's k_hi
envelope is unknown at build time, so the channel's k=0 row (which IS
the wrap's k_hi via view aliasing) gets filled when the wrap finishes
marching.  :meth:`shore.mesh.Mesh.march` calls
:meth:`fill_channel_blocks` after the wrap march completes; users who
drive marching manually must call :meth:`fill_channel_blocks` before
using the channel block geometry.

Default face BCs:

- Wrap k_lo = WALL (STL surface).
- Wrap k_hi = SHARED (with channel k_lo).
- Channel k_lo = SHARED (with wrap k_hi).
- Channel k_hi = FREE (far-field box face).
- Wrap and channel lateral seams (j_hi↔j_lo, i_hi↔cap, i_lo↔cap)
  mirror the cubed-sphere topology on both sides.

### `OHTopology.build` <Badge type="info" text="method" /> {#OHTopology-build}

```python
build(
    surface: np.ndarray | None,
    nk: int,
    ni: int | None = None,
    nj: int | None = None,
    spacing_k: Spacing1D | None = None,
    mesh: object | None = None,
    projector: object | None = None,
    theta_cap: float = 0.05 * np.pi,
    i_spacing: str = 'uniform',
    i_beta: float = 3.0,
    j_spacing: str = 'uniform',
    j_beta: float = 3.0,
    **kwargs,
) -> list[HexBlock]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/topology.py#L1894">source</a>

Build the 12-block OH topology.

`nk` is the wrap depth (number of layers including k=0).  The
channel depth comes from the OHTopology constructor.  Returns
``[sub0..sub3, cap_north, cap_south, chan_sub0..chan_sub3,
chan_north, chan_south]``.

### `OHTopology.fill_channel_blocks` <Badge type="info" text="method" /> {#OHTopology-fill_channel_blocks}

```python
fill_channel_blocks(blocks: list[HexBlock]) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/topology.py#L2039">source</a>

Build channel block interiors after the wrap has been marched.

Walks the 6 wrap-channel pairs, computes each channel's outer face
on the user-supplied bounding box, and fills the channel's
`k=1..nk_channel-1` interior by ruled-surface frustum interpolation.

The wrap-channel seam (channel k=0 = wrap k_hi) is already populated
because the two faces share memory.  Only k=1..nk_channel-1 needs
filling.

When `self.c1_seam` is true (default), the channel's first step
at every (i, j) is pinned to the wrap's last marching step at the
same (i, j) — bit-exact C1 spacing match across the seam.  Pass
`c1_seam=False` to the constructor to disable and use the raw
spacing law instead (cell-size jump at the seam).

## `CHTopology` <Badge type="tip" text="class" /> {#CHTopology}

```python
CHTopology(
    *,
    chord_min: float,
    chord_max: float,
    normal_min: float,
    normal_max: float,
    nk_channel: int,
    freestream_world: np.ndarray | None = None,
    span_axis: str = 'z',
    k_law_channel: str = 'uniform',
    k_beta_channel: float = 3.0,
    c1_seam: bool = True,
    i_split_body_upper: int | None = None,
    i_split_body_lower: int | None = None,
    ni_outflow: int = 8,
)
```
Bases: `Topology`<br>
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/topology.py#L2179">source</a>

Hybrid C-grid wrap + H-channel far-field for sharp- or blunt-TE bodies.

Mirrors :class:`OHTopology` but with a C-grid wrap (3 or 4 blocks)
instead of a cubed-sphere wrap.  The user-supplied bounding "box"
is a rectangle in the chord/normal plane (extruded across the wrap's
span extents), defined by `chord_min`, `chord_max`, `normal_min`,
`normal_max` in world coords.  `freestream_world` and
`span_axis` select which world axes carry chord and span.

**Block layout** returned by :meth:`build`:

CH-3 (sharp TE) — 12 blocks:

```
[upper, body, lower,                              # wrap
 chan_upper, chan_body_top, chan_body_front,
 chan_body_bot, chan_lower,                       # k_hi-side channels
 chan_outflow_top_upper, chan_outflow_wake_upper,
 chan_outflow_wake_lower, chan_outflow_top_lower] # outflow (4 sub-blocks)

```
CH-4 (blunt TE) — 14 blocks:

```
[upper, body, lower, wake_base,                   # wrap
 chan_upper, chan_body_top, chan_body_front,
 chan_body_bot, chan_lower,                       # k_hi-side channels
 chan_outflow_top_upper, chan_outflow_wake_upper,
 chan_outflow_wake_lower, chan_outflow_top_lower,
 chan_outflow_wake_base]                          # outflow (5 sub-blocks)

```
The outflow region is partitioned per side into a *top* sub-block
(conformal SHARED to `chan_upper.i_lo` / `chan_lower.i_hi` —
same `nk_channel`) and a *wake* sub-block (conformal SHARED to
the wrap's `upper.i_lo` / `lower.i_hi` — same `nk_wrap`).
Adjacent sub-blocks meet at envelope-far-tip seams (SHARED) and
at the wake-cut seam (SHARED across upper and lower wake halves;
via `chan_outflow_wake_base` for blunt TE).

The 5 `k_hi`-side channels (`chan_upper`, `chan_body_top`,
`chan_body_front`, `chan_body_bot`, `chan_lower`) are
SHARED-by-view with their wrap partners — wrap.k_hi and chan.k_lo
bit-identical via numpy view aliasing on a single `nodes_all`
allocation per wrap segment.  The body wrap's k_hi is split into
3 i-segments at user-tunable indices `i_split_body_upper` and
`i_split_body_lower` (defaults: `ni_body // 3` and
`2 * ni_body // 3`), each sub-channel sharing only its own
chunk of body.k_hi via numpy view aliasing.  Adjacent
chan_body_* sub-channels share their TE-direction seam by
overlapping on the split column (identity AxisMap).

The 2 (or 3) outflow channels are NOT view-aliased: their k=0
inner face has shape `(nk_wrap, nj, 3)` derived from the wrap's
`i_lo` / `i_hi` faces (a cross-axis remap from the wrap
convention), so :meth:`fill_channel_blocks` copies the face data
explicitly each time.  The SHARED face-graph entries are still
emitted so cc.par writers see the connection.

**Default face BCs**:

- Wrap k_lo = WALL on body, SHARED on upper/lower (wake-cut /
  wake-base).
- Wrap k_hi = SHARED with the corresponding chan.k_lo.
- Wrap i_lo / i_hi at the wake far-tips = SHARED with chan_outflow_*.
- Channel k_lo = SHARED.  Channel k_hi = FREE (frame edge).
- Channel j_lo / j_hi = FREE (wing tips, Phase 2 will optionally
  add spanwise tip channels).

Channels do NOT march; :meth:`shore.mesh.Mesh.march` runs only
the C-grid wrap marcher and then calls :meth:`fill_channel_blocks`.

**Parameters**

- **`freestream_world`** — Length-3 world-frame freestream vector.  `None` defaults to the +chord_axis (chord_axis = first non-span world axis).
- **`span_axis`** — `"x"`, `"y"`, or `"z"`.  Default `"z"`.
- **`chord_min, chord_max, normal_min, normal_max`** — Frame extents in world coordinates.  Must contain the marched wrap envelope (validated in :meth:`fill_channel_blocks`).
- **`nk_channel`** — Total channel-block layers (incl. shared-seam slot at k=0). Must be &gt;= 2.
- **`k_law_channel`** — Spacing law along channel k: `'uniform'` / `'tanh'` / `'tanh2'`.  Default `'uniform'`.
- **`k_beta_channel`** — Clustering strength for tanh / tanh2.
- **`c1_seam`** — When True (default), the channel's first step at every (i, j) is pinned to the wrap's last marching step at the same node — bit-exact C1 spacing match across the wrap-channel seam.  Only applies to the 5 `k_hi`-side channels (the outflow channels bridge wrap *side* faces and have no wrap "last step" along their own k-axis to inherit from).
- **`i_split_body_upper, i_split_body_lower`** — Body-block i-indices at which `body.k_hi` is split into the 3 chan_body_* sub-channels.  `None` (default) = `ni_body // 3` and `2 * ni_body // 3`.

### `CHTopology.is_blunt` <Badge type="info" text="property" /> {#CHTopology-is_blunt}

```python
is_blunt: bool
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/topology.py#L2333">source</a>

True when build() detected a blunt TE (CH-4).

### `CHTopology.build` <Badge type="info" text="method" /> {#CHTopology-build}

```python
build(
    surface: np.ndarray | None,
    nk: int,
    *,
    ni_body: int | None = None,
    ni_wake: int | None = None,
    n_stations: int | None = None,
    wake_length: float | None = None,
    wake_growth: float = 1.15,
    te_dihedral_deg: float = 30.0,
    spacing_k: Spacing1D | None = None,
    mesh: object | None = None,
    nk_w: int = 3,
    **kwargs,
) -> list[HexBlock]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/topology.py#L2337">source</a>

Build the CH topology block list.

`nk` is the wrap depth.  The channel depth is taken from
`CHTopology.nk_channel` set at construction.

### `CHTopology.fill_channel_blocks` <Badge type="info" text="method" /> {#CHTopology-fill_channel_blocks}

```python
fill_channel_blocks(blocks: list[HexBlock]) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/topology.py#L3026">source</a>

Populate every channel block's interior after wrap marching.

Walks the channel blocks in order, computes each channel's
outer quad on the user-specified frame, and fills the
`k=1..nk_channel-1` interior by ruled-surface frustum
interpolation.  C1 seam pinning (per-(i, j) first step =
wrap's last marching step) is applied to the 5 k_hi-side
channels when `self.c1_seam` is True.
