# C-grid topology

A C-grid is a structured grid topology that **wraps a body with a
sharp or blunt trailing edge** (airfoil, hydrofoil, wing) and
**extends two wake branches downstream**. Looking at a span station,
the grid lines form the letter C: the outer envelope opens
downstream of the TE, leaving a wake-cut where the upper and lower
surfaces meet.

The C-grid is a **first-class topology** in SHORE: build it with
`shore mesh --topology cgrid` (or `topology = "cgrid"` in a TOML
config). It coexists with the [O-grid](o-grid) and
[cubed-sphere](cubed-sphere) topologies as a peer.

## Block layouts

`shore mesh --topology cgrid` supports three layouts via
`--cgrid-blocks`:

| Layout | Blocks | TE | Use case |
|---|---|---|---|
| **3-block** (default) | upper / body / lower | sharp | The Xall-correct sharp-TE C-grid; SHARED seams at the upper-TE column, lower-TE column, and wake cut |
| **4-block** | upper / body / lower / wake_base | blunt | 3-block + a wake-base block W bridging the TE thickness gap; W's `i_hi` face is the bluff-base WALL |
| **1-block** (legacy) | one concatenated block | sharp only | Single hex block with the wake cut encoded as an i-direction overlap; **Xall cannot apply BCs cleanly** on this overlap (the j-counter "flips direction" along the self-overlapping edge).  Kept for opt-in compatibility with downstream tools that don't have the j-flip issue. |

::: warning Migration from v0.1.0
v0.1.0 shipped the 1-block layout as the default.  v0.2.0 makes
**3-block the default**.  If you have a v0.1.0 workflow that relies
on the single-block layout, pass `--cgrid-blocks 1` to opt in.  See
[the v0.2.0 changelog](/guide/changelog) for the migration note.
:::

## 3-block layout (sharp TE, default)

The 3-block layout splits the C-grid into three hex blocks:

```
                      BLOCK upper (upper wake)
upper-wake-far-tip ──────────────────────── upper-TE
                                                │
                                                │ BLOCK body (suction → LE → pressure)
                                                │
                                            lower-TE
lower-wake-far-tip ──────────────────────── lower-TE
                      BLOCK lower (lower wake)
```

Per-block dimensions:

- **`upper`** `(ni_wake, n_stations, nk, 3)` — i sweeps `upper-wake-far-tip → upper-TE`.
- **`body`** `(ni_body, n_stations, nk, 3)` — i sweeps `upper-TE → suction → LE → pressure → lower-TE`.
- **`lower`** `(ni_wake, n_stations, nk, 3)` — i sweeps `lower-TE → lower-wake-far-tip`.

j is the span direction; k is wall-normal extrusion.

### SHARED seams

Three SHARED face pairs:

| Seam | Faces | AxisMap (this side) | Notes |
|---|---|---|---|
| Upper-TE column | `upper.i_hi ↔ body.i_lo` | identity on `(j, k)` | Same-axis i↔i; node coincidence via numpy view aliasing. |
| Lower-TE column | `body.i_hi ↔ lower.i_lo` | identity on `(j, k)` | Same as above. |
| Wake cut       | `upper.k_lo ↔ lower.k_lo` | `((i, True), (j, False))` | Same-axis k↔k; **i flipped** because upper sweeps far-tip→TE while lower sweeps TE→far-tip.  This i-reversal is the manifestation of what was the "j-counter flips direction" pathology in the 1-block layout's self-overlap. |

For sharp TE the wake-cut nodes coincide bit-exactly in space at
k=0 — the SHARED face is a true memory-aliased view.  Xall sees a
clean face-to-face adjacency on all three seams and applies BCs
correctly.

### Face boundary conditions (defaults)

| Block | i_lo | i_hi | j_lo | j_hi | k_lo | k_hi |
|---|---|---|---|---|---|---|
| `upper` | `FREE` (far-tip) | `SHARED → body.i_lo` | `FREE` | `FREE` | `SHARED → lower.k_lo` | `FREE` |
| `body`  | `SHARED → upper.i_hi` | `SHARED → lower.i_lo` | `FREE` | `FREE` | `WALL` (airfoil surface) | `FREE` |
| `lower` | `SHARED → body.i_hi` | `FREE` (far-tip) | `FREE` | `FREE` | `SHARED → upper.k_lo` | `FREE` |

For a chimera assembly, the `FREE` faces become BC code `40`
(chimera fringe) or `-9` (extrapolation) in a downstream `cc.par`.

## 4-block layout (blunt TE)

For a body with a blunt (finite-thickness) TE, the upper and lower
wake k=0 strips are **parallel but separated by the TE thickness**.
The 3-block layout's wake-cut SHARED face would have a non-zero gap
between U.k_lo and L.k_lo — invalid as SHARED — and the TE base
wall of the airfoil would not be meshed.

The 4-block layout adds a **wake-base block W** that fills the gap:

```
y axis up
   ▲    upper-wake-far-tip ────── upper-TE     (Block upper, k_lo at y = +te/2)
   │                                  │
   │                                  │ Block body
   │                                  │
   │                                  │
   │       far-tip W ────── TE base wall      (Block wake_base)
   │                                  │
   │    lower-wake-far-tip ────── lower-TE    (Block lower, k_lo at y = -te/2)
   │
```

- **`wake_base`** `(ni_wake, n_stations, nk_w, 3)` — i along the
  wake (far-tip → TE), j is span, k spans the TE thickness from
  `y = -te/2` to `y = +te/2`.

The wake_base block is meaningful only when the upper and lower
**wake k=0 strips** sit at `y = ±te/2` over their **full chord-axis
extent** (not just at the TE-side endpoint).  When the body is
already blunt (`naca0012_blunt.stl` and similar), the wake builder
in [`build_cgrid_k0`](/reference/api-cap-equator) anchors each wake
branch at the body's TE node and runs along the freestream, so the
wake strips are at the right offset by construction.  When the body
is **sharp** but the topology is upgraded to 4-block (typically by
the conformal-cap path — see [Tip caps](#tip-caps-closed-finite-wing)),
[`shore.surface.cgrid.round_body_te_in_surface`](/reference/api-cap-equator)
runs after k=0 build to shift the body's TE columns *and* the entire
wake k=0 strips by `±te/2` along the airfoil-thickness direction, so
the wake_base block has uniform `te_thickness` y-spread end to end.

Two extra SHARED seams (replacing the 3-block wake-cut):

| Seam | Faces | AxisMap (W side) |
|---|---|---|
| W upper face ↔ upper k_lo | `wake_base.k_hi ↔ upper.k_lo` | identity on `(i, j)` (W's i and upper's i match) |
| W lower face ↔ lower k_lo | `wake_base.k_lo ↔ lower.k_lo` | `((i, True), (j, False))` (W's i and lower's i are reversed) |

W's `i_hi` face — the cross-section at the TE — is **`WALL`**: this
is the bluff-base wall of the airfoil that the 3-block layout
cannot represent.  The 1-block legacy layout simply ignores it.

W is not marched; its geometry is fully determined at build time by
linear interpolation in k between the lower- and upper-wake k=0
strips.  The marcher's per-block loop skips it.

## Building from the CLI

```bash
# Sharp-TE 3-block (default):
shore mesh naca0012.stl --topology cgrid \
    --ni-body 120 --ni-wake 30 --n-stations 5 --wake-length 15 \
    --nk 30 --ds 5e-4 --growth 1.15 \
    -o naca0012
# → naca0012_upper.geo + naca0012_body.geo + naca0012_lower.geo

# Blunt-TE 4-block:
shore mesh naca0012_blunt.stl --topology cgrid --cgrid-blocks 4 \
    --ni-body 120 --ni-wake 30 --n-stations 5 --wake-length 15 \
    --nk 30 --ds 5e-4 --growth 1.15 \
    -o naca0012_blunt
# → naca0012_blunt_{upper,body,lower,wake_base}.geo

# Conformal tip caps (sharp body + caps): the C-grid is auto-promoted
# to 4 blocks, the body's TE is rounded over its full span, and 6 cap
# blocks per tip are emitted with airfoil-wall nodes bit-coincident
# with the body wrap — no chimera fringe at the cap-body seam.
shore mesh naca0012.stl --topology cgrid \
    --cap-tips both --cap-te-thickness-frac 0.002 \
    --ni-body 120 --ni-wake 30 --n-stations 5 --wake-length 15 \
    --nk 20 --ds 5e-4 --growth 1.15 \
    -o naca0012
# → 16 blocks total: 4 C-grid + 6 cap_jlo + 6 cap_jhi.

# Legacy 1-block (sharp only):
shore mesh naca0012.stl --topology cgrid --cgrid-blocks 1 \
    --ni-body 120 --ni-wake 30 --n-stations 5 --wake-length 15 \
    --nk 30 --ds 5e-4 --growth 1.15 \
    -o naca0012
# → naca0012.geo
```

Required flags (any layout):

- `--ni-body` — per-station body wrap node count.
- `--ni-wake` — per-branch wake node count (TE included).
- `--n-stations` — number of span-axis stations (becomes block `nj`).
- `--wake-length` — physical length of each wake branch from the TE,
  in body-coordinate units. A typical value is 15-20 chord lengths.

Optional flags:

- `--cgrid-blocks {1,3,4}` — block layout (default `3`).
- `--span-axis x|y|z` — defaults to `z`. The body's longest non-span
  axis is taken as the chord direction; the freestream defaults to
  `+x_world`-projected-onto-the-station-plane.
- `--wake-growth` — geometric growth ratio for wake i-spacing
  (default `1.15`). Wake cells start at the body's TE-cell length and
  grow toward the far field.
- `--te-dihedral-deg` — sharpness threshold for TE detection (default
  `30°`). Lower for blunter TEs (e.g. `20°`); higher to reject
  spurious tessellation seams.

Conformal-cap flags (require `--topology cgrid`):

- `--cap-tips {none,jlo,jhi,both}` — emit conformal tip caps (default
  `none`).  See [Tip caps — Conformal cap (recommended)](#conformal-cap-recommended).
- `--cap-te-thickness-frac` — synthetic TE thickness as a fraction of
  chord (default `0.002`).  Used only on sharp-TE bodies.
- `--cap-te-thickness` — absolute synthetic TE thickness; mutually
  exclusive with `--cap-te-thickness-frac`.
- `--cap-le-frac`, `--cap-te-frac` — chord fractions for the cap's
  LE-cut and TE-cut (default `0.05` and `0.95`).
- `--cap-tip-inset` — cap spanwise wall thickness (default `1e-3 *
  chord`).
- `--cap-n-cap-j` — cap j-axis (airfoil-thickness) node count (default
  `12`).
- `--cap-j-spacing-cgrid`, `--cap-j-beta-cgrid` — cap j-axis spacing
  law and clustering strength.
- `--cap-nk` — cap spanwise resolution (default = `--nk`).

The output is one `<stem>_<label>.geo` file per block (or one
`<stem>.geo` for `--cgrid-blocks 1`).  When `--cap-tips != none`, cap
blocks are written alongside the C-grid blocks as
`<stem>_cap_{jlo,jhi}_<label>.geo`.

## Python API

```python
from shore.mesh import Body, Mesh, CGridTopology, Spacing1D

body = Body.from_stl("naca0012.stl", skip_hygiene=True)
topo = CGridTopology()
blocks = topo.build(
    surface=None,
    nk=30,
    mesh=body.mesh,
    spacing_k=Spacing1D(law="geometric", ds=5e-4, growth=1.15),
    ni_body=120, ni_wake=30, n_stations=5,
    wake_length=15.0,
    span_axis="z",
    blocks=3,    # 3 (default), 4 (blunt-TE), or 1 (legacy)
)
m = Mesh(body=body, blocks=blocks, topology=topo)
m.march(ds=5e-4, growth=1.15, smooth_strength=0.1, smooth_iters=1)
m.write_geo("naca0012")
# → naca0012_upper.geo, naca0012_body.geo, naca0012_lower.geo
```

To produce a **conformal capped** wing in one call, pass the
`cap_tips` family of kwargs to `topo.build(...)`.  Sharp-TE bodies
get auto-promoted to the 4-block layout and the body's TE is rounded
over its full span:

```python
blocks = topo.build(
    surface=None, nk=20, mesh=body.mesh,
    spacing_k=Spacing1D(law="geometric", ds=5e-4, growth=1.15),
    ni_body=120, ni_wake=30, n_stations=5, wake_length=15.0,
    span_axis="z",
    cap_tips="both",                  # <- conformal caps
    cap_te_thickness_frac=0.002,      # 0.2% of chord (default)
    cap_le_frac=0.05, cap_te_frac=0.95,
    cap_tip_inset=0.02, cap_n_cap_j=12,
)
# → 16 blocks: 4 C-grid + 6 cap_jlo + 6 cap_jhi.
m = Mesh(body=body, blocks=blocks, topology=topo)
m.march(ds=5e-4, growth=1.15)         # caps are not marched
m.write_geo("naca0012_capped")
```

`Mesh.from_stl()` does **not** support `CGridTopology` because the
C-grid takes different parameters (`ni_body` / `ni_wake` /
`n_stations` / `wake_length`) that don't fit the legacy
`(ni, nj, nk)` signature. Build the block list directly with
`CGridTopology().build(...)`, then construct the `Mesh` from the
resulting blocks as shown above.

## Pipeline (algorithm summary)

The C-grid k=0 layer is built in four stages — see
[Algorithm — C-grid construction](/algorithm/c-grid) for derivations.

1. **TE detection** (`detect_te_edges`) — find the TE polyline(s)
   via dihedral classification.  Sharp TEs produce one strand;
   blunt TEs produce two (upper + lower).
2. **Span slicing** (`slice_by_span`) — slice the STL with
   `n_stations` planes perpendicular to the span axis.  Each slice
   is reordered to start at the upper-TE and resampled by
   arclength fraction.  For blunt TEs the body wrap terminates at
   the lower-TE node, excluding the TE-base segment.
3. **C-curve construction** (`build_c_curve`) — append upper and
   lower wake branches per station with geometric i-spacing.
4. **Topology assembly** + **hyperbolic march** — split the
   per-station polyline into 1, 3, or 4 blocks (`per_block_layers`
   on `CGridK0Result`); wire SHARED seams; march each block layer
   by layer with seam-aware ghost rows at the upper-TE and
   lower-TE seams.

## Body requirements

- **STL must be watertight enough that
  `trimesh.intersections.mesh_plane` can stitch a closed cross
  section at every span station.** Multi-component slices
  (interior holes, internal cavities) are rejected.
- **TE must be approximately span-aligned**, within ~30° (controls
  the candidate-edge filter inside `detect_te_edges`). Highly swept
  wings need a coarser threshold or a manually specified TE.
- **The body's chord axis must dominate** — the longer of the two
  non-span axes is auto-picked as the chord direction; the
  freestream defaults to chord-axis +1.

Two end-to-end examples run on a NACA 0012 fixture without external
data:

- [`examples/06-c-grid/visualize_naca0012.py`](https://github.com/szaghi/shore/blob/main/examples/06-c-grid/visualize_naca0012.py)
  — Python script with optional PyVista viewer.
- [`examples/06-c-grid/visualize_naca0012.sh`](https://github.com/szaghi/shore/blob/main/examples/06-c-grid/visualize_naca0012.sh)
  — bash equivalent using only `shore` CLI calls (`shore profile naca`
  → `shore body extrude` → `shore mesh --topology cgrid` → `shore export`).

## Tip caps (closed finite wing)

The C-grid wraps the wing's outer surface but leaves the spanwise
boundaries (`j=0`, `j=n_stations-1`) FREE — i.e. the wing is implicitly
infinite-span or terminates at a Chimera-fringe boundary.  For a
finite wing with a closed tip face, SHORE provides **two paths**:

1. **Conformal cap** (recommended) — ``shore mesh --topology cgrid
   --cap-tips both``.  The C-grid body's TE is auto-rounded over its
   full span by ``cap_te_thickness_frac`` (default 0.2% chord) so the
   cap and body meshes match bit-exactly at the airfoil-wall interface.
   The cap blocks are emitted as part of the same ``shore mesh``
   invocation; no separate command is needed.  Block count: 4 C-grid
   + 6 cap blocks per requested tip.  See [Conformal cap](#conformal-cap-recommended)
   below.

2. **Chimera-fringe cap** (independent overset) — ``shore cap-tip``
   command on an existing C-grid body block.  The body wrap stays
   sharp; the cap auto-rounds *locally* and overlays as an overset
   chimera component.  Use this when you want to add caps to an
   already-meshed C-grid without rebuilding it, or when chimera
   interpolation at the cap-body seam is acceptable.  Block count:
   5 cap blocks per tip.  See [Chimera-fringe cap](#chimera-fringe-cap)
   below.

Both paths produce **5 hex blocks** in an H-style "T-fill":

- ``cap_le``       — wraps the LE from upper-LE-cut to lower-LE-cut.
- ``cap_upper``    — chord band on the suction side, between LE-cut and
  TE-cut, above the chord cut at y=0.
- ``cap_lower``    — pressure-side mirror of ``cap_upper``.
- ``cap_te_upper`` — TE-region suction side, between TE-cut and the
  airfoil's TE-base segment.
- ``cap_te_lower`` — pressure-side mirror of ``cap_te_upper``.

A 6th block ``cap_wake_base`` is added when the user's body wrap is
blunt (i.e. has a real ``wake_base`` block) and ``is_blunt=True`` is
passed to ``CGridTipCap`` together with the wake_base cross-section.

The 5-block layout has **no singular axis** as long as the airfoil
cross-section has a non-zero TE thickness.  Sharp-TE inputs are
auto-rounded by ``CGridTipCap`` — see [Auto-round on sharp TE](#auto-round-on-sharp-te) below.

Default face BCs:

- **k_lo**: `WALL` — the actual wing-tip wall closure.
- **k_hi**: `FREE` — chimera fringe receiving from the C-grid body wrap.
- **j_hi**: `WALL` — the airfoil suction / pressure / LE / TE-base
  surface, mirrored into the cap.
- **j_lo**: `SHARED` between ``cap_upper`` and ``cap_lower`` (chord cut),
  and between ``cap_te_upper`` and ``cap_te_lower`` (chord cut continuation
  into the TE region).
- **i_lo / i_hi**: `SHARED` at the LE-cuts (cap_le ↔ cap_upper / cap_lower)
  and at the TE-cuts (cap_upper / cap_lower ↔ cap_te_upper / cap_te_lower).
  Outer ``i_lo`` of cap_le and outer ``i_hi`` of the TE blocks are `FREE`.

### Conformal cap (recommended)

Driven by ``shore mesh --topology cgrid --cap-tips {jlo|jhi|both}`` (or
the equivalent ``cap_tips`` kwarg on ``CGridTopology.build``).

When ``cap_tips != "none"``:

1. **Body TE rounding.**  For a sharp-TE body the C-grid topology
   auto-rounds the body's TE **over its full span** by ``cap_te_thickness``
   (or ``cap_te_thickness_frac * chord``, default 0.2%).  The body's
   i=0 and i=ni-1 columns are split apart by the chosen thickness in
   the airfoil-thickness direction at every spanwise station.  The
   upper / lower wake k=0 strips are shifted in bulk to ``y = ±te/2``
   along their full chord-axis extent so the ``wake_base`` block has
   uniform TE-thickness y-spread end to end (not just at the TE-side
   endpoint).
2. **Auto-promote 3 → 4 blocks.**  The 3-block sharp-TE C-grid layout
   is incompatible with non-zero TE thickness, so ``cap_tips != "none"``
   forces the 4-block layout (``upper``, ``body``, ``lower``,
   ``wake_base``).  A one-line stderr note announces this.
3. **Cap-block construction.**  Each requested tip emits 6 cap blocks
   (``cap_<tip>_le``, ``cap_<tip>_upper``, ``cap_<tip>_lower``,
   ``cap_<tip>_te_upper``, ``cap_<tip>_te_lower``, ``cap_<tip>_wake_base``).
   The cap's airfoil-wall nodes are **bit-coincident** with the body's
   k=0 wall nodes at the spanwise tip (built from the same rounded
   body cross-section).  The ``cap_<tip>_wake_base`` block's k_hi face
   is bit-coincident with the body's wake_base block's spanwise-extreme
   face.

The result is a **fully conformal mesh**: cap and body share their
interface nodes exactly, no Chimera fringe needed at the cap-body
seam.  The downside is that the body's TE is no longer geometrically
sharp (it has the synthetic ``te_thickness`` over its full span); this
is a small aerodynamic artifact at the TE / wing-tip junction,
typically below the discretisation noise floor.

```bash
shore mesh naca0012.stl --topology cgrid \
    --cap-tips both \
    --cap-te-thickness-frac 0.002 \
    --cap-le-frac 0.05 --cap-te-frac 0.95 \
    --cap-tip-inset 0.02 --cap-n-cap-j 12 \
    --ni-body 120 --ni-wake 30 --n-stations 5 --wake-length 15 \
    --nk 20 --ds 5e-4 --growth 1.15 \
    -o wing
```

Outputs: ``wing_upper.geo``, ``wing_body.geo``, ``wing_lower.geo``,
``wing_wake_base.geo`` (4 C-grid) plus 6 cap blocks per requested tip
(``wing_cap_jlo_*.geo`` and/or ``wing_cap_jhi_*.geo``).  At
``--cap-tips both`` the total is **16 blocks**.  No separate
``shore cap-tip`` invocation is needed.

The single ``shore mesh`` call produces a perfectly-matched mesh; downstream
``shore cc-par`` sees the cap blocks as part of the same connected
component (no chimera-fringe BCs at the cap-body seam).

Worked example:
[`examples/06-c-grid/visualize_naca0012_cap_tip.py`](https://github.com/szaghi/shore/blob/main/examples/06-c-grid/visualize_naca0012_cap_tip.py)
(or its bash sibling) builds a NACA 0012 wing with closed tips on
both sides via the conformal path and emits a 16-block ParaView ``.vtm``.

### Chimera-fringe cap

The standalone ``shore cap-tip`` command produces an overlapping cap
component for a body wrap that **stays sharp**.  The cap-body seam
relies on Chimera fringe interpolation at solver runtime.  Use this
when:

- You want to add caps to an existing C-grid `.geo` output without
  re-running ``shore mesh``.
- You want to keep the body's geometry sharp at the TE everywhere
  (only the cap is locally rounded).
- You have an STL body wrap that's already sharp and you don't want
  to alter it.

The cap-body interface is **not bit-coincident** in this path; cap and
body have a small geometric offset where the cap's synthetic TE
thickness doesn't match the body's sharp TE.  Xall handles this via
its standard chimera-fringe machinery.

CLI:

```bash
shore cap-tip naca0012_body.geo \
    --tip jlo \
    --le-frac 0.05 --te-frac 0.95 \
    --tip-inset 0.02 --n-cap-j 12 \
    -o cap_jlo --adjacency-json cap_jlo.adjacency.json
```

See the [`shore cap-tip` reference](/reference/cap-tip) for the full
CLI surface.

### Auto-round on sharp TE

Sharp-TE airfoils have i=0 and i=ni-1 of the cross-section coincident
at a single point.  An H-fill quad block whose corner sits on this
sharp point inevitably collapses to a singular axis (zero internal
angle).  Both cap paths handle this by **auto-rounding**:

1. Detect that the cross-section is sharp.
2. Synthesise a blunt cross-section by **splitting the TE node** into
   two points spaced ``te_thickness`` apart along the airfoil-thickness
   direction (anti-symmetric about the original TE point).
3. Build the 5-block cap from the synthesised blunt cross-section.

The user's **C-grid body wrap is left untouched** — only the cap is
rounded, locally at the wing tip, for its own structural soundness.
Aerodynamically this introduces a small artifact at the TE on the
wing-tip wall (a tiny TE-base segment, default 0.2% of chord), but it's
the only way to close a sharp-TE wing tip with a singular-axis-free
hex topology.

The synthetic TE thickness is controlled via:

- ``te_thickness_frac`` (default ``0.002`` = 0.2% chord) — fraction of
  the airfoil's chord extent.
- ``te_thickness`` — absolute value, mutually exclusive with
  ``te_thickness_frac``.
- ``auto_round_te=False`` — refuses sharp inputs with a clear error.

When ``CGridTipCap`` auto-rounds it prints a one-line stderr note:

```
CGridTipCap: sharp TE detected; auto-rounding cap with te_thickness=0.002 (the body wrap is unchanged).
```

For blunt-TE inputs (i=0 and i=ni-1 already distinct), no rounding is
applied — the cap uses the body's existing TE thickness.

### Conformal vs chimera-fringe: how the two paths differ

The cap fills the airfoil-interior cross-section at the wing-tip
plane.  The body wrap fills the fluid extruded outward from the
airfoil surface.  These two 3D regions meet only on a **1D edge**
(the airfoil profile line at the tip plane), not on a 2D face.  Two
approaches:

**Conformal path** (the recommended `--cap-tips` flag on `shore
mesh`).  The body wrap is rounded so its TE has a finite thickness
matching the cap's, and the cap blocks are constructed reading the
*same* cross-section nodes as the body wrap.  Cap and body share
their interface nodes **bit-coincidentally** by construction — they
hold equal float values, even though they live in independent
arrays (no view aliasing).  Operationally this is identical to a
SHARED face seam from Xall's point of view: chimera search lands
trivially with donor-receiver mapping = identity at every shared
node.  The formal SHARED-face graph in the adjacency JSON does *not*
declare a cap-body face seam (the cap and body remain formally
independent components), but the geometry guarantees zero
interpolation error at the interface.

**Chimera-fringe path** (the standalone `shore cap-tip` command).
The body wrap stays sharp-TE; the cap is built independently with
its own auto-rounded local cross-section.  Cap and body do **not**
have coincident nodes at the interface (the cap's synthetic TE
thickness differs from the body's sharp TE).  The cap-body coupling
relies on chimera-fringe interpolation at solver runtime — well-
suited to assemblies where the cap is added to a pre-existing
C-grid `.geo` output.

A true SHARED 2D face between cap and body is not expressible in
SHORE's current connection model because the geometric overlap is
1D only.  Adding 1D-edge SHARED would require a refactor of `Face`,
the adjacency JSON schema, the `cc.par` writer, the splitter, and
several debug viewers — a substantial effort that the conformal
path's bit-coincident-node approach avoids.

### Cap construction details (both paths)

The cap is **fully volumetric at build time** — its nodes are
determined parametrically from the body's tip cross-section plus the
cap parameters.  `Mesh.march()` is a silent no-op for cap blocks.

The standalone Python API for the chimera-fringe path:

```python
from shore.io.geo import read_geo
from shore.mesh import CGridTipCap, Mesh
from shore.topology.cgrid_cap import find_symmetric_le_cuts

body = read_geo("naca0012_body.geo")           # (ni, n_stations, nk, 3)
cross = body[:, 0, 0, :]                        # jlo tip, k=0 (the wall)
le_up, le_lo = find_symmetric_le_cuts(cross, le_frac=0.05, span_axis_world=2)

topo = CGridTipCap(
    le_cut_upper_idx=le_up,
    le_cut_lower_idx=le_lo,
    span_dir=(0.0, 0.0, -1.0),     # outward at the jlo tip
    tip_inset=0.02,
    n_cap_j=12,
    is_blunt=False,                 # set True + wake_base_cross_section for blunt
)
cap = Mesh.from_array(surface=cross, topology=topo, nk=body.shape[2])
cap.write_geo("cap_jlo")            # 5 .geo files
```

For the conformal path use `CGridTopology.build(cap_tips=...)`
directly (see [Python API](#python-api) above).

Worked examples on a NACA 0012 wing:

- [`examples/06-c-grid/visualize_naca0012_cap_tip.py`](https://github.com/szaghi/shore/blob/main/examples/06-c-grid/visualize_naca0012_cap_tip.py)
  — Python pipeline using the conformal path.
- [`examples/06-c-grid/visualize_naca0012_cap_tip.sh`](https://github.com/szaghi/shore/blob/main/examples/06-c-grid/visualize_naca0012_cap_tip.sh)
  — bash equivalent using `shore mesh --cap-tips both`.

Both produce 16-block ParaView `.vtm` files; the cap-body interface
is bit-coincident.

### Future: STL/STEP cap surfaces

The current ``CGridTipCap.build`` accepts a 1D cross-section
polyline as ``surface``.  A future extension will accept a 2D
sliced STL/STEP cross-section ``(ni, nj_thin, 3)`` so the cap can
be driven directly by a CAD/STL surface input — the same way the
C-grid body wrap is.  The signature stays the same.

## See also

- [Algorithm — C-grid construction](/algorithm/c-grid) — full derivation
  of TE detection, span slicing, C-curve build, and the SHARED-seam
  wiring for 3-block / 4-block topologies.
- [CH hybrid topology](/topologies/ch-hybrid) — the C-grid wrap
  composed with a far-field H-channel (12-/14-block contiguous
  body-to-frame mesh).
- [`shore mesh` reference](/reference/mesh) — full CLI command reference.
- [Pipeline overview](/guide/pipeline) — where the C-grid fits in the
  body-fitted workflow.
