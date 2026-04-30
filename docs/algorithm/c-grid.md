# C-grid construction

The C-grid topology turns a sharp-TE body STL (airfoil, hydrofoil,
wing) into a single-block hexahedral mesh that wraps the body and
extends two wake branches downstream of the trailing edge. This page
describes the **k=0 surface layer** — the algorithmic input to the
shared [hyperbolic marcher](hyperbolic-marching). The
[topology page](/topologies/c-grid) covers the consumer-facing API.

The pipeline runs in four stages, all in
[`shore.surface.cgrid`](https://github.com/szaghi/shore/blob/main/src/shore/surface/cgrid.py):

```text
trimesh.Trimesh
       │
       ▼  detect_te_edges          → TEDetectionResult
       │
       ▼  slice_by_span (x N)      → list[SpanSliceResult]
       │
       ▼  build_c_curve (x N)      → list[CCurveResult]
       │
       ▼  lift to 3D
              │
              ▼
       (ni_total, n_stations, 3) k=0 layer
```

Each stage has its own dataclass result type so a downstream user can
plug into the pipeline at any level — for instance, supply a manual
TE polyline by skipping `detect_te_edges` and feeding a
`TEDetectionResult` directly into `slice_by_span`.

## Trailing-edge detection

`detect_te_edges(mesh, *, dihedral_deg=30.0, span_axis="z")` finds the
TE polyline by inspecting **dihedral angles** at every shared edge of
the STL. The classifier runs three filters in sequence; the
intersection survives:

1. **Dihedral filter** — keep edges whose face-normal angle exceeds
   `dihedral_deg` (default `30°`). Sharp body features
   (TE, possibly LE on a Gurney flap, the cap-fan ridges of an
   extruded STL) all pass this filter.
2. **Span-axis filter** — keep edges whose direction is within ~30°
   of the span axis. The TE runs along the span; LE-projection
   ridges and end-cap fans run perpendicular to it. This filter
   alone rejects most spurious features.
3. **Chord-extreme filter** — among the surviving connected
   components ("groups"), keep only those whose maximum chord
   coordinate is within `5e-4` of the global maximum chord. The
   literal TE sits at chord-max; cosine-spacing artefacts on
   lat-lon-style sections sit slightly inside (≈ 0.2% of chord)
   and get filtered out.

Connected components are built via union-find on the candidate-edge
graph. Components are then sorted by descending span extent; only
those spanning ≥ 80% of the longest are retained, so a long TE strand
beats short tessellation seams.

**Sharp vs blunt**. A sharp TE produces a single component spanning
the body. A blunt (finite-thickness) TE produces two components
(the upper and lower edge of the TE strip), both spanning the body
and separated by the TE thickness. `TEDetectionResult.is_blunt` is
`True` iff at least two components survive and their average
inter-strand distance is below 5% of the body's bounding-box
diagonal; `te_thickness` reports the measured spacing.

The result's `polyline` property returns the longest strand sorted
along `span_axis`; consumers that ignore blunt TEs can use it
directly.

## Span slicing

`slice_by_span(mesh, *, n_stations, ni_body, span_axis="z", te_result=None)`
slices the STL with `n_stations` planes perpendicular to the span
axis. Stations are placed uniformly between the body's span min and
max, with a 2% margin at each end so the planes don't graze
end-cap triangles.

For each plane:

1. **Cut** with `trimesh.intersections.mesh_plane`, returning a list
   of 3D segments.
2. **Drop the span axis** to get 2D segments in the local
   `(chord, thickness)` frame. The chord axis is the longer of the
   two non-span body extents; thickness is the remaining axis.
3. **Stitch** the segments into a closed polyline. The stitcher
   builds a vertex-rounded adjacency graph and walks it; T-junctions
   (vertices appearing in > 2 segments) and multi-component cuts
   are rejected. Duplicate segments (which `mesh_plane` emits when
   the cutting plane intersects a triangulation diagonal exactly)
   are deduped before adjacency check.
4. **Reorder around the TE.** With `te_result`, interpolate the
   TE polyline at this station's span coordinate to get a 2D anchor;
   without it, fall back to the polyline endpoint at maximum chord.
   Rotate the closed loop so the closest node is at index 0, then
   walk one step forward — if thickness drops, flip the loop so
   the walk climbs to the suction surface first.
5. **Resample by arclength fraction.** Cumulative arclength is
   computed across the polyline, divided into `ni_body - 1` equal
   intervals, and node positions are linearly interpolated along
   each segment. This **common parametrisation** is what makes
   tapered / swept wings line up across stations: `i = 0` is the
   upper-TE on every station, `i ≈ ni_body / 2` is the LE region
   on every station regardless of how the chord varies along the
   span.

The result is a `SpanSliceResult` per station with:
`polyline_2d` (`(ni_body, 2)`), `te_position_2d` (`(2,)`), and
`span_coord` (scalar).

## C-curve

`build_c_curve(slice_result, *, ni_wake, wake_length, freestream_2d, wake_growth=1.15)`
adds two wake branches to a station's body wrap, producing the
`(ni_total, 2)` C-curve in the local 2D frame.

The body wrap is a sequence ordered upper-TE → LE → lower-TE. The
function appends:

- **Upper wake**: `ni_wake` nodes from the upper-wake far-field tip
  to the upper-TE. Wake direction is the supplied `freestream_2d`
  unit vector. i-spacing is **geometric** with growth ratio
  `wake_growth`, rooted at the body's TE-cell length
  (`||body[1] - body[0]||`), then **globally rescaled** so the
  cumulative arclength sums to exactly `wake_length`. The branch
  is reversed (far-field first, TE last) before concatenation so
  i = 0 sits at the upper-wake far-field tip.
- **Lower wake**: same construction at the lower-TE, oriented
  outward (TE first, far-field last) and concatenated after the
  body wrap.

The shared TE node is dropped on each side during concatenation, so
the final C-curve has `ni_total = 2 * (ni_wake - 1) + ni_body`
nodes.

For **sharp TEs** the upper and lower wake branches start at the
same point (the TE), so the C-curve is geometrically self-coincident
along the wake-cut at k=0. For **blunt TEs** the two wake roots are
separated by the TE thickness; the upper branch starts at the
upper-TE and the lower at the lower-TE, leaving a finite TE-base
patch on the body.

## Lifting to 3D

`build_cgrid_k0(mesh, *, ni_body, ni_wake, n_stations, wake_length, ...)`
chains the three stages together: detect the TE, slice every
station, build a C-curve per station, and lift each 2D C-curve
back into 3D by re-inserting the span coordinate. The lifting is
straightforward — the chord and thickness coordinates of the local
2D frame map back to fixed world axes (the chord axis and thickness
axis selected from the mesh bounds), and the span axis carries the
station's `span_coord`.

The result is a single `(ni_total, n_stations, 3)` k=0 layer.
`CGridK0Result.per_block_layers()` returns three numpy *views* into
that layer — `(upper, body, lower)` segments — that share the upper-TE
and lower-TE columns by view aliasing.  These views are what
`CGridTopology.build` carves into `HexBlock`s for the 3-block (and
4-block) layouts.

## Topology assembly

`CGridTopology.build(..., blocks=N)` dispatches the per-block layer
triple into 1, 3, or 4 hex blocks, each with its own SHARED-seam
wiring.

### 3-block layout (sharp TE, default)

Three hex blocks `upper / body / lower`, with three SHARED face
pairs:

| Seam | Faces | AxisMap (this side) | AxisMap rationale |
|---|---|---|---|
| Upper-TE column | `upper.i_hi ↔ body.i_lo` | identity on `(j, k)` | Same i-direction on both sides; node coincidence via numpy view aliasing of the shared TE column. |
| Lower-TE column | `body.i_hi ↔ lower.i_lo` | identity on `(j, k)` | Same. |
| Wake cut | `upper.k_lo ↔ lower.k_lo` | `((i, True), (j, False))` | Same-axis k↔k.  Upper sweeps far-tip→TE while lower sweeps TE→far-tip — same physical wake-cut line, opposite parametric directions, hence the i-flip. |

The wake-cut SHARED face is the structurally-honest replacement for
what was the **single-block C-grid's i-direction self-overlap**:

```
1-block layout (legacy):
  i = 0..ni_wake-1            → upper wake (k=0 at y=0 for sharp TE)
  i = ni_wake-1..ni_wake-1+ni_body-1
                              → body wrap (TE→LE→TE)
  i = ni_wake-1+ni_body-1..ni_total-1
                              → lower wake (k=0 also at y=0)

The upper-wake and lower-wake k=0 lines are bit-coincident in space
but topologically distinct in i.  At the wake-cut edge, the j-counter
"flips direction" — because index ``i = 5`` on upper-wake's k=0 is
the same physical point as index ``i = ni_total - 6`` on lower-wake's
k=0, but they sit on opposite sides of the same parametric line.
Xall's BC dispatcher cannot apply BCs cleanly at this self-overlap.
```

The 3-block layout makes the same geometry into a clean SHARED face
between two distinct blocks, with the i-flip recorded explicitly on
the AxisMap.  Xall sees a face-to-face adjacency and applies BCs
correctly.

### 4-block layout (blunt TE)

For a body with non-zero TE thickness, the upper and lower wake k=0
strips are separated in space by `te_thickness`.  The 3-block
layout's wake-cut SHARED face would have a non-zero gap — invalid
as SHARED — and the TE base wall of the airfoil would not be
meshed at all.

The 4-block layout adds a **wake-base block W** that fills the gap.
W's shape is `(ni_wake, n_stations, nk_w, 3)` with i along the wake
(matching upper's direction) and k spanning the TE thickness gap
from `y = -te/2` (= lower's k=0 line, with i reversed) to
`y = +te/2` (= upper's k=0 line, identity).  W's `i_hi` face is the
**TE base WALL**.

The wake-cut SHARED face from the 3-block layout becomes two
separate seams in the 4-block layout:

| Seam | Faces | AxisMap (W side) | AxisMap rationale |
|---|---|---|---|
| W upper face | `wake_base.k_hi ↔ upper.k_lo` | identity on `(i, j)` | W's i and upper's i both sweep far-tip→TE; node coincidence via direct copy at build time. |
| W lower face | `wake_base.k_lo ↔ lower.k_lo` | `((i, True), (j, False))` | W's i sweeps far-tip→TE, lower's i sweeps TE→far-tip; same i-flip as in the 3-block wake-cut. |

W is **not marched** — its geometry is fully determined at build time
by linear interpolation in k between the lower- and upper-wake k=0
strips.  `Mesh._march_cgrid_three`'s per-block loop iterates only
over `[upper, body, lower]` even when there are 4 blocks.

For blunt TE, `slice_by_span` also terminates the body wrap at the
**lower-TE node** rather than continuing along the TE base — this
keeps `body.i_hi = lower-TE` and lets W mesh the TE base separately.

### 1-block layout (legacy)

`blocks=1` returns the single concatenated block described in the
pathology callout above.  The wake-cut self-overlap is silent.
Kept for opt-in compatibility with downstream tools that don't need
clean Xall BCs.

## See also

- [C-grid topology](/topologies/c-grid) — consumer-facing API and
  CLI flags.
- [Hyperbolic marching](hyperbolic-marching) — the shared marcher
  the C-grid feeds into after k=0 construction.
- [Spacing laws](spacing-laws) — the geometric / tanh distributions
  used for wake i-spacing and wall-normal k-spacing.
