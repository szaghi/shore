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

The result is a single `(ni_total, n_stations, 3)` k=0 layer ready
for the standard hyperbolic marcher to consume — `CGridTopology.build`
allocates the volume array, fills k=0 from this layer, and returns
a single `HexBlock` with `periodic_i = periodic_j = False`. The marcher
in `Mesh._march_cgrid` runs the same hyperbolic step + Laplacian
smoothing pipeline as `_march_ogrid`, with both axes non-periodic
and the open Jacobian guard.

## See also

- [C-grid topology](/topologies/c-grid) — consumer-facing API and
  CLI flags.
- [Hyperbolic marching](hyperbolic-marching) — the shared marcher
  the C-grid feeds into after k=0 construction.
- [Spacing laws](spacing-laws) — the geometric / tanh distributions
  used for wake i-spacing and wall-normal k-spacing.
