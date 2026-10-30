---
title: shore.surface.cgrid
outline: [2, 3]
---

# `shore.surface.cgrid`

C-grid surface construction: TE detection, span slicing, C-curve build.

A C-grid is a single-block structured mesh that wraps around a body
with a sharp trailing edge — typically an airfoil or hydrofoil — by
cutting the topological cylinder along the wake.  Three pieces of
geometry processing turn an STL into the C-grid k=0 layer:

1. :func:`detect_te_edges` — find the TE polyline via dihedral-angle
   classification on triangle-pair edges.
2. :func:`slice_by_span` — slice the STL with planes perpendicular
   to the span axis, producing one cross-section polyline per
   station, ordered upper-TE → LE → lower-TE.
3. :func:`build_c_curve` — add upper + lower wake branches to the
   per-station body polyline, returning the (ni_total, 2) C-curve.

The orchestrator :func:`build_cgrid_k0` chains them together to
produce the (ni, n_stations, 3) k=0 layer that
:class:`shore.mesh.topology.CGridTopology` consumes.

Conventions
-----------
- `span_axis` defaults to `"z"`; the body extends from
  `span_axis_min` to `span_axis_max` along this axis.
- The C-curve i-direction starts at the **upper TE** (i=0), sweeps
  forward over the suction surface, around the LE, back over the
  pressure surface, ends at the **lower TE** (i=ni_body-1).  The
  upper wake branch then extends i from `-ni_wake` to `0`
  (concatenated below the body wrap in i-space, so total i runs
  `[upper-wake-far, upper-TE, ..., LE, ..., lower-TE, lower-wake-far]`).
- Wake branches extend in the **global freestream direction**
  (default `+x_world`), regardless of local TE orientation.
- Per-station polylines are re-parameterised by **arclength fraction**
  so i = 0.5 sits at the LE on every station, regardless of chord
  variation along the span (tapered / swept wings).

## `TEDetectionResult` <Badge type="tip" text="dataclass" /> {#TEDetectionResult}

```python
TEDetectionResult(
    polylines: list[np.ndarray],
    is_blunt: bool,
    te_thickness: float,
    span_axis: str,
    max_dihedral_deg: float,
    n_candidate_groups: int,
)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/surface/cgrid.py#L57">source</a>

Outcome of :func:`detect_te_edges`.

**Attributes**

- **`polylines`** — List of ordered `(n_pts, 3)` arrays of TE vertex coordinates, each sorted by `span_axis`.  A **sharp** TE produces one polyline; a **blunt** TE produces two (the upper-edge and lower-edge feature lines of the finite-thickness TE strip).  Order: by descending span extent.
- **`is_blunt`** — `True` iff `len(polylines) >= 2` and the two strands are within `blunt_te_max_thickness * chord` of each other (i.e. they're the upper / lower of the same TE, not unrelated features).  Consumers use this flag to decide whether to treat the section as sharp- or blunt-TE.
- **`te_thickness`** — For blunt TEs (`is_blunt == True`), the average distance between the two strands at matched span coordinates.  `0.0` for sharp TEs.
- **`span_axis`** — The axis label (`"x"` / `"y"` / `"z"`) the polylines are sorted on.
- **`max_dihedral_deg`** — Largest dihedral angle observed in the input mesh, in degrees. Useful for diagnostics when no TE is found.
- **`n_candidate_groups`** — Number of disconnected span-aligned candidate groups discovered.

### `TEDetectionResult.polyline` <Badge type="info" text="property" /> {#TEDetectionResult-polyline}

```python
polyline: np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/surface/cgrid.py#L98">source</a>

Convenience: the longest TE polyline.

For sharp TEs this is the only polyline; for blunt TEs it's
either the upper or lower strand (consumers that don't care
about blunt TEs can use this directly).

## `detect_te_edges` <Badge type="info" text="function" /> {#detect_te_edges}

```python
detect_te_edges(
    mesh: trimesh.Trimesh,
    *,
    dihedral_deg: float = 30.0,
    span_axis: str = 'z',
) -> TEDetectionResult
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/surface/cgrid.py#L113">source</a>

Find the trailing-edge polyline of a body STL by dihedral classification.

For every edge shared by exactly two triangles, compute the
dihedral angle between the two face normals.  An edge whose
dihedral exceeds `dihedral_deg` (i.e. the body folds by more
than `180° - dihedral_deg` flat-line at that edge) is a TE
candidate.  Candidates are clustered by vertex connectivity; the
largest connected group is taken as the trailing edge.  The
resulting polyline is ordered by `span_axis` coordinate.

**Parameters**

- **`mesh`** — Watertight `trimesh.Trimesh` of the body.
- **`dihedral_deg`** — Sharpness threshold in degrees.  An edge is a TE candidate when `angle(face_normals) >= dihedral_deg`.  Default `30°` (≈ feature line at body curvature).  Lower (e.g. `20°`) for blunter TEs; higher (e.g. `45°`) to reject spurious tessellation seams.
- **`span_axis`** — `"x"`, `"y"`, or `"z"`.  The polyline is sorted along this axis at the end.  Default `"z"`.

**Returns**

- `TEDetectionResult` — `polyline` is the ordered `(n_te_pts, 3)` array of TE vertices.

**Raises**

- `MeshInputError` — If no edge exceeds the threshold (body is too smooth — try a smaller `dihedral_deg` or supply the TE manually), or if the largest candidate group is too small (&lt; 2 points) to form a polyline, or if the resulting polyline isn't approximately aligned with `span_axis`.
- `ParameterError` — On invalid arguments.

## `SpanSliceResult` <Badge type="tip" text="dataclass" /> {#SpanSliceResult}

```python
SpanSliceResult(polyline_2d: np.ndarray, te_position_2d: np.ndarray, span_coord: float)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/surface/cgrid.py#L389">source</a>

Per-station cross-section polyline plus its TE anchor.

**Attributes**

- **`polyline_2d`** — `(ni_body, 2)` arclength-resampled cross-section polyline in the local in-plane coordinate frame (chord axis x thickness axis).  Ordered upper-TE -&gt; LE -&gt; lower-TE.
- **`te_position_2d`** — `(2,)` location of the upper-TE anchor on this station, in the same local 2D frame.  For sharp TE this coincides with `polyline_2d[0]` and `polyline_2d[-1]`; for blunt TE it sits between them.
- **`span_coord`** — Span-axis coordinate of this station.

## `slice_by_span` <Badge type="info" text="function" /> {#slice_by_span}

```python
slice_by_span(
    mesh: trimesh.Trimesh,
    *,
    n_stations: int,
    ni_body: int,
    span_axis: str = 'z',
    te_result: TEDetectionResult | None = None,
) -> list[SpanSliceResult]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/surface/cgrid.py#L412">source</a>

Slice the STL at `n_stations` planes perpendicular to span.

Each slice produces one closed cross-section polyline.  The
polyline is rotated into a local 2D frame (the span axis dropped),
re-ordered so it starts at the upper TE and walks around the LE,
and resampled by **arclength fraction** to a common `ni_body`
node count.  This common parametrisation lets stations with
different chord lengths (tapered / swept wings) line up at the
same parametric i-position — i = 0 is upper-TE on every station,
i = ni_body//2 is approximately the LE on every station,
i = ni_body-1 is lower-TE.

**Parameters**

- **`mesh`** — Watertight body STL.
- **`n_stations`** — Number of span stations to slice.  Stations are placed uniformly between the body's span-axis min and max, with the first and last stations strictly inside the body (offset by 2% of the span extent) to avoid degenerate slices on the end-cap faces.
- **`ni_body`** — Common point count of the resampled cross-section polyline. Must be `>= 8`.  Even values give symmetric LE/TE placement; odd values place a single node exactly at the LE.
- **`span_axis`** — `"x"`, `"y"`, or `"z"`.  Default `"z"`.
- **`te_result`** — Optional `TEDetectionResult` from :func:`detect_te_edges`.  When supplied, the per-station polyline ordering anchors at the TE; without it, the function picks the polyline endpoint closest to the chord-axis maximum as the TE — which works for symmetric airfoils but is fragile for cambered sections.

**Returns**

- `list[SpanSliceResult]` — One entry per station, sorted by span coordinate.

**Raises**

- `MeshInputError` — If a slice produces no polyline (station outside the body), a non-closed polyline (body has a hole at this station), or multiple disconnected polylines (body is multi-component).
- `ParameterError` — On invalid arguments.

## `CCurveResult` <Badge type="tip" text="dataclass" /> {#CCurveResult}

```python
CCurveResult(
    polyline_2d: np.ndarray,
    ni_wake_upper: int,
    ni_wake_lower: int,
    span_coord: float,
)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/surface/cgrid.py#L758">source</a>

Per-station C-curve: body wrap with upper + lower wake branches.

**Attributes**

- **`polyline_2d`** — `(ni_total, 2)` array in the local 2D frame (chord axis x thickness axis).  Layout, in i-order: - `[0 : ni_wake_upper]` — upper wake, far-field → upper-TE - `[ni_wake_upper : ni_wake_upper + ni_body]` — body, upper-TE → LE → lower-TE - `[ni_wake_upper + ni_body : ni_total]` — lower wake, lower-TE → far-field Equivalently: `i = 0` is the upper-wake far-field tip; `i = ni_total - 1` is the lower-wake far-field tip.  The TE sits at `i = ni_wake_upper` (upper) and `i = ni_wake_upper + ni_body - 1` (lower).
- **`ni_wake_upper, ni_wake_lower`** — Node counts of each wake branch (including the TE node, which is shared with the body).  When the C-curve is flattened, the shared TE node appears twice — once as the upper-wake terminus, once as the body i=0 — but the consumer typically uses `polyline_2d` directly, where each TE is a single node.
- **`span_coord`** — Span-axis coordinate of this station.

## `build_c_curve` <Badge type="info" text="function" /> {#build_c_curve}

```python
build_c_curve(
    slice_result: SpanSliceResult,
    *,
    ni_wake: int,
    wake_length: float,
    freestream_2d: np.ndarray,
    wake_growth: float = 1.15,
) -> CCurveResult
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/surface/cgrid.py#L791">source</a>

Add upper + lower wake branches to a station body polyline.

The body polyline (output of :func:`slice_by_span`) is a wrap from
upper-TE around the LE to lower-TE.  This function appends two
wake branches in the global freestream direction:

- **Upper branch**: starts at the upper-wake far-field tip,
  contracts toward the upper-TE; concatenated **before** the
  body wrap.
- **Lower branch**: starts at the lower-TE, expands toward the
  lower-wake far-field tip; concatenated **after** the body wrap.

Wake i-spacing is **geometric**: starting from the body's TE-cell
edge length and growing by `wake_growth` per cell, until the
wake reaches `wake_length` (truncated to `ni_wake` nodes).

**Parameters**

- **`slice_result`** — Output of :func:`slice_by_span` for one station — a body wrap polyline ordered upper-TE → LE → lower-TE.
- **`ni_wake`** — Number of nodes per wake branch.  Each branch has this many nodes, **including** the TE node it shares with the body. Must be `>= 2`.  Total C-curve length is `2 * (ni_wake - 1) + ni_body`.
- **`wake_length`** — Length of each wake branch from the TE, in body-coordinate units.  Typical value: 20 \* chord.
- **`freestream_2d`** — Length-2 unit vector `(u_chord, u_thick)` giving the freestream direction in the station's local frame.  Should approximate `[+1, 0]` for unswept wings; for swept wings, the projection of the global freestream onto the station plane.
- **`wake_growth`** — Geometric growth ratio for wake i-spacing.  `1.0` gives uniform spacing; `1.15` (default) gives mild stretching.

**Returns**

- `CCurveResult` — `polyline_2d` is the full C-curve in the local 2D frame.

**Raises**

- `ParameterError` — On invalid arguments.

## `CGridK0Result` <Badge type="tip" text="dataclass" /> {#CGridK0Result}

```python
CGridK0Result(
    nodes: np.ndarray,
    ni_body: int,
    ni_wake: int,
    span_axis: str,
    te_result: TEDetectionResult,
)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/surface/cgrid.py#L949">source</a>

3D k=0 surface layer of a C-grid topology.

**Attributes**

- **`nodes`** — `(ni_total, n_stations, 3)` array of world-frame node coordinates.  `ni_total = 2 * (ni_wake - 1) + ni_body`.
- **`ni_body, ni_wake`** — Counts used to build the curve.  `ni_body` is the body wrap node count; `ni_wake` is the per-branch wake node count (each wake branch shares its TE node with the body — see :class:`CCurveResult` for the exact slice layout).
- **`span_axis`** — `"x"`, `"y"`, or `"z"`.
- **`te_result`** — The :class:`TEDetectionResult` produced during k=0 build.

### `CGridK0Result.ni_total` <Badge type="info" text="property" /> {#CGridK0Result-ni_total}

```python
ni_total: int
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/surface/cgrid.py#L975">source</a>

Total i-direction node count of the concatenated layer.

### `CGridK0Result.n_stations` <Badge type="info" text="property" /> {#CGridK0Result-n_stations}

```python
n_stations: int
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/surface/cgrid.py#L980">source</a>

Span-direction station count.

### `CGridK0Result.per_block_layers` <Badge type="info" text="method" /> {#CGridK0Result-per_block_layers}

```python
per_block_layers() -> tuple[np.ndarray, np.ndarray, np.ndarray]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/surface/cgrid.py#L984">source</a>

Return per-block k=0 layers for the 3-block C-grid topology.

Splits the concatenated `nodes` array along i into three
contiguous segments, one per :class:`CGridTopology` block:

- **upper wake**: `nodes[0 : ni_wake]` — sweeps
  `upper-wake-far-tip → upper-TE`.
- **body**: `nodes[ni_wake - 1 : ni_wake - 1 + ni_body]` —
  sweeps `upper-TE → suction → LE → pressure → lower-TE`.
- **lower wake**: `nodes[ni_wake - 1 + ni_body - 1 : ni_total]` —
  sweeps `lower-TE → lower-wake-far-tip`.

The slices overlap by one node at each TE column (upper wake
last = body first, body last = lower wake first), so the
returned arrays share TE columns by **numpy view aliasing**
— modifying a TE node through one view updates the partner
view automatically.  This is the bit-exact node coincidence
the 3-block topology's SHARED seams rely on.

**Returns**

- **`(upper, body, lower)`** (tuple of np.ndarray) — Each `(ni_block, n_stations, 3)` view into `self.nodes`. Shapes: `(ni_wake, n_stations, 3)`, `(ni_body, n_stations, 3)`, `(ni_wake, n_stations, 3)`.

## `build_cgrid_k0` <Badge type="info" text="function" /> {#build_cgrid_k0}

```python
build_cgrid_k0(
    mesh: trimesh.Trimesh,
    *,
    ni_body: int,
    ni_wake: int,
    n_stations: int,
    wake_length: float,
    span_axis: str = 'z',
    freestream_world: np.ndarray | None = None,
    te_dihedral_deg: float = 30.0,
    wake_growth: float = 1.15,
) -> CGridK0Result
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/surface/cgrid.py#L1021">source</a>

Build the C-grid k=0 surface layer from a body STL.

Pipeline:

1. :func:`detect_te_edges` — find the TE polyline.
2. :func:`slice_by_span` — slice the STL into `n_stations`
   cross-section polylines, each resampled to `ni_body` nodes
   and ordered upper-TE → LE → lower-TE.
3. :func:`build_c_curve` — for each station, append upper +
   lower wake branches in the (per-station projected) freestream
   direction.
4. Lift the resulting 2D C-curves back into 3D by re-inserting
   the span coordinate, producing the world-frame
   `(ni_total, n_stations, 3)` k=0 layer.

**Parameters**

- **`mesh`** — Watertight body STL.
- **`ni_body`** — Per-station body wrap node count (upper-TE → LE → lower-TE).
- **`ni_wake`** — Per-branch wake node count (shared TE node included).
- **`n_stations`** — Number of span-axis stations.  `>= 2`.
- **`wake_length`** — Length of each wake branch from the TE, in body-coordinate units.
- **`span_axis`** — `"x"` / `"y"` / `"z"`.  Default `"z"`.
- **`freestream_world`** — Length-3 vector giving the freestream direction in world coordinates.  If `None`, defaults to the unit vector along the body's dominant chord axis (the largest non-span body extent), pointing **toward** the TE.
- **`te_dihedral_deg`** — Dihedral threshold for TE detection.  Default `30°`.
- **`wake_growth`** — Geometric growth ratio for wake i-spacing.

**Returns**

- `CGridK0Result` — `nodes` is the `(ni_total, n_stations, 3)` world-frame array, ready to be marched outward by a :class:`HexBlock`.

**Raises**

- `MeshInputError, ParameterError` — Propagated from the pipeline stages.

## `round_body_te_in_surface` <Badge type="info" text="function" /> {#round_body_te_in_surface}

```python
round_body_te_in_surface(
    surface: np.ndarray,
    ni_wake: int,
    ni_body: int,
    te_thickness: float,
    span_axis: str,
) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/surface/cgrid.py#L1211">source</a>

Round the body's TE in-place over the full span of a C-grid k=0 surface.

Splits the body's i=0 (upper-TE) and i=ni_body-1 (lower-TE) columns
apart by `te_thickness` along the airfoil-thickness direction at
every spanwise station.  Also shifts each wake block's k=0 strip
in bulk so the upper wake sits at `+te_thickness/2` and the lower
wake at `-te_thickness/2` along their full chord-axis extent —
not just at the TE-side endpoints.

Without the wake-strip shift the wake blocks' k=0 strips converge
back to `y=0` immediately downstream of the TE (the wake builder
in :func:`build_cgrid_k0` runs each branch from a sharp anchor at
y=0), producing a degenerate near-zero-thickness `wake_base`
block in the 4-block C-grid layout — a spurious singular line just
downstream of the TE.

The body sub-region of the concatenated k=0 array is
`surface[ni_wake-1 : ni_wake-1 + ni_body, :, :]`.  The upper-wake
region is `surface[: ni_wake, :, :]` (i=0 far-tip → i=ni_wake-1
upper-TE; the latter aliases the body's i=0 column).  The lower-
wake region is `surface[ni_wake-1 + ni_body - 1 :, :, :]`
(i=0 lower-TE → i=ni_wake-1 far-tip; the former aliases the body's
i=ni_body-1 column).

**Parameters**

- **`surface`** — `(ni_total, n_stations, 3)` k=0 surface array as produced by :func:`build_cgrid_k0`.  Must have a sharp TE at every station.
- **`ni_wake, ni_body`** — Block-count integers from the C-grid construction.
- **`te_thickness`** — World-frame TE thickness to apply.  Must be &gt; 0.
- **`span_axis`** — `'x'` / `'y'` / `'z'`.

**Raises**

- `ParameterError` — On invalid arguments or non-sharp input.
