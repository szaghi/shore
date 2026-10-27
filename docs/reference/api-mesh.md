# shore.mesh.Mesh

The OOP entry point for SHORE: a multi-block structured mesh.
`Mesh.from_stl(...)` builds the blocks (surface projection + topology
wiring) and `Mesh.march()` extrudes them into a volume by hyperbolic
marching. Every `shore mesh` CLI invocation is a thin wrapper around
this class.

```python
from shore.mesh import CubedSphereTopology, Mesh

mesh = Mesh.from_stl(
    "sphere.stl",
    topology=CubedSphereTopology(),
    ni=40, nj=60, nk=30,
    theta_cap_deg=30.0,
)
mesh.march(ds=1e-3, growth=1.15)
mesh.write_geo("sphere")          # 6 sphere_<label>.geo files
```

## Class

```python
class Mesh:
    body:     Body
    blocks:   list[HexBlock]
    topology: Topology
```

| Attribute | Description |
|---|---|
| `body` | The triangulated body surface and its cached BVH (only set when built via `from_stl`; `None` when built via [`from_array`](#meshfrom_array)). |
| `blocks` | The list of [`HexBlock`](#hexblock-attribute) instances produced by the topology. Order matches the topology's canonical block list (sub0..sub3, cap_north, cap_south for `CubedSphereTopology`). |
| `topology` | The [`Topology`](./api-topology) instance that produced `blocks`. Determines which marcher `march()` dispatches to (cubed-sphere, O-grid, C-grid, or OH hybrid). |

## `Mesh.from_stl`

```python
@classmethod
def from_stl(
    cls,
    path: str | Path,
    topology: Topology,
    ni: int,
    nj: int,
    nk: int,
    theta_cap: float | None = None,
    theta_cap_deg: float | None = None,
    i_spacing: str = "uniform",
    i_beta: float = 3.0,
    j_spacing: str = "uniform",
    j_beta: float = 3.0,
    j_cluster_lon: float | None = None,
    spacing_k: Spacing1D | None = None,
) -> Mesh
```

Load an STL, project the surface, and wire the blocks. **Does not
march** — call `Mesh.march()` afterwards to fill the wall-normal
layers.

**Parameters**

| Name | Type | Default | Description |
|---|---|---|---|
| `path` | `str \| Path` | — | Path to the body's STL file. Loaded via [`load_stl`](./api-surface-io). |
| `topology` | `Topology` | — | A wired [`Topology`](./api-topology) builder — typically `CubedSphereTopology()`, `OGridTopology()`, `CGridTopology(...)`, or `OHTopology(...)`. |
| `ni`, `nj` | `int` | — | Surface **node counts** (the `(ni, nj)` shape of each block's k=0 layer). For `CubedSphereTopology`, `nj` must be a multiple of 4 (rounded down with a warning if not). |
| `nk` | `int` | — | Wall-normal **node count** including the wall row at k=0. With `nk = 30` the marcher extrudes 29 cells from the wall. |
| `theta_cap` | `float \| None` | `0.05*pi` | Polar exclusion angle in radians. Mutually exclusive with `theta_cap_deg`. |
| `theta_cap_deg` | `float \| None` | — | Same as `theta_cap` in degrees. |
| `i_spacing` | `str` | `"uniform"` | Surface i-direction spacing law: `"uniform"` / `"tanh"` / `"tanh2"`. |
| `i_beta` | `float` | `3.0` | Clustering strength for tanh / tanh2 on i. |
| `j_spacing` | `str` | `"uniform"` | Surface j-direction spacing law: `"uniform"` or `"tanh2"`. |
| `j_beta` | `float` | `3.0` | Clustering strength for tanh2 on j. |
| `j_cluster_lon` | `float \| None` | `None` | (legacy single-block / O-grid only) longitude in radians to cluster the j-distribution toward. Ignored by `CubedSphereTopology`. |
| `spacing_k` | `Spacing1D \| None` | `Spacing1D(law="geometric", ds=1e-3, growth=1.15)` | Wall-normal spacing law. See [`Spacing1D`](./api-spacing). |

**Returns**

`Mesh` — ready to call `.march()`.

**Raises**

- `MeshInputError` — STL malformed / not watertight / not star-shaped from anchor.
- `ParameterError` — invalid spacing law, conflicting angle args, `theta_cap` outside `(0, π/2)`.

---

## `Mesh.from_array`

```python
@classmethod
def from_array(
    cls,
    surface: np.ndarray,
    topology: Topology,
    nk: int,
    spacing_k: Spacing1D | None = None,
    **kwargs,
) -> Mesh
```

Build a `Mesh` from a pre-built surface array (no STL).  Used by
topologies whose input is a numeric array — typically a cross-section
polyline from another mesh — rather than a triangulated body.

The resulting `Mesh` has `body = None`; methods that depend on
`self.body` (e.g. `blend_normals_k > 0` STL-normal blending in
[`Mesh.march`](#meshmarch)) will fail loudly if called.

For [`CGridTipCap`](./api-topology#cgridtipcap), `surface` is a
`(ni_body, 3)` cross-section polyline and `nk` is the cap's
spanwise (k-axis) node count.

```python
from shore.io.geo import read_geo
from shore.mesh import CGridTipCap, Mesh

body = read_geo("naca0012_body.geo")
cross = body[:, 0, 0, :]                       # jlo tip, k=0 (the wall)
topo = CGridTipCap(
    le_cut_upper_idx=39, le_cut_lower_idx=80,
    span_dir=(0.0, 0.0, -1.0), tip_inset=0.02, n_cap_j=12,
)
cap = Mesh.from_array(surface=cross, topology=topo, nk=20)
cap.write_geo("cap_jlo")
```

---

## `Mesh.march`

```python
def march(
    self,
    ds: float = 1e-3,
    growth: float = 1.15,
    smooth_strength: float = 0.2,
    smooth_iters: int = 2,
    steps: np.ndarray | None = None,
    layer_callback: Callable[[int, float], None] | None = None,
    blend_normals_k: int = 0,
) -> None
```

Extrude every block in `self.blocks` from layer `k=0` outward, in
place. Dispatches to the per-topology marcher
(`CubedSphereTopology` → `_march_caps`, `OGridTopology` →
`_march_ogrid`, `CGridTopology` → `_march_cgrid_*` per block count,
`OHTopology` → `_march_caps_blocks` for the wrap then
`OHTopology.fill_channel_blocks` for the deferred channel
interior).  For [`CGridTipCap`](./api-topology#cgridtipcap), the cap
is fully volumetric at build time and `march()` is a **silent no-op**.

**Parameters**

| Name | Type | Default | Description |
|---|---|---|---|
| `ds` | `float` | `1e-3` | First-layer thickness for the geometric spacing law. Ignored when `steps` is supplied or when `spacing_k` overrides via per-edge data. |
| `growth` | `float` | `1.15` | Layer-to-layer thickness ratio for the geometric law. |
| `smooth_strength` | `float` | `0.2` | Laplacian blend factor per layer in `[0, 1]`. |
| `smooth_iters` | `int` | `2` | Number of Laplacian sweeps applied after each extrusion step. |
| `steps` | `np.ndarray \| None` | `None` | Pre-computed per-layer step sizes, shape `(nk-1,)`. Highest precedence — overrides `ds`/`growth` and any per-edge `spacing_k`. |
| `layer_callback` | `Callable[[int, float], None] \| None` | `None` | Called as `layer_callback(k, jmin_global)` after each accepted layer. The CLI's Rich progress bar uses this hook. |
| `blend_normals_k` | `int` | `0` | Number of layers over which to blend STL surface normals into the grid-computed normals. `0` disables. At layer `k < blend_normals_k`, blend weight is `alpha = 1 - k / blend_normals_k`, so `k=0` extrudes along pure STL normals and the blend reaches zero at `k = blend_normals_k`. **Deprecated** (ADR 0019): it showed no benefit over the grid normals; `DeprecationWarning` for `> 0`, `ParameterError` on C-grid topologies. Requires `self.body.mesh`. |

**k-step precedence**

`Mesh.march` resolves the wall-normal step schedule by precedence
(highest first):

1. The `steps` argument, when supplied.
2. Block 0's k-edges: if all 4 carry a user-overridden
   [`Spacing1D`](./api-spacing) (replaced after construction), and
   every other block agrees, use that schedule. Disagreement →
   `UserWarning`, fall back.
3. Block 0's `spacing_k` attribute set at `from_stl` time.
4. The `(ds, growth)` geometric defaults.

All blocks share one schedule because k-coupled SHARED / DIRICHLET
seams require coincident node positions layer by layer. See
[Per-edge spacing](/guide/per-edge-spacing) for how to use
precedence (2).

**Raises**

- `MarchingError` — negative-Jacobian cell detected (with `tol = -1e-6` to absorb the corner-junction noise floor).
- `NotImplementedError` — `self.topology` has no registered marcher.

---

## I/O methods

### `Mesh.write_geo`

```python
def write_geo(self, stem: str | Path) -> None
```

Write each block to `<stem>_<label>.geo`. For
`CubedSphereTopology`: 6 files (`<stem>_sub0.geo` ...
`<stem>_cap_south.geo`). For `OGridTopology`: one file
(`<stem>.geo`). For `OHTopology`: 12 files (6 wrap +
6 `<stem>_chan_<label>.geo`).

### `Mesh.write_vtk` / `Mesh.export`

```python
def write_vtk(self, stem: str | Path, with_jacobian: bool = True) -> None
def export   (self, stem: str | Path, with_jacobian: bool = True) -> None
```

Write all blocks as a single PyVista MultiBlock `.vtm` for
ParaView. `export` is an alias for `write_vtk` matching the CLI
verb. Per-cell Jacobian is attached as `CellData` by default. See
[`shore.io.vtk`](./api-vtk).

## Inspection methods

### `Mesh.info`

```python
def info(self) -> dict
```

Return mesh metadata: `n_blocks`, ordered `labels`, and a
per-block dict with `shape`, `n_points`, `n_cells`, `xyz_min`,
`xyz_max`. The `shore info` CLI is the one-line wrapper.

### `Mesh.check`

```python
def check(self, warn_threshold: float = 1e-6) -> dict
```

Per-block layer-by-layer Jacobian check. Returns
`{"ok": bool, "blocks": {label: {"jmins": [...], "any_negative": ..., "any_warn": ...}}}`.
The marcher's per-layer guard catches inversions during generation;
`check` is the independent post-hoc cross-check.

### `Mesh.view`

```python
def view(
    self, *,
    jupyter: bool = False,
    off_screen: bool = False,
    show_edges: bool = True,
    opacity: float = 0.6,
    colorby: str = "Jacobian",
) -> pv.Plotter
```

Open an interactive PyVista window with all blocks rendered. The
`shore view` CLI is the wrapper. Requires the `[vis]` extra.

## HexBlock attribute

Each entry in `Mesh.blocks` is a [`HexBlock`](./api-topology#hexblock):

| Attribute | Description |
|---|---|
| `label` | Block name (`"sub0"`, `"cap_north"`, ...). |
| `nodes` | Coordinate array, shape `(ni, nj, nk, 3)`. After `march()`, layers `1..nk-1` are filled. |
| `faces` | Six [`Face`](./api-topology#face) instances, one per face name. |
| `edges` | Twelve [`Edge`](./api-spacing#edge) instances, one per axis-parallel edge. |

## See also

- [`shore.mesh.topology`](./api-topology) — `Topology`, `CubedSphereTopology`, `OGridTopology`, `CGridTopology`, `OHTopology`, `HexBlock`, `Face`, `AxisMap`.
- [`shore.mesh.spacing`](./api-spacing) — `Spacing1D` and the per-edge data model.
- [`shore mesh` CLI](./mesh) — the command that wraps this class.
- [Hyperbolic marching](/algorithm/hyperbolic-marching) — the algorithm.
- [Cubed-sphere topology](/topologies/cubed-sphere) — block layout, seam wiring.
