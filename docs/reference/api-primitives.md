# shore.primitive

Parametric structured-grid primitives — `BoxMesh`,
`AnnularCylinderMesh`, `FlatCapsCylinder`. Each is a `PrimitiveMesh`
subclass: it owns a list of [`HexBlock`](./api-topology#hexblock)
instances and inherits `write_geo` / `write_vtk` from the base.
Construction is driven by a single `build(...)` classmethod per
primitive; no STL projection, no marching.

For the conceptual overview, see
[Primitives](/topologies/primitives). For the CLI wrappers, see
[`shore primitive`](./primitive).

## `PrimitiveMesh` (base)

```python
class PrimitiveMesh(ABC):
    blocks:        list[HexBlock]
    label_default: str = "primitive"
```

| Method | Description |
|---|---|
| `mesh.write_geo(output)` | Single-block primitives write `output` verbatim (one `.geo`). Multi-block primitives write `output_<label>.geo` per block. |
| `mesh.write_vtk(output, with_jacobian=True)` | Write a single PyVista MultiBlock `.vtm` for ParaView, regardless of block count. |

## `BoxMesh`

Single-block axis-aligned Cartesian box.

```python
class BoxMesh(PrimitiveMesh):
    label_default = "box"

    @classmethod
    def build(
        cls,
        vmin: tuple[float, float, float],
        vmax: tuple[float, float, float],
        ni: int, nj: int, nk: int,
        i_law: str = "uniform",
        j_law: str = "uniform",
        k_law: str = "uniform",
        i_beta: float = 3.0,
        j_beta: float = 3.0,
        k_beta: float = 3.0,
        edge_pins: dict[str, Spacing1D] | None = None,
        label: str = "box",
    ) -> BoxMesh
```

| Parameter | Description |
|---|---|
| `vmin`, `vmax` | World-space corners. Each component of `vmax` must strictly exceed `vmin`. |
| `ni`, `nj`, `nk` | Node counts along the three axes. All `>= 2`. |
| `i_law` / `j_law` / `k_law` | Default per-axis spacing law: `"uniform"` / `"tanh"` / `"tanh2"`. |
| `i_beta` / `j_beta` / `k_beta` | Tanh / tanh2 clustering strength per axis. |
| `edge_pins` | Optional sparse override `{edge_name: Spacing1D}`. When any edge of an axis is pinned, the axis switches to **transfinite interpolation** between the 4 parallel edges' distributions; the other axes stay tensor-product. See [Per-edge spacing](/guide/per-edge-spacing). |
| `label` | Block label (default `"box"`). |

Tensor-product example:

```python
from shore.primitive import BoxMesh

box = BoxMesh.build(
    vmin=(0.0, 0.0, 0.0),
    vmax=(10.0, 5.0, 3.0),
    ni=40, nj=30, nk=20,
    i_law="tanh", i_beta=4.0,    # cluster cells at the inflow plane
)
box.write_geo("channel.geo")
```

Per-edge override:

```python
from shore.mesh.spacing import Spacing1D

box = BoxMesh.build(
    vmin=(0.0, 0.0, 0.0), vmax=(10.0, 5.0, 3.0),
    ni=40, nj=30, nk=20,
    edge_pins={"j_lo_k_lo": Spacing1D(law="tanh", beta=4.0)},
    # i-axis interior is now a TFI blend of the 4 parallel i-edges,
    # only one of which is tanh-clustered.
)
```

## `AnnularCylinderMesh`

Single-block annular cylinder, optionally a partial sector.

```python
class AnnularCylinderMesh(PrimitiveMesh):
    label_default = "annulus"

    @classmethod
    def build(
        cls,
        r_in: float,
        r_out: float,
        z0: float,
        z1: float,
        ni: int, nj: int, nk: int,
        theta_start: float = 0.0,
        theta_end:   float | None = None,
        axis:   str | tuple[float, float, float] = "z",
        center: tuple[float, float, float] = (0.0, 0.0, 0.0),
        i_law: str = "uniform",
        j_law: str = "uniform",
        k_law: str = "uniform",
        i_beta: float = 3.0,
        j_beta: float = 3.0,
        k_beta: float = 3.0,
        edge_pins: dict[str, Spacing1D] | None = None,
        label: str = "annulus",
    ) -> AnnularCylinderMesh
```

**Index meaning:**

- `i` (axis 0): radial, `r_in → r_out`
- `j` (axis 1): azimuthal, `theta_start → theta_end`
- `k` (axis 2): axial, along `axis`

| Parameter | Description |
|---|---|
| `r_in`, `r_out` | Inner and outer radii (`r_out > r_in >= 0`). |
| `z0`, `z1` | Axial extent along `axis` (`z1 > z0`). |
| `theta_start` | Sector start in radians. |
| `theta_end` | Sector end in radians. **`None`** → full annulus (j-periodic). |
| `axis` | `"x"` / `"y"` / `"z"` (case-insensitive) or a 3-vector. |
| `center` | World-space anchor of the axis line. |
| `i_law` / `j_law` / `k_law`, `*_beta` | Per-axis spacing laws and strengths. |
| `edge_pins` | Per-edge override (Box-style). For a full annulus, j-edges (`i_lo_k_lo`, `i_hi_k_lo`, `i_lo_k_hi`, `i_hi_k_hi`) are **rejected** because j is periodic and per-edge endpoint pinning is meaningless on a periodic axis. |

Full annulus → uniform j is enforced; tanh / tanh2 on `j_law` is
rejected with `ParameterError`. Use a partial sector when you need
azimuthal clustering.

```python
from shore.primitive import AnnularCylinderMesh

# Full annulus around a near-body sphere, default uniform spacing
ann = AnnularCylinderMesh.build(
    r_in=1.5, r_out=5.0, z0=-2.0, z1=2.0,
    ni=30, nj=64, nk=40,
)
ann.blocks[0].periodic_j   # True
```

## `FlatCapsCylinder`

5-block flat-caps cylinder. Disk interior at every axial slice is
decomposed into a central square + 4 trapezoidal sub-blocks
(`center`, `sub_e`, `sub_n`, `sub_w`, `sub_s`). No central hole, no
periodic j — Chimera-fringe-free background.

```python
class FlatCapsCylinder(PrimitiveMesh):
    label_default = "flat_caps"

    @classmethod
    def build(
        cls,
        r_out: float,
        z0: float,
        z1: float,
        ni: int,
        nc: int,
        nk: int,
        inner_frac: float = 0.5,
        axis:   str | tuple[float, float, float] = "z",
        center: tuple[float, float, float] = (0.0, 0.0, 0.0),
        i_law: str = "uniform",
        c_law: str = "uniform",
        k_law: str = "uniform",
        i_beta: float = 3.0,
        c_beta: float = 3.0,
        k_beta: float = 3.0,
        edge_pins: dict[str, Spacing1D] | None = None,
    ) -> FlatCapsCylinder
```

| Parameter | Description |
|---|---|
| `r_out` | Outer disk radius. |
| `z0`, `z1` | Axial extent (`z1 > z0`). |
| `ni` | Radial node count for each sub-block (`>= 2`). |
| `nc` | Tangential node count along each side of the central square (`>= 2`). |
| `nk` | Axial node count (`>= 2`). |
| `inner_frac` | Ratio `r_inner / r_out` of the central square's corner-radius to the outer disk radius (`0 < inner_frac < 1`, default `0.5`). |
| `axis`, `center` | Cylinder axis and world-space anchor. |
| `i_law`, `c_law`, `k_law`, `*_beta` | Radial / tangential / axial spacing laws and strengths. |
| `edge_pins` | **Only k-axis edges accepted.** Pinning any i- or c-axis edge raises `ParameterError`: the radial and tangential distributions are pinned by the geometric construction (per-j first-step matching at the central-square seam) and overriding them would require rewriting that geometry. K-edge pins switch the axial coordinate to TFI between the 4 parallel k-edges of every block, allowing different axial spacing on the inner radius vs the outer circumference at the same k-layer. |

```python
from shore.primitive import FlatCapsCylinder

fc = FlatCapsCylinder.build(
    r_out=10.0, z0=-10.0, z1=50.0,
    ni=20, nc=24, nk=40,
)
[b.label for b in fc.blocks]
# ['center', 'sub_e', 'sub_n', 'sub_w', 'sub_s']

fc.write_geo("background")    # background_center.geo, background_sub_e.geo, ...
fc.write_vtk("background")    # background.vtm
```

## See also

- [`shore primitive` CLI](./primitive) — the CLI wrappers.
- [Primitives overview](/topologies/primitives) — when to use each.
- [`shore.mesh.spacing`](./api-spacing) — `Spacing1D`.
- [`shore.mesh.topology`](./api-topology) — `HexBlock`.
- [Per-edge spacing](/guide/per-edge-spacing) — TFI walkthrough.
