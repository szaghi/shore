# shore.mesh.topology

Block-layout builders. A `Topology` takes a body surface and produces
the wired list of [`HexBlock`](#hexblock) instances that
[`Mesh.from_stl`](./api-mesh) hands to `Mesh.march()`. SHORE ships
several topologies; new ones plug into the same interface.

## Public exports

| Name | Role |
|---|---|
| `Topology` | ABC; subclass to add a new topology. |
| `CubedSphereTopology` | 6-block cubed-sphere topology — production default. |
| `OGridTopology` | Single-block lat-lon O-grid. |
| `CGridTopology` | C-grid for sharp- / blunt-TE airfoils (3-block / 4-block / legacy 1-block). |
| `OHTopology` | 12-block hybrid: cubed-sphere wrap + H-channel out to an axis-aligned bounding box. |
| `HexBlock` | Structured hexahedral block (`(ni, nj, nk, 3)` nodes + 6 faces + 12 edges + 8 vertices). |
| `Face`, `FaceBC`, `AxisMap`, `FACE_ALONG_AXES` | Face metadata used for multi-block wiring (re-exported from `shore.mesh.face`). |

## `Topology` (ABC)

```python
class Topology(ABC):
    @abstractmethod
    def build(self, surface: np.ndarray, nk: int, **kwargs) -> list[HexBlock]:
        ...
```

A topology is responsible for:

1. Allocating the per-block coordinate arrays (typically slices of a
   shared ring buffer to make SHARED seams memory-aliased).
2. Filling layer `k=0` (either by ray-cast through the body or by
   accepting a pre-projected surface).
3. Wiring [`Face`](#face) objects with the correct `bc`, `partner`,
   and [`AxisMap`](#axismap).
4. Returning the blocks in a canonical order (layer 0 of every block
   is the body; layers 1..nk-1 are filled later by `Mesh.march()`).

## `CubedSphereTopology`

```python
class CubedSphereTopology(Topology):
    def build(
        self,
        surface: np.ndarray | None,    # ignored; kept for ABC compat
        nk: int,
        ni: int | None = None,
        nj: int | None = None,
        spacing_k: Spacing1D | None = None,
        projector: object | None = None,
        theta_cap: float = 0.05 * np.pi,
        i_spacing: str = "uniform",
        i_beta: float = 3.0,
        j_spacing: str = "uniform",
        j_beta: float = 3.0,
    ) -> list[HexBlock]
```

Builds the 6-block cubed-sphere decomposition: 4 equatorial
sub-blocks (`sub0..sub3`) + 2 polar caps (`cap_north`, `cap_south`),
each cap k=0 face built by gnomonic flat-square projection (a
tangent-plane square radially projected onto the body), the equator
by slerp-meridian interpolation between the two cap seams.

Block list returned: `[sub0, sub1, sub2, sub3, cap_north, cap_south]`.

`nj` is rounded down to a multiple of 4 with a `UserWarning` if not
already so (each sub-block spans one 90° azimuthal sector).

See [Cubed-sphere topology](/topologies/cubed-sphere) for the full
algorithm, [k=0 mesh construction](/algorithm/k0-construction) for
the gnomonic seed + equator pinning, and
[Seam-aware normals](/algorithm/seam-aware-normals) for the
ghost-row mechanism.

## `OGridTopology`

```python
class OGridTopology(Topology):
    def build(
        self,
        surface: np.ndarray | None,
        nk: int,
        ni: int | None = None,
        nj: int | None = None,
        spacing_k: Spacing1D | None = None,
        projector: object | None = None,
        theta_cap: float = 0.05 * np.pi,
        i_spacing: str = "uniform",
        i_beta: float = 3.0,
        j_spacing: str = "uniform",
        j_beta: float = 3.0,
    ) -> list[HexBlock]
```

Single-block lat-lon O-grid topology. The j-direction is periodic;
i-boundaries are polar singular axis rows. Returns a single-element
list `[block]`.

`OGridTopology` and `CubedSphereTopology` are peers — the cubed
sphere is *not* a wrapped O-grid. They share only the underlying
spacing primitives. Pick the O-grid via `--topology ogrid` when you
need a single-block periodic-j sphere covering, typically as a
Chimera background.

See [O-grid topology](/topologies/o-grid).

## `OHTopology`

```python
class OHTopology(Topology):
    def __init__(
        self,
        box_vmin: tuple[float, float, float],
        box_vmax: tuple[float, float, float],
        nk_channel: int,
        k_law_channel: str = "uniform",
        k_beta_channel: float = 3.0,
        c1_seam: bool = True,
    ): ...

    def build(
        self,
        surface: np.ndarray | None,    # ignored; kept for ABC compat
        nk: int,
        ni: int | None = None,
        nj: int | None = None,
        spacing_k: Spacing1D | None = None,
        projector: object | None = None,
        i_spacing: str = "uniform",
        i_beta: float = 3.0,
        j_spacing: str = "uniform",
        j_beta: float = 3.0,
        **kwargs,
    ) -> list[HexBlock]
```

12-block conforming hybrid: a 6-block cubed-sphere wrap (marched
outward from the STL) plus a 6-block H-channel (frustum blocks
bridging the wrap envelope to the user-supplied axis-aligned
bounding box). Wrap and channel are SHARED-by-view at the seam —
one contiguous numpy buffer per pair, bit-exact with no copy.

Block list returned: `[sub0, sub1, sub2, sub3, cap_north, cap_south,
chan_sub0, chan_sub1, chan_sub2, chan_sub3, chan_north, chan_south]`.

Each wrap sub-block faces exactly one box face. OH pins
`theta_cap = pi/4` and `phi_offset = pi/4` internally regardless of
the user's `theta_cap` argument; this is what makes the 1:1
wrap-to-channel mapping geometrically possible (cap-equator seams
fall on cube edges, not inside cube faces).

Channel block interiors are filled **after marching**, by
`OHTopology.fill_channel_blocks(blocks)`. `Mesh.march()` calls this
automatically when the topology is `OHTopology`; users who drive
marching manually must call it themselves before using the channel
geometry.

When `c1_seam=True` (default), the channel's first cell at every
`(i, j)` is pinned to the wrap's last marching step at the same
`(i, j)`, and the remaining channel cells grow geometrically out to
the box face — bit-exact C1 spacing across the seam, no second-cell
jump. The `k_law_channel` / `k_beta_channel` arguments apply only
when `c1_seam=False`.

See [OH hybrid topology](/topologies/oh-hybrid) for the geometry,
[`shore.topology.oh_channel`](https://github.com/szaghi/shore/blob/main/src/shore/topology/oh_channel.py)
for the channel-block builder, and
[`examples/08-oh-hybrid/`](https://github.com/szaghi/shore/tree/main/examples/08-oh-hybrid)
for a runnable demo.

## `HexBlock`

```python
class HexBlock(Block):
    label:       str
    nodes:       np.ndarray  # (ni, nj, nk, 3)
    periodic_i:  bool        # default False
    periodic_j:  bool        # default True (changed by topologies that wire seams)
```

| Method / property | Returns |
|---|---|
| `block.face(name)` | The [`Face`](#face) named e.g. `"i_lo"` / `"k_hi"`. |
| `block.faces` | All 6 faces, in canonical order. |
| `block.edge(name)` | The [`Edge`](./api-spacing#edge) named e.g. `"j_lo_k_hi"`. |
| `block.edges` | All 12 edges. |
| `block.vertices` | All 8 vertices. |

The `nodes` array is owned externally — typically a slice of the
topology's ring allocation, so neighbouring blocks alias their
shared faces by numpy view.

## `Face`

```python
@dataclass
class Face:
    index:           int
    name:            FaceName        # 'i_lo' | 'i_hi' | 'j_lo' | 'j_hi' | 'k_lo' | 'k_hi'
    bc:              FaceBC
    partner:         Face | None     # back-reference, None for WALL/FREE
    axis_map:        AxisMap | None  # required when partner is set
    body_projected:  bool            # True iff bc == WALL
```

`Face.nodes(block)` returns a 2D numpy view of the face's node slab
(read/write). `Face.partner` carries a back-reference to the
partnered face on another block; for SHARED seams the two `nodes`
views alias the same underlying memory.

## `FaceBC`

```python
class FaceBC(enum.Enum):
    WALL      = "wall"        # k=0 body surface, Dirichlet from projection
    FREE      = "free"        # outer boundary, marched freely
    SHARED    = "shared"      # shared with another block via numpy view
    DIRICHLET = "dirichlet"   # pinned from a partner face (cap boundaries)
    PERIODIC  = "periodic"    # wraps to another face on the same block
```

The downstream consumers (cc.par writer, adjacency JSON, debug
viewers) dispatch on this enum. See
[`shore cc-par` BC mapping](./cc-par#shore--cc-par-bc-mapping) for
the Xall translation.

## `AxisMap`

```python
@dataclass(frozen=True)
class AxisMap:
    along: tuple[
        tuple[AxisLabel, bool],   # (partner_axis, flipped)
        tuple[AxisLabel, bool],
    ]
```

Stamped on every connection face (`SHARED` / `DIRICHLET` /
`PERIODIC`) at wiring time. Describes how this face's two
along-seam axes (defined by `FACE_ALONG_AXES[face_name]`) map onto
the partner face's along-seam axes — by axis label and per-axis
flip flag.

| Constructor | Use |
|---|---|
| `AxisMap.identity(face_name)` | Same-axis, no-flip mapping for `face_name`. The most common case (same-axis SHARED seams). |
| `AxisMap.from_pairs(first, second)` | Construct from two `(axis, flipped)` pairs. |

| Method | Returns |
|---|---|
| `axis_map.partner_axes()` | The two partner along-axis labels. |
| `axis_map.flips()` | The two flip flags. |
| `axis_map.inverse(this_face, partner_face)` | The partner-side `AxisMap` for the same seam. |
| `axis_map.validate_against_partner(this_face, partner_face)` | Raises `ValueError` if the partner-axis labels are not exactly `partner_face`'s along-axes. |

The cc.par writer derives the `cc.par` orientation digit (e.g.
`135` / `145` / `235` / `415`) from `axis_map` + `partner.name`,
without any per-topology lookup table. Any topology that stamps
valid `AxisMap` values on its seams is fully covered by the
existing writer.

## `FACE_ALONG_AXES`

```python
FACE_ALONG_AXES: dict[FaceName, tuple[AxisLabel, AxisLabel]] = {
    "i_lo": ("j", "k"),
    "i_hi": ("j", "k"),
    "j_lo": ("i", "k"),
    "j_hi": ("i", "k"),
    "k_lo": ("i", "j"),
    "k_hi": ("i", "j"),
}
```

The canonical along-axis order for each face name. Used by
`AxisMap.identity` and as the order the two `AxisMap.along` pairs
are interpreted in. New topologies should respect this ordering when
constructing their seam mappings.

## See also

- [`shore.mesh.Mesh`](./api-mesh) — the consumer of the block list.
- [`shore.mesh.spacing`](./api-spacing) — `Edge`, `Spacing1D`.
- [Cubed-sphere topology](/topologies/cubed-sphere) — block layout.
- [O-grid topology](/topologies/o-grid) — single-block alternative.
- [OH hybrid topology](/topologies/oh-hybrid) — wrap + channel hybrid.
- [Topology overview](/topologies/) — taxonomy and how to choose.
