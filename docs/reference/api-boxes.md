# shore.boxes

`boxes_from_stl` — generate marker [`Box`](./api-boxes-io)
instances from a body STL. Three methods: `voxel-fill` (default,
guaranteed body cover), `aabb` and `obb` (single-box, best-effort).

For the algorithm derivation, see
[Algorithm — marker boxes](/algorithm/boxes). For the CLI and
ParaView visualisation, see [Marker boxes](./boxes).

## `boxes_from_stl`

```python
def boxes_from_stl(
    stl_path: str | Path,
    associated_block: str,
    *,
    max_outward_gap: float,
    method: Literal["aabb", "obb", "voxel-fill"] = "voxel-fill",
    face_samples_per_axis: int = 4,
) -> list[Box]
```

Build a list of [`Box`](./api-boxes-io) instances covering the
body. All boxes share `associated_block`.

| Parameter | Description |
|---|---|
| `stl_path` | Path to the body's STL file. Loaded via [`load_stl`](./api-surface-io). |
| `associated_block` | SHORE block label every returned box references. Resolved to a 1-based index by [`write_cc_par`](./api-cc-par). |
| `max_outward_gap` | Maximum allowed distance, in STL world units, from any point on a box's surface to the nearest STL surface point — measured **outward only** (points on the box surface that lie inside the body don't count). Strictly positive. Should be smaller than the body's near-mesh wall-normal thickness so the boxes fit inside the wall. |
| `method` | `"voxel-fill"` (default), `"aabb"`, or `"obb"`. See [Methods](#methods). |
| `face_samples_per_axis` | `aabb` / `obb` post-check sampling density. Each face is sampled on a `k × k` grid (8 corners + `6 * k²` samples per box) for the outward-gap measurement. Default 4. Ignored by `voxel-fill`. |

**Returns** `list[Box]`. Length is 1 for `aabb` / `obb`; many for
`voxel-fill` (post-merge count tracks body surface area).

**Raises**

- `ParameterError` — invalid parameters or unknown `method`.
- `MeshInputError` — STL malformed; or `voxel-fill` requested on a
  non-watertight STL (single-box methods don't need watertightness).

## Methods

### `"voxel-fill"` (default)

Voxelise the body's AABB at side `max_outward_gap / sqrt(3)`. Include
every voxel that overlaps the body
(`signed_distance(centre) > -half_diag`) — both centre-inside and
body-shell voxels. Greedy-merge face-adjacent grid blocks into
slab-shaped AABBs whenever the merged block's 8 corners still
satisfy `outward_gap <= max_outward_gap`.

**Coverage guaranteed** — the body interior and surface are tiled
by the box union by construction. Box vertex outward gap to the
body is bounded by `max_outward_gap`. Watertight STL required.

### `"aabb"` (best effort)

Single world-axis-aligned bounding box of the body. Tight, but the
bounding cube's corners typically stick out by `r * (sqrt(3) - 1)`
past a circumscribed sphere of radius `r` — much further than
typical near-mesh thickness budgets. Emits a `UserWarning` when
the measured outward gap exceeds the budget.

Useful for quick previews on body-fitted-cube bodies or when the
budget is very loose.

### `"obb"` (best effort)

Single oriented bounding box from area-weighted PCA. Tight, with
the same best-effort warning behaviour as `"aabb"`. Better than
`"aabb"` for elongated bodies (the OBB aligns with the body's
principal axis).

## Example

```python
from shore.boxes import boxes_from_stl
from shore.io.cc_par import write_cc_par

boxes = boxes_from_stl(
    "sphere.stl",
    associated_block="sub0",
    max_outward_gap=0.30,
    method="voxel-fill",         # default
)
print(len(boxes))   # several dozen for a unit-radius icosphere

write_cc_par(
    "wall.adjacency.json",
    "wall.cc.par",
    grd_filename="wall.grd",
    boxes=boxes,
)
```

For sharply concave geometries with `aabb` / `obb`, raise
`face_samples_per_axis` to 6–8 for a stricter outward-gap
post-check.

## See also

- [Marker boxes](./boxes) — full reference (algorithm, CLI,
  ParaView).
- [`shore boxes-from-stl` CLI](./boxes#generating-boxes-from-an-stl-shore-boxes-from-stl) —
  command-line wrapper.
- [`shore.io.boxes`](./api-boxes-io) — the `Box` dataclass.
- [`shore.io.cc_par.write_cc_par`](./api-cc-par) — bake boxes into
  `cc.par`.
- [Algorithm — marker boxes](/algorithm/boxes) — voxel-fill +
  greedy-merge derivation.
