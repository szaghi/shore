# shore.io.boxes

The `Box` data class — one Xall `cc.par` type-9 marker box — plus
JSON serialisation helpers. Used by
[`shore.io.cc_par.write_cc_par(boxes=...)`](./api-cc-par) and the
`shore cc-par --boxes` CLI.

For the cc.par boxes section format, the algorithm that generates
boxes from an STL, and the ParaView visualiser, see
[Marker boxes](./boxes).

## `Box`

```python
@dataclass(frozen=True)
class Box:
    vertices:          np.ndarray      # (8, 3) float64
    associated_block:  str
```

One type-9 hexahedron with 8 corner vertices in the canonical
"i-fastest" Fortran-column-major order:

| Index | `(i, j, k)` |
|---|---|
| 0 (`v1`) | `(lo, lo, lo)` |
| 1 (`v2`) | `(hi, lo, lo)` |
| 2 (`v3`) | `(lo, hi, lo)` |
| 3 (`v4`) | `(hi, hi, lo)` |
| 4 (`v5`) | `(lo, lo, hi)` |
| 5 (`v6`) | `(hi, lo, hi)` |
| 6 (`v7`) | `(lo, hi, hi)` |
| 7 (`v8`) | `(hi, hi, hi)` |

The `(i, j, k)` directions need not be world-aligned — any rotated
or sheared parallelepiped is a valid `Box` as long as the corner
correspondence above is preserved.

`Box.__post_init__` validates the hexahedron at construction:
non-finite values, wrong shape, empty `associated_block`,
degenerate / coplanar / self-intersecting (bow-tie) / inverted
geometry all raise `ParameterError`.

`associated_block` is a SHORE block label; at write time the
[`cc.par` writer](./api-cc-par) resolves it to a 1-based index
into the blocks list. The box then inherits the block's
`level / group / priority` for chimera precedence.

### `Box.from_aabb`

```python
@classmethod
def from_aabb(
    cls,
    vmin: tuple[float, float, float] | np.ndarray,
    vmax: tuple[float, float, float] | np.ndarray,
    associated_block: str,
) -> Box
```

Convenience constructor for an axis-aligned `Box`. Builds the 8
canonical-order vertices from the two diagonal corners of an AABB.

```python
from shore.io.boxes import Box

box = Box.from_aabb(
    vmin=(-0.6, -0.6, -0.6),
    vmax=(+0.6, +0.6, +0.6),
    associated_block="sub0",
)
```

## JSON helpers

```python
def box_to_json(box: Box) -> dict[str, Any]
def box_from_json(payload: dict[str, Any]) -> Box
```

Round-trip a `Box` through the JSON schema that `shore cc-par
--boxes` consumes:

```json
{
  "associated_block": "sub0",
  "vertices": [[-0.6, -0.6, -0.6], [+0.6, -0.6, -0.6], ...]
}
```

Lists of boxes are written as a JSON array of these objects; that
is the exact format `shore boxes-from-stl` emits.

## See also

- [Marker boxes](./boxes) — full reference (algorithm, CLI, ParaView).
- [`shore.boxes`](./api-boxes) — `boxes_from_stl` (STL → list of
  `Box`).
- [`shore.io.cc_par`](./api-cc-par) — `write_cc_par(boxes=...)`.
- [Algorithm — marker boxes](/algorithm/boxes) — voxel-fill +
  greedy-merge derivation.
