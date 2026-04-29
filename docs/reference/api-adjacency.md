# shore.io.adjacency

Read / write the JSON sidecar that captures a multi-block adjacency
graph. SHORE's `.geo` format is purely geometric (node coordinates,
no connectivity); the adjacency JSON serialises every block's face
BCs, partner references, and `axis_map` orientation so that the
downstream chimera writer ([`shore cc-par`](./cc-par)) can produce
the right connection patches.

For the JSON schema and version semantics, see the
[`shore split` reference — Adjacency JSON](./split#adjacency-json).

## Public exports

| Name | Role |
|---|---|
| `write_adjacency_json` | Serialise a list of `HexBlock` to JSON. |
| `read_adjacency_json` | Parse + version-check; returns the raw payload dict. |
| `summarise_adjacency` | Format a one-screen summary of an adjacency payload (used by `shore split -v`). |
| `ADJACENCY_VERSION` | Current schema version (`2`). |

## `write_adjacency_json`

```python
def write_adjacency_json(
    blocks: list[HexBlock],
    path: str | Path,
) -> Path
```

Walk every block's 6 faces and emit a JSON file capturing each
face's `bc`, `partner` (label + face), and `axis_map` (canonical
2-pair representation per [`AxisMap`](./api-topology#axismap)).

The output schema:

```json
{
  "version": 2,
  "blocks": [
    {
      "label": "sub0__k0",
      "shape": [20, 9, 5],
      "faces": {
        "i_lo": {
          "bc": "DIRICHLET",
          "partner": {"label": "cap_south", "face": "j_lo"},
          "axis_map": [["i", false], ["k", false]]
        },
        "k_lo": {"bc": "WALL", "partner": null, "axis_map": null},
        ...
      }
    },
    ...
  ]
}
```

`axis_map` is a 2-element list of `[partner_axis, flipped]` pairs
in the canonical along-axis order for *this face's* name (see
[`FACE_ALONG_AXES`](./api-topology#face-along-axes)). Faces without
a partner (WALL, FREE) carry `null` for both `partner` and
`axis_map`.

## `read_adjacency_json`

```python
def read_adjacency_json(path: str | Path) -> dict[str, Any]
```

Parse and version-check an adjacency JSON file. Returns the raw
parsed dict (suitable for hand-inspection or for passing into
`write_cc_par`).

**Raises**

- `MeshInputError` — missing file, malformed JSON, missing
  required keys, or version mismatch.

**Versioning.** This build emits and requires version 2. Version 1
sidecars are rejected with a clear regenerate-the-sidecar message —
not a soft fallback, because the v1 schema lacked `axis_map` and
the orientation cannot be recovered from labels alone. Sidecars
are derived artefacts (re-run the topology / split that produced
them).

## `summarise_adjacency`

```python
def summarise_adjacency(payload: dict[str, Any]) -> str
```

Format a one-screen human-readable summary: per-block face BC
counts and a sample of the SHARED / DIRICHLET partners. Used by
`shore split --verbose`.

## Example

```python
from shore.io.adjacency import write_adjacency_json, read_adjacency_json
from shore.split import BlockSplit, split_blocks

# Subdivide and write the sidecar
new_blocks = split_blocks(mesh.blocks, [BlockSplit("sub0", "k", at=[10])])
write_adjacency_json(new_blocks, "mesh.adjacency.json")

# Re-read and feed downstream
payload = read_adjacency_json("mesh.adjacency.json")
print(payload["version"], len(payload["blocks"]))
```

## See also

- [`shore split` CLI](./split) — produces the sidecar.
- [`shore.split`](./api-split) — Python API for the split.
- [`shore cc-par`](./cc-par) — the downstream consumer.
- [`shore.mesh.face.AxisMap`](./api-topology#axismap) — what the
  `axis_map` JSON pair represents.
