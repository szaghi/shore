---
title: shore.io.adjacency
outline: [2, 3]
---

# `shore.io.adjacency`

Adjacency JSON sidecar for split block lists.

The .geo format SHORE writes is purely geometric — node coordinates,
no connectivity beyond a single block.  Once a topology has been split
into chunks, the downstream solver needs to know which face of which
block is SHARED with which face of which other block, which faces are
WALL / FREE / DIRICHLET, etc.  This module reads and writes that
graph as a JSON sidecar.

Format
------

The file describes a list of named blocks; each block carries its
shape and its 6 face boundary conditions.  Faces with a partner
(SHARED, DIRICHLET, PERIODIC) record the partner's `label`,
`face` name, and an `axis_map` describing how this face's
along-seam axes map onto the partner face's along-seam axes (see
:class:`shore.mesh.face.AxisMap`).  Faces without a partner (WALL,
FREE) carry `null` for both `partner` and `axis_map`.

The `axis_map` field is a 2-element list of `[partner_axis, flipped]`
pairs in the canonical along-axis order for *this face's* name (see
:data:`shore.mesh.face.FACE_ALONG_AXES`).

Example:

```
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
        "i_hi": {
          "bc": "DIRICHLET",
          "partner": {"label": "cap_north", "face": "i_lo"},
          "axis_map": [["j", false], ["k", false]]
        },
        "j_lo": {
          "bc": "SHARED",
          "partner": {"label": "sub3__k0", "face": "j_hi"},
          "axis_map": [["i", false], ["k", false]]
        },
        "j_hi": {
          "bc": "SHARED",
          "partner": {"label": "sub1__k0", "face": "j_lo"},
          "axis_map": [["i", false], ["k", false]]
        },
        "k_lo": {"bc": "WALL", "partner": null, "axis_map": null},
        "k_hi": {
          "bc": "SHARED",
          "partner": {"label": "sub0__k1", "face": "k_lo"},
          "axis_map": [["i", false], ["j", false]]
        }
      }
    },
    ...
  ]
}

```
Versioning
----------

Version 3 adds block-to-geometry associations (ADR 0003). A top-level
`associations` table lists each distinct association once (`geometry`
key, `kind` and entity `keys`); each block may carry an `associated`
map from its faces, edges and corner indices to entries of that table:

```
"associations": [{"geometry": "body", "kind": "body", "keys": []}],
"blocks": [{..., "associated": {
    "faces": {"k_lo": 0},
    "edges": {"j_lo_k_lo": 0, "j_hi_k_lo": 0, "i_lo_k_lo": 0, "i_hi_k_lo": 0},
    "vertices": {"0": 0, "2": 0, "4": 0, "6": 0}}}]

```
Associations are optional: a sidecar without any differs from version 2
only in its version number. Readers accept versions 2 and 3.

Version 1 (deprecated) stored only `{label, face}` partners with no
axis_map.  This proved insufficient for any consumer that needs to
write a connection orientation (cc.par, overset connectivity tools)
because the axis correspondence between two faces of a topology is not
recoverable from labels alone — see the bug fixed in v2.  Files with
`version: 1` are rejected with a clear regenerate-the-sidecar
message; this is acceptable because adjacency sidecars are derived
artefacts, not source files.

## `write_adjacency_json` <Badge type="info" text="function" /> {#write_adjacency_json}

```python
write_adjacency_json(blocks: list[HexBlock], path: str | Path) -> Path
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/io/adjacency.py#L112">source</a>

Write the adjacency graph of *blocks* to a JSON sidecar.

**Parameters**

- **`blocks`** — List of HexBlocks with face BCs already wired (typically the output of :func:`shore.split.split_blocks` or a topology builder).
- **`path`** — Output path; `.json` is appended if absent.

**Returns**

- `Path` — The actual path written.

## `read_adjacency_json` <Badge type="info" text="function" /> {#read_adjacency_json}

```python
read_adjacency_json(path: str | Path) -> dict[str, Any]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/io/adjacency.py#L274">source</a>

Read an adjacency JSON sidecar and return its raw dict.

Validates only the version field and the top-level `blocks` key.
The caller is responsible for matching the sidecar against actual
`.geo` files.

**Raises**

- `MeshInputError` — On missing file, malformed JSON, missing `version` / `blocks` keys, or unsupported version.

## `summarise_adjacency` <Badge type="info" text="function" /> {#summarise_adjacency}

```python
summarise_adjacency(payload: dict[str, Any]) -> str
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/io/adjacency.py#L361">source</a>

Return a one-block-per-line human summary string for printing.

Useful for the `shore split` CLI's verbose output.  Walks the
parsed dict (does not need real HexBlocks).

## `wire_blocks_from_adjacency` <Badge type="info" text="function" /> {#wire_blocks_from_adjacency}

```python
wire_blocks_from_adjacency(
    blocks: list[HexBlock],
    payload: dict[str, Any],
    *,
    source: str = '<adjacency>',
) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/io/adjacency.py#L416">source</a>

Wire the face BCs of *blocks* in place from a parsed adjacency *payload*.

`.geo` files carry only geometry, so :class:`HexBlock`s built from
them default every face to WALL (k_lo) / FREE.  This restores the real
boundary graph (SHARED / DIRICHLET / PERIODIC seams with their partner
pointers and axis_maps, plus any per-face `bc_code` overrides) onto
those blocks, so a subsequent :func:`shore.split.split_blocks` can
preserve the pre-existing inter-block seams instead of losing them.

Every block in *payload* must resolve to exactly one block in *blocks*
via :func:`_resolve_label` (exact label, else unique `_<label>`
suffix).  A sidecar block that matches zero or several input blocks is a
hard :class:`MeshInputError` — silently skipping a mismatch is what let
a cubed-sphere's DIRICHLET/SHARED seams silently vanish in a combined
multi-component split.

Raises :class:`MeshInputError` on an unknown FaceBC name, a malformed
axis_map, or an unresolved / ambiguous block or partner label.

## `apply_associations` <Badge type="info" text="function" /> {#apply_associations}

```python
apply_associations(
    block: HexBlock,
    block_payload: dict[str, Any],
    table: list[Association],
) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/io/adjacency.py#L471">source</a>

Restore a block's associations from its sidecar entry (version 3).
