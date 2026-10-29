"""Adjacency JSON sidecar for split block lists.

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
(SHARED, DIRICHLET, PERIODIC) record the partner's ``label``,
``face`` name, and an ``axis_map`` describing how this face's
along-seam axes map onto the partner face's along-seam axes (see
:class:`shore.mesh.face.AxisMap`).  Faces without a partner (WALL,
FREE) carry ``null`` for both ``partner`` and ``axis_map``.

The ``axis_map`` field is a 2-element list of ``[partner_axis, flipped]``
pairs in the canonical along-axis order for *this face's* name (see
:data:`shore.mesh.face.FACE_ALONG_AXES`).

Example::

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

Versioning
----------

Version 3 adds block-to-geometry associations (ADR 0003). A top-level
``associations`` table lists each distinct association once (``geometry``
key, ``kind`` and entity ``keys``); each block may carry an ``associated``
map from its faces, edges and corner indices to entries of that table::

    "associations": [{"geometry": "body", "kind": "body", "keys": []}],
    "blocks": [{..., "associated": {
        "faces": {"k_lo": 0},
        "edges": {"j_lo_k_lo": 0, "j_hi_k_lo": 0, "i_lo_k_lo": 0, "i_hi_k_lo": 0},
        "vertices": {"0": 0, "2": 0, "4": 0, "6": 0}}}]

Associations are optional: a sidecar without any differs from version 2
only in its version number. Readers accept versions 2 and 3.

Version 1 (deprecated) stored only ``{label, face}`` partners with no
axis_map.  This proved insufficient for any consumer that needs to
write a connection orientation (cc.par, overset connectivity tools)
because the axis correspondence between two faces of a topology is not
recoverable from labels alone — see the bug fixed in v2.  Files with
``version: 1`` are rejected with a clear regenerate-the-sidecar
message; this is acceptable because adjacency sidecars are derived
artefacts, not source files.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from shore.exceptions import MeshInputError
from shore.mesh.association import Association
from shore.mesh.block import HexBlock
from shore.mesh.face import FACE_ALONG_AXES, AxisMap, FaceBC, FaceName

ADJACENCY_VERSION = 3
_READABLE_VERSIONS = (2, 3)

_FACE_NAMES: tuple[FaceName, ...] = ("i_lo", "i_hi", "j_lo", "j_hi", "k_lo", "k_hi")


def write_adjacency_json(blocks: list[HexBlock], path: str | Path) -> Path:
    """Write the adjacency graph of *blocks* to a JSON sidecar.

    Parameters
    ----------
    blocks :
        List of HexBlocks with face BCs already wired (typically the
        output of :func:`shore.split.split_blocks` or a topology
        builder).
    path :
        Output path; ``.json`` is appended if absent.

    Returns
    -------
    Path
        The actual path written.
    """
    out = Path(path)
    if out.suffix.lower() != ".json":
        out = out.with_suffix(".json")

    # Build a label-keyed lookup so partner-face references can carry
    # both labels (the partner Face object holds no back-reference to
    # its parent block in the runtime objects, so we walk every block).
    by_face_id: dict[int, tuple[str, FaceName]] = {}
    for blk in blocks:
        for f in blk.faces():
            by_face_id[id(f)] = (blk.label, f.name)

    table: list[Association] = []
    index: dict[Association, int] = {}

    def assoc_id(a: Association) -> int:
        if a not in index:
            index[a] = len(table)
            table.append(a)
        return index[a]

    patches = _shared_face_patches(blocks)
    payload: dict[str, Any] = {"version": ADJACENCY_VERSION, "blocks": []}
    for blk in blocks:
        ni, nj, nk, _ = blk.nodes.shape
        face_payload: dict[str, dict[str, Any]] = {}
        for f in blk.faces():
            partner_obj: dict[str, str] | None = None
            if f.partner is not None:
                p_lookup = by_face_id.get(id(f.partner))
                if p_lookup is not None:
                    partner_obj = {"label": p_lookup[0], "face": p_lookup[1]}
            axis_map_obj = _axis_map_to_json(f.axis_map)
            face_obj: dict[str, Any] = {
                "bc": f.bc.name,
                "partner": partner_obj,
                "axis_map": axis_map_obj,
            }
            # Optional per-face natural-BC override (FREE faces only).  Emit
            # only when set so existing sidecars round-trip byte-for-byte.
            if f.bc_code is not None:
                face_obj["bc_code"] = f.bc_code
            # A face shared with several partners (ADR 0036): one sub-range per
            # partner. Emitted only for such faces, so other sidecars are unchanged.
            if (blk.label, f.name) in patches:
                face_obj["patches"] = patches[(blk.label, f.name)]
            face_payload[f.name] = face_obj
        block_obj: dict[str, Any] = {
            "label": blk.label,
            "shape": [ni, nj, nk],
            "faces": face_payload,
        }
        associated = {
            "faces": {f.name: assoc_id(f.association) for f in blk.faces() if f.association},
            "edges": {e.name: assoc_id(e.association) for e in blk.edges() if e.association},
            "vertices": {
                str(v.index): assoc_id(v.association) for v in blk.vertices() if v.association
            },
        }
        if any(associated.values()):
            block_obj["associated"] = {k: v for k, v in associated.items() if v}
        payload["blocks"].append(block_obj)

    if table:
        payload = {
            "version": ADJACENCY_VERSION,
            "associations": [a.to_json() for a in table],
            "blocks": payload["blocks"],
        }
    out.write_text(json.dumps(payload, indent=2) + "\n")
    return out


def _shared_face_patches(blocks: list[HexBlock]) -> dict[tuple[str, str], list[dict[str, Any]]]:
    """Sub-range patches of the faces with several partners (ADR 0036).

    Each patch: ``range`` (0-based node indices ``[[i0, j0, k0], [i1, j1,
    k1]]``), ``partner`` and ``axis_map`` (this face's along-axes onto the
    partner face's), located by the faces' nodes.
    """
    from shore.mesh.seams import AXES, seams, shared_faces

    shared = shared_faces(blocks)
    if not shared:
        return {}
    points = np.concatenate([b.nodes.reshape(-1, 3) for b in blocks])
    scale = float(np.linalg.norm(points.max(axis=0) - points.min(axis=0))) or 1.0
    names = "ijk"
    out: dict[tuple[str, str], list[dict[str, Any]]] = {}
    for label, sides in seams(blocks, 1e-9 * scale, only=shared).items():
        for seam in sides:
            if (label, seam.face) not in shared:
                continue
            along = [
                [names[abs(seam.transform[AXES[a]]) - 1], seam.transform[AXES[a]] < 0]
                for a in FACE_ALONG_AXES[seam.face]  # type: ignore[index]
            ]
            out.setdefault((label, seam.face), []).append(
                {
                    "range": (seam.rng - 1).tolist(),
                    "partner": {"label": seam.donor, "face": seam.donor_face},
                    "axis_map": along,
                }
            )
    for entries in out.values():
        entries.sort(key=lambda e: e["range"][0])
    return out


def _axis_map_to_json(axis_map: AxisMap | None) -> list[list[Any]] | None:
    """Serialise an :class:`AxisMap` to its JSON form, or ``None``."""
    if axis_map is None:
        return None
    return [
        [axis_map.along[0][0], axis_map.along[0][1]],
        [axis_map.along[1][0], axis_map.along[1][1]],
    ]


def _axis_map_from_json(payload: Any, *, source: str) -> AxisMap | None:
    """Parse the JSON form of an :class:`AxisMap`.

    Accepts ``None`` (no map) or a 2-element list of
    ``[axis_label, flipped_bool]`` pairs.  Raises
    :class:`MeshInputError` on malformed input — the reader is the last
    line of defence between a hand-edited or stale sidecar and the
    consumers (cc.par writer etc.) that depend on the field.
    """
    if payload is None:
        return None
    if not isinstance(payload, list) or len(payload) != 2:
        raise MeshInputError(f"{source}: axis_map must be a 2-element list, got {payload!r}")
    pairs: list[tuple[str, bool]] = []
    for entry in payload:
        if not isinstance(entry, list) or len(entry) != 2:
            raise MeshInputError(f"{source}: axis_map entry must be [axis, flipped], got {entry!r}")
        axis, flipped = entry
        if axis not in ("i", "j", "k"):
            raise MeshInputError(f"{source}: axis_map axis label must be 'i'/'j'/'k', got {axis!r}")
        if not isinstance(flipped, bool):
            raise MeshInputError(f"{source}: axis_map flipped flag must be a bool, got {flipped!r}")
        pairs.append((axis, flipped))
    return AxisMap.from_pairs(pairs[0], pairs[1])


def read_adjacency_json(path: str | Path) -> dict[str, Any]:
    """Read an adjacency JSON sidecar and return its raw dict.

    Validates only the version field and the top-level ``blocks`` key.
    The caller is responsible for matching the sidecar against actual
    ``.geo`` files.

    Raises
    ------
    MeshInputError
        On missing file, malformed JSON, missing ``version`` /
        ``blocks`` keys, or unsupported version.
    """
    p = Path(path)
    if not p.exists():
        raise MeshInputError(f"Adjacency sidecar not found: {p}")
    try:
        data = json.loads(p.read_text())
    except json.JSONDecodeError as exc:
        raise MeshInputError(f"Invalid JSON in {p}: {exc}") from exc
    if not isinstance(data, dict):
        raise MeshInputError(f"{p}: expected a JSON object at top level")
    if "version" not in data:
        raise MeshInputError(f"{p}: missing required key 'version'")
    if data["version"] == 1:
        raise MeshInputError(
            f"{p}: adjacency sidecar is version 1; this build emits and "
            f"requires version {ADJACENCY_VERSION}.  Version 1 stored "
            "partner labels but no axis_map, which made cc.par-style "
            "consumers ambiguous on cross-axis seams.  Regenerate the "
            "sidecar by re-running the topology / split that produced it."
        )
    if data["version"] not in _READABLE_VERSIONS:
        raise MeshInputError(
            f"{p}: unsupported adjacency version {data['version']!r} "
            f"(this build reads {', '.join(map(str, _READABLE_VERSIONS))})"
        )
    if "blocks" not in data:
        raise MeshInputError(f"{p}: missing required key 'blocks'")

    # Validate axis_map shapes eagerly so a malformed sidecar fails here
    # rather than deep inside the cc.par writer.
    for blk in data["blocks"]:
        label = blk.get("label", "<unknown>")
        for face_name, face in blk.get("faces", {}).items():
            _axis_map_from_json(
                face.get("axis_map"),
                source=f"{p} ({label}.{face_name})",
            )
            _validate_bc_code(face, source=f"{p} ({label}.{face_name})")
    # Parse associations eagerly too; keep them as Association objects
    table = [
        Association.from_json(a, source=f"{p} (associations[{i}])")
        for i, a in enumerate(data.get("associations", []))
    ]
    data["associations"] = table
    for blk in data["blocks"]:
        for group in blk.get("associated", {}).values():
            for name, idx in group.items():
                if not isinstance(idx, int) or not 0 <= idx < len(table):
                    raise MeshInputError(
                        f"{p} ({blk.get('label')}.{name}): association index {idx!r}"
                        f" out of range (table has {len(table)})"
                    )
    return data


def _validate_bc_code(face: Any, *, source: str) -> None:
    """Validate an optional per-face ``bc_code`` override.

    A ``bc_code`` is only meaningful on a ``FREE`` face and must be a
    positive Xall natural-BC code in ``[1, 19]`` (``nccnat``).  Absent is
    fine (the cc.par writer falls back to the per-block default).
    """
    code = face.get("bc_code")
    if code is None:
        return
    if face.get("bc") != "FREE":
        raise MeshInputError(
            f"{source}: bc_code is only valid on FREE faces, got bc={face.get('bc')!r}"
        )
    if not isinstance(code, int) or isinstance(code, bool):
        raise MeshInputError(f"{source}: bc_code must be an integer, got {code!r}")
    if not 1 <= code <= 19:
        raise MeshInputError(f"{source}: bc_code must be a natural-BC code in [1, 19], got {code}")


def summarise_adjacency(payload: dict[str, Any]) -> str:
    """Return a one-block-per-line human summary string for printing.

    Useful for the ``shore split`` CLI's verbose output.  Walks the
    parsed dict (does not need real HexBlocks).
    """
    lines: list[str] = []
    for blk in payload.get("blocks", []):
        label = blk["label"]
        shape = blk["shape"]
        bc_summary: list[str] = []
        for face_name in _FACE_NAMES:
            face = blk["faces"][face_name]
            bc = face["bc"]
            partner = face["partner"]
            if partner is None:
                bc_summary.append(f"{face_name}={bc}")
            else:
                bc_summary.append(f"{face_name}={bc}->{partner['label']}.{partner['face']}")
        lines.append(f"{label:24s} {tuple(shape)!s:>20s}   " + "  ".join(bc_summary))
    return "\n".join(lines)


def _resolve_label(sidecar_label: str, blocks: list[HexBlock], *, source: str) -> HexBlock:
    """Map a sidecar block label onto exactly one input :class:`HexBlock`.

    The sidecar carries the *topology* labels (``sub0``, ``cap_north``),
    but ``shore split`` derives a block's label from its ``.geo`` filename
    stem, stripping only the common prefix shared by *all* inputs.  In a
    multi-component split (e.g. sphere + background) there is no common
    prefix, so the sphere blocks keep their full stem ``sphere_o_sub0``
    while the sidecar still says ``sub0``.  Match exact first, then the
    unique block whose label ends with ``_<sidecar_label>`` (the unstripped
    prefix).  Zero or multiple matches is a hard error — a silent skip here
    is exactly what let the cubed-sphere's seams vanish in a combined split.
    """
    exact = [b for b in blocks if b.label == sidecar_label]
    if len(exact) == 1:
        return exact[0]
    suffix = [b for b in blocks if b.label.endswith("_" + sidecar_label)]
    if len(suffix) == 1:
        return suffix[0]
    candidates = exact + suffix
    if not candidates:
        raise MeshInputError(
            f"{source}: sidecar block {sidecar_label!r} matches no input block "
            f"(have {sorted(b.label for b in blocks)!r}).  The sidecar labels must "
            f"match the .geo stems shore split derives; check the input set."
        )
    raise MeshInputError(
        f"{source}: sidecar block {sidecar_label!r} is ambiguous — matches "
        f"{sorted(b.label for b in candidates)!r}.  Disambiguate the .geo filenames."
    )


def wire_blocks_from_adjacency(
    blocks: list[HexBlock],
    payload: dict[str, Any],
    *,
    source: str = "<adjacency>",
) -> None:
    """Wire the face BCs of *blocks* in place from a parsed adjacency *payload*.

    ``.geo`` files carry only geometry, so :class:`HexBlock`\\ s built from
    them default every face to WALL (k_lo) / FREE.  This restores the real
    boundary graph (SHARED / DIRICHLET / PERIODIC seams with their partner
    pointers and axis_maps, plus any per-face ``bc_code`` overrides) onto
    those blocks, so a subsequent :func:`shore.split.split_blocks` can
    preserve the pre-existing inter-block seams instead of losing them.

    Every block in *payload* must resolve to exactly one block in *blocks*
    via :func:`_resolve_label` (exact label, else unique ``_<label>``
    suffix).  A sidecar block that matches zero or several input blocks is a
    hard :class:`MeshInputError` — silently skipping a mismatch is what let
    a cubed-sphere's DIRICHLET/SHARED seams silently vanish in a combined
    multi-component split.

    Raises :class:`MeshInputError` on an unknown FaceBC name, a malformed
    axis_map, or an unresolved / ambiguous block or partner label.
    """
    for blk_payload in payload.get("blocks", []):
        label = blk_payload.get("label")
        blk = _resolve_label(label, blocks, source=source)
        for face_name, face in blk_payload.get("faces", {}).items():
            bc_name = face["bc"]
            try:
                bc = FaceBC[bc_name]
            except KeyError as exc:
                raise MeshInputError(
                    f"{source} ({label}.{face_name}): unknown FaceBC {bc_name!r}"
                ) from exc
            partner_ref = face.get("partner")
            if bc in (FaceBC.SHARED, FaceBC.DIRICHLET, FaceBC.PERIODIC):
                if partner_ref is None:
                    raise MeshInputError(
                        f"{source} ({label}.{face_name}): bc={bc_name} requires a partner"
                    )
                partner_blk = _resolve_label(
                    partner_ref["label"], blocks, source=f"{source} ({label}.{face_name} partner)"
                )
                partner_face = partner_blk.face(partner_ref["face"])
                axis_map = _axis_map_from_json(
                    face.get("axis_map"), source=f"{source} ({label}.{face_name})"
                )
                blk._set_face_bc(face_name, bc, partner=partner_face, axis_map=axis_map)
            else:
                blk._set_face_bc(face_name, bc, bc_code=face.get("bc_code"))
        apply_associations(blk, blk_payload, payload.get("associations", []))


def apply_associations(
    block: HexBlock, block_payload: dict[str, Any], table: list[Association]
) -> None:
    """Restore a block's associations from its sidecar entry (version 3)."""
    associated = block_payload.get("associated", {})
    for name, idx in associated.get("faces", {}).items():
        block.face(name).association = table[idx]
    for name, idx in associated.get("edges", {}).items():
        block.edge(name).association = table[idx]
    vertices = block.vertices()
    for name, idx in associated.get("vertices", {}).items():
        vertices[int(name)].association = table[idx]


# ---------------------------------------------------------------------------
# Re-export the FaceBC enum names for callers building or verifying payloads
# ---------------------------------------------------------------------------

VALID_BC_NAMES = frozenset(bc.name for bc in FaceBC)
