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

from shore.exceptions import MeshInputError
from shore.mesh.block import HexBlock
from shore.mesh.face import AxisMap, FaceBC, FaceName

ADJACENCY_VERSION = 2

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
            face_payload[f.name] = {
                "bc": f.bc.name,
                "partner": partner_obj,
                "axis_map": axis_map_obj,
            }
        payload["blocks"].append(
            {
                "label": blk.label,
                "shape": [ni, nj, nk],
                "faces": face_payload,
            }
        )

    out.write_text(json.dumps(payload, indent=2) + "\n")
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
    if data["version"] != ADJACENCY_VERSION:
        raise MeshInputError(
            f"{p}: unsupported adjacency version {data['version']!r} "
            f"(this build supports {ADJACENCY_VERSION})"
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
    return data


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


# ---------------------------------------------------------------------------
# Re-export the FaceBC enum names for callers building or verifying payloads
# ---------------------------------------------------------------------------

VALID_BC_NAMES = frozenset(bc.name for bc in FaceBC)
