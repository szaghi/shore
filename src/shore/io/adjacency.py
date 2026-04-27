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
(SHARED, DIRICHLET, PERIODIC) record the partner's ``label`` and
``face`` name; faces without a partner (WALL, FREE) carry ``null``.

Example::

    {
      "version": 1,
      "blocks": [
        {
          "label": "sub0__k0",
          "shape": [20, 9, 5],
          "faces": {
            "i_lo": {"bc": "DIRICHLET", "partner": {"label": "cap_south", "face": "k_lo"}},
            "i_hi": {"bc": "DIRICHLET", "partner": {"label": "cap_north", "face": "k_lo"}},
            "j_lo": {"bc": "SHARED",    "partner": {"label": "sub3__k0", "face": "j_hi"}},
            "j_hi": {"bc": "SHARED",    "partner": {"label": "sub1__k0", "face": "j_lo"}},
            "k_lo": {"bc": "WALL",      "partner": null},
            "k_hi": {"bc": "SHARED",    "partner": {"label": "sub0__k1", "face": "k_lo"}}
          }
        },
        ...
      ]
    }
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from shore.exceptions import MeshInputError
from shore.mesh.block import HexBlock
from shore.mesh.face import FaceBC, FaceName

ADJACENCY_VERSION = 1

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
            partner_obj = None
            if f.partner is not None:
                p_lookup = by_face_id.get(id(f.partner))
                if p_lookup is not None:
                    partner_obj = {"label": p_lookup[0], "face": p_lookup[1]}
            face_payload[f.name] = {"bc": f.bc.name, "partner": partner_obj}
        payload["blocks"].append(
            {
                "label": blk.label,
                "shape": [ni, nj, nk],
                "faces": face_payload,
            }
        )

    out.write_text(json.dumps(payload, indent=2) + "\n")
    return out


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
    if data["version"] != ADJACENCY_VERSION:
        raise MeshInputError(
            f"{p}: unsupported adjacency version {data['version']!r} "
            f"(this build supports {ADJACENCY_VERSION})"
        )
    if "blocks" not in data:
        raise MeshInputError(f"{p}: missing required key 'blocks'")
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
