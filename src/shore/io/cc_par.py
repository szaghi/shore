"""Xnavis ``cc.par`` writer.

``cc.par`` is the chimera pre-processor input file consumed by the
Xnavis / Xall tooling.  It pairs with a ``.grd`` mesh file to declare
per-block multigrid level / group / priority and to enumerate every
patch (face boundary) on every block.  Some patches are simple BCs
(walls, inflow); others are *connections* between adjacent blocks.

This module is a one-way bridge: it consumes the SHORE
``.adjacency.json`` sidecar (see :mod:`shore.io.adjacency`) and emits a
fully populated ``cc.par`` ASCII file ready for Xnavis ingestion.  It
does **not** read ``cc.par`` files — the Xnavis side has its own
``cc_par_object`` reader.

Format reference
----------------
``cc.par`` has 5 sections:

1. **Header** (10 lines): grd filename, output basename, 3 booleans,
   blank, 2 multigrid levels, blank, boundary-layer thickness,
   numerical-beach width, blank.
2. **Blocks**: ``N`` followed by ``N`` lines of
   ``level group priority ! comment``.
3. **Patches**: ``N`` followed by ``N`` lines of
   ``block face BC family Imin Imax Jmin Jmax Kmin Kmax ! patch_index``.
   ``face`` is 1..6 for Imin/Imax/Jmin/Jmax/Kmin/Kmax.  Extents are in
   **cell counts**, i.e. ``n_nodes - 1`` per axis.  ``BC > 99`` with
   ``family > 0`` denotes a connection; ``family`` is the partner
   patch's 1-based index, and the BC value's three digits encode the
   partner-axis orientation (see :func:`_orientation_code`).
4. **Edges**: always ``0`` for SHORE-generated meshes.
5. **Boxes**: hole-cutting boxes (user-supplied for chimera donor
   zones).  SHORE has no boxes, so we always emit ``0``.

Reference implementation: ``overset_exploded_cc_par_object.F90`` in
``Xnavis-projects/overset-exploded`` (the Fortran writer this Python
output is meant to round-trip with).

BC code table
-------------
Mapping from SHORE's ``FaceBC`` enum to ``cc.par`` BC codes:

- ``WALL``       → ``1`` (``parete``, viscous wall) + ``family = 0``.
- ``FREE``       → caller-chosen via :class:`BlockMeta.free_face_bc`,
  defaulting to ``40`` (``offgen``, chimera face cell-center).  Body-
  fitted and buffer meshes use ``40``; far-field background patches
  use a natural BC (e.g. ``-9``, extrapolation).
- ``SHARED``     → connection orientation code (``>= 100``) +
  ``family = partner patch index``.  For same-axis seams (i↔i, j↔j,
  k↔k, no flip — the only kind SHORE produces between primitive splits
  and within the cubed-sphere equatorial ring), the canonical code is
  ``135`` on both sides.  See :func:`_orientation_code`.
- ``DIRICHLET``  → connection orientation code, same as ``SHARED``.
  In the cubed-sphere topology, DIRICHLET is used for cap-equator
  seams which are *cross-axis* (the sub-block's i-axis maps to one of
  the cap's lateral faces); the orientation code differs per sub_q
  and per pole.  See :data:`_CAP_EQUATOR_ORIENT`.
- ``PERIODIC``   → ``cc.par`` has no first-class periodic BC; emits
  ``40`` with a ``UserWarning``.

Cap-equator partner reconciliation
----------------------------------
SHORE's adjacency JSON records the sub-block side of cap-equator
seams (``sub_q.i_lo / i_hi → cap_*.k_lo``) but the cap side carries
``partner: null`` because the topology builder leaves it unset.  This
writer reconstructs the symmetric cap-side patch from the
:func:`shore.topology.sphere_cap._overwrite_cap_boundaries` mapping
(see :data:`_CAP_EQUATOR_PARTNER`); the cap-side ``i_lo / i_hi /
j_lo / j_hi`` faces are emitted as connections back to the sub-block,
overriding any DIRICHLET-without-partner placeholder in the JSON.

Limitations
-----------
- Only the cubed-sphere cap-equator cross-axis layout is hard-coded.
  Future topologies introducing additional cross-axis seams must
  extend :data:`_CAP_EQUATOR_PARTNER` and :data:`_CAP_EQUATOR_ORIENT`.
- No ``boxes`` section — user-supplied hole-cutting boxes for chimera
  donor zones are out of scope; the writer emits ``0 ! boxes number``.
- No ``edges`` section — likewise unsupported by ``cc.par`` parsers
  for SHORE-style adjacency.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from shore.exceptions import MeshInputError, ParameterError
from shore.io.adjacency import read_adjacency_json

# ---------------------------------------------------------------------------
# Face-name ↔ cc.par face-index translation
# ---------------------------------------------------------------------------

_FACE_TO_INDEX: dict[str, int] = {
    "i_lo": 1,  # Imin
    "i_hi": 2,  # Imax
    "j_lo": 3,  # Jmin
    "j_hi": 4,  # Jmax
    "k_lo": 5,  # Kmin
    "k_hi": 6,  # Kmax
}
# Iteration order matters: cc.par patch indices are assigned in this order
# per block.  Matches the canonical _FACE_NAMES order in adjacency.py.
_FACE_NAMES: tuple[str, ...] = ("i_lo", "i_hi", "j_lo", "j_hi", "k_lo", "k_hi")

# ---------------------------------------------------------------------------
# Default BC codes (from Xnavis modpar.f90)
# ---------------------------------------------------------------------------

_BC_WALL = 1  # parete (viscous wall)
_BC_OFFGEN = 40  # offset for chimera cell-center BC; default for FREE faces

# Same-axis SHARED / DIRICHLET seam code (i→Imin-inc, j→Jmin-inc, k→Kmin-inc).
# Used when adjacent blocks share the same coordinate frame and the seam is
# bit-exact axis-aligned.  See module docstring.
_SAME_AXIS_CODE = 135

# ---------------------------------------------------------------------------
# Cubed-sphere cap-equator cross-axis partner / orientation tables
# ---------------------------------------------------------------------------
#
# Derived from
# :func:`shore.topology.sphere_cap._overwrite_cap_boundaries`, which assigns:
#
#   North cap (sub.i_hi pole → cap_north[..., 0, :]):
#     cap_north[ 0, :, :] = sub0[ni-1, :, :]              # cap.i_lo edge
#     cap_north[ :,-1, :] = sub1[ni-1, :, :]              # cap.j_hi edge
#     cap_north[-1, :, :] = sub2[ni-1, ::-1, :]           # cap.i_hi edge (reversed)
#     cap_north[ :, 0, :] = sub3[ni-1, ::-1, :]           # cap.j_lo edge (reversed)
#
#   South cap (sub.i_lo pole → cap_south[..., 0, :]):
#     cap_south[ :, 0, :] = sub0[0, :, :]                 # cap.j_lo edge
#     cap_south[-1, :, :] = sub1[0, :, :]                 # cap.i_hi edge
#     cap_south[ :,-1, :] = sub2[0, ::-1, :]              # cap.j_hi edge (reversed)
#     cap_south[ 0, :, :] = sub3[0, ::-1, :]              # cap.i_lo edge (reversed)
#
# Orientation digits per Xnavis convention (overset/calcorn.f90 lines 119-125,
# 392-411):
#
#   l=1: partner's I-axis at Imin, increasing
#   l=2: partner's I-axis at Imax, decreasing (flip)
#   l=3: partner's J-axis at Jmin, increasing
#   l=4: partner's J-axis at Jmax, decreasing (flip)
#   l=5: partner's K-axis at Kmin, increasing
#   l=6: partner's K-axis at Kmax, decreasing (flip)
#
# Code = li*100 + lj*10 + lk where (li, lj, lk) describes how this patch's
# (i, j, k) maps onto the partner's axes.
#
# All cap-equator seams in the cubed sphere have the wall-normal axis (k)
# preserved between sub and cap (both grow outward), so on every line below
# either lk=5 (preserved, increasing) or — when the corresponding partner
# face perpendicular is k — lk encodes that.  The partner-face perpendicular
# is the cap's k (face 5/6, but the cap-side seam face is one of i_lo/i_hi/
# j_lo/j_hi, perpendicular = i or j); k is therefore an *along-seam* axis
# on the cap side too, and remains lk=5.  For the sub side, the patch face
# is i_lo or i_hi → perpendicular axis = i.  See per-pair derivation below.

# (sub_q, pole) → (cap_face_name, sub_face_name, sub_orient, cap_orient)
_CAP_EQUATOR_PARTNER: dict[tuple[int, str], tuple[str, str, int, int]] = {
    # ----- North pole: sub.i_hi (face 2) → cap_north -----
    # sub0: cap_north[0, :, :] ← sub0[ni-1, :, :]
    #   - sub's i (perpendicular, Imax) → cap's i (Imin, increasing inward) → li=1
    #   - sub's j (along seam) → cap's j (preserved, increasing)            → lj=3
    #   - sub's k (along seam) → cap's k (preserved, increasing)            → lk=5
    #   sub-side code: 135
    #   cap-side (cap_north.i_lo, perpendicular = cap's i, Imin):
    #   - cap's i (perp Imin) → sub's i (Imax, decreasing inward)           → li=2
    #   - cap's j (along)     → sub's j (preserved)                          → lj=3
    #   - cap's k (along)     → sub's k (preserved)                          → lk=5
    #   cap-side code: 235
    (0, "north"): ("i_lo", "i_hi", 135, 235),
    # sub1: cap_north[:, -1, :] ← sub1[ni-1, :, :]
    #   sub's pole row's j-index runs along cap's i-index (cap's j is fixed at -1).
    #   - sub's i (perp Imax) → cap's j (Jmax, decreasing inward)            → li=4
    #   - sub's j (along)     → cap's i (Imin, increasing)                   → lj=1
    #   - sub's k (along)     → cap's k                                      → lk=5
    #   sub-side code: 415
    #   cap-side (cap_north.j_hi, perpendicular = cap's j, Jmax):
    #   - cap's i (along)     → sub's j (Jmin, increasing)                   → li=3
    #   - cap's j (perp Jmax) → sub's i (Imax, decreasing inward)            → lj=2
    #   - cap's k (along)     → sub's k                                      → lk=5
    #   cap-side code: 325
    (1, "north"): ("j_hi", "i_hi", 415, 325),
    # sub2: cap_north[-1, :, :] ← sub2[ni-1, ::-1, :]  (j-axis reversed)
    #   - sub's i (perp Imax) → cap's i (Imax, decreasing inward)            → li=2
    #   - sub's j (along, reversed) → cap's j (Jmax, decreasing)             → lj=4
    #   - sub's k (along)     → cap's k                                      → lk=5
    #   sub-side code: 245
    #   cap-side (cap_north.i_hi):
    #   - cap's i (perp Imax) → sub's i (Imax, decreasing)                   → li=2
    #   - cap's j (along)     → sub's j (Jmax, decreasing — reversed pole)   → lj=4
    #   - cap's k (along)     → sub's k                                      → lk=5
    #   cap-side code: 245
    (2, "north"): ("i_hi", "i_hi", 245, 245),
    # sub3: cap_north[:, 0, :] ← sub3[ni-1, ::-1, :]
    #   - sub's i (perp Imax) → cap's j (Jmin, increasing inward)            → li=3
    #   - sub's j (along, reversed) → cap's i (Imax, decreasing)             → lj=2
    #   - sub's k (along)     → cap's k                                      → lk=5
    #   sub-side code: 325
    #   cap-side (cap_north.j_lo):
    #   - cap's i (along)     → sub's j (Jmax, decreasing — reversed)       → li=4
    #   - cap's j (perp Jmin) → sub's i (Imax, decreasing inward)            → lj=2
    #   - cap's k (along)     → sub's k                                      → lk=5
    #   cap-side code: 425
    (3, "north"): ("j_lo", "i_hi", 325, 425),
    # ----- South pole: sub.i_lo (face 1) → cap_south -----
    # sub0: cap_south[:, 0, :] ← sub0[0, :, :]
    #   - sub's i (perp Imin) → cap's j (Jmin, increasing inward)            → li=3
    #   - sub's j (along)     → cap's i (Imin, increasing)                   → lj=1
    #   - sub's k (along)     → cap's k                                      → lk=5
    #   sub-side code: 315
    #   cap-side (cap_south.j_lo):
    #   - cap's i (along)     → sub's j (Jmin, increasing)                   → li=3
    #   - cap's j (perp Jmin) → sub's i (Imin, increasing inward)            → lj=1
    #   - cap's k (along)     → sub's k                                      → lk=5
    #   cap-side code: 315
    (0, "south"): ("j_lo", "i_lo", 315, 315),
    # sub1: cap_south[-1, :, :] ← sub1[0, :, :]
    #   - sub's i (perp Imin) → cap's i (Imax, decreasing inward)            → li=2
    #   - sub's j (along)     → cap's j (Jmin, increasing)                   → lj=3
    #   - sub's k (along)     → cap's k                                      → lk=5
    #   sub-side code: 235
    #   cap-side (cap_south.i_hi):
    #   - cap's i (perp Imax) → sub's i (Imin, increasing inward)            → li=1
    #   - cap's j (along)     → sub's j (Jmin, increasing)                   → lj=3
    #   - cap's k (along)     → sub's k                                      → lk=5
    #   cap-side code: 135
    (1, "south"): ("i_hi", "i_lo", 235, 135),
    # sub2: cap_south[:, -1, :] ← sub2[0, ::-1, :]
    #   - sub's i (perp Imin) → cap's j (Jmax, decreasing inward)            → li=4
    #   - sub's j (along, reversed) → cap's i (Imax, decreasing)             → lj=2
    #   - sub's k (along)     → cap's k                                      → lk=5
    #   sub-side code: 425
    #   cap-side (cap_south.j_hi):
    #   - cap's i (along)     → sub's j (Jmax, decreasing — reversed)        → li=4
    #   - cap's j (perp Jmax) → sub's i (Imin, increasing inward)            → lj=1
    #   - cap's k (along)     → sub's k                                      → lk=5
    #   cap-side code: 415
    (2, "south"): ("j_hi", "i_lo", 425, 415),
    # sub3: cap_south[0, :, :] ← sub3[0, ::-1, :]
    #   - sub's i (perp Imin) → cap's i (Imin, increasing inward)            → li=1
    #   - sub's j (along, reversed) → cap's j (Jmax, decreasing)             → lj=4
    #   - sub's k (along)     → cap's k                                      → lk=5
    #   sub-side code: 145
    #   cap-side (cap_south.i_lo):
    #   - cap's i (perp Imin) → sub's i (Imin, increasing inward)            → li=1
    #   - cap's j (along)     → sub's j (Jmax, decreasing — reversed)        → lj=4
    #   - cap's k (along)     → sub's k                                      → lk=5
    #   cap-side code: 145
    (3, "south"): ("i_lo", "i_lo", 145, 145),
}


# ---------------------------------------------------------------------------
# Public dataclasses
# ---------------------------------------------------------------------------


@dataclass
class BlockMeta:
    """Per-block metadata for a ``cc.par`` blocks-section line.

    Parameters
    ----------
    level :
        Overset level (priority among overlapping bodies).  Default 1.
    group :
        Group index.  Default 0.
    priority :
        Within-group priority.  Default 1.
    comment :
        Trailing ``! comment`` for the line.  Defaults to the block label.
    free_face_bc :
        BC code emitted for ``FREE`` faces of this block.  Default ``40``
        (``offgen``, chimera-cell-center) — appropriate for body-fitted /
        buffer meshes.  Use a negative natural-BC code (e.g. ``-9`` for
        extrapolation, ``-3`` for inflow) on far-field background blocks.
    """

    level: int = 1
    group: int = 0
    priority: int = 1
    comment: str = ""
    free_face_bc: int = _BC_OFFGEN


@dataclass
class _Patch:
    """Internal: one resolved patch line about to be written."""

    block_index: int  # 1-based
    face_index: int  # 1..6
    bc: int
    family: int  # partner patch index (1-based) or 0
    extents: tuple[int, int, int, int, int, int]  # Imin Imax Jmin Jmax Kmin Kmax in cells
    patch_index: int = 0  # filled in after global ordering


@dataclass
class _BlockEntry:
    """Internal: one resolved block line + its 6 patches."""

    label: str
    shape_cells: tuple[int, int, int]  # (ni-1, nj-1, nk-1)
    meta: BlockMeta
    patches: list[_Patch] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Public writer
# ---------------------------------------------------------------------------


def write_cc_par(
    adjacency: dict[str, Any] | str | Path,
    path: str | Path,
    *,
    grd_filename: str,
    output_basename: str = "cc.",
    save_ghost_cells: bool = True,
    increase_overlap: bool = False,
    extend_internal_wall: bool = False,
    multigrid_levels: tuple[int, int] = (4, 1),
    boundary_layer_thickness: float = -0.1,
    numerical_beach: float = -1.0,
    block_metadata: dict[str, BlockMeta] | None = None,
) -> Path:
    """Write a ``cc.par`` file from a SHORE adjacency payload.

    Parameters
    ----------
    adjacency :
        Either a parsed adjacency dict (as produced by
        :func:`shore.io.adjacency.write_adjacency_json` and read back via
        :func:`shore.io.adjacency.read_adjacency_json`) or a path to a
        ``.adjacency.json`` file on disk.
    path :
        Output ``cc.par`` path.
    grd_filename :
        Name of the ``.grd`` mesh file Xnavis should pair with this
        ``cc.par`` (header line 1).
    output_basename :
        Output basename Xnavis will use for its derived files (header
        line 2).  Default ``"cc."``.
    save_ghost_cells, increase_overlap, extend_internal_wall :
        Header booleans (lines 3-5).
    multigrid_levels :
        Two-element tuple ``(N_levels, finest_level_active)`` (line 7).
    boundary_layer_thickness, numerical_beach :
        Header reals (lines 9-10).
    block_metadata :
        Optional ``{label: BlockMeta}`` dict.  Any block not listed gets
        a default ``BlockMeta()``.

    Returns
    -------
    Path
        The actual path written.
    """
    payload = _resolve_payload(adjacency)
    blocks_meta = block_metadata or {}

    # Pass 1: build _BlockEntry list with shape and meta.
    label_to_index: dict[str, int] = {}
    entries: list[_BlockEntry] = []
    for idx, blk in enumerate(payload["blocks"], start=1):
        label = blk["label"]
        ni, nj, nk = blk["shape"]
        meta = blocks_meta.get(label, BlockMeta())
        if not meta.comment:
            meta = BlockMeta(
                level=meta.level,
                group=meta.group,
                priority=meta.priority,
                comment=label,
                free_face_bc=meta.free_face_bc,
            )
        entries.append(
            _BlockEntry(
                label=label,
                shape_cells=(ni - 1, nj - 1, nk - 1),
                meta=meta,
            )
        )
        label_to_index[label] = idx

    # Pass 2: build patches per block in canonical face order.
    for entry in entries:
        block_idx = label_to_index[entry.label]
        block_payload = next(b for b in payload["blocks"] if b["label"] == entry.label)
        for face_name in _FACE_NAMES:
            face = block_payload["faces"][face_name]
            patch = _resolve_patch(
                block_idx=block_idx,
                block_label=entry.label,
                face_name=face_name,
                face_payload=face,
                shape_cells=entry.shape_cells,
                free_face_bc=entry.meta.free_face_bc,
            )
            entry.patches.append(patch)

    # Pass 3: assign global 1-based patch indices.
    patch_counter = 0
    by_label_face: dict[tuple[str, str], int] = {}
    for entry in entries:
        for patch, face_name in zip(entry.patches, _FACE_NAMES, strict=True):
            patch_counter += 1
            patch.patch_index = patch_counter
            by_label_face[(entry.label, face_name)] = patch_counter

    # Pass 4: resolve `family` references for connection patches.
    # The first pass left family = -1 for connection patches whose partner
    # patch index wasn't yet known; we now have the (label, face) → index map.
    for entry in entries:
        block_payload = next(b for b in payload["blocks"] if b["label"] == entry.label)
        for patch, face_name in zip(entry.patches, _FACE_NAMES, strict=True):
            if patch.family != -1:
                continue
            # Re-fetch the partner from the resolved cap-equator table or the
            # adjacency JSON.
            partner_label, partner_face = _partner_for(
                entry.label, face_name, block_payload["faces"][face_name]
            )
            if partner_label is None:
                # No partner found: degrade to non-connection.
                patch.bc = _BC_OFFGEN
                patch.family = 0
                continue
            partner_idx = by_label_face.get((partner_label, partner_face))
            if partner_idx is None:
                raise ParameterError(
                    f"cc.par: block {entry.label!r} face {face_name!r} points to "
                    f"unknown partner ({partner_label!r}, {partner_face!r})"
                )
            patch.family = partner_idx

    # Pass 5: validate symmetry for connections.
    _validate_symmetry(entries)

    # Pass 6: serialise.
    out = _format_cc_par(
        entries=entries,
        grd_filename=grd_filename,
        output_basename=output_basename,
        save_ghost_cells=save_ghost_cells,
        increase_overlap=increase_overlap,
        extend_internal_wall=extend_internal_wall,
        multigrid_levels=multigrid_levels,
        boundary_layer_thickness=boundary_layer_thickness,
        numerical_beach=numerical_beach,
    )
    out_path = Path(path)
    out_path.write_text(out)
    return out_path


# ---------------------------------------------------------------------------
# Internals
# ---------------------------------------------------------------------------


def _resolve_payload(adjacency: dict[str, Any] | str | Path) -> dict[str, Any]:
    """Accept either a parsed dict or a path to a JSON file."""
    if isinstance(adjacency, (str, Path)):
        return read_adjacency_json(adjacency)
    if not isinstance(adjacency, dict):
        raise MeshInputError(
            "cc.par: adjacency must be a parsed payload dict or a path to a "
            f"JSON file, got {type(adjacency).__name__}"
        )
    return adjacency


def _resolve_patch(
    *,
    block_idx: int,
    block_label: str,
    face_name: str,
    face_payload: dict[str, Any],
    shape_cells: tuple[int, int, int],
    free_face_bc: int,
) -> _Patch:
    """Translate one face entry into a ``_Patch``.

    Connection patches get ``family = -1`` here; pass 4 fills it in from
    the global label-face index.
    """
    bc_name = face_payload["bc"]
    face_idx = _FACE_TO_INDEX[face_name]
    extents = _face_extents(face_name, shape_cells)

    if bc_name == "WALL":
        return _Patch(block_idx, face_idx, _BC_WALL, 0, extents)
    if bc_name == "FREE":
        return _Patch(block_idx, face_idx, free_face_bc, 0, extents)
    if bc_name == "PERIODIC":
        warnings.warn(
            f"cc.par: PERIODIC BC on {block_label!r}.{face_name} not "
            "natively supported by Xnavis; emitting offgen=40 instead.  "
            "Use Xnavis chimera boxes for periodic stitching.",
            UserWarning,
            stacklevel=3,
        )
        return _Patch(block_idx, face_idx, _BC_OFFGEN, 0, extents)
    if bc_name in ("SHARED", "DIRICHLET"):
        # Connection: BC = orientation code, family = -1 for now.
        bc_code = _orientation_code(block_label, face_name, face_payload)
        return _Patch(block_idx, face_idx, bc_code, -1, extents)
    raise ParameterError(f"cc.par: unknown FaceBC name {bc_name!r} on {block_label!r}.{face_name}")


def _face_extents(
    face_name: str, shape_cells: tuple[int, int, int]
) -> tuple[int, int, int, int, int, int]:
    """Return ``(Imin, Imax, Jmin, Jmax, Kmin, Kmax)`` cell extents for *face*.

    Cell indices are 0-based in cc.par's convention: a face perpendicular to
    axis ``a`` pins both extents on that axis to the same boundary value, and
    spans the other two axes from 0 to ``n_cells``.
    """
    ni_c, nj_c, nk_c = shape_cells
    if face_name == "i_lo":
        return (0, 0, 0, nj_c, 0, nk_c)
    if face_name == "i_hi":
        return (ni_c, ni_c, 0, nj_c, 0, nk_c)
    if face_name == "j_lo":
        return (0, ni_c, 0, 0, 0, nk_c)
    if face_name == "j_hi":
        return (0, ni_c, nj_c, nj_c, 0, nk_c)
    if face_name == "k_lo":
        return (0, ni_c, 0, nj_c, 0, 0)
    if face_name == "k_hi":
        return (0, ni_c, 0, nj_c, nk_c, nk_c)
    raise ParameterError(f"cc.par: unknown face name {face_name!r}")


def _orientation_code(block_label: str, face_name: str, face_payload: dict[str, Any]) -> int:
    """Choose the cc.par orientation code for a connection patch.

    Same-axis seams (i↔i, j↔j, k↔k) use ``_SAME_AXIS_CODE = 135``.
    Cross-axis cubed-sphere cap-equator seams look up the per-pair code
    from :data:`_CAP_EQUATOR_PARTNER`.
    """
    cap_pair = _cap_equator_lookup(block_label, face_name)
    if cap_pair is not None:
        _partner_label, _partner_face, sub_code, cap_code = cap_pair
        return cap_code if block_label.startswith("cap_") else sub_code

    # Default: same-axis SHARED / DIRICHLET seam.
    return _SAME_AXIS_CODE


def _cap_equator_lookup(block_label: str, face_name: str) -> tuple[str, str, int, int] | None:
    """Return ``(partner_label, partner_face, sub_code, cap_code)`` for a
    cubed-sphere cap-equator seam, or ``None`` when *block.face* isn't one.
    """
    # Sub-block side: label like 'sub0', 'sub1'... and face is i_lo or i_hi.
    if block_label.startswith("sub") and face_name in ("i_lo", "i_hi"):
        try:
            q = int(block_label[3:].split("_", 1)[0])
        except ValueError:
            return None
        if not (0 <= q <= 3):
            return None
        pole = "south" if face_name == "i_lo" else "north"
        entry = _CAP_EQUATOR_PARTNER.get((q, pole))
        if entry is None:
            return None
        cap_face, _sub_face, sub_code, cap_code = entry
        cap_label = f"cap_{pole}"
        return cap_label, cap_face, sub_code, cap_code

    # Cap side: label like 'cap_north' / 'cap_south', face is i_lo/i_hi/j_lo/j_hi.
    if block_label in ("cap_north", "cap_south") and face_name in (
        "i_lo",
        "i_hi",
        "j_lo",
        "j_hi",
    ):
        pole = "north" if block_label == "cap_north" else "south"
        for (q, p), (cap_face_t, sub_face_t, sub_code, cap_code) in _CAP_EQUATOR_PARTNER.items():
            if p == pole and cap_face_t == face_name:
                return f"sub{q}", sub_face_t, sub_code, cap_code
        return None

    return None


def _partner_for(
    block_label: str, face_name: str, face_payload: dict[str, Any]
) -> tuple[str | None, str | None]:
    """Resolve ``(partner_label, partner_face)`` for a connection face.

    For cap-equator seams the JSON's ``partner`` is unreliable (cap side
    has ``null``), so consult the cap-equator table first.
    """
    cap_pair = _cap_equator_lookup(block_label, face_name)
    if cap_pair is not None:
        partner_label, partner_face, _sub_code, _cap_code = cap_pair
        return partner_label, partner_face

    p = face_payload.get("partner")
    if p is None:
        return None, None
    return p["label"], p["face"]


def _validate_symmetry(entries: list[_BlockEntry]) -> None:
    """Verify that every connection patch is paired symmetrically.

    Mirrors the consistency check in ``overset/readpar.f90`` (lines 219-228):
    if patch ``p`` has ``BC > 99`` and ``family == q``, then patch ``q``'s
    family must equal ``p``.
    """
    by_index: dict[int, _Patch] = {}
    for entry in entries:
        for patch in entry.patches:
            by_index[patch.patch_index] = patch
    for patch in by_index.values():
        if patch.bc <= 99:
            continue
        if patch.family <= 0:
            raise ParameterError(
                f"cc.par: connection patch {patch.patch_index} (block "
                f"{patch.block_index}, face {patch.face_index}) has BC={patch.bc} "
                "but no family — every connection must point at a partner."
            )
        partner = by_index.get(patch.family)
        if partner is None:
            raise ParameterError(
                f"cc.par: patch {patch.patch_index} family={patch.family} is out of range"
            )
        if partner.family != patch.patch_index:
            raise ParameterError(
                f"cc.par: asymmetric connection — patch {patch.patch_index} "
                f"points at {patch.family}, but patch {patch.family} points at "
                f"{partner.family} instead.  Both sides of a SHARED/DIRICHLET "
                "seam must reference each other."
            )


def _format_cc_par(
    *,
    entries: list[_BlockEntry],
    grd_filename: str,
    output_basename: str,
    save_ghost_cells: bool,
    increase_overlap: bool,
    extend_internal_wall: bool,
    multigrid_levels: tuple[int, int],
    boundary_layer_thickness: float,
    numerical_beach: float,
) -> str:
    """Serialise *entries* into the cc.par textual format."""
    lines: list[str] = []

    # ----- Header -----
    lines.append(f"'{grd_filename}'")
    lines.append(f"'{output_basename}'")
    lines.append(".true. " if save_ghost_cells else ".false.")
    lines.append(".true. " if increase_overlap else ".false.")
    lines.append(".true. " if extend_internal_wall else ".false.")
    lines.append("")
    lines.append(f"{multigrid_levels[0]:3d} {multigrid_levels[1]:3d}")
    lines.append("")
    lines.append(f"{boundary_layer_thickness: .4E}")
    lines.append(f"{numerical_beach: .4E}")
    lines.append("")

    # ----- Blocks -----
    lines.append(f"{len(entries):d} ! blocks number")
    lines.append("")
    for entry in entries:
        lines.append(
            f"{entry.meta.level:9d} {entry.meta.group:9d} {entry.meta.priority:9d} "
            f"! {entry.meta.comment}"
        )
    lines.append("")

    # ----- Patches -----
    total_patches = sum(len(e.patches) for e in entries)
    lines.append(f"{total_patches:d} ! patches number")
    lines.append("")
    for entry in entries:
        for patch in entry.patches:
            ext = patch.extents
            lines.append(
                f"{patch.block_index:9d} {patch.face_index:9d} {patch.bc:9d} "
                f"{patch.family:9d} {ext[0]:9d} {ext[1]:9d} {ext[2]:9d} "
                f"{ext[3]:9d} {ext[4]:9d} {ext[5]:9d} ! {patch.patch_index:9d}"
            )
    lines.append("")

    # ----- Edges (always 0 for SHORE) -----
    lines.append("0 ! edges number")
    lines.append("")

    # ----- Boxes (always 0 for SHORE) -----
    lines.append("0 ! boxes number")
    lines.append("")

    return "\n".join(lines)
