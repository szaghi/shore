"""Xall ``cc.par`` writer.

``cc.par`` is the chimera pre-processor input file consumed by the
Xall / Xall tooling.  It pairs with a ``.grd`` mesh file to declare
per-block multigrid level / group / priority and to enumerate every
patch (face boundary) on every block.  Some patches are simple BCs
(walls, inflow); others are *connections* between adjacent blocks.

This module is a one-way bridge: it consumes the SHORE
``.adjacency.json`` sidecar (see :mod:`shore.io.adjacency`) and emits a
fully populated ``cc.par`` ASCII file ready for Xall ingestion.  It
does **not** read ``cc.par`` files — the Xall side has its own
``cc_par_object`` reader.

Orientation codes are derived **from the adjacency JSON's
``axis_map`` field** — there are no per-topology lookup tables in
this writer.  Any topology that produces a fully-wired adjacency
(every SHARED/DIRICHLET face carries ``axis_map`` and ``partner``)
generates correct ``cc.par`` output.

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
- ``SHARED`` / ``DIRICHLET`` → connection orientation code (``>= 100``)
  + ``family = partner patch index``.  The 3-digit orientation code
  ``lᵢ lⱼ lₖ`` is derived from this face's ``axis_map`` plus the
  partner face's name (see :func:`_orientation_code_from_axis_map`).
- ``PERIODIC``   → ``cc.par`` has no first-class periodic BC; emits
  ``40`` with a ``UserWarning``.

Limitations
-----------
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
from shore.io.adjacency import _axis_map_from_json, read_adjacency_json
from shore.mesh.face import FACE_ALONG_AXES, AxisMap

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
# Default BC codes (from Xall modpar.f90)
# ---------------------------------------------------------------------------

_BC_WALL = 1  # parete (viscous wall)
_BC_OFFGEN = 40  # offset for chimera cell-center BC; default for FREE faces

# ---------------------------------------------------------------------------
# Orientation-code derivation (Xall convention from overset/calcorn.f90 lines
# 119-125, 392-411).  A connection patch's BC is a 3-digit code lᵢ lⱼ lₖ where
# each digit describes how *this side's* axis (i, j, k) maps onto the
# partner's axes:
#
#   l=1: partner's I-axis at Imin, increasing
#   l=2: partner's I-axis at Imax, decreasing (flipped)
#   l=3: partner's J-axis at Jmin, increasing
#   l=4: partner's J-axis at Jmax, decreasing (flipped)
#   l=5: partner's K-axis at Kmin, increasing
#   l=6: partner's K-axis at Kmax, decreasing (flipped)
#
# Code = li * 100 + lj * 10 + lk.
#
# This writer derives the code purely from the adjacency JSON's
# ``axis_map`` field plus the partner face name — no per-topology
# tables are needed (they used to live here for the cubed-sphere
# cap-equator seams; they were removed once topologies started
# carrying axis_map directly).
# ---------------------------------------------------------------------------

_AXIS_TO_INDEX: dict[str, int] = {"i": 0, "j": 1, "k": 2}

# Face name → (perpendicular axis label, partner-side at-max flag).
# E.g. "i_lo" → ("i", False), "j_hi" → ("j", True).
_PERP_AXIS_AND_MAX: dict[str, tuple[str, bool]] = {
    "i_lo": ("i", False),
    "i_hi": ("i", True),
    "j_lo": ("j", False),
    "j_hi": ("j", True),
    "k_lo": ("k", False),
    "k_hi": ("k", True),
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
        Name of the ``.grd`` mesh file Xall should pair with this
        ``cc.par`` (header line 1).
    output_basename :
        Output basename Xall will use for its derived files (header
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
            "natively supported by Xall; emitting offgen=40 instead.  "
            "Use Xall chimera boxes for periodic stitching.",
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
    """Compute the cc.par 3-digit orientation code from the face's axis_map.

    The code ``lᵢ lⱼ lₖ`` describes how *this side's* axes (i, j, k) map
    onto the partner block's axes; see the module-level docstring for
    digit semantics.  We need only:

    - This face's name (gives the perpendicular axis on this side).
    - The partner face's name (gives the partner-side perpendicular
      axis and whether it is at min/max — which determines the
      perpendicular digit).
    - The :class:`AxisMap` (gives the two along-axis correspondences,
      including their flip flags).

    Raises
    ------
    ParameterError
        On a connection face missing ``partner`` or ``axis_map`` (the
        adjacency sidecar reader normally catches this earlier).
    """
    partner_payload = face_payload.get("partner")
    if partner_payload is None:
        raise ParameterError(
            f"cc.par: connection face {block_label!r}.{face_name} has no partner "
            "in the adjacency JSON — every SHARED/DIRICHLET face must record "
            "its partner block + face name."
        )
    partner_face_name = partner_payload["face"]

    axis_map = _axis_map_from_json(
        face_payload.get("axis_map"), source=f"cc.par ({block_label}.{face_name})"
    )
    if axis_map is None:
        raise ParameterError(
            f"cc.par: connection face {block_label!r}.{face_name} has no "
            "axis_map in the adjacency JSON — every SHARED/DIRICHLET face "
            "must carry an axis_map (regenerate the sidecar with the "
            "current SHORE build)."
        )
    return _orientation_code_from_axis_map(face_name, partner_face_name, axis_map)


def _orientation_code_from_axis_map(
    face_name: str,
    partner_face_name: str,
    axis_map: AxisMap,
) -> int:
    """Pure function: derive ``lᵢ lⱼ lₖ`` from face names + axis_map.

    Independent of any topology table; testable in isolation.
    """
    perp_axis_self, _self_at_max = _PERP_AXIS_AND_MAX[face_name]
    perp_axis_partner, partner_perp_at_max = _PERP_AXIS_AND_MAX[partner_face_name]

    # Build the (partner_axis, flipped) tuple per this side's axis label.
    # Perpendicular axis on this side maps to the partner's perpendicular
    # axis: partner is *at* its boundary, so direction is fixed by which
    # boundary (Imin/Imax/...).
    self_along: tuple[str, str] = FACE_ALONG_AXES[face_name]  # ("j", "k") for i_lo/i_hi
    by_axis: dict[str, tuple[str, bool]] = {
        perp_axis_self: (perp_axis_partner, partner_perp_at_max),
        self_along[0]: axis_map.along[0],
        self_along[1]: axis_map.along[1],
    }

    # Compose lᵢ, lⱼ, lₖ for this side's i, j, k respectively.
    digits: list[int] = []
    for self_axis in ("i", "j", "k"):
        partner_axis, flipped = by_axis[self_axis]
        digits.append(_digit_for(partner_axis, flipped))
    return digits[0] * 100 + digits[1] * 10 + digits[2]


def _digit_for(partner_axis: str, flipped_or_at_max: bool) -> int:
    """Return the cc.par ``l`` digit for one (partner_axis, direction) pair.

    ``flipped_or_at_max`` overloads two semantics that map to the same
    digit choice:

    - For an *along-seam* mapping it is the AxisMap flip flag.
    - For the *perpendicular* mapping it is whether the partner face is
      at its axis maximum (``*_hi``).

    Both reduce to: ``True`` ⇒ partner-axis at max, decreasing
    (digits 2/4/6); ``False`` ⇒ partner-axis at min, increasing
    (digits 1/3/5).
    """
    base = _AXIS_TO_INDEX[partner_axis] * 2 + 1  # 1, 3, 5 for i, j, k
    return base + (1 if flipped_or_at_max else 0)


def _partner_for(
    block_label: str, face_name: str, face_payload: dict[str, Any]
) -> tuple[str | None, str | None]:
    """Resolve ``(partner_label, partner_face)`` from the adjacency JSON.

    With v2 sidecars every connection face carries its partner — no
    topology lookup tables.  Returns ``(None, None)`` when the face has
    no partner (which is treated as a degraded non-connection patch).
    """
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
    for block_idx, entry in enumerate(entries, start=1):
        lines.append(
            f"{entry.meta.level:9d} {entry.meta.group:9d} {entry.meta.priority:9d} "
            f"! {block_idx} : {entry.meta.comment}"
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
