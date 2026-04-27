"""Split structured blocks along ijk vertex indices, preserving adjacency.

Why splitting matters
---------------------

CFD solvers and parallel decomposers often need *smaller* blocks than
the topology builder naturally produces — for load balancing across
ranks, GPU streaming-multiprocessor occupancy, halo-exchange
granularity, or just to fit a single block within an L2 cache.  SHORE
emits its blocks at the natural topological scale (one cubed-sphere
sub-block, one flat-caps sub-block, etc.); this module is the
post-processing pass that subdivides those blocks into chunks while
recording every internal and external SHARED adjacency that results.

What this is, what it isn't
---------------------------

This is a *vertex-index* splitter.  Given a list of
:class:`~shore.mesh.block.HexBlock` instances and a list of
:class:`BlockSplit` directives — each saying "split block X along
axis A at vertex indices [i1, i2, ...]" — return a fresh list of new
blocks plus all the face BC plumbing.  The original block list is
untouched (deep-copied for safety).

The original starting topology is preserved in spirit: if you had a
6-block cubed sphere, splitting ``sub0`` into 2 chunks gives you 7
blocks (2 sub0 chunks + sub1 + sub2 + sub3 + cap_north + cap_south).
What changes is only the internal block count and the adjacency
graph; world coordinates of every node are bit-exactly preserved.

Out of scope for v1 (deliberate):

- Splitting a single block along multiple axes in one pass.  Compose
  by calling :func:`split_blocks` twice on the result of the first
  call (each call is single-axis-per-block).
- Non-vertex-aligned splits (interpolation between existing nodes).
- Auto-balancing splits across SHARED partners.  When you split block
  A at index 32 and A is SHARED with block B at the affected face,
  B must also be split at the matching index — otherwise the splitter
  raises :class:`ParameterError`.
- Splitting periodic faces (cubed-sphere sub-block ↔ sub-block
  j-seams).  The sub-blocks are aliased numpy views into a contiguous
  ring; splitting fights with that aliasing and is rejected.

Adjacency rules
---------------

For a block split along axis A at one or more indices, the new chunks
expose two kinds of SHARED faces:

- **Internal seams**: the cut planes themselves.  Chunk *k*'s ``A_hi``
  face is SHARED with chunk *(k+1)*'s ``A_lo`` face; the partner is
  set on both sides bit-exactly (no node duplication: the chunks
  aliase the same column from the source block via numpy views).
- **External faces** (every face *not* perpendicular to the split
  axis): each chunk gets the original block's BC for that face.  When
  the original BC was SHARED with another block, the partner pointer
  on the partnered block is updated to point at the right chunk
  (whichever chunk now exposes the face the partner used to share).

  When the original BC was SHARED with a partner that *itself* gets
  split, both blocks must be split at compatible indices on the
  shared axis; otherwise the seam can't be reconstructed.
"""

from __future__ import annotations

from dataclasses import dataclass

from shore.exceptions import ParameterError
from shore.mesh.block import HexBlock
from shore.mesh.face import Face, FaceBC, FaceName

# Axes label → numpy index axis on a (ni, nj, nk, 3) HexBlock node array.
_AXIS_INDEX = {"i": 0, "j": 1, "k": 2}

# Pair of perpendicular face names along an axis: (low-side, high-side).
_AXIS_FACES: dict[str, tuple[FaceName, FaceName]] = {
    "i": ("i_lo", "i_hi"),
    "j": ("j_lo", "j_hi"),
    "k": ("k_lo", "k_hi"),
}


@dataclass(frozen=True)
class BlockSplit:
    """One split directive: which block, which axis, at which indices.

    Parameters
    ----------
    label :
        Label of the block to split (must match a block's ``label`` in the
        input list).
    axis :
        Axis to split along: ``"i"``, ``"j"``, or ``"k"``.
    at :
        Sequence of *vertex* indices at which to cut.  Each index must
        satisfy ``1 <= idx <= n - 1`` along the chosen axis (so each
        resulting chunk has at least 2 vertices in that direction).
        Indices are sorted and deduplicated; ``[32]`` and ``[16, 32]``
        produce 2 and 3 chunks respectively.

    Examples
    --------
    >>> BlockSplit(label="sub0", axis="i", at=[32])
    BlockSplit(label='sub0', axis='i', at=(32,))
    >>> BlockSplit(label="background_center", axis="k", at=[10, 20, 30])
    BlockSplit(label='background_center', axis='k', at=(10, 20, 30))
    """

    label: str
    axis: str
    at: tuple[int, ...]

    def __init__(self, label: str, axis: str, at) -> None:
        if axis not in _AXIS_INDEX:
            raise ParameterError(f"BlockSplit: axis must be 'i', 'j', or 'k'; got {axis!r}")
        cuts = tuple(sorted({int(x) for x in at}))
        if not cuts:
            raise ParameterError(f"BlockSplit({label!r}): at=[] is not a split")
        if cuts[0] < 1:
            raise ParameterError(
                f"BlockSplit({label!r}): cut index {cuts[0]} must be >= 1 "
                "(each chunk needs >= 2 vertices on the split axis)"
            )
        # ``label`` and ``axis`` are validated against actual blocks at
        # apply time (split_blocks); upper-bound on cuts is checked there.
        object.__setattr__(self, "label", label)
        object.__setattr__(self, "axis", axis)
        object.__setattr__(self, "at", cuts)


def split_blocks(
    blocks: list[HexBlock],
    splits: list[BlockSplit],
) -> list[HexBlock]:
    """Apply every split directive and return a fresh list of new blocks.

    Blocks not named in *splits* pass through unchanged (deep-copied).
    Each split block is replaced by ``len(at) + 1`` chunks; internal
    seams between consecutive chunks are wired as ``FaceBC.SHARED``
    with bit-exact node coincidence.

    External SHARED faces are preserved: when the original block had a
    SHARED face whose partner was *another* block in *blocks*, the
    partner's face BC is updated to point at the right new chunk.  When
    both blocks are split at a shared face, the splits must be at
    compatible indices (this implementation requires *identical*
    indices on the shared axis).

    Parameters
    ----------
    blocks :
        The input block list (typically ``Mesh.blocks`` or the output
        of a primitive's ``.build``).  Not mutated.
    splits :
        List of :class:`BlockSplit` directives.  Each directive names
        one block; a single directive may not list the same block
        twice (split along multiple axes in one pass is not supported).

    Returns
    -------
    list[HexBlock]
        Fresh list of new blocks.  Order: split blocks are replaced
        in-place by their chunks (in axis-ascending order); unsplit
        blocks keep their original order.

    Raises
    ------
    ParameterError
        On unknown block labels, out-of-range cut indices, multiple
        splits of the same block in one call, or incompatible splits
        across a SHARED seam.
    """
    # Validate: every split's label exists, no duplicates.
    by_label = {blk.label: blk for blk in blocks}
    seen_labels = set()
    for sp in splits:
        if sp.label not in by_label:
            raise ParameterError(
                f"split_blocks: no block named {sp.label!r} in the input "
                f"(have {sorted(by_label)!r})"
            )
        if sp.label in seen_labels:
            raise ParameterError(
                f"split_blocks: block {sp.label!r} appears in multiple "
                "BlockSplit directives; v1 supports only one axis per block "
                "per call.  Compose by calling split_blocks twice."
            )
        seen_labels.add(sp.label)
        # Bound check on cut indices.
        n = by_label[sp.label].nodes.shape[_AXIS_INDEX[sp.axis]]
        if sp.at[-1] >= n:
            raise ParameterError(
                f"split_blocks({sp.label!r}, axis={sp.axis!r}): cut index "
                f"{sp.at[-1]} out of range (axis has {n} vertices, max cut = {n - 1})"
            )

    splits_by_label = {sp.label: sp for sp in splits}

    # First pass: build the per-block chunk lists (no inter-block wiring yet).
    # Layout: ``chunks[label]`` = list of new HexBlocks (one per chunk) for a
    # split block, OR a single-element list ``[deep_copy(block)]`` for
    # unsplit blocks.
    chunks: dict[str, list[HexBlock]] = {}
    for blk in blocks:
        if blk.label in splits_by_label:
            chunks[blk.label] = _split_one(blk, splits_by_label[blk.label])
        else:
            # Unsplit blocks pass through with a deep-copied face graph
            # so external partner pointers can be safely rewired.
            chunks[blk.label] = [_clone_block(blk)]

    # Second pass: rewire external SHARED faces.  For each *original*
    # block, walk its 6 faces; when a face was SHARED in the original
    # block list, find the new partner chunk(s) and re-wire on both sides.
    _rewire_external_seams(blocks, chunks, splits_by_label)

    # Flatten preserving original order.
    out: list[HexBlock] = []
    for blk in blocks:
        out.extend(chunks[blk.label])
    return out


# ---------------------------------------------------------------------------
# Per-block split
# ---------------------------------------------------------------------------


def _split_one(block: HexBlock, sp: BlockSplit) -> list[HexBlock]:
    """Split *block* along ``sp.axis`` at ``sp.at``; wire internal seams.

    Returns the list of chunks in axis-ascending order.  Each chunk's
    nodes are a contiguous slice of ``block.nodes`` (no copy — they
    aliase the original allocation); copying happens only inside
    :func:`split_blocks` when a block is unsplit.
    """
    axis_idx = _AXIS_INDEX[sp.axis]
    lo_face, hi_face = _AXIS_FACES[sp.axis]
    nodes = block.nodes
    n = nodes.shape[axis_idx]

    # Endpoints of each chunk on the split axis.
    starts = (0, *sp.at)
    stops = (*sp.at, n - 1)
    # Each chunk spans [start, stop] inclusive (vertex indices); chunk
    # k's "stop" vertex is shared with chunk (k+1)'s "start" vertex.

    new_blocks: list[HexBlock] = []
    for k, (s, e) in enumerate(zip(starts, stops, strict=True)):
        sub_nodes = _slice_axis(nodes, axis_idx, s, e + 1)  # +1 for slice end
        new_label = f"{block.label}__{sp.axis}{k}"
        new_blk = HexBlock(
            label=new_label,
            nodes=sub_nodes,
            periodic_i=block.periodic_i if sp.axis != "i" else False,
            periodic_j=block.periodic_j if sp.axis != "j" else False,
        )
        # Inherit the original block's *external* face BCs along the
        # axes perpendicular to the split.  The split axis's two faces
        # are filled in below: lo_face inherits the source's lo_face
        # only on the first chunk, hi_face inherits only on the last.
        for f in block.faces():
            if f.name == lo_face and k != 0:
                continue
            if f.name == hi_face and k != len(starts) - 1:
                continue
            new_blk._set_face_bc(f.name, f.bc, partner=f.partner)
        new_blocks.append(new_blk)

    # Wire internal seams: chunk k's hi_face SHARED with chunk (k+1)'s lo_face.
    for k in range(len(new_blocks) - 1):
        a = new_blocks[k]
        b = new_blocks[k + 1]
        a._set_face_bc(hi_face, FaceBC.SHARED, partner=b.face(lo_face))
        b._set_face_bc(lo_face, FaceBC.SHARED, partner=a.face(hi_face))

    return new_blocks


def _slice_axis(nodes, axis_idx: int, start: int, stop: int):
    """Return ``nodes[..., start:stop, ...]`` along axis *axis_idx*."""
    slicer = [slice(None)] * nodes.ndim
    slicer[axis_idx] = slice(start, stop)
    return nodes[tuple(slicer)]


def _clone_block(block: HexBlock) -> HexBlock:
    """Deep-clone a HexBlock's face graph; node array is *not* copied.

    The face graph (``_faces`` / partner pointers) needs to be
    independent so the splitter can rewire it without mutating the
    original input.  Node data stays aliased — we never need to mutate
    it during splitting.
    """
    clone = HexBlock(
        label=block.label,
        nodes=block.nodes,
        periodic_i=block.periodic_i,
        periodic_j=block.periodic_j,
    )
    for f in block.faces():
        clone._set_face_bc(f.name, f.bc, partner=f.partner)
    return clone


# ---------------------------------------------------------------------------
# External seam rewiring
# ---------------------------------------------------------------------------


def _rewire_external_seams(
    original_blocks: list[HexBlock],
    chunks: dict[str, list[HexBlock]],
    splits: dict[str, BlockSplit],
) -> None:
    """Rewire SHARED partner pointers across split boundaries.

    For each original block, walk its 6 faces; when a face was SHARED
    with another block, redirect the partner pointer to the right chunk
    on each side.  Modifies the ``chunks`` dict in place.
    """
    # Build a fresh by-label map of the *new* blocks for label lookups.
    new_by_label: dict[str, list[HexBlock]] = chunks

    # Track which (orig_label, face_name) pairs we've already wired so we
    # don't double-process from both sides.
    done: set[tuple[str, str]] = set()

    for orig in original_blocks:
        for f in orig.faces():
            if f.bc != FaceBC.SHARED or f.partner is None:
                continue
            # Identify the partner block by walking ``original_blocks``.
            partner_blk = _find_block_for_face(original_blocks, f.partner)
            if partner_blk is None:
                continue
            partner_face_name = f.partner.name

            key1 = (orig.label, f.name)
            key2 = (partner_blk.label, partner_face_name)
            if key1 in done or key2 in done:
                continue
            done.add(key1)
            done.add(key2)

            _wire_external_pair(
                orig=orig,
                orig_face=f.name,
                partner=partner_blk,
                partner_face=partner_face_name,
                chunks=new_by_label,
                splits=splits,
            )


def _find_block_for_face(blocks: list[HexBlock], face: Face) -> HexBlock | None:
    """Return the block whose ``faces()`` list contains *face* (by identity)."""
    for blk in blocks:
        for f in blk.faces():
            if f is face:
                return blk
    return None


def _wire_external_pair(
    orig: HexBlock,
    orig_face: FaceName,
    partner: HexBlock,
    partner_face: FaceName,
    chunks: dict[str, list[HexBlock]],
    splits: dict[str, BlockSplit],
) -> None:
    """Wire one SHARED seam between (orig, orig_face) and (partner, partner_face).

    Four cases:
    - Neither block split: unsplit pass-through (both chunks have a single
      element); just point the partner refs at each other's faces.
    - One block split, perpendicular face on the split: the split axis
      doesn't run along the SHARED face → the face appears on every
      chunk (at the same position).  This is a *broken* seam: the
      partner's single face cannot point at multiple chunks.  Reject.
    - One block split, the split-axis face is the SHARED face: only the
      first or last chunk exposes that face.  Wire those two.
    - Both blocks split: the split axes must produce matching index
      sequences on the shared face; pair chunks 1:1.
    """
    orig_split = splits.get(orig.label)
    partner_split = splits.get(partner.label)
    orig_chunks = chunks[orig.label]
    partner_chunks = chunks[partner.label]

    # The split axis "runs along" a face when the face is one of the two
    # endpoints of that axis: e.g. axis 'i' runs along faces 'i_lo' and
    # 'i_hi'.  Otherwise the split axis is *perpendicular* to the face,
    # meaning every chunk produced by the split still exposes the face
    # at full extent (just sliced into pieces along the perpendicular
    # axis).
    orig_axis_face = orig_split is not None and orig_face.startswith(orig_split.axis)
    partner_axis_face = partner_split is not None and partner_face.startswith(partner_split.axis)

    # Case 1: neither block was split.
    if orig_split is None and partner_split is None:
        a = orig_chunks[0]
        b = partner_chunks[0]
        a._set_face_bc(orig_face, FaceBC.SHARED, partner=b.face(partner_face))
        b._set_face_bc(partner_face, FaceBC.SHARED, partner=a.face(orig_face))
        return

    # Case 2: split runs along the SHARED face on one side but the other
    # side is unsplit, OR split axes differ between the two sides.  The
    # face would multiply on one side without a matching partition on
    # the other — not reconstructible.
    one_unsplit = orig_split is None or partner_split is None
    perp_split_unmatched = not orig_axis_face or not partner_axis_face
    axes_differ = (
        orig_split is not None
        and partner_split is not None
        and orig_split.axis != partner_split.axis
    )
    if (one_unsplit and perp_split_unmatched) or (axes_differ and perp_split_unmatched):
        raise ParameterError(
            f"split_blocks: SHARED seam between {orig.label!r}({orig_face}) "
            f"and {partner.label!r}({partner_face}) crosses a split axis on "
            "one side but not the other — the partner block would also need "
            "to be split with matching indices.  Add the corresponding "
            "BlockSplit on the partner block, or skip the split."
        )

    # Case 3: both split along the same axis perpendicular to the SHARED
    # face, with matching cuts → pair chunks 1:1.
    if (
        orig_split is not None
        and partner_split is not None
        and not orig_axis_face
        and not partner_axis_face
    ):
        if orig_split.axis != partner_split.axis or orig_split.at != partner_split.at:
            raise ParameterError(
                f"split_blocks: SHARED seam between {orig.label!r}({orig_face}) "
                f"and {partner.label!r}({partner_face}): both blocks are split "
                f"on a perpendicular axis but with different (axis, indices) "
                f"({orig_split.axis}{list(orig_split.at)} vs "
                f"{partner_split.axis}{list(partner_split.at)}).  Use the same "
                "split on both sides."
            )
        for a, b in zip(orig_chunks, partner_chunks, strict=True):
            a._set_face_bc(orig_face, FaceBC.SHARED, partner=b.face(partner_face))
            b._set_face_bc(partner_face, FaceBC.SHARED, partner=a.face(orig_face))
        return

    # Cases 4 + 5: split axis aligns with the SHARED face on every split
    # side.  Only the first/last chunk exposes ``axis_lo`` / ``axis_hi``.
    a_chunks = _exposed_chunks(orig_chunks, orig_face, orig_split)
    b_chunks = _exposed_chunks(partner_chunks, partner_face, partner_split)
    if len(a_chunks) != len(b_chunks):
        raise ParameterError(
            f"split_blocks: SHARED seam between {orig.label!r}({orig_face}) "
            f"and {partner.label!r}({partner_face}) split into {len(a_chunks)} "
            f"vs {len(b_chunks)} pieces — split indices must match on both "
            "sides of a SHARED seam."
        )
    for a, b in zip(a_chunks, b_chunks, strict=True):
        a._set_face_bc(orig_face, FaceBC.SHARED, partner=b.face(partner_face))
        b._set_face_bc(partner_face, FaceBC.SHARED, partner=a.face(orig_face))


def _exposed_chunks(
    chunk_list: list[HexBlock],
    face_name: FaceName,
    split: BlockSplit | None,
) -> list[HexBlock]:
    """Return the chunks that still expose *face_name* after splitting.

    For an unsplit block (``split is None``), the face is on the only
    chunk.  For a block split along the axis of *face_name*, only the
    first chunk exposes ``<axis>_lo`` and only the last chunk exposes
    ``<axis>_hi``.  When *split* operates on a different axis, every
    chunk exposes the face — but :func:`_wire_external_pair` rejects
    that case before we get here.
    """
    if split is None:
        return [chunk_list[0]]
    axis = split.axis
    if face_name == f"{axis}_lo":
        return [chunk_list[0]]
    if face_name == f"{axis}_hi":
        return [chunk_list[-1]]
    # Face perpendicular to the split: every chunk exposes it.
    return list(chunk_list)
