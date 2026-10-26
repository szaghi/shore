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

Capabilities:

- **Multi-axis in one call.**  A block may be cut on i, j and k at once
  (one :class:`BlockSplit` per axis, same label); it becomes an
  ``ni x nj x nk`` chunk grid.
- **Automatic paired co-splits.**  A cut that subdivides a SHARED /
  DIRICHLET seam is propagated through the seam's AxisMap onto the
  partner block (to a fixpoint), so the caller never authors paired
  directives and the two sides stay conformal.  Cut *indices* on the
  two sides need not be equal — they are mapped (with flip) through the
  AxisMap.

Out of scope (deliberate):

- Non-vertex-aligned splits (interpolation between existing nodes).
- Splitting periodic faces (cubed-sphere sub-block ↔ sub-block
  j-seams).  The sub-blocks are aliased numpy views into a contiguous
  ring; splitting fights with that aliasing and is rejected.

Adjacency rules
---------------

A block split into a chunk grid exposes two kinds of seams:

- **Internal seams**: the cut planes themselves.  Adjacent chunks are
  wired SHARED with identity AxisMaps and bit-exact node coincidence
  (the chunks alias the same column from the source via numpy views).
- **External seams**: each original face is inherited only by the chunks
  that still expose it; an original inter-block SHARED/DIRICHLET seam is
  re-paired chunk-sheet to chunk-sheet onto the matching partner chunks
  (each side keeps its own BC — DIRICHLET stays DIRICHLET).
"""

from __future__ import annotations

from dataclasses import dataclass, field

from shore.exceptions import ParameterError
from shore.mesh.block import HexBlock
from shore.mesh.face import FACE_ALONG_AXES, AxisMap, Face, FaceBC, FaceName

# Axes label → numpy index axis on a (ni, nj, nk, 3) HexBlock node array.
_AXIS_INDEX = {"i": 0, "j": 1, "k": 2}

# Pair of perpendicular face names along an axis: (low-side, high-side).
_AXIS_FACES: dict[str, tuple[FaceName, FaceName]] = {
    "i": ("i_lo", "i_hi"),
    "j": ("j_lo", "j_hi"),
    "k": ("k_lo", "k_hi"),
}


def _inherit_face(target: HexBlock, source_face: Face) -> None:
    """Copy a face's BC + partner + axis_map from *source_face* onto *target*.

    Used by both the per-chunk inheritance pass in :func:`_split_one`
    and by :func:`_clone_block` for unsplit pass-through.  Connection
    faces (SHARED / DIRICHLET / PERIODIC) carry an axis_map that
    survives the partner-pointer rewiring untouched: slicing a face
    along a perpendicular axis does not change the along-axis
    correspondence at the seam.
    """
    target._set_face_bc(
        source_face.name,
        source_face.bc,
        partner=source_face.partner,
        axis_map=source_face.axis_map,
        bc_code=source_face.bc_code,
    )


def _inherit_associations(
    chunk: HexBlock, block: HexBlock, span: dict[str, tuple[int, int]]
) -> None:
    """Give *chunk*'s vertices, edges and faces the association of the
    lowest-dimensional *block* entity they lie on (ADR 0003).

    A chunk entity is fixed at ``lo``/``hi`` along some axes; each such side
    lies on the source boundary only where the chunk touches it. Three
    boundary constraints name a source vertex, two an edge, one a face; none
    (an internal cut) leaves the entity unassociated.
    """
    n = {ax: block.nodes.shape[_AXIS_INDEX[ax]] for ax in "ijk"}
    vertex_by_ijk = {v.ijk: v for v in block.vertices()}

    def source_association(fixed: dict[str, str]):
        on = {}
        for ax, side in fixed.items():
            s, e = span[ax]
            if (side == "lo" and s == 0) or (side == "hi" and e == n[ax] - 1):
                on[ax] = side
        if len(on) == 3:
            ijk = tuple(0 if on[ax] == "lo" else -1 for ax in "ijk")
            return vertex_by_ijk[ijk].association
        if len(on) == 2:
            a1, a2 = sorted(on)
            return block.edge(f"{a1}_{on[a1]}_{a2}_{on[a2]}").association
        if len(on) == 1:
            ((ax, side),) = on.items()
            return block.face(f"{ax}_{side}").association
        return None

    for f in chunk.faces():
        f.association = source_association({f.name[0]: f.name[2:]})
    for e in chunk.edges():
        p = e.name.split("_")
        e.association = source_association({p[0]: p[1], p[2]: p[3]})
    for v in chunk.vertices():
        fixed = {ax: ("lo" if idx == 0 else "hi") for ax, idx in zip("ijk", v.ijk, strict=True)}
        v.association = source_association(fixed)


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
    *,
    propagate_seams: bool = True,
) -> list[HexBlock]:
    """Apply split directives and return a fresh, fully-wired block list.

    A block may be cut on **several axes in one call**, and cuts that
    subdivide a SHARED / DIRICHLET seam are **auto-propagated to the
    partner block** so the caller never has to author paired co-splits.
    The result is a conformal multi-block decomposition with every
    internal cut-plane and every original inter-block seam wired (partner
    pointers + AxisMaps) onto the right chunks.

    Parameters
    ----------
    blocks :
        Input block list (e.g. ``Mesh.blocks`` or a primitive's
        ``.build``).  Not mutated (chunks alias node memory read-only).
    splits :
        :class:`BlockSplit` directives.  Multiple directives may name the
        same block (one per axis); their cuts are unioned per axis.
    propagate_seams :
        When True (default), a cut on an in-plane axis of a SHARED /
        DIRICHLET face is mirrored through the seam's AxisMap onto the
        partner block (to a fixpoint), so both sides of every seam are
        cut consistently.  Set False only to reproduce the legacy
        "caller must pre-pair every co-split" behaviour.

    Returns
    -------
    list[HexBlock]
        Fresh blocks.  An original block split into an ``ni x nj x nk``
        chunk grid is replaced (in row-major i-major-then-j-then-k order)
        by its chunks, labelled ``<label>__i<a>__j<b>__k<c>`` (only the
        cut axes appear in the suffix); unsplit blocks keep their label
        and original position.

    Raises
    ------
    ParameterError
        Unknown label, out-of-range cut, or a seam whose two sides cannot
        be made conformal (e.g. a partner with an incompatible extent).
    """
    by_label = {blk.label: blk for blk in blocks}

    # 1. Collect requested cuts: cuts[label][axis] = sorted set of vertex idx.
    cuts: dict[str, dict[str, set[int]]] = {b.label: {} for b in blocks}
    for sp in splits:
        if sp.label not in by_label:
            raise ParameterError(
                f"split_blocks: no block named {sp.label!r} in the input "
                f"(have {sorted(by_label)!r})"
            )
        n = by_label[sp.label].nodes.shape[_AXIS_INDEX[sp.axis]]
        if sp.at[-1] >= n:
            raise ParameterError(
                f"split_blocks({sp.label!r}, axis={sp.axis!r}): cut index "
                f"{sp.at[-1]} out of range (axis has {n} vertices, max cut = {n - 1})"
            )
        cuts[sp.label].setdefault(sp.axis, set()).update(sp.at)

    # 2. Propagate cuts across SHARED/DIRICHLET seams to a fixpoint, so both
    #    sides of every subdivided seam carry matching cuts.
    if propagate_seams:
        _propagate_seam_cuts(blocks, by_label, cuts)

    # 3. Build the per-block chunk grids (node-aliased slices, internal
    #    cut-plane seams wired).  chunks[label] = ChunkGrid.
    grids: dict[str, _ChunkGrid] = {}
    for blk in blocks:
        grids[blk.label] = _split_one_multi(blk, cuts[blk.label])

    # 4. Rewire the original inter-block seams onto the matching chunk pairs.
    _rewire_external_seams(blocks, grids, cuts)

    # 5. Flatten preserving original block order.
    out: list[HexBlock] = []
    for blk in blocks:
        out.extend(grids[blk.label].flat)
    return out


# ---------------------------------------------------------------------------
# Seam cut propagation
# ---------------------------------------------------------------------------


def _propagate_seam_cuts(
    blocks: list[HexBlock],
    by_label: dict[str, HexBlock],
    cuts: dict[str, dict[str, set[int]]],
) -> None:
    """Mirror in-plane cuts across SHARED/DIRICHLET seams until stable.

    A face perpendicular to axis A has two in-plane axes.  When a block is
    cut on an in-plane axis of a connection face, the partner block must
    receive the corresponding cut (mapped through the seam's AxisMap, with
    index reversal when the seam is flipped) so the two sides stay
    conformal after splitting.  Iterates to a fixpoint because a
    propagated cut can itself induce further cuts on a third block.

    Cuts on the axis *perpendicular* to a face do not affect that seam's
    in-plane partition, so they are never propagated through it (a k-cut
    of a block whose i_lo face is SHARED only changes how that i_lo face
    is sliced along k — handled by the partner's own k-cut if any, via the
    same propagation on the k-perpendicular faces).
    """
    # Identity lookup: which (label, face) each Face object belongs to.
    owner: dict[int, tuple[str, FaceName]] = {}
    for blk in blocks:
        for f in blk.faces():
            owner[id(f)] = (blk.label, f.name)

    connection = (FaceBC.SHARED, FaceBC.DIRICHLET, FaceBC.PERIODIC)
    changed = True
    while changed:
        changed = False
        for blk in blocks:
            for f in blk.faces():
                if f.bc not in connection or f.partner is None or f.axis_map is None:
                    continue
                p_owner = owner.get(id(f.partner))
                if p_owner is None:
                    continue
                p_label, _p_face = p_owner
                # This face's two in-plane axes, in canonical order.
                my_axes = FACE_ALONG_AXES[f.name]
                partner_axes = f.axis_map.partner_axes()
                flips = f.axis_map.flips()
                partner_blk = by_label[p_label]
                for slot, my_ax in enumerate(my_axes):
                    src = cuts[blk.label].get(my_ax)
                    if not src:
                        continue
                    p_ax = partner_axes[slot]
                    flip = flips[slot]
                    # Partner extent (vertex count) along the mapped axis.
                    p_n = partner_blk.nodes.shape[_AXIS_INDEX[p_ax]]
                    dst = cuts[p_label].setdefault(p_ax, set())
                    for c in src:
                        mapped = (p_n - 1 - c) if flip else c
                        if 0 < mapped < p_n - 1 and mapped not in dst:
                            dst.add(mapped)
                            changed = True


def propagate_cuts(
    blocks: list[HexBlock],
    requested: dict[str, dict[str, set[int]]],
) -> dict[str, dict[str, set[int]]]:
    """Return *requested* cuts after seam-propagation to a fixpoint.

    Pure helper for the balance planner: given per-block per-axis cut sets
    (vertex indices), return a deep copy extended with every cut induced
    across SHARED/DIRICHLET seams — exactly the set :func:`split_blocks`
    would execute.  Lets the planner reason about the *actual* chunking
    (and enforce the multigrid cell floor) on the propagated result, so its
    plan matches the split output block-for-block.
    """
    by_label = {b.label: b for b in blocks}
    cuts: dict[str, dict[str, set[int]]] = {
        b.label: {ax: set(requested.get(b.label, {}).get(ax, ())) for ax in ("i", "j", "k")}
        for b in blocks
    }
    _propagate_seam_cuts(blocks, by_label, cuts)
    return cuts


def chunk_axis_sizes(n_vertices: int, at: set[int]) -> list[int]:
    """Cell counts of each chunk produced by cutting an axis at *at*."""
    return [e - s for s, e in _chunk_bounds(n_vertices, at)]


# ---------------------------------------------------------------------------
# Per-block multi-axis split
# ---------------------------------------------------------------------------


@dataclass
class _ChunkGrid:
    """The chunk decomposition of one original block.

    ``chunks[(ia, ja, ka)]`` is the HexBlock for the chunk at chunk-grid
    coordinate ``(ia, ja, ka)`` (0-based per axis).  ``n`` is the chunk
    count per axis ``(ni, nj, nk)``.  An unsplit axis has count 1.
    """

    label: str
    n: tuple[int, int, int]
    chunks: dict[tuple[int, int, int], HexBlock] = field(default_factory=dict)

    @property
    def flat(self) -> list[HexBlock]:
        """Chunks in row-major (i-major, then j, then k) order."""
        ni, nj, nk = self.n
        return [self.chunks[(a, b, c)] for a in range(ni) for b in range(nj) for c in range(nk)]

    def exposed(self, face_name: FaceName) -> dict[tuple[int, int], HexBlock]:
        """Chunks exposing *face_name*, keyed by their two in-plane chunk coords.

        For ``i_lo`` the exposing chunks are those at ia=0 (the whole
        j x k chunk sheet); the key is each chunk's (ja, ka) in the
        :data:`FACE_ALONG_AXES` order for *face_name*.
        """
        ni, nj, nk = self.n
        perp = face_name[0]  # 'i' / 'j' / 'k'
        lo = face_name.endswith("_lo")
        extent = {"i": ni, "j": nj, "k": nk}
        out: dict[tuple[int, int], HexBlock] = {}
        for (a, b, c), blk in self.chunks.items():
            coord = {"i": a, "j": b, "k": c}
            if coord[perp] != (0 if lo else extent[perp] - 1):
                continue
            ax0, ax1 = FACE_ALONG_AXES[face_name]
            out[(coord[ax0], coord[ax1])] = blk
        return out


def _chunk_bounds(n_vertices: int, at: set[int]) -> list[tuple[int, int]]:
    """Vertex [start, stop_inclusive] spans for an axis cut at *at*."""
    cut = sorted(at)
    starts = [0, *cut]
    stops = [*cut, n_vertices - 1]
    return list(zip(starts, stops, strict=True))


def _split_one_multi(block: HexBlock, axis_cuts: dict[str, set[int]]) -> _ChunkGrid:
    """Split *block* on every requested axis at once into a chunk grid.

    *axis_cuts* maps ``'i'/'j'/'k'`` to a set of vertex cut indices (may
    be empty/missing for an uncut axis).  Each chunk's nodes are a
    contiguous numpy slice of ``block.nodes`` (aliased, not copied).
    External faces inherit the source BC only on the chunks that expose
    them; internal cut-plane seams are wired SHARED with identity maps.
    An entirely unsplit block yields a 1x1x1 grid (its sole chunk is a
    face-graph clone, node array aliased — safe to rewire).
    """
    bounds = {
        ax: _chunk_bounds(block.nodes.shape[_AXIS_INDEX[ax]], axis_cuts.get(ax) or set())
        for ax in ("i", "j", "k")
    }
    ni, nj, nk = (len(bounds["i"]), len(bounds["j"]), len(bounds["k"]))
    grid = _ChunkGrid(label=block.label, n=(ni, nj, nk))
    last = {"i": ni - 1, "j": nj - 1, "k": nk - 1}

    for a in range(ni):
        for b in range(nj):
            for c in range(nk):
                (si, ei), (sj, ej), (sk, ek) = bounds["i"][a], bounds["j"][b], bounds["k"][c]
                sub = block.nodes[si : ei + 1, sj : ej + 1, sk : ek + 1, :]
                suffix = "".join(
                    f"__{ax}{idx}"
                    for ax, idx, cnt in (("i", a, ni), ("j", b, nj), ("k", c, nk))
                    if cnt > 1
                )
                chunk = HexBlock(
                    label=block.label + suffix,
                    nodes=sub,
                    periodic_i=block.periodic_i if ni == 1 else False,
                    periodic_j=block.periodic_j if nj == 1 else False,
                )
                # Inherit each source external face only on the chunks that
                # still expose it (e.g. i_lo only on the ia=0 sheet).
                coord = {"i": a, "j": b, "k": c}
                for f in block.faces():
                    perp = f.name[0]
                    want = 0 if f.name.endswith("_lo") else last[perp]
                    if coord[perp] == want:
                        _inherit_face(chunk, f)
                _inherit_associations(chunk, block, {"i": (si, ei), "j": (sj, ej), "k": (sk, ek)})
                grid.chunks[(a, b, c)] = chunk

    # Wire internal cut-plane seams between adjacent chunks on each cut axis.
    for ax in ("i", "j", "k"):
        lo_face, hi_face = _AXIS_FACES[ax]
        cnt = {"i": ni, "j": nj, "k": nk}[ax]
        for key, chunk in grid.chunks.items():
            a, b, c = key
            idx = {"i": a, "j": b, "k": c}[ax]
            if idx == cnt - 1:
                continue
            nxt_key = {"i": (a + 1, b, c), "j": (a, b + 1, c), "k": (a, b, c + 1)}[ax]
            nxt = grid.chunks[nxt_key]
            chunk._set_face_bc(
                hi_face,
                FaceBC.SHARED,
                partner=nxt.face(lo_face),
                axis_map=AxisMap.identity(hi_face),
            )
            nxt._set_face_bc(
                lo_face,
                FaceBC.SHARED,
                partner=chunk.face(hi_face),
                axis_map=AxisMap.identity(lo_face),
            )

    return grid


# ---------------------------------------------------------------------------
# External seam rewiring
# ---------------------------------------------------------------------------


def _rewire_external_seams(
    original_blocks: list[HexBlock],
    grids: dict[str, _ChunkGrid],
    cuts: dict[str, dict[str, set[int]]],
) -> None:
    """Re-pair every original inter-block seam onto the matching chunk sheets.

    Cut propagation (:func:`_propagate_seam_cuts`) has already guaranteed
    both sides of every SHARED/DIRICHLET seam carry conformal in-plane
    cuts, so each exposing chunk on one side has exactly one partner chunk
    on the other.  This walks each original seam once, maps in-plane chunk
    coordinates through the seam's AxisMap, and wires the partner pointers
    (preserving each side's original BC — DIRICHLET stays DIRICHLET).
    """
    owner: dict[int, str] = {}
    for blk in original_blocks:
        for f in blk.faces():
            owner[id(f)] = blk.label

    connection = (FaceBC.SHARED, FaceBC.DIRICHLET, FaceBC.PERIODIC)
    done: set[tuple[str, str]] = set()
    for orig in original_blocks:
        for f in orig.faces():
            if f.bc not in connection or f.partner is None:
                continue
            p_label = owner.get(id(f.partner))
            if p_label is None:
                continue
            p_face = f.partner.name
            if (orig.label, f.name) in done or (p_label, p_face) in done:
                continue
            done.add((orig.label, f.name))
            done.add((p_label, p_face))
            if f.axis_map is None:
                raise ParameterError(
                    f"split_blocks: seam {orig.label!r}({f.name}) ↔ "
                    f"{p_label!r}({p_face}) has no axis_map — input topology "
                    f"must be fully wired (run a topology builder / sidecar first)."
                )
            partner_bc = next(
                pf.bc for pf in original_blocks_face(original_blocks, p_label, p_face)
            )
            _wire_seam_sheets(
                a_grid=grids[orig.label],
                a_face=f.name,
                a_map=f.axis_map,
                a_bc=f.bc,
                b_grid=grids[p_label],
                b_face=p_face,
                b_bc=partner_bc,
            )


def original_blocks_face(blocks: list[HexBlock], label: str, face: FaceName):
    """Yield the Face object for *label*.*face* among *blocks* (0 or 1 items)."""
    for blk in blocks:
        if blk.label == label:
            yield blk.face(face)


def _wire_seam_sheets(
    *,
    a_grid: _ChunkGrid,
    a_face: FaceName,
    a_map: AxisMap,
    a_bc: FaceBC,
    b_grid: _ChunkGrid,
    b_face: FaceName,
    b_bc: FaceBC,
) -> None:
    """Pair the exposing chunk sheets of one seam side to the other and wire.

    *a_map* maps *a_face*'s two in-plane axes onto *b_face*'s along-axes,
    each with a flip flag.  For each exposing chunk on side A at in-plane
    coord ``(ka0, ka1)`` (in ``FACE_ALONG_AXES[a_face]`` order), the
    matching B chunk's coord on each B along-axis is the same chunk index,
    reversed when that axis is flipped.
    """
    a_sheet = a_grid.exposed(a_face)
    b_sheet = b_grid.exposed(b_face)
    if len(a_sheet) != len(b_sheet):
        raise ParameterError(
            f"split_blocks: seam {a_grid.label!r}({a_face}) ↔ {b_grid.label!r}({b_face}) "
            f"has {len(a_sheet)} vs {len(b_sheet)} exposing chunks after split — "
            f"seam cut propagation failed to make the two sides conformal."
        )

    b_along = FACE_ALONG_AXES[b_face]
    partner_axes = a_map.partner_axes()  # b-axis for each a along-axis slot
    flips = a_map.flips()
    # chunk counts on b per axis, for flip arithmetic
    b_counts = {"i": b_grid.n[0], "j": b_grid.n[1], "k": b_grid.n[2]}
    # position of each b along-axis within b's (along0, along1) key tuple
    b_slot_of = {b_along[0]: 0, b_along[1]: 1}

    a_inv = a_map  # this side's map
    b_inv = a_map.inverse(a_face, b_face)  # partner side's map (for its axis_map)

    for a_coord, a_chunk in a_sheet.items():
        # Map a_coord (ordered by a_along) onto a b-key (ordered by b_along).
        b_key = [0, 0]
        for slot, a_idx in enumerate(a_coord):
            b_ax = partner_axes[slot]
            cnt = b_counts[b_ax]
            mapped = (cnt - 1 - a_idx) if flips[slot] else a_idx
            b_key[b_slot_of[b_ax]] = mapped
        b_chunk = b_sheet[(b_key[0], b_key[1])]
        a_chunk._set_face_bc(a_face, a_bc, partner=b_chunk.face(b_face), axis_map=a_inv)
        b_chunk._set_face_bc(b_face, b_bc, partner=a_chunk.face(a_face), axis_map=b_inv)
