"""Blocking: a block skeleton edited by commands, premeshed into blocks (ADR 0042).

A blocking is the topology of a multi-block grid without its nodes: skeleton
vertices (points), hexahedral blocks over them (8 vertex ids each, in
``HexBlock.vertices()`` order: index ``4 i + 2 j + k``), a node count on
every edge, and associations of vertices, edges and faces with geometry.
Blocks touching each other share vertex ids, so a face or an edge is the
same entity in every block that has it.

Editing keeps the blocking conforming:
- ``split`` cuts a block across one axis and carries the cut through every
  block sharing a cut edge (the whole sheet), at the same node;
- ``delete`` removes blocks (the hole they leave is a boundary).

``premesh`` makes ordinary ``HexBlock``\\ s:
- vertices are projected onto their associations;
- each edge is a straight line of equally spaced nodes, projected onto its
  association and re-equidistributed along the projected polyline;
- each face is a TFI of its four edges, projected onto its association;
- each block is a TFI of its six faces.
Edges and faces are made once and copied into every block that has them,
so seams match bit for bit. Faces shared by two blocks become SHARED
seams; boundary faces with an association are WALL, the others FREE.

An entity without its own association takes the one of the entities above
it, lowest dimension first, when they agree: an edge the association of its
associated faces, a vertex that of its associated edges, else of its faces.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

import numpy as np

from shore.exceptions import ParameterError
from shore.mesh.face import FACE_ALONG_AXES, AxisMap, FaceBC

if TYPE_CHECKING:
    from shore.mesh.association import Association, AssociationResolver
    from shore.mesh.block import HexBlock

AXES = ("i", "j", "k")
_BIT = {"i": 4, "j": 2, "k": 1}

# HexBlock edge names -> (axis, local vertex at the axis' low end, at its high end)
EDGES: dict[str, tuple[str, int, int]] = {}
for _axis in AXES:
    _a, _b = (x for x in AXES if x != _axis)
    for _hb in (0, 1):
        for _ha in (0, 1):
            _name = f"{_a}_{'hi' if _ha else 'lo'}_{_b}_{'hi' if _hb else 'lo'}"
            _lo = _ha * _BIT[_a] + _hb * _BIT[_b]
            EDGES[_name] = (_axis, _lo, _lo + _BIT[_axis])
del _axis, _a, _b, _ha, _hb, _name, _lo


def face_corners(face: str) -> np.ndarray:
    """Local vertex indices (2, 2) of *face*, over its along-axes (FACE_ALONG_AXES)."""
    axis, side = face[0], face[2:]
    p, q = FACE_ALONG_AXES[face]  # type: ignore[index]
    base = _BIT[axis] if side == "hi" else 0
    return np.array([[base + a * _BIT[p] + b * _BIT[q] for b in (0, 1)] for a in (0, 1)])


FACES = tuple(FACE_ALONG_AXES)

EntityKey = tuple  # ("v", id) | ("e", a, b) a < b | ("f", a, b, c, d) sorted


@dataclass
class Blocking:
    """A block skeleton (see the module docstring)."""

    points: np.ndarray  # (V, 3) vertex positions
    blocks: list[tuple[str, tuple[int, ...]]]  # (label, 8 vertex ids)
    nodes: dict[tuple[int, int], int]  # edge (a < b) -> node count
    associations: dict[EntityKey, Association] = field(default_factory=dict)
    next_label: int = 1

    def copy(self) -> Blocking:
        return dataclasses.replace(
            self,
            points=self.points.copy(),
            blocks=list(self.blocks),
            nodes=dict(self.nodes),
            associations=dict(self.associations),
        )

    # -- lookup ---------------------------------------------------------------

    def labels(self) -> list[str]:
        return [label for label, _ in self.blocks]

    def block(self, label: str) -> tuple[int, ...]:
        for name, verts in self.blocks:
            if name == label:
                return verts
        raise ParameterError(f"blocking: no block {label!r}; blocks: {self.labels()}.")

    def edge_count(self, a: int, b: int) -> int:
        return self.nodes[(min(a, b), max(a, b))]

    def entity(self, label: str, name: str) -> tuple[EntityKey, str]:
        """The key of face, edge or vertex *name* of block *label*, and which it is."""
        verts = self.block(label)
        if name in FACE_ALONG_AXES:
            return ("f", *sorted(verts[c] for c in face_corners(name).ravel())), "face"
        if name in EDGES:
            _, lo, hi = EDGES[name]
            a, b = sorted((verts[lo], verts[hi]))
            return ("e", a, b), "edge"
        if name.startswith("v") and name[1:].isdigit() and int(name[1:]) < 8:
            return ("v", verts[int(name[1:])]), "vertex"
        raise ParameterError(
            f"blocking: {name!r} is not a face, edge or vertex (v0..v7) of block {label!r}."
        )

    def shape(self, label: str) -> tuple[int, int, int]:
        verts = self.block(label)
        return tuple(  # type: ignore[return-value]
            self.edge_count(verts[0], verts[_BIT[a]]) for a in AXES
        )

    def edges(self) -> dict[tuple[int, int], list[tuple[int, str]]]:
        """Every edge (a < b) -> the (block index, edge name) that have it."""
        out: dict[tuple[int, int], list[tuple[int, str]]] = {}
        for n, (_, verts) in enumerate(self.blocks):
            for name, (_, lo, hi) in EDGES.items():
                key = (min(verts[lo], verts[hi]), max(verts[lo], verts[hi]))
                out.setdefault(key, []).append((n, name))
        return out

    def polylines(self) -> list[tuple[str, str, np.ndarray]]:
        """Every edge once, as ``(block, edge name, nodes (n, 3))``: its nodes
        equally spaced on the straight line between its skeleton vertices, in
        the direction of the edge in the first block that has it."""
        out = []
        for key, owners in self.edges().items():
            b, name = owners[0]
            label, verts = self.blocks[b]
            _, lo, hi = EDGES[name]
            t = np.linspace(0.0, 1.0, self.nodes[key])[:, None]
            out.append(
                (label, name, (1.0 - t) * self.points[verts[lo]] + t * self.points[verts[hi]])
            )
        return out

    def faces(self) -> dict[tuple[int, ...], list[tuple[int, str]]]:
        """Every face (sorted vertex ids) -> the (block index, face name) that have it."""
        out: dict[tuple[int, ...], list[tuple[int, str]]] = {}
        for n, (_, verts) in enumerate(self.blocks):
            for name in FACES:
                key = tuple(sorted(verts[c] for c in face_corners(name).ravel()))
                out.setdefault(key, []).append((n, name))
        return out


# -- construction and editing ---------------------------------------------------


def box(vmin: Any, vmax: Any, ni: int, nj: int, nk: int, label: str = "b0") -> Blocking:
    """One block between the corners *vmin* and *vmax*, with *ni* x *nj* x *nk* nodes."""
    lo, hi = np.asarray(vmin, dtype=np.float64), np.asarray(vmax, dtype=np.float64)
    if lo.shape != (3,) or hi.shape != (3,) or not np.all(hi > lo):
        raise ParameterError(f"blocking: need vmin < vmax on every axis, got {lo} and {hi}.")
    counts = {"i": ni, "j": nj, "k": nk}
    if any(int(n) != n or n < 2 for n in counts.values()):
        raise ParameterError(f"blocking: node counts must be integers >= 2, got {counts}.")
    points = np.array(
        [[(hi if v & _BIT[a] else lo)[d] for d, a in enumerate(AXES)] for v in range(8)]
    )
    nodes = {(lo_v, hi_v): int(counts[axis]) for axis, lo_v, hi_v in EDGES.values()}
    return Blocking(points, [(label, tuple(range(8)))], nodes)


def split(blocking: Blocking, label: str, axis: str, fraction: float) -> Blocking:
    """Cut block *label* across *axis* at the node nearest *fraction* of it.

    Every block sharing a cut edge is cut too, at the same node, so the
    blocking stays conforming. The low part keeps the label, the high part
    gets a new one (``b<n>``).
    """
    if axis not in AXES:
        raise ParameterError(f"split_blocking: axis must be one of {AXES}, got {axis!r}.")
    if not 0.0 < fraction < 1.0:
        raise ParameterError(f"split_blocking: fraction must be in (0, 1), got {fraction}.")
    verts = blocking.block(label)
    edges = blocking.edges()
    index = {name: n for n, (name, _) in enumerate(blocking.blocks)}
    lo0, hi0 = verts[0], verts[_BIT[axis]]
    n = blocking.edge_count(lo0, hi0)
    m = round(fraction * (n - 1))
    if not 0 < m < n - 1:
        raise ParameterError(
            f"split_blocking: {label} has {n} nodes along {axis}; fraction {fraction} "
            "does not fall on an interior node."
        )
    cuts: dict[tuple[int, int], int] = {}  # edge (a < b) -> node index counted from a
    split_blocks: dict[int, tuple[str, int]] = {}  # block index -> (axis, node from its lo end)
    todo = [(index[label], axis, m)]
    while todo:
        b, ax, mm = todo.pop()
        if b in split_blocks:
            if split_blocks[b] != (ax, mm):
                raise ParameterError(
                    f"split_blocking: the cut reaches {blocking.blocks[b][0]} twice, "
                    "along different axes or nodes."
                )
            continue
        split_blocks[b] = (ax, mm)
        bv = blocking.blocks[b][1]
        for eax, lo, hi in EDGES.values():
            if eax != ax:
                continue
            a, c = bv[lo], bv[hi]
            key, from_a = ((a, c), mm) if a < c else ((c, a), blocking.edge_count(a, c) - 1 - mm)
            if cuts.setdefault(key, from_a) != from_a:
                raise ParameterError(f"split_blocking: edge {key} would be cut twice.")
            for other, oname in edges[key]:
                oax, olo, _ = EDGES[oname]
                ov = blocking.blocks[other][1]
                om = from_a if ov[olo] == key[0] else blocking.edge_count(*key) - 1 - from_a
                todo.append((other, oax, om))
    out = blocking.copy()
    new_vertex: dict[tuple[int, int], int] = {}
    points = [out.points]
    for nv, key in enumerate(sorted(cuts), start=len(out.points)):
        a, c = key
        count = blocking.nodes[key]
        t = cuts[key] / (count - 1)
        points.append(((1.0 - t) * blocking.points[a] + t * blocking.points[c])[None])
        new_vertex[key] = nv
        del out.nodes[key]
        out.nodes[(a, nv)] = cuts[key] + 1
        out.nodes[(c, nv)] = count - cuts[key]
        assoc = out.associations.pop(("e", a, c), None)
        if assoc is not None:
            out.associations[("e", a, nv)] = assoc
            out.associations[("e", c, nv)] = assoc
    out.points = np.concatenate(points)

    def cut_vertex(a: int, c: int) -> int:
        return new_vertex[(min(a, c), max(a, c))]

    blocks: list[tuple[str, tuple[int, ...]]] = []
    for b, (name, bv) in enumerate(blocking.blocks):
        if b not in split_blocks:
            blocks.append((name, bv))
            continue
        ax = split_blocks[b][0]
        bit = _BIT[ax]
        lo_part, hi_part = list(bv), list(bv)
        for v in range(8):
            if v & bit:
                lo_part[v] = cut_vertex(bv[v - bit], bv[v])
            else:
                hi_part[v] = cut_vertex(bv[v], bv[v + bit])
        hi_label = f"b{out.next_label}"
        out.next_label += 1
        blocks += [(name, tuple(lo_part)), (hi_label, tuple(hi_part))]
        # the cut face's edges take the counts of the edges they run beside
        for oax, lo, hi in EDGES.values():
            if oax != ax and not lo & bit:
                count = blocking.edge_count(bv[lo], bv[hi])
                p, q = cut_vertex(bv[lo], bv[lo + bit]), cut_vertex(bv[hi], bv[hi + bit])
                out.nodes[(min(p, q), max(p, q))] = count
        # faces containing the cut axis: both halves keep the face's association
        for face in FACES:
            if face[0] == ax:
                continue
            corners = [bv[c] for c in face_corners(face).ravel()]
            assoc = out.associations.pop(("f", *sorted(corners)), None)
            if assoc is not None:
                for part in (lo_part, hi_part):
                    key = ("f", *sorted(part[c] for c in face_corners(face).ravel()))
                    out.associations[key] = assoc
    out.blocks = blocks
    return out


def delete(blocking: Blocking, labels: list[str] | tuple[str, ...]) -> Blocking:
    """Remove blocks *labels*; their faces left alone become boundary faces."""
    for label in labels:
        blocking.block(label)
    remaining = [(name, v) for name, v in blocking.blocks if name not in labels]
    if not remaining:
        raise ParameterError("delete_blocks: cannot delete every block of the blocking.")
    out = blocking.copy()
    out.blocks = remaining
    return out


def associate(
    blocking: Blocking, label: str, entity: str, association: Association | None
) -> Blocking:
    """Set (or with None clear) the association of *entity* of block *label*."""
    key, _ = blocking.entity(label, entity)
    out = blocking.copy()
    if association is None:
        out.associations.pop(key, None)
    else:
        out.associations[key] = association
    return out


# -- premesh ------------------------------------------------------------------------


def effective_associations(blocking: Blocking) -> dict[EntityKey, Association]:
    """Every entity's association: its own, else inherited (module docstring)."""
    own = blocking.associations
    faces = blocking.faces()
    edges = blocking.edges()
    out: dict[EntityKey, Association] = {}
    face_of_edge: dict[tuple[int, int], set] = {e: set() for e in edges}
    for fkey in faces:
        assoc = own.get(("f", *fkey))
        if assoc is None:
            continue
        out[("f", *fkey)] = assoc
        for a in fkey:
            for b in fkey:
                if (a, b) in face_of_edge:
                    face_of_edge[(a, b)].add(assoc)
    edge_of_vertex: dict[int, set] = {}
    face_of_vertex: dict[int, set] = {}
    for (a, b), from_faces in face_of_edge.items():
        assoc = own.get(("e", a, b))
        if assoc is not None:
            edge_of_vertex.setdefault(a, set()).add(assoc)
            edge_of_vertex.setdefault(b, set()).add(assoc)
        elif len(from_faces) == 1:
            assoc = next(iter(from_faces))
        if assoc is not None:
            out[("e", a, b)] = assoc
        for v in (a, b):
            face_of_vertex.setdefault(v, set()).update(from_faces)
    for v in {v for _, verts in blocking.blocks for v in verts}:
        assoc = own.get(("v", v))
        if assoc is None:
            candidates = edge_of_vertex.get(v) or face_of_vertex.get(v) or set()
            assoc = next(iter(candidates)) if len(candidates) == 1 else None
        if assoc is not None:
            out[("v", v)] = assoc
    return out


def _orient(grid: np.ndarray, src: np.ndarray, dst: np.ndarray) -> np.ndarray:
    """*grid* (na, nb, 3) with corners *src* (2, 2), re-indexed to corners *dst*."""
    for transpose in (False, True):
        g, s = (grid.transpose(1, 0, 2), src.T) if transpose else (grid, src)
        for fa in (False, True):
            for fb in (False, True):
                ss = s[::-1] if fa else s
                ss = ss[:, ::-1] if fb else ss
                if np.array_equal(ss, dst):
                    gg = g[::-1] if fa else g
                    return gg[:, ::-1] if fb else gg
    raise ParameterError(f"blocking: face corners {src.tolist()} and {dst.tolist()} differ.")


def _edge_grid(points: np.ndarray, n: int, assoc: Any, resolver: Any) -> np.ndarray:
    """Equally spaced nodes from points[0] to points[1], projected onto *assoc*."""
    from shore.mesh.redistribute import resample

    t = np.linspace(0.0, 1.0, n)
    line = (1.0 - t)[:, None] * points[0] + t[:, None] * points[1]
    if assoc is not None and n > 2:
        for _ in range(3):
            line[1:-1] = resolver.project(assoc, line[1:-1])
            line = resample(line, t)
        line[1:-1] = resolver.project(assoc, line[1:-1])
    return line


def tfi_face(i_lo: np.ndarray, i_hi: np.ndarray, j_lo: np.ndarray, j_hi: np.ndarray) -> np.ndarray:
    """Coons patch (na, nb, 3) from its four boundary rows.

    *i_lo* / *i_hi* are the rows at a = 0 / a = na - 1 (length nb), *j_lo* /
    *j_hi* those at b = 0 / b = nb - 1 (length na). The boundary is copied.
    """
    na, nb = len(j_lo), len(i_lo)
    u = np.linspace(0.0, 1.0, na)[:, None, None]
    v = np.linspace(0.0, 1.0, nb)[None, :, None]
    grid = (
        (1 - u) * i_lo[None] + u * i_hi[None] + (1 - v) * j_lo[:, None] + v * j_hi[:, None]
        - (1 - u) * (1 - v) * j_lo[0] - u * (1 - v) * j_lo[-1]
        - (1 - u) * v * j_hi[0] - u * v * j_hi[-1]
    )  # fmt: skip
    grid[0], grid[-1], grid[:, 0], grid[:, -1] = i_lo, i_hi, j_lo, j_hi
    return grid


def tfi_volume(x: np.ndarray) -> None:
    """Fill the interior of *x* (ni, nj, nk, 3) by TFI of its six faces, in place."""
    ni, nj, nk, _ = x.shape
    u = np.linspace(0.0, 1.0, ni)[:, None, None, None]
    v = np.linspace(0.0, 1.0, nj)[None, :, None, None]
    w = np.linspace(0.0, 1.0, nk)[None, None, :, None]
    uu, vv, ww = (1 - u, u), (1 - v, v), (1 - w, w)
    ends = (0, -1)
    faces = (
        sum(uu[a] * x[ends[a]][None] for a in (0, 1))
        + sum(vv[a] * x[:, ends[a]][:, None] for a in (0, 1))
        + sum(ww[a] * x[:, :, ends[a]][:, :, None] for a in (0, 1))
    )
    edges = (
        sum(uu[a] * vv[b] * x[ends[a], ends[b]][None, None] for a in (0, 1) for b in (0, 1))
        + sum(uu[a] * ww[c] * x[ends[a], :, ends[c]][None, :, None] for a in (0, 1) for c in (0, 1))
        + sum(vv[b] * ww[c] * x[:, ends[b], ends[c]][:, None, None] for b in (0, 1) for c in (0, 1))
    )
    corners = sum(
        uu[a] * vv[b] * ww[c] * x[ends[a], ends[b], ends[c]]
        for a in (0, 1)
        for b in (0, 1)
        for c in (0, 1)
    )
    x[1:-1, 1:-1, 1:-1] = (faces - edges + corners)[1:-1, 1:-1, 1:-1]


def premesh(blocking: Blocking, resolver: AssociationResolver | None) -> list[HexBlock]:
    """The blocking's blocks as ``HexBlock``\\ s (module docstring)."""
    from shore.mesh.block import HexBlock

    assoc = effective_associations(blocking)
    if assoc and resolver is None:
        raise ParameterError("premesh: the blocking has associations but no geometry resolver.")

    def project(key: EntityKey, pts: np.ndarray) -> np.ndarray:
        a = assoc.get(key)
        return pts if a is None else resolver.project(a, pts)  # type: ignore[union-attr]

    used = sorted({v for _, verts in blocking.blocks for v in verts})
    points = blocking.points.copy()
    for v in used:
        points[v] = project(("v", v), blocking.points[v])
    edge_grids = {
        (a, b): _edge_grid(points[[a, b]], blocking.nodes[(a, b)], assoc.get(("e", a, b)), resolver)
        for a, b in blocking.edges()
    }

    def edge(a: int, b: int) -> np.ndarray:
        return edge_grids[(a, b)] if a < b else edge_grids[(b, a)][::-1]

    face_grids: dict[tuple[int, ...], tuple[np.ndarray, np.ndarray]] = {}
    for fkey, owners in blocking.faces().items():
        b, name = owners[0]
        c = np.array(blocking.blocks[b][1])[face_corners(name)]
        grid = tfi_face(edge(c[0, 0], c[0, 1]), edge(c[1, 0], c[1, 1]),
                        edge(c[0, 0], c[1, 0]), edge(c[0, 1], c[1, 1]))  # fmt: skip
        face_assoc = assoc.get(("f", *fkey))
        if face_assoc is not None:
            grid[1:-1, 1:-1] = resolver.project(face_assoc, grid[1:-1, 1:-1])  # type: ignore[union-attr]
        face_grids[fkey] = (grid, c)

    out: list[HexBlock] = []
    for label, verts in blocking.blocks:
        ni, nj, nk = blocking.shape(label)
        x = np.zeros((ni, nj, nk, 3))
        v = np.array(verts)
        for name in FACES:
            fkey = tuple(sorted(v[face_corners(name)].ravel()))
            grid, src = face_grids[fkey]
            _face_view(x, name)[...] = _orient(grid, src, v[face_corners(name)])
        tfi_volume(x)
        out.append(HexBlock(label, x, periodic_i=False, periodic_j=False))
    _wire(blocking, out, assoc)
    return out


def _face_view(x: np.ndarray, name: str) -> np.ndarray:
    axis, side = "ijk".index(name[0]), (0 if name[2:] == "lo" else -1)
    return x[(slice(None),) * axis + (side,)]


def _wire(blocking: Blocking, blocks: list[HexBlock], assoc: dict[EntityKey, Any]) -> None:
    """Face BCs, seams with their axis maps, and associations of the premeshed blocks."""
    for fkey, owners in blocking.faces().items():
        if len(owners) > 2:
            raise ParameterError(f"premesh: face {fkey} belongs to {len(owners)} blocks.")
        if len(owners) == 1:
            ((b, name),) = owners
            bc = FaceBC.WALL if ("f", *fkey) in assoc else FaceBC.FREE
            blocks[b]._set_face_bc(name, bc)  # type: ignore[arg-type]
            continue
        (b0, f0), (b1, f1) = owners
        c0 = np.array(blocking.blocks[b0][1])[face_corners(f0)]
        c1 = np.array(blocking.blocks[b1][1])[face_corners(f1)]
        blocks[b0]._set_face_bc(f0, FaceBC.SHARED, blocks[b1].face(f1), _axis_map(c0, c1, f1))  # type: ignore[arg-type]
        blocks[b1]._set_face_bc(f1, FaceBC.SHARED, blocks[b0].face(f0), _axis_map(c1, c0, f0))  # type: ignore[arg-type]
    for (label, verts), block in zip(blocking.blocks, blocks, strict=True):
        for face in block.faces():
            key, _ = blocking.entity(label, face.name)
            face.association = assoc.get(key)
        for e in block.edges():
            _, lo, hi = EDGES[e.name]
            a, b = sorted((verts[lo], verts[hi]))
            e.association = assoc.get(("e", a, b))
        for vertex in block.vertices():
            vertex.association = assoc.get(("v", verts[vertex.index]))


def _axis_map(mine: np.ndarray, theirs: np.ndarray, partner_face: str) -> AxisMap:
    """How a face with corners *mine* maps onto its partner's, corners *theirs*."""
    partner_axes = FACE_ALONG_AXES[partner_face]  # type: ignore[index]
    pos = {int(theirs[a, b]): (a, b) for a in (0, 1) for b in (0, 1)}
    pairs = []
    for end in (mine[1, 0], mine[0, 1]):
        p0, p1 = pos[int(mine[0, 0])], pos[int(end)]
        d = 0 if p0[0] != p1[0] else 1
        pairs.append((partner_axes[d], bool(p1[d] < p0[d])))
    return AxisMap.from_pairs(pairs[0], pairs[1])
