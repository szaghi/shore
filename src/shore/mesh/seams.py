"""Point-matched seams between block faces (ADR 0035, 0036).

The face graph names one partner per face, but a face may be shared with
several smaller faces: CH's ``body.k_hi`` with its three channel blocks,
whose ``k_lo`` faces all name ``body.k_hi``. Every connected pair of faces is
taken once here, and the smaller face is located inside the larger one by
its nodes, which gives each seam its exact index ranges on both sides.

Indices follow the CGNS SIDS: 1-based, ranges ``[[imin, jmin, kmin], [imax,
jmax, kmax]]``, and a ``Transform`` ``T`` with ``Index2 = T.(Index1 - Begin1)
+ Begin2`` (entry ``a`` is ``+-(b + 1)``: index ``a`` of this block runs along
index ``b`` of the donor, reversed when negative).
"""

from __future__ import annotations

from collections.abc import Collection
from typing import TYPE_CHECKING, NamedTuple

import numpy as np

from shore.exceptions import ParameterError
from shore.mesh.face import FACE_ALONG_AXES, FaceBC

if TYPE_CHECKING:
    from shore.mesh.block import HexBlock

AXES = {"i": 0, "j": 1, "k": 2}
CONNECTED = (FaceBC.SHARED, FaceBC.DIRICHLET, FaceBC.PERIODIC)


class Seam(NamedTuple):
    """One side of a point-matched seam: ``face`` of this block over ``rng``,
    onto ``donor``'s ``donor_face`` over ``donor_rng`` (begin = the image of
    ``rng``'s begin), with ``transform``."""

    face: str
    rng: np.ndarray
    donor: str
    donor_face: str
    donor_rng: np.ndarray
    transform: tuple[int, int, int]


def seams(
    blocks: list[HexBlock], tol: float, *, only: Collection[tuple[str, str]] | None = None
) -> dict[str, list[Seam]]:
    """Per block label, its seams; every connected pair of faces once on each side.

    *tol* is the largest node distance accepted (absolute). With *only*,
    just the pairs with a side in it (``(label, face)``) are located.

    Raises
    ------
    ParameterError
        A face connected to a block outside *blocks*, or a pair whose nodes
        do not coincide.
    """
    owner = {id(face): b for b in blocks for face in b.faces()}
    out: dict[str, list[Seam]] = {b.label: [] for b in blocks}
    seen: set[frozenset[int]] = set()
    for block in blocks:
        for face in block.faces():
            if face.bc not in CONNECTED:
                continue
            donor = owner.get(id(face.partner))
            if donor is None:
                raise ParameterError(
                    f"{block.label}.{face.name} is connected to a block outside the given blocks."
                )
            key = frozenset((id(face), id(face.partner)))
            if key in seen:
                continue
            seen.add(key)
            if only is not None and not (
                (block.label, face.name) in only or (donor.label, face.partner.name) in only
            ):
                continue
            t = transform_vector(face.name, face.partner.name, face.axis_map.along)
            sides = [(block, face.name), (donor, face.partner.name)]
            if face_size(donor, face.partner.name) < face_size(block, face.name):
                sides.reverse()
                t = inverse_transform(t)
            (small, s_face), (large, l_face) = sides
            rng = point_range(small, s_face)
            begin = _locate(small, s_face, rng, large, l_face, t, tol)
            end = affine(t, rng[0], begin, rng[1])
            out[small.label].append(
                Seam(s_face, rng, large.label, l_face, np.array([begin, end], "i4"), t)
            )
            t2 = inverse_transform(t)
            rng2 = np.array([np.minimum(begin, end), np.maximum(begin, end)], "i4")
            donor2 = np.array([affine(t2, begin, rng[0], j) for j in rng2], "i4")
            out[large.label].append(Seam(l_face, rng2, small.label, s_face, donor2, t2))
    return out


def shared_faces(blocks: list[HexBlock]) -> set[tuple[str, str]]:
    """``(label, face)`` of the faces with more than one partner.

    A face is shared with several faces when more than one face names it as
    partner, or when a face naming it is not its own partner.
    """
    label = {id(f): b.label for b in blocks for f in b.faces()}
    named: dict[int, set[int]] = {}
    out: set[tuple[str, str]] = set()
    by_id = {id(f): f for b in blocks for f in b.faces()}
    for b in blocks:
        for f in b.faces():
            if f.bc in CONNECTED and f.partner is not None and id(f.partner) in label:
                named.setdefault(id(f.partner), set()).add(id(f))
    for fid, namers in named.items():
        face = by_id[fid]
        if len(namers) > 1 or any(face.partner is not by_id[n] for n in namers):
            out.add((label[fid], face.name))
    return out


def transform_vector(
    face: str, partner: str, along: tuple[tuple[str, bool], tuple[str, bool]]
) -> tuple[int, int, int]:
    """``Transform`` of a face connection from its AxisMap.

    Along the faces this is the AxisMap; across them a step into this block
    is a step out of the donor.
    """
    out = [0, 0, 0]
    for axis, (partner_axis, flipped) in zip(FACE_ALONG_AXES[face], along, strict=True):
        out[AXES[axis]] = (-1 if flipped else 1) * (AXES[partner_axis] + 1)
    normal, donor_normal = AXES[face[0]], AXES[partner[0]]
    inward = 1 if face.endswith("_lo") else -1  # +1 index step into this block
    donor_outward = 1 if partner.endswith("_hi") else -1
    out[normal] = inward * donor_outward * (donor_normal + 1)
    return out[0], out[1], out[2]


def inverse_transform(t: tuple[int, int, int]) -> tuple[int, int, int]:
    out = [0, 0, 0]
    for a, v in enumerate(t):
        out[abs(v) - 1] = (1 if v > 0 else -1) * (a + 1)
    return out[0], out[1], out[2]


def affine(
    t: tuple[int, int, int], origin: np.ndarray, image: np.ndarray, index: np.ndarray
) -> np.ndarray:
    """``T (index - origin) + image``."""
    out = np.array(image, dtype=np.int64).copy()
    for a, v in enumerate(t):
        out[abs(v) - 1] += (1 if v > 0 else -1) * (int(index[a]) - int(origin[a]))
    return out


def point_range(block: HexBlock, face: str) -> np.ndarray:
    """1-based ``[[imin, jmin, kmin], [imax, jmax, kmax]]`` of *face*."""
    shape = block.nodes.shape[:3]
    lo = [1, 1, 1]
    hi = list(shape)
    a = AXES[face[0]]
    lo[a] = hi[a] = 1 if face.endswith("_lo") else shape[a]
    return np.array([lo, hi], "i4")


def face_size(block: HexBlock, face: str) -> int:
    shape = block.nodes.shape[:3]
    return int(np.prod([n for a, n in enumerate(shape) if a != AXES[face[0]]]))


def _face_nodes(nodes: np.ndarray, rng: np.ndarray) -> np.ndarray:
    return nodes[tuple(slice(lo - 1, hi) for lo, hi in zip(rng[0], rng[1], strict=True))]


def _map_indices(
    begin: np.ndarray, transform: tuple[int, int, int], shape: tuple
) -> tuple[np.ndarray, ...]:
    """0-based donor indices of a face block of *shape* whose first node maps to *begin*."""
    grids = np.meshgrid(*(np.arange(n) for n in shape), indexing="ij")  # I - Begin
    out = [np.full(shape, int(begin[b]) - 1, dtype=np.int64) for b in range(3)]
    for a, v in enumerate(transform):
        out[abs(v) - 1] = out[abs(v) - 1] + (1 if v > 0 else -1) * grids[a]
    return tuple(out)


def _locate(
    small: HexBlock,
    s_face: str,
    rng: np.ndarray,
    large: HexBlock,
    l_face: str,
    t: tuple[int, int, int],
    tol: float,
) -> np.ndarray:
    """Donor index (1-based) of *rng*'s begin, where its nodes coincide with *large*'s."""
    here = _face_nodes(small.nodes, rng)
    extent = rng[1] - rng[0]
    face = point_range(large, l_face)
    choices = []
    for b in range(3):
        a = next(a for a in range(3) if abs(t[a]) == b + 1)
        n = int(extent[a])
        lo, hi = int(face[0][b]), int(face[1][b])
        choices.append(range(lo, hi - n + 1) if t[a] > 0 else range(lo + n, hi + 1))
    shape = large.nodes.shape[:3]
    for bi in choices[0]:
        for bj in choices[1]:
            for bk in choices[2]:
                begin = np.array([bi, bj, bk])
                idx = _map_indices(begin, t, here.shape[:3])
                if any(
                    int(i.min()) < 0 or int(i.max()) >= n for i, n in zip(idx, shape, strict=True)
                ):
                    continue
                if float(np.abs(here - large.nodes[idx]).max()) <= tol:
                    return begin
    raise ParameterError(
        f"{small.label}.{s_face} and {large.label}.{l_face} are not point-matched: "
        "their nodes do not coincide."
    )
