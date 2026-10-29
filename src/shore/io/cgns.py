"""CGNS export of structured block lists (ADR 0035).

One CGNS Base per component, one structured Zone per block, written in the
CGNS/HDF5 file mapping with ``h5py`` (optional dependency, ``shore[cgns]``).

- ``GridCoordinates``: ``CoordinateX/Y/Z``, double precision.
- ``ZoneGridConnectivity``: one ``GridConnectivity1to1_t`` per connected
  face (SHARED, DIRICHLET, PERIODIC), with ``PointRange``,
  ``PointRangeDonor`` and ``Transform`` derived from the face's
  :class:`~shore.mesh.face.AxisMap`. CGNS 1-to-1 connections are point
  matched: a face whose nodes do not coincide with its partner's is refused.
- ``ZoneBC``: one ``BC_t`` per other face: WALL is ``BCWall``, SYMMETRY is
  ``BCSymmetryPlane`` and FREE (an overset fringe or a natural boundary,
  decided downstream) is ``BCTypeUserDefined`` with the descriptor
  ``shore:FREE``.

Indices are 1-based and Fortran-ordered as the SIDS require; an array of
SIDS dimensions ``(d1, d2, ...)`` is stored as an HDF5 dataset of shape
``(..., d2, d1)``. The O-grid's periodic j seam is not a face of the block
graph (as in ``cc.par``), so it is not written either.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np

from shore.exceptions import ParameterError
from shore.mesh.face import FACE_ALONG_AXES, FaceBC

if TYPE_CHECKING:
    from shore.mesh.block import HexBlock

CGNS_VERSION = 3.4  # nothing beyond the 3.4 SIDS is written
NAME_LENGTH = 32
_AXES = {"i": 0, "j": 1, "k": 2}
_CONNECTED = (FaceBC.SHARED, FaceBC.DIRICHLET, FaceBC.PERIODIC)
_BC_TYPE = {
    FaceBC.WALL: "BCWall",
    FaceBC.SYMMETRY: "BCSymmetryPlane",
    FaceBC.FREE: "BCTypeUserDefined",
}


def write_cgns(
    components: Mapping[str, Sequence[HexBlock]],
    path: str | Path,
    *,
    tolerance: float = 1e-9,
) -> Path:
    """Write *components* (``{base name: blocks}``) as one CGNS/HDF5 file.

    Parameters
    ----------
    components :
        Base name -> its blocks; the bases are written in this order and
        the zones in block order. Connections must stay inside a base.
    path :
        Output file (conventionally ``.cgns``).
    tolerance :
        Largest distance between a connected face's node and its partner's,
        relative to the base's bounding-box diagonal.

    Raises
    ------
    ParameterError
        A base or zone name longer than 32 characters or repeated, a
        connection to a block outside its base, or a connection whose
        nodes do not coincide.
    """
    h5py = _h5py()
    path = Path(path)
    _check_names(list(components), "base")
    with h5py.File(path, "w", track_order=True, libver="earliest") as f:
        _root(f)
        _node(f, "CGNSLibraryVersion", "CGNSLibraryVersion_t", np.array([CGNS_VERSION], "f4"))
        for name, blocks in components.items():
            _write_base(f, name, list(blocks), tolerance)
    return path


def _write_base(f: Any, name: str, blocks: list[HexBlock], tolerance: float) -> None:
    labels = [b.label for b in blocks]
    _check_names(labels, "zone")
    owner = {id(face): b for b in blocks for face in b.faces()}
    points = np.concatenate([b.nodes.reshape(-1, 3) for b in blocks])
    scale = float(np.linalg.norm(points.max(axis=0) - points.min(axis=0))) or 1.0
    links = _links(blocks, owner, name, tolerance * scale)
    base = _node(f, name, "CGNSBase_t", np.array([3, 3], "i4"))
    for block in blocks:
        ni, nj, nk = block.nodes.shape[:3]
        sizes = np.array([[ni, nj, nk], [ni - 1, nj - 1, nk - 1], [0, 0, 0]], "i4")
        zone = _node(base, block.label, "Zone_t", sizes)
        _node(zone, "ZoneType", "ZoneType_t", _chars("Structured"))
        coords = _node(zone, "GridCoordinates", "GridCoordinates_t")
        for c, axis in enumerate("XYZ"):
            data = np.ascontiguousarray(block.nodes[..., c].transpose(2, 1, 0), dtype="f8")
            _node(coords, f"Coordinate{axis}", "DataArray_t", data)
        connections = links[block.label]
        if connections:
            zgc = _node(zone, "ZoneGridConnectivity", "ZoneGridConnectivity_t")
            names: set[str] = set()
            for face_name, rng, donor, rng_donor, transform in connections:
                conn_name = _unique(f"{face_name}_{donor}", names)
                conn = _node(zgc, conn_name, "GridConnectivity1to1_t", _chars(donor))
                _node(conn, "Transform", '"int[IndexDimension]"', np.array(transform, "i4"))
                _node(conn, "PointRange", "IndexRange_t", rng)
                _node(conn, "PointRangeDonor", "IndexRange_t", rng_donor)
        others = [f for f in block.faces() if f.bc in _BC_TYPE]
        if others:
            zbc = _node(zone, "ZoneBC", "ZoneBC_t")
            for face in others:
                bc = _node(zbc, face.name, "BC_t", _chars(_BC_TYPE[face.bc]))
                _node(bc, "PointRange", "IndexRange_t", _point_range(block, face.name))
                if face.bc is FaceBC.FREE:
                    _node(bc, "SHORE", "Descriptor_t", _chars("shore:FREE"))


_Link = tuple[str, np.ndarray, str, np.ndarray, tuple[int, int, int]]


def _links(
    blocks: list[HexBlock], owner: dict[int, HexBlock], base: str, tol: float
) -> dict[str, list[_Link]]:
    """Per zone, its 1-to-1 connections ``(face, range, donor, donor range, transform)``.

    Every connected pair of faces is written once in each of its two zones.
    The smaller face is located inside the larger one by its nodes: a face
    may be shared with several smaller faces (CH's ``body.k_hi`` with its
    three channel blocks), while the face graph names one partner per face.
    """
    out: dict[str, list[_Link]] = {b.label: [] for b in blocks}
    seen: set[frozenset[int]] = set()
    for block in blocks:
        for face in block.faces():
            if face.bc not in _CONNECTED:
                continue
            donor = owner.get(id(face.partner))
            if donor is None:
                raise ParameterError(
                    f"cgns: {block.label}.{face.name} is connected to a block outside "
                    f"base {base!r}."
                )
            key = frozenset((id(face), id(face.partner)))
            if key in seen:
                continue
            seen.add(key)
            t = transform_vector(face.name, face.partner.name, face.axis_map.along)
            sides = [(block, face.name), (donor, face.partner.name)]
            if _face_size(donor, face.partner.name) < _face_size(block, face.name):
                sides.reverse()
                t = _inverse(t)
            (small, s_face), (large, l_face) = sides
            rng = _point_range(small, s_face)
            begin = _locate(small, s_face, rng, large, l_face, t, tol)
            end = _affine(t, rng[0], begin, rng[1])
            out[small.label].append((s_face, rng, large.label, np.array([begin, end], "i4"), t))
            t2 = _inverse(t)
            rng2 = np.array([np.minimum(begin, end), np.maximum(begin, end)], "i4")
            donor2 = [_affine(t2, begin, rng[0], j) for j in rng2]
            out[large.label].append((l_face, rng2, small.label, np.array(donor2, "i4"), t2))
    return out


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
    face = _point_range(large, l_face)
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
        f"cgns: {small.label}.{s_face} and {large.label}.{l_face} are not point-matched; "
        "CGNS 1-to-1 connections need coincident nodes."
    )


def _face_size(block: HexBlock, face: str) -> int:
    shape = block.nodes.shape[:3]
    return int(np.prod([n for a, n in enumerate(shape) if a != _AXES[face[0]]]))


def _inverse(t: tuple[int, int, int]) -> tuple[int, int, int]:
    out = [0, 0, 0]
    for a, v in enumerate(t):
        out[abs(v) - 1] = (1 if v > 0 else -1) * (a + 1)
    return out[0], out[1], out[2]


def _affine(
    t: tuple[int, int, int], origin: np.ndarray, image: np.ndarray, index: np.ndarray
) -> np.ndarray:
    """``T (index - origin) + image`` (SIDS ``Index2 = T.(Index1 - Begin1) + Begin2``)."""
    out = np.array(image, dtype=np.int64).copy()
    for a, v in enumerate(t):
        out[abs(v) - 1] += (1 if v > 0 else -1) * (int(index[a]) - int(origin[a]))
    return out


def _unique(name: str, used: set[str]) -> str:
    candidate, n = name[:NAME_LENGTH], 1
    while candidate in used:
        n += 1
        suffix = f"_{n}"
        candidate = name[: NAME_LENGTH - len(suffix)] + suffix
    used.add(candidate)
    return candidate


def transform_vector(
    face: str, partner: str, along: tuple[tuple[str, bool], tuple[str, bool]]
) -> tuple[int, int, int]:
    """CGNS ``Transform`` of a 1-to-1 face connection from its AxisMap.

    Entry ``a`` is ``+-(b + 1)``: index ``a`` of this zone runs along index
    ``b`` of the donor, reversed when negative. Along the faces this is the
    AxisMap; across them a step into this zone is a step out of the donor.
    """
    out = [0, 0, 0]
    for axis, (partner_axis, flipped) in zip(FACE_ALONG_AXES[face], along, strict=True):
        out[_AXES[axis]] = (-1 if flipped else 1) * (_AXES[partner_axis] + 1)
    normal, donor_normal = _AXES[face[0]], _AXES[partner[0]]
    inward = 1 if face.endswith("_lo") else -1  # +1 index step into this zone
    donor_outward = 1 if partner.endswith("_hi") else -1
    out[normal] = inward * donor_outward * (donor_normal + 1)
    return out[0], out[1], out[2]


def _point_range(block: HexBlock, face: str) -> np.ndarray:
    """1-based ``[[imin, jmin, kmin], [imax, jmax, kmax]]`` of *face*."""
    shape = block.nodes.shape[:3]
    lo = [1, 1, 1]
    hi = list(shape)
    a = _AXES[face[0]]
    lo[a] = hi[a] = 1 if face.endswith("_lo") else shape[a]
    return np.array([lo, hi], "i4")


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


# ---------------------------------------------------------------------------
# CGNS/HDF5 file mapping
# ---------------------------------------------------------------------------


def _h5py() -> Any:
    try:
        import h5py
    except ImportError as exc:  # pragma: no cover - exercised without h5py
        raise ParameterError("CGNS export needs h5py: pip install 'shore-mesh[cgns]'.") from exc
    return h5py


def _check_names(names: list[str], kind: str) -> None:
    for n in names:
        if not n or len(n) > NAME_LENGTH or "/" in n:
            raise ParameterError(f"cgns: {kind} name {n!r} must be 1-32 characters without '/'.")
    if len(set(names)) != len(names):
        raise ParameterError(f"cgns: {kind} names must be unique, got {names}.")


def _chars(text: str) -> np.ndarray:
    return np.frombuffer(text.encode("ascii"), dtype="i1").copy()


def _string(text: str, size: int) -> np.bytes_:
    return np.bytes_(text.encode("ascii").ljust(size, b"\0"))


def _attrs(group: Any, name: str, label: str, dtype: str) -> None:
    group.attrs.create("name", _string(name, 33), dtype="S33")
    group.attrs.create("label", _string(label, 33), dtype="S33")
    group.attrs.create("type", _string(dtype, 3), dtype="S3")
    group.attrs.create("flags", np.array([1], "i4"))


def _root(f: Any) -> None:
    _attrs(f, "HDF5 MotherNode", "Root Node of HDF5 File", "MT")
    f.create_dataset(" format", data=_chars("IEEE_LITTLE_32"))
    version = np.frombuffer(b"HDF5 Version 1.8.21".ljust(33, b"\0"), "i1")
    f.create_dataset(" hdf5version", data=version)


_TYPES = {"i1": "C1", "i4": "I4", "i8": "I8", "f4": "R4", "f8": "R8"}


def _node(parent: Any, name: str, label: str, data: np.ndarray | None = None) -> Any:
    group = parent.create_group(name, track_order=True)
    if data is None:
        _attrs(group, name, label, "MT")
        return group
    data = np.asarray(data)
    _attrs(group, name, label, _TYPES[data.dtype.str[1:]])
    group.create_dataset(" data", data=data)
    return group
