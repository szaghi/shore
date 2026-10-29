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
from shore.mesh.face import FaceBC
from shore.mesh.seams import inverse_transform as _inverse  # noqa: F401  (re-export)
from shore.mesh.seams import point_range, seams, transform_vector  # noqa: F401

if TYPE_CHECKING:
    from shore.mesh.block import HexBlock

CGNS_VERSION = 3.4  # nothing beyond the 3.4 SIDS is written
NAME_LENGTH = 32
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
    points = np.concatenate([b.nodes.reshape(-1, 3) for b in blocks])
    scale = float(np.linalg.norm(points.max(axis=0) - points.min(axis=0))) or 1.0
    try:
        links = seams(blocks, tolerance * scale)
    except ParameterError as exc:
        raise ParameterError(f"cgns: base {name!r}: {exc}") from exc
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
            for seam in connections:
                conn_name = _unique(f"{seam.face}_{seam.donor}", names)
                conn = _node(zgc, conn_name, "GridConnectivity1to1_t", _chars(seam.donor))
                transform = np.array(seam.transform, "i4")
                _node(conn, "Transform", '"int[IndexDimension]"', transform)
                _node(conn, "PointRange", "IndexRange_t", seam.rng)
                _node(conn, "PointRangeDonor", "IndexRange_t", seam.donor_rng)
        others = [f for f in block.faces() if f.bc in _BC_TYPE]
        if others:
            zbc = _node(zone, "ZoneBC", "ZoneBC_t")
            for face in others:
                bc = _node(zbc, face.name, "BC_t", _chars(_BC_TYPE[face.bc]))
                _node(bc, "PointRange", "IndexRange_t", point_range(block, face.name))
                if face.bc is FaceBC.FREE:
                    _node(bc, "SHORE", "Descriptor_t", _chars("shore:FREE"))


def _unique(name: str, used: set[str]) -> str:
    candidate, n = name[:NAME_LENGTH], 1
    while candidate in used:
        n += 1
        suffix = f"_{n}"
        candidate = name[: NAME_LENGTH - len(suffix)] + suffix
    used.add(candidate)
    return candidate


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
