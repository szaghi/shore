"""Plot3D export of structured block lists (ADR 0035).

Multi-block, whole (no IBLANK), 3-D grid file, Fortran unformatted
sequential, double precision, little endian (the layout ParaView's and
VTK's PLOT3D readers and most CFD tools read):

  Record 1 : nblocks                       - int32
  Record 2 : (ni, nj, nk) per block        - 3 * nblocks int32, *node* counts
  Record 2 + n : X, Y, Z of block n        - 3 * ni * nj * nk float64,
                                             each i-fastest

Plot3D carries no block names, connectivity or boundary conditions; CGNS
(:mod:`shore.io.cgns`) does.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import TYPE_CHECKING, BinaryIO

import numpy as np

from shore.exceptions import ParameterError

if TYPE_CHECKING:
    from shore.mesh.block import HexBlock


def write_plot3d(blocks: Sequence[HexBlock], path: str | Path) -> Path:
    """Write *blocks* as one multi-block Plot3D grid file at *path*."""
    if not blocks:
        raise ParameterError("plot3d: no blocks to write.")
    path = Path(path)
    dims = np.array([b.nodes.shape[:3] for b in blocks], dtype="<i4").ravel()
    with path.open("wb") as f:
        _record(f, np.array([len(blocks)], dtype="<i4"))
        _record(f, dims)
        for b in blocks:
            xyz = np.concatenate([b.nodes[..., c].ravel(order="F") for c in range(3)])
            _record(f, xyz.astype("<f8"))
    return path


def _record(f: BinaryIO, data: np.ndarray) -> None:
    marker = np.array([data.nbytes], dtype="<i4").tobytes()
    f.write(marker)
    f.write(data.tobytes())
    f.write(marker)
