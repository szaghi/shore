"""Xall binary ``.grd`` writer.

Format
------
Fortran unformatted sequential binary (each ``write`` statement wraps its
payload in 4-byte record-length markers):

  Record 1   : nb                          - int32, number of blocks
  Records 2..nb+1 : ni_n, nj_n, nk_n      - 3x int32 per block,
                                             stored as *cell* counts
                                             (= node count - 1)
  For each block n = 1..nb:
    Record nb+1+3(n-1)+1 : X flat          - ni_n*nj_n*nk_n float64,
                                             i-fastest (Fortran column-major)
    Record nb+1+3(n-1)+2 : Y flat
    Record nb+1+3(n-1)+3 : Z flat

Reference implementation: ``geogrd/src/geogrd.f90`` in
``Xnavis-projects/geogrd``.

Note on node vs. cell counts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
``geogrd`` reads ``.geo`` files whose first line is the *node* count
``(ni, nj, nk)`` (the SHORE convention) and immediately subtracts 1
before writing to ``.grd`` (``ni = ni - 1``).  SHORE ``(ni, nj, nk, 3)``
arrays therefore map directly: pass the array as-is and this writer
performs the subtraction.
"""

from __future__ import annotations

import struct
from pathlib import Path

import numpy as np


def write_grd(
    path: str | Path,
    geo_paths: list[str | Path],
) -> Path:
    """Write multiple structured blocks to an Xall binary ``.grd`` file.

    Reads each block from an ASCII ``.geo`` file and packs all blocks into
    one Fortran unformatted binary ``.grd``.  The input list determines both
    the block count and block ordering, which must match the ordering used
    when generating the companion ``cc.par``.

    Parameters
    ----------
    path :
        Output file path.  ``.grd`` extension is appended if absent.
    geo_paths :
        Ordered list of ``.geo`` file paths, one per block.

    Returns
    -------
    Path
        The path actually written.
    """
    from shore.io.geo import read_geo

    out = Path(path)
    if out.suffix.lower() != ".grd":
        out = out.with_suffix(".grd")

    nb = len(geo_paths)
    if nb == 0:
        raise ValueError("write_grd: geo_paths list must not be empty")

    grids = [read_geo(p) for p in geo_paths]

    with out.open("wb") as f:
        _write_record(f, np.array([nb], dtype=np.int32))

        for grid in grids:
            ni, nj, nk, _ = grid.shape
            # Store cell counts (= node counts - 1), matching geogrd convention.
            _write_record(f, np.array([ni - 1, nj - 1, nk - 1], dtype=np.int32))

        for grid in grids:
            coords = grid.astype(np.float64, copy=False)
            for axis in range(3):  # x, y, z
                flat = coords[:, :, :, axis].ravel(order="F")
                _write_record(f, flat)

    return out


def _write_record(f, data: np.ndarray) -> None:
    """Write one Fortran unformatted record (4-byte length markers)."""
    nbytes = data.nbytes
    marker = struct.pack("=i", nbytes)  # native-endian int32
    f.write(marker)
    f.write(data.tobytes())
    f.write(marker)
