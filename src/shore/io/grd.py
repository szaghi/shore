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


def merge_grd(inputs: list[str | Path], out: str | Path) -> Path:
    """Concatenate multiple ``.grd`` files into one.

    Used when an Xall run combines several independent Chimera
    components (e.g. wall + background): the merged ``.grd`` lists all
    blocks of *inputs[0]* first, then *inputs[1]*, etc., with the header
    block count rewritten to the sum.  Block payloads are passed
    through byte-for-byte; no node values are touched.

    The companion :func:`shore.io.cc_par.merge_cc_par` shifts patch
    block-indices and family pointers to match this concatenation.
    Invoke them as a pair, with the **same** input ordering.

    Parameters
    ----------
    inputs :
        Two-or-more input ``.grd`` paths, in the desired output order.
    out :
        Output path; ``.grd`` extension is appended if absent.

    Returns
    -------
    Path
        The path actually written.

    Raises
    ------
    ValueError
        If ``inputs`` has fewer than 2 entries, or any input is empty
        / malformed.
    """
    if len(inputs) < 2:
        raise ValueError(f"merge_grd: need at least 2 input files, got {len(inputs)}")

    out_path = Path(out)
    if out_path.suffix.lower() != ".grd":
        out_path = out_path.with_suffix(".grd")

    parsed = [_parse_grd_for_merge(Path(p)) for p in inputs]
    total_nb = sum(p[0] for p in parsed)

    # ``.grd`` layout is [header][all shape records][all coord records].
    # Naive byte-level concatenation would interleave shapes and coords
    # of consecutive files; we must re-interleave instead: emit every
    # input's shape block first, then every input's coord block.
    with out_path.open("wb") as f:
        _write_record(f, np.array([total_nb], dtype=np.int32))
        for _nb, shapes, _coords in parsed:
            f.write(shapes)
        for _nb, _shapes, coords in parsed:
            f.write(coords)

    return out_path


def _parse_grd_for_merge(path: Path) -> tuple[int, bytes, bytes]:
    """Return ``(n_blocks, shape_records_bytes, coord_records_bytes)``.

    ``.grd`` layout: header record (nb int32), then ``nb`` shape records
    of 3 int32 each, then ``nb`` blocks of 3 coord records (X, Y, Z) of
    ``ni*nj*nk`` float64 each.  We split at the shape/coord boundary so
    :func:`merge_grd` can interleave shapes-then-coords across multiple
    inputs without decoding any node data.

    Each Fortran unformatted record is wrapped in a 4-byte length
    marker on each side; we use the marker to skip the payload byte
    range.
    """
    if not path.exists():
        raise FileNotFoundError(f"merge_grd: {path} not found")
    data = path.read_bytes()
    if len(data) < 12:
        raise ValueError(f"merge_grd: {path} is too short to be a .grd file")

    header_len = struct.unpack("=i", data[0:4])[0]
    if header_len != 4:
        raise ValueError(
            f"merge_grd: {path} header record length is {header_len} "
            "(expected 4 for a single int32 n_blocks); not a .grd file?"
        )
    nb = struct.unpack("=i", data[4:8])[0]
    trailer = struct.unpack("=i", data[8:12])[0]
    if trailer != 4:
        raise ValueError(f"merge_grd: {path} header trailer is {trailer} (expected 4)")
    if nb <= 0:
        raise ValueError(f"merge_grd: {path} reports {nb} blocks (must be > 0)")

    # nb shape records, each (4 + 12 + 4) = 20 bytes.
    shapes_start = 12
    shapes_end = shapes_start + nb * 20
    if shapes_end > len(data):
        raise ValueError(
            f"merge_grd: {path} truncated in shape section "
            f"(need {shapes_end} bytes, have {len(data)})"
        )

    # The remaining bytes are exactly the coord records (X/Y/Z per block).
    return nb, data[shapes_start:shapes_end], data[shapes_end:]
