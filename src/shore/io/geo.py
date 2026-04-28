"""Xall_GPU geo format writer/reader.

Format (ASCII, one block per file)
-----------------------------------
ni nj nk
x(0,0,0)  x(1,0,0)  ...   (all X, i-major, then j, then k)
y(0,0,0)  ...
z(0,0,0)  ...

Coordinate ordering matches Fortran column-major: i is the fastest index.
The grid array passed in must have shape (ni, nj, nk, 3) with the last
axis being (x, y, z).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from shore.exceptions import MeshInputError


def write_geo(path: str | Path, grid: np.ndarray) -> None:
    """Write a single structured block to an ASCII geo file.

    Parameters
    ----------
    path :
        Output file path.
    grid :
        Coordinates array, shape (ni, nj, nk, 3).
    """
    path = Path(path)
    ni, nj, nk, _ = grid.shape

    with path.open("w") as f:
        f.write(f"{ni} {nj} {nk}\n")
        for coord in range(3):  # x, then y, then z
            flat = grid[:, :, :, coord].ravel(order="F")  # Fortran (i-major)
            # 6 values per line for readability
            for start in range(0, len(flat), 6):
                line = "  ".join(f"{v:.10e}" for v in flat[start : start + 6])
                f.write(line + "\n")


def read_geo(path: str | Path) -> np.ndarray:
    """Read a single structured block from an ASCII geo file.

    Parameters
    ----------
    path :
        Path to the .geo file.

    Returns
    -------
    numpy.ndarray, shape (ni, nj, nk, 3)
    """
    path = Path(path)
    with path.open("r") as f:
        ni, nj, nk = map(int, f.readline().split())
        data = np.fromstring(f.read(), sep=" ")

    n = ni * nj * nk
    if data.size != 3 * n:
        raise MeshInputError(f"Expected {3 * n} values for ({ni},{nj},{nk}) grid, got {data.size}.")

    grid = np.empty((ni, nj, nk, 3), dtype=np.float64)
    for coord in range(3):
        grid[:, :, :, coord] = data[coord * n : (coord + 1) * n].reshape((ni, nj, nk), order="F")

    return grid
