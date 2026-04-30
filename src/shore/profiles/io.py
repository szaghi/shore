"""Selig ``.dat`` reader and writer for 2D airfoil sections.

The Selig format is the de-facto standard used by the UIUC airfoil
database and most CFD pre-processors:

- An optional first line with the airfoil name (any non-numeric text).
- ``N`` rows of two whitespace-separated floats: ``x y``.
- Coordinates are normalised to chord = 1, ordered upper-TE -> LE
  -> lower-TE.

This module is the single entry / exit point between SHORE's in-memory
``(N, 2)`` numpy arrays and on-disk ``.dat`` files.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from shore.exceptions import MeshInputError, ParameterError


def read_section(path: str | Path) -> np.ndarray:
    """Load a Selig-format ``.dat`` airfoil section file.

    Parameters
    ----------
    path :
        Path to a ``.dat`` file.  May start with a non-numeric name
        line (UIUC convention); rows are then ``x y`` whitespace-
        separated floats.

    Returns
    -------
    numpy.ndarray, shape ``(N, 2)``
        Coordinates as read from the file, no reordering.

    Raises
    ------
    MeshInputError
        If the file is missing, empty, or contains malformed rows.
    """
    p = Path(path)
    if not p.exists():
        raise MeshInputError(f"read_section: file not found: {p}")
    coords: list[tuple[float, float]] = []
    with p.open() as f:
        for lineno, raw in enumerate(f, start=1):
            line = raw.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 2:
                # Likely a name line — skip if it's the first data
                # line, error otherwise (a stray non-numeric mid-file
                # is a malformed file, not a name).
                if not coords:
                    continue
                raise MeshInputError(f"read_section: malformed line {lineno} in {p}: {raw!r}")
            try:
                x = float(parts[0])
                y = float(parts[1])
            except ValueError:
                if not coords:
                    # Name line — skip.
                    continue
                raise MeshInputError(
                    f"read_section: non-numeric values on line {lineno} in {p}: {raw!r}"
                ) from None
            coords.append((x, y))
    if not coords:
        raise MeshInputError(f"read_section: no coordinate rows found in {p}")
    return np.asarray(coords, dtype=np.float64)


def write_section(
    path: str | Path,
    coords: np.ndarray,
    *,
    name: str | None = None,
) -> None:
    """Write a 2D section to a Selig-format ``.dat`` file.

    Parameters
    ----------
    path :
        Output path.  Parent directory must exist.
    coords :
        ``(N, 2)`` array of section coordinates in Selig order.
    name :
        Optional first-line name comment.  Default ``None`` (no
        name line written; only coordinate rows).

    Raises
    ------
    ParameterError
        If ``coords`` is not ``(N, 2)`` with ``N >= 4``.
    """
    arr = np.asarray(coords, dtype=np.float64)
    if arr.ndim != 2 or arr.shape[1] != 2:
        raise ParameterError(f"write_section: coords must have shape (N, 2), got {arr.shape}")
    if arr.shape[0] < 4:
        raise ParameterError(f"write_section: coords must have at least 4 rows, got {arr.shape[0]}")
    p = Path(path)
    with p.open("w") as f:
        if name:
            f.write(f"{name}\n")
        for x, y in arr:
            f.write(f"{x:>12.7f} {y:>12.7f}\n")
