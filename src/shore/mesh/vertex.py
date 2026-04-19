"""Vertex: corner descriptor for a structured block."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from shore.mesh.block import Block


@dataclass
class Vertex:
    """One corner of a structured block.

    Parameters
    ----------
    index :
        0-based index within the block's vertex list.
    ijk :
        ``(i, j, k)`` corner indices into ``block.nodes``.  Use ``-1`` for
        the last index along an axis (standard Python convention).
    body_projected :
        True when this vertex lies on the k=0 body wall.
    """

    index: int
    ijk: tuple[int, int, int]
    body_projected: bool = False

    def coords(self, block: Block) -> np.ndarray:
        """Return the (3,) coordinate vector for this vertex."""
        i, j, k = self.ijk
        return block.nodes[i, j, k, :]
