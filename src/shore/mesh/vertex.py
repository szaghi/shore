"""Vertex: corner descriptor for a structured block."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from shore.mesh.association import Association
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
    association :
        Geometry this corner is projected onto (ADR 0003), or ``None``.
    """

    index: int
    ijk: tuple[int, int, int]
    association: Association | None = None

    @property
    def body_projected(self) -> bool:
        """True when associated with a whole body (legacy flag, read-only)."""
        return self.association is not None and self.association.kind.value == "body"

    def coords(self, block: Block) -> np.ndarray:
        """Return the (3,) coordinate vector for this vertex."""
        i, j, k = self.ijk
        return block.nodes[i, j, k, :]
