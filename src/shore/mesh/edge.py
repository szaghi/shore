"""Edge: arc descriptor for a structured block edge."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal

import numpy as np

from shore.mesh.spacing import Spacing1D

if TYPE_CHECKING:
    from shore.mesh.block import Block

# Canonical edge names: the two axes that are *constant* along this edge.
# E.g. "j_lo_k_lo" runs along i with j=0, k=0.
EdgeName = str


@dataclass
class Edge:
    """One edge of a structured block.

    Parameters
    ----------
    index :
        0-based index within the block's edge list.
    name :
        Canonical name encoding the two constant axes, e.g. ``"j_lo_k_lo"``.
    axis :
        The free axis along which this edge runs (``'i'``, ``'j'``, or
        ``'k'``), or ``None`` for non-axis-aligned edges (future use).
    vertex_indices :
        Indices of the two endpoint vertices in the block's vertex list.
    n_nodes :
        Number of nodes along this edge.
    spacing :
        1D spacing law for node distribution along this edge.
    body_projected :
        True when all nodes of this edge lie on the k=0 body wall.
    """

    index: int
    name: EdgeName
    axis: Literal["i", "j", "k"] | None
    vertex_indices: tuple[int, int]
    n_nodes: int
    spacing: Spacing1D
    body_projected: bool = False

    def nodes(self, block: Block) -> np.ndarray:
        """Return a (n_nodes, 3) view of this edge's node coordinates."""
        return _edge_slice(block, self.name)


def _edge_slice(block: Block, name: EdgeName) -> np.ndarray:
    """Extract the node arc for a named edge from block.nodes."""
    nd = block.nodes   # (ni, nj, nk, 3)
    match name:
        # along-i edges (j and k fixed)
        case "j_lo_k_lo":
            return nd[:, 0,  0,  :]
        case "j_hi_k_lo":
            return nd[:, -1, 0,  :]
        case "j_lo_k_hi":
            return nd[:, 0,  -1, :]
        case "j_hi_k_hi":
            return nd[:, -1, -1, :]
        # along-j edges (i and k fixed)
        case "i_lo_k_lo":
            return nd[0,  :, 0,  :]
        case "i_hi_k_lo":
            return nd[-1, :, 0,  :]
        case "i_lo_k_hi":
            return nd[0,  :, -1, :]
        case "i_hi_k_hi":
            return nd[-1, :, -1, :]
        # along-k edges (i and j fixed)
        case "i_lo_j_lo":
            return nd[0,  0,  :, :]
        case "i_hi_j_lo":
            return nd[-1, 0,  :, :]
        case "i_lo_j_hi":
            return nd[0,  -1, :, :]
        case "i_hi_j_hi":
            return nd[-1, -1, :, :]
        case _:
            raise ValueError(f"Unknown edge name '{name}'.")
