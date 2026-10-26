"""Edge: arc descriptor for a structured block edge."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal

import numpy as np

from shore.mesh.spacing import Spacing1D

if TYPE_CHECKING:
    from shore.mesh.association import Association
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
    association :
        Geometry this edge's nodes are projected onto (ADR 0003), or ``None``.
    first_step, last_step :
        Optional per-vertex first / last interval lengths.  When set, the
        topology / primitive builder treats them as a C1 matching request
        at that endpoint of the edge: the first (or last) cell along the
        edge will have this physical length, with the interior shaped by
        ``spacing``.  When unset (``None``), the edge is built purely from
        ``spacing`` with no endpoint pinning.

        Setting these does *not* require ``spacing.law == 'pinned'``: the
        consumer is free to read the steps and switch laws as needed.  At
        corners where multiple incident edges request different first/last
        steps, the builder applies a precedence rule (cap / body geometry
        wins) and emits a ``UserWarning`` naming the overridden edges.
    """

    index: int
    name: EdgeName
    axis: Literal["i", "j", "k"] | None
    vertex_indices: tuple[int, int]
    n_nodes: int
    spacing: Spacing1D
    association: Association | None = None
    first_step: float | None = None
    last_step: float | None = None

    @property
    def body_projected(self) -> bool:
        """True when associated with a whole body (legacy flag, read-only)."""
        return self.association is not None and self.association.kind.value == "body"

    def nodes(self, block: Block) -> np.ndarray:
        """Return a (n_nodes, 3) view of this edge's node coordinates."""
        return _edge_slice(block, self.name)


def _edge_slice(block: Block, name: EdgeName) -> np.ndarray:
    """Extract the node arc for a named edge from block.nodes."""
    nd = block.nodes  # (ni, nj, nk, 3)
    match name:
        # along-i edges (j and k fixed)
        case "j_lo_k_lo":
            return nd[:, 0, 0, :]
        case "j_hi_k_lo":
            return nd[:, -1, 0, :]
        case "j_lo_k_hi":
            return nd[:, 0, -1, :]
        case "j_hi_k_hi":
            return nd[:, -1, -1, :]
        # along-j edges (i and k fixed)
        case "i_lo_k_lo":
            return nd[0, :, 0, :]
        case "i_hi_k_lo":
            return nd[-1, :, 0, :]
        case "i_lo_k_hi":
            return nd[0, :, -1, :]
        case "i_hi_k_hi":
            return nd[-1, :, -1, :]
        # along-k edges (i and j fixed)
        case "i_lo_j_lo":
            return nd[0, 0, :, :]
        case "i_hi_j_lo":
            return nd[-1, 0, :, :]
        case "i_lo_j_hi":
            return nd[0, -1, :, :]
        case "i_hi_j_hi":
            return nd[-1, -1, :, :]
        case _:
            raise ValueError(f"Unknown edge name '{name}'.")
