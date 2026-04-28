"""Face and FaceBC: face descriptors for structured blocks."""

from __future__ import annotations

import enum
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Literal

import numpy as np

if TYPE_CHECKING:
    from shore.mesh.block import Block

FaceName = Literal["i_lo", "i_hi", "j_lo", "j_hi", "k_lo", "k_hi"]


class FaceBC(enum.Enum):
    """Boundary condition type for a block face."""

    WALL = "wall"  # k=0 body surface, Dirichlet from projection
    FREE = "free"  # outer boundary, marched freely
    SHARED = "shared"  # shared with another block via numpy view
    DIRICHLET = "dirichlet"  # pinned from a partner face (cap boundaries)
    PERIODIC = "periodic"  # wraps to another face on the same block


@dataclass
class Face:
    """One face of a structured block.

    Parameters
    ----------
    index :
        0-based index within the block's face list.
    name :
        Canonical face name: ``'i_lo'``, ``'i_hi'``, ``'j_lo'``, ``'j_hi'``,
        ``'k_lo'`` (wall), ``'k_hi'`` (outer).
    bc :
        Boundary condition type.
    partner :
        The partnered Face for SHARED / DIRICHLET / PERIODIC BCs.  None
        for WALL and FREE.
    body_projected :
        True when all nodes of this face lie on the k=0 body wall (i.e.
        bc == WALL).
    """

    index: int
    name: FaceName
    bc: FaceBC
    partner: Face | None = field(default=None, repr=False)
    body_projected: bool = False

    def nodes(self, block: Block) -> np.ndarray:
        """Return a view of this face's node coordinates in block.nodes."""
        return _face_slice(block, self.name)


def _face_slice(block: Block, name: FaceName) -> np.ndarray:
    """Extract the node sub-array for a named face from block.nodes."""
    nd = block.nodes  # (ni, nj, nk, 3)
    match name:
        case "i_lo":
            return nd[0, :, :, :]
        case "i_hi":
            return nd[-1, :, :, :]
        case "j_lo":
            return nd[:, 0, :, :]
        case "j_hi":
            return nd[:, -1, :, :]
        case "k_lo":
            return nd[:, :, 0, :]
        case "k_hi":
            return nd[:, :, -1, :]
        case _:
            raise ValueError(f"Unknown face name '{name}'.")
