"""Face and FaceBC: face descriptors for structured blocks."""

from __future__ import annotations

import enum
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Literal

import numpy as np

if TYPE_CHECKING:
    from shore.mesh.association import Association
    from shore.mesh.block import Block

FaceName = Literal["i_lo", "i_hi", "j_lo", "j_hi", "k_lo", "k_hi"]
AxisLabel = Literal["i", "j", "k"]


# The two along-seam parametric axes for each face name, in canonical order.
# A face perpendicular to axis ``a`` has the other two axes lying *in* the
# face plane.  The first entry of each tuple is treated as the face's
# "first along-axis", the second as the "second along-axis", for the
# purpose of AxisMap.
FACE_ALONG_AXES: dict[FaceName, tuple[AxisLabel, AxisLabel]] = {
    "i_lo": ("j", "k"),
    "i_hi": ("j", "k"),
    "j_lo": ("i", "k"),
    "j_hi": ("i", "k"),
    "k_lo": ("i", "j"),
    "k_hi": ("i", "j"),
}


class FaceBC(enum.Enum):
    """Boundary condition type for a block face."""

    WALL = "wall"  # k=0 body surface, Dirichlet from projection
    FREE = "free"  # outer boundary, marched freely
    SHARED = "shared"  # shared with another block via numpy view
    DIRICHLET = "dirichlet"  # pinned from a partner face (cap boundaries)
    PERIODIC = "periodic"  # wraps to another face on the same block
    SYMMETRY = "symmetry"  # lies on a symmetry plane (half models, ADR 0013)


@dataclass(frozen=True)
class AxisMap:
    """Mapping from a face's along-seam axes onto its partner's along-seam axes.

    A face perpendicular to one block axis has two *along-seam* axes
    lying in the face plane (e.g. an ``i_lo`` face has along-axes
    ``(j, k)``).  When the face is SHARED / DIRICHLET / PERIODIC with
    another face on a partner block, every consumer (cc.par writer,
    overset connectivity, debug viewers, etc.) needs to know how the
    parametric (i, j, k) of one side maps onto the partner's (i, j, k):
    same axis or different axis, increasing or reversed.

    ``along`` is a 2-tuple of ``(partner_axis, flipped)`` pairs, in the
    order given by :data:`FACE_ALONG_AXES` for *this face's* name.
    For an ``i_lo`` face whose ``j``-row coincides with the partner's
    ``i``-row in the *same* direction, and whose ``k``-row coincides
    with the partner's ``k``-row in the same direction, the AxisMap is
    ``along=(("i", False), ("k", False))``.

    Two consistency requirements are enforced at construction:

    - The two ``partner_axis`` entries must differ (a 2D face cannot
      degenerate onto a single 1D partner axis).
    - The two ``partner_axis`` entries must be drawn from the partner
      face's own along-axes — but that is verified at wiring time by
      :meth:`AxisMap.validate_against_partner`, since :class:`AxisMap`
      itself is decoupled from any specific face name.

    Attributes
    ----------
    along :
        ``((partner_axis_for_first_along, flipped),
        (partner_axis_for_second_along, flipped))``.

    Examples
    --------
    Same-axis seam, no flips (the most common case)::

        AxisMap.identity("i_lo")  # -> along=(("j", False), ("k", False))

    Cross-axis seam with one reversal (e.g. an ``i_lo`` face whose
    ``j``-row maps to the partner's reversed ``i``-row, ``k`` preserved)::

        AxisMap.from_pairs(("i", True), ("k", False))
    """

    along: tuple[tuple[AxisLabel, bool], tuple[AxisLabel, bool]]

    def __post_init__(self) -> None:
        a0, a1 = self.along[0][0], self.along[1][0]
        if a0 == a1:
            raise ValueError(f"AxisMap: partner along-axes must differ, got {a0!r} twice")
        for entry in self.along:
            ax, flip = entry
            if ax not in ("i", "j", "k"):
                raise ValueError(f"AxisMap: invalid axis label {ax!r}")
            if not isinstance(flip, bool):
                raise TypeError(f"AxisMap: 'flipped' must be bool, got {type(flip).__name__}")

    # ------------------------------------------------------------------
    # Constructors
    # ------------------------------------------------------------------

    @classmethod
    def from_pairs(
        cls,
        first: tuple[AxisLabel, bool],
        second: tuple[AxisLabel, bool],
    ) -> AxisMap:
        """Construct an :class:`AxisMap` from two ``(axis, flipped)`` pairs."""
        return cls(along=(first, second))

    @classmethod
    def identity(cls, face_name: FaceName) -> AxisMap:
        """Same-axis, no-flip mapping for *face_name*.

        Maps each of *face_name*'s along-axes to the partner's
        same-named axis with no reversal — appropriate when both blocks
        share node arrays in identical (i, j, k) orientation (the most
        common case: same-axis SHARED seams).
        """
        a, b = FACE_ALONG_AXES[face_name]
        return cls(along=((a, False), (b, False)))

    # ------------------------------------------------------------------
    # Inspection helpers
    # ------------------------------------------------------------------

    def partner_axes(self) -> tuple[AxisLabel, AxisLabel]:
        """Return the two partner along-axis labels (without the flip flags)."""
        return self.along[0][0], self.along[1][0]

    def flips(self) -> tuple[bool, bool]:
        """Return the two flip flags (True = reversed indexing)."""
        return self.along[0][1], self.along[1][1]

    def inverse(
        self,
        this_face: FaceName,
        partner_face: FaceName,
    ) -> AxisMap:
        """Return the AxisMap for the *partner* side of the same seam.

        If ``self`` describes how *this_face*'s along-axes map onto
        *partner_face*'s along-axes, the inverse describes how
        *partner_face*'s along-axes map onto *this_face*'s along-axes.
        Both sides agree on flip parity (a flip is symmetric).
        """
        my_axes = FACE_ALONG_AXES[this_face]
        partner_axes_canon = FACE_ALONG_AXES[partner_face]
        # Build reverse mapping: partner_axis -> (my_axis, flipped)
        rev: dict[AxisLabel, tuple[AxisLabel, bool]] = {}
        for my_ax, (p_ax, flip) in zip(my_axes, self.along, strict=True):
            rev[p_ax] = (my_ax, flip)
        first = rev[partner_axes_canon[0]]
        second = rev[partner_axes_canon[1]]
        return AxisMap(along=(first, second))

    def validate_against_partner(
        self,
        this_face: FaceName,
        partner_face: FaceName,
    ) -> None:
        """Check that ``self.partner_axes()`` are exactly *partner_face*'s along-axes.

        Raises
        ------
        ValueError
            If any partner-axis label in *self* is not an along-axis of
            *partner_face*, or if the two label slots are not a
            permutation of *partner_face*'s along-axes.
        """
        partner_axes_canon = set(FACE_ALONG_AXES[partner_face])
        my_partner_axes = set(self.partner_axes())
        if my_partner_axes != partner_axes_canon:
            raise ValueError(
                f"AxisMap on {this_face!r} ↔ {partner_face!r}: partner axes "
                f"{tuple(sorted(my_partner_axes))} do not match "
                f"{partner_face!r}'s along-axes "
                f"{tuple(sorted(partner_axes_canon))}"
            )


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
    axis_map :
        Mapping from this face's along-axes onto the partner's
        along-axes.  Required when ``partner`` is set; ``None``
        otherwise.  See :class:`AxisMap`.
    bc_code :
        Optional explicit downstream natural-BC code for a ``FREE``
        face, overriding the per-block default the cc.par writer would
        otherwise apply (``BlockMeta.free_face_bc`` / ``--free-bc``).
        ``None`` means "use the per-block default".  Only meaningful on
        ``FREE`` faces; must be a positive Xall natural-BC code in
        ``[1, nccnat=19]`` (e.g. ``9`` free/far-field, ``3``/``4``
        inflow/outflow).  Lets a single block carry heterogeneous outer
        BCs — inflow on one face, extrapolation on the others.
    association :
        Geometry this face's nodes are projected onto (ADR 0003), or ``None``.
    """

    index: int
    name: FaceName
    bc: FaceBC
    # `partner` carries a back-reference to a Face on another block.  Two
    # SHARED faces point at each other, so dataclass-generated __eq__ would
    # recurse forever — exclude from equality (and identity is the right
    # comparison anyway: faces are owned by blocks, not values).
    partner: Face | None = field(default=None, repr=False, compare=False)
    axis_map: AxisMap | None = field(default=None, repr=False)
    bc_code: int | None = None
    association: Association | None = field(default=None, repr=False)

    @property
    def body_projected(self) -> bool:
        """True when associated with a whole body (legacy flag, read-only)."""
        return self.association is not None and self.association.kind.value == "body"

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
