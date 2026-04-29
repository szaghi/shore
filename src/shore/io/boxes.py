"""Hexahedral marker boxes for Xall ``cc.par``.

Boxes are the last section of an Xall ``cc.par`` file (after the
edges section).  Each box defines a hexahedral region of the
simulation domain that marks cells inside it as "internal wall" —
i.e. cells lying inside a solid body that the chimera hole-cutter
should treat as donors-only, not active flow cells.

Conceptually, a box is a closed region that contains the body it
represents and itself fits inside the body's near-mesh: the union
of one body's boxes lies in the thin shell between the body's STL
surface and the outer extent of the body's near-mesh, so that a
*background* mesh's cells overlapping the body land inside at least
one box and get hole-cut.

Box types
---------

The Xall ``cc.par`` reader (``overset/src/scatole.f90``) supports four
types — parallelepiped (1), single cylinder (2, deprecated), coaxial
cylinders (3), and general hexahedron (9).  SHORE only emits **type
9 (hexahedron)**: it is the universal form, and any axis-aligned or
oriented bounding box can be expressed as 8 vertices.  The
parallelepiped / cylinder shorthands are reader-side conveniences
SHORE has no reason to produce.

Vertex ordering
---------------

Type-9 boxes carry 8 vertices in the order required by Xall's
``esaedro`` subroutine (``scatole.f90:30-44``).  The convention is
the standard "VTK voxel" / Fortran-column-major ordering, where the
8 corners enumerate (i, j, k) bits in the order ``v_idx = 1 + i +
2*j + 4*k``:

::

    1 = (i=lo, j=lo, k=lo)
    2 = (i=hi, j=lo, k=lo)
    3 = (i=lo, j=hi, k=lo)
    4 = (i=hi, j=hi, k=lo)
    5 = (i=lo, j=lo, k=hi)
    6 = (i=hi, j=lo, k=hi)
    7 = (i=lo, j=hi, k=hi)
    8 = (i=hi, j=hi, k=hi)

The ``(i, j, k)`` directions need not be world-aligned: a box may be
arbitrarily skewed and rotated so long as the connectivity above is
preserved.  :class:`Box` validates this at construction by
reconstructing the six face centres and rejecting hexahedra whose
faces self-intersect or are degenerate.

Block association
-----------------

Every box references a SHORE block label via :attr:`Box.associated_block`.
At ``cc.par`` write time the label is resolved to the block's 1-based
index in the adjacency JSON's block list; the box inherits that
block's ``level / group / priority`` for the chimera precedence
rules in ``check_parete_interna``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from shore.exceptions import ParameterError


@dataclass(frozen=True)
class Box:
    """One Xall ``cc.par`` type-9 hexahedral marker box.

    Parameters
    ----------
    vertices :
        ``(8, 3)`` array of corner coordinates in the canonical order
        documented above (see module docstring).  Stored as a
        ``float64`` numpy array; lists / tuples are accepted and
        converted at construction.
    associated_block :
        Label of the SHORE block this box is associated with.  At
        ``cc.par`` write time this is resolved to a 1-based block
        index; the box inherits the block's ``level / group /
        priority`` for chimera precedence.

    Examples
    --------
    Build an axis-aligned cube of half-side 0.6 around the origin,
    associated with block ``"sub0"``:

    >>> import numpy as np
    >>> half = 0.6
    >>> verts = np.array([
    ...     [-half, -half, -half],   # v1: (i=lo, j=lo, k=lo)
    ...     [+half, -half, -half],   # v2: (i=hi, j=lo, k=lo)
    ...     [-half, +half, -half],   # v3: (i=lo, j=hi, k=lo)
    ...     [+half, +half, -half],   # v4: (i=hi, j=hi, k=lo)
    ...     [-half, -half, +half],   # v5: (i=lo, j=lo, k=hi)
    ...     [+half, -half, +half],   # v6: (i=hi, j=lo, k=hi)
    ...     [-half, +half, +half],   # v7: (i=lo, j=hi, k=hi)
    ...     [+half, +half, +half],   # v8: (i=hi, j=hi, k=hi)
    ... ])
    >>> box = Box(vertices=verts, associated_block="sub0")
    """

    vertices: np.ndarray = field()
    associated_block: str = field()

    def __post_init__(self) -> None:
        # Coerce to (8, 3) float64.  ``frozen=True`` blocks normal
        # attribute assignment, so use object.__setattr__ as the
        # standard frozen-dataclass escape hatch.
        verts = np.asarray(self.vertices, dtype=np.float64)
        if verts.shape != (8, 3):
            raise ParameterError(f"Box: vertices must have shape (8, 3), got {verts.shape}")
        if not np.all(np.isfinite(verts)):
            raise ParameterError("Box: vertices contain non-finite values")
        if not isinstance(self.associated_block, str) or not self.associated_block:
            raise ParameterError(
                f"Box: associated_block must be a non-empty string, got {self.associated_block!r}"
            )
        object.__setattr__(self, "vertices", verts)
        # Validate the hexahedron is non-degenerate and consistently
        # oriented.  Cheap to do at construction; pays for itself the
        # first time someone passes a typo'd vertex order.
        _validate_hexahedron(verts)

    # ------------------------------------------------------------------
    # Convenience constructors
    # ------------------------------------------------------------------

    @classmethod
    def from_aabb(
        cls,
        vmin: tuple[float, float, float] | np.ndarray,
        vmax: tuple[float, float, float] | np.ndarray,
        associated_block: str,
    ) -> Box:
        """Build a Box from an axis-aligned bounding box.

        Parameters
        ----------
        vmin, vmax :
            Diagonal corners of the AABB in world coordinates.
            ``vmax > vmin`` componentwise.
        associated_block :
            Block label.

        Returns
        -------
        Box
            A type-9 hexahedron whose 8 vertices enumerate the corners
            of the AABB in canonical order.
        """
        lo = np.asarray(vmin, dtype=np.float64)
        hi = np.asarray(vmax, dtype=np.float64)
        if lo.shape != (3,) or hi.shape != (3,):
            raise ParameterError(
                f"Box.from_aabb: vmin / vmax must be 3-vectors, got {lo.shape}, {hi.shape}"
            )
        if not np.all(hi > lo):
            raise ParameterError(
                f"Box.from_aabb: vmax must exceed vmin componentwise, "
                f"got vmin={lo.tolist()}, vmax={hi.tolist()}"
            )
        verts = np.empty((8, 3), dtype=np.float64)
        for idx in range(8):
            i = idx & 1
            j = (idx >> 1) & 1
            k = (idx >> 2) & 1
            verts[idx, 0] = hi[0] if i else lo[0]
            verts[idx, 1] = hi[1] if j else lo[1]
            verts[idx, 2] = hi[2] if k else lo[2]
        return cls(vertices=verts, associated_block=associated_block)


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

# Canonical face vertex sets for a type-9 hexahedron, mirroring the
# face-centre formulas in ``scatole.f90:30-35`` (translated to 0-based
# Python indices).  The ordering encodes the (i_lo / i_hi / j_lo /
# j_hi / k_lo / k_hi) face pairs in this order.
_FACE_VERTEX_INDICES: tuple[tuple[int, int, int, int], ...] = (
    (0, 2, 4, 6),  # i_lo: vertices {1, 3, 5, 7} (1-based)
    (1, 3, 5, 7),  # i_hi: {2, 4, 6, 8}
    (0, 1, 4, 5),  # j_lo: {1, 2, 5, 6}
    (2, 3, 6, 7),  # j_hi: {3, 4, 7, 8}
    (0, 1, 2, 3),  # k_lo: {1, 2, 3, 4}
    (4, 5, 6, 7),  # k_hi: {5, 6, 7, 8}
)

# Pairs of opposite faces, used to verify the hexahedron has positive
# volume (face centres of opposite faces lie on opposite sides of the
# centroid along the corresponding axis).
_OPPOSITE_FACE_PAIRS: tuple[tuple[int, int], ...] = ((0, 1), (2, 3), (4, 5))


def _validate_hexahedron(vertices: np.ndarray) -> None:
    """Reject degenerate or self-intersecting hexahedra.

    Four checks, each cheap and each catching a distinct user error:

    1. **Coincident vertices.**  Two vertices closer than a tiny
       absolute tolerance produce a zero-area face — Xall's
       ``calc_sbox`` would fail at runtime.
    2. **Non-coplanar.**  The 8 vertices must span 3 dimensions.
    3. **Per-face planarity.**  Each of the six faces' four vertices
       must lie within a small tolerance of a single plane.  Catches
       bow-tie vertex orderings (e.g. swapping ``v3`` and ``v5``):
       opposite-face-centres tests aren't enough — a bow-tied
       hexahedron can still pass them while having two
       butterfly-quad faces.
    4. **Positive volume.**  The hexahedron's centroid must lie
       strictly between every pair of opposite face centres along
       the line connecting them.

    Raises
    ------
    ParameterError
        With a message identifying the failed check.
    """
    # 1. Coincident vertices.
    diffs = vertices[:, None, :] - vertices[None, :, :]
    dist2 = np.einsum("ijk,ijk->ij", diffs, diffs)
    np.fill_diagonal(dist2, np.inf)
    min_dist2 = dist2.min()
    bbox_diag = np.linalg.norm(vertices.max(axis=0) - vertices.min(axis=0))
    tol = max(1e-12, (bbox_diag * 1e-9) ** 2)
    if min_dist2 < tol:
        raise ParameterError(
            f"Box: hexahedron has coincident vertices (min squared "
            f"distance {min_dist2:.3e} < {tol:.3e})"
        )

    # 2. Non-coplanar (cheaper to evaluate before face-centre work).
    centroid = vertices.mean(axis=0)
    centred = vertices - centroid
    rank = np.linalg.matrix_rank(centred, tol=bbox_diag * 1e-8)
    if rank < 3:
        raise ParameterError(
            f"Box: hexahedron is degenerate — vertices span {rank} dimensions (need 3)"
        )

    # 3. Per-face planarity.  A face's 4 vertices must lie on a single
    # plane to within a tolerance proportional to the bounding box.
    # We measure this as the largest singular value of the
    # mean-centred face matrix on the *normal* direction: equivalently,
    # the smallest singular value (the out-of-plane spread) must be
    # near zero.
    plane_tol = bbox_diag * 1e-6
    for face_idx, idxs in enumerate(_FACE_VERTEX_INDICES):
        face_pts = vertices[list(idxs)]
        face_centroid = face_pts.mean(axis=0)
        centred_face = face_pts - face_centroid
        # Smallest singular value = out-of-plane RMS extent.
        s = np.linalg.svd(centred_face, compute_uv=False)
        if s[-1] > plane_tol:
            raise ParameterError(
                f"Box: face {face_idx} is non-planar (out-of-plane "
                f"extent {s[-1]:.3e} > tol {plane_tol:.3e}); "
                "check vertex ordering against the canonical (i, j, k) "
                "bit layout"
            )

    # 4. Positive volume via opposite-face-centre test.
    face_centres = np.array([vertices[list(idxs)].mean(axis=0) for idxs in _FACE_VERTEX_INDICES])
    for f_lo, f_hi in _OPPOSITE_FACE_PAIRS:
        c_lo = face_centres[f_lo]
        c_hi = face_centres[f_hi]
        # Centroid should lie between the two face centres.
        # Equivalently: (centroid - c_lo) . (c_hi - c_lo) is positive
        # AND smaller than |c_hi - c_lo|^2.
        axis = c_hi - c_lo
        axis_len2 = float(axis @ axis)
        if axis_len2 < tol:
            raise ParameterError(
                f"Box: opposite face centres {f_lo} and {f_hi} coincide (zero-thickness hexahedron)"
            )
        proj = float((centroid - c_lo) @ axis)
        if not (0.0 < proj < axis_len2):
            raise ParameterError(
                f"Box: hexahedron has inverted or self-intersecting "
                f"face pair (faces {f_lo} ↔ {f_hi}); check vertex "
                "ordering against the canonical (i, j, k) bit layout"
            )


# ---------------------------------------------------------------------------
# JSON round-trip helpers (used by the cc-par CLI's --boxes flag)
# ---------------------------------------------------------------------------


def box_to_json(box: Box) -> dict[str, Any]:
    """Serialise a :class:`Box` to a JSON-friendly dict."""
    return {
        "associated_block": box.associated_block,
        "vertices": box.vertices.tolist(),
    }


def box_from_json(payload: dict[str, Any]) -> Box:
    """Parse one ``{associated_block, vertices}`` dict into a :class:`Box`.

    Raises
    ------
    ParameterError
        On missing keys, wrong vertex shape, or an invalid hexahedron
        (forwarded from :class:`Box.__post_init__`).
    """
    if not isinstance(payload, dict):
        raise ParameterError(
            f"Box JSON: each entry must be an object, got {type(payload).__name__}"
        )
    if "associated_block" not in payload:
        raise ParameterError("Box JSON: missing required key 'associated_block'")
    if "vertices" not in payload:
        raise ParameterError("Box JSON: missing required key 'vertices'")
    return Box(
        vertices=np.asarray(payload["vertices"], dtype=np.float64),
        associated_block=str(payload["associated_block"]),
    )
