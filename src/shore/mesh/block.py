"""Block ABC and HexBlock: structured grid block classes."""

from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np

from shore.mesh.edge import Edge, EdgeName
from shore.mesh.face import AxisMap, Face, FaceBC, FaceName
from shore.mesh.spacing import Spacing1D
from shore.mesh.vertex import Vertex


class Block(ABC):
    """Abstract base for a structured multi-dimensional node block.

    Subclasses define the concrete topology (number of vertices, edges,
    faces) and provide named accessors.  The only shared state is ``nodes``
    (the coordinate array) and ``label``.
    """

    label: str
    nodes: np.ndarray  # (*dims, 3)

    @property
    @abstractmethod
    def n_vertices(self) -> int: ...

    @property
    @abstractmethod
    def n_edges(self) -> int: ...

    @property
    @abstractmethod
    def n_faces(self) -> int: ...

    @abstractmethod
    def faces(self) -> list[Face]: ...

    @abstractmethod
    def edges(self) -> list[Edge]: ...

    @abstractmethod
    def vertices(self) -> list[Vertex]: ...

    @abstractmethod
    def wall_face(self) -> Face:
        """The k=0 body-wall face."""
        ...

    @abstractmethod
    def outer_face(self) -> Face:
        """The k=nk-1 outer boundary face."""
        ...

    def face(self, name: FaceName) -> Face:
        """Look up a face by canonical name."""
        for f in self.faces():
            if f.name == name:
                return f
        raise KeyError(f"No face named '{name}' on block '{self.label}'.")

    def edge(self, name: EdgeName) -> Edge:
        """Look up an edge by canonical name."""
        for e in self.edges():
            if e.name == name:
                return e
        raise KeyError(f"No edge named '{name}' on block '{self.label}'.")


class HexBlock(Block):
    """Structured hexahedral block: (ni, nj, nk, 3) node array.

    Topology
    --------
    - 8 vertices at the 8 corners.
    - 12 edges: 4 along each axis direction.
    - 6 faces: i_lo, i_hi, j_lo, j_hi, k_lo (wall), k_hi (outer).

    The ``nodes`` array is owned externally (e.g. a slice of a topology's
    ring allocation) so that adjacent blocks can share numpy views.

    Parameters
    ----------
    label :
        Human-readable block identifier (e.g. ``'sub0'``, ``'cap_north'``).
    nodes :
        Coordinate array, shape ``(ni, nj, nk, 3)``.  May be a view.
    periodic_i, periodic_j :
        Periodicity flags forwarded to marching and smoothing.
    spacing_i, spacing_j, spacing_k :
        Default spacing laws for the three axis directions.  Individual
        edges may override by carrying a different ``Spacing1D`` instance.
    """

    _n_vertices = 8
    _n_edges = 12
    _n_faces = 6

    def __init__(
        self,
        label: str,
        nodes: np.ndarray,
        periodic_i: bool = False,
        periodic_j: bool = True,
        spacing_i: Spacing1D | None = None,
        spacing_j: Spacing1D | None = None,
        spacing_k: Spacing1D | None = None,
    ) -> None:
        self.label = label
        self.nodes = nodes
        self.periodic_i = periodic_i
        self.periodic_j = periodic_j

        ni, nj, nk, _ = nodes.shape

        _sp_i = spacing_i or Spacing1D(law="uniform")
        _sp_j = spacing_j or Spacing1D(law="uniform")
        _sp_k = spacing_k or Spacing1D(law="geometric", ds=1e-3)

        # Stored for the per-edge spacing resolver: the constructor assigns
        # the same Spacing1D instance to all 12 edges of each axis, so an
        # ``edge.spacing is block.<default>`` identity check tells the
        # resolver whether a given edge has been overridden post-construction.
        self._default_spacing_i = _sp_i
        self._default_spacing_j = _sp_j
        self._default_spacing_k = _sp_k

        # --- build vertex list (8 corners) ---
        # vertex index order: (i_lo/hi) x (j_lo/hi) x (k_lo/hi) = 000..111
        self._vertices: list[Vertex] = [
            Vertex(0, (0, 0, 0), body_projected=True),
            Vertex(1, (0, 0, -1), body_projected=False),
            Vertex(2, (0, -1, 0), body_projected=True),
            Vertex(3, (0, -1, -1), body_projected=False),
            Vertex(4, (-1, 0, 0), body_projected=True),
            Vertex(5, (-1, 0, -1), body_projected=False),
            Vertex(6, (-1, -1, 0), body_projected=True),
            Vertex(7, (-1, -1, -1), body_projected=False),
        ]

        # --- build edge list (12 edges) ---
        # Along-i edges (axis='i', j and k fixed)
        self._edges: list[Edge] = [
            Edge(0, "j_lo_k_lo", "i", (0, 4), ni, _sp_i, body_projected=True),
            Edge(1, "j_hi_k_lo", "i", (2, 6), ni, _sp_i, body_projected=True),
            Edge(2, "j_lo_k_hi", "i", (1, 5), ni, _sp_i, body_projected=False),
            Edge(3, "j_hi_k_hi", "i", (3, 7), ni, _sp_i, body_projected=False),
            # Along-j edges (axis='j', i and k fixed)
            Edge(4, "i_lo_k_lo", "j", (0, 2), nj, _sp_j, body_projected=True),
            Edge(5, "i_hi_k_lo", "j", (4, 6), nj, _sp_j, body_projected=True),
            Edge(6, "i_lo_k_hi", "j", (1, 3), nj, _sp_j, body_projected=False),
            Edge(7, "i_hi_k_hi", "j", (5, 7), nj, _sp_j, body_projected=False),
            # Along-k edges (axis='k', i and j fixed)
            Edge(8, "i_lo_j_lo", "k", (0, 1), nk, _sp_k, body_projected=False),
            Edge(9, "i_hi_j_lo", "k", (4, 5), nk, _sp_k, body_projected=False),
            Edge(10, "i_lo_j_hi", "k", (2, 3), nk, _sp_k, body_projected=False),
            Edge(11, "i_hi_j_hi", "k", (6, 7), nk, _sp_k, body_projected=False),
        ]

        # --- build face list (6 faces) ---
        # BCs are set to FREE by default; topology wires the real BCs after
        # construction via _set_face_bc().
        self._faces: list[Face] = [
            Face(0, "i_lo", FaceBC.FREE, body_projected=False),
            Face(1, "i_hi", FaceBC.FREE, body_projected=False),
            Face(2, "j_lo", FaceBC.FREE, body_projected=False),
            Face(3, "j_hi", FaceBC.FREE, body_projected=False),
            Face(4, "k_lo", FaceBC.WALL, body_projected=True),
            Face(5, "k_hi", FaceBC.FREE, body_projected=False),
        ]

    # ------------------------------------------------------------------
    # Abstract interface implementation
    # ------------------------------------------------------------------

    @property
    def n_vertices(self) -> int:
        return self._n_vertices

    @property
    def n_edges(self) -> int:
        return self._n_edges

    @property
    def n_faces(self) -> int:
        return self._n_faces

    def faces(self) -> list[Face]:
        return self._faces

    def edges(self) -> list[Edge]:
        return self._edges

    def vertices(self) -> list[Vertex]:
        return self._vertices

    def wall_face(self) -> Face:
        return self._faces[4]  # k_lo

    def outer_face(self) -> Face:
        return self._faces[5]  # k_hi

    # ------------------------------------------------------------------
    # Mutation helpers (called by Topology during wiring)
    # ------------------------------------------------------------------

    def _set_face_bc(
        self,
        name: FaceName,
        bc: FaceBC,
        partner: Face | None = None,
        axis_map: AxisMap | None = None,
        bc_code: int | None = None,
    ) -> None:
        """Set this block's face BC and (for connection BCs) partner + axis_map.

        Parameters
        ----------
        name :
            Face to mutate.
        bc :
            New :class:`FaceBC`.
        partner :
            Partner face when ``bc`` is SHARED / DIRICHLET / PERIODIC;
            must be ``None`` for WALL / FREE.
        axis_map :
            Required when ``partner`` is supplied: describes how this
            face's along-axes map onto the partner face's along-axes
            (see :class:`AxisMap`).  Validated against
            ``partner.name`` here so any topology builder that wires a
            seam without an axis_map fails loudly at construction
            rather than silently producing a wrong cc.par downstream.
        bc_code :
            Optional explicit natural-BC code for a ``FREE`` face (see
            :class:`shore.mesh.face.Face`).  Only valid when ``bc`` is
            FREE; rejected otherwise so a mis-wired seam can't silently
            carry a stray code.
        """
        connection = bc in (FaceBC.SHARED, FaceBC.DIRICHLET, FaceBC.PERIODIC)
        if connection:
            if partner is None:
                raise ValueError(
                    f"Block {self.label!r}.{name}: bc={bc.name} requires a partner Face"
                )
            if axis_map is None:
                raise ValueError(
                    f"Block {self.label!r}.{name}: bc={bc.name} requires an axis_map "
                    f"(partner={partner.name!r}); pass AxisMap.identity({name!r}) "
                    f"for a same-axis no-flip seam"
                )
            axis_map.validate_against_partner(name, partner.name)
        else:
            if partner is not None:
                raise ValueError(
                    f"Block {self.label!r}.{name}: bc={bc.name} must not have a partner"
                )
            if axis_map is not None:
                raise ValueError(
                    f"Block {self.label!r}.{name}: bc={bc.name} must not have an axis_map"
                )
        if bc_code is not None and bc is not FaceBC.FREE:
            raise ValueError(
                f"Block {self.label!r}.{name}: bc_code is only valid on FREE faces, "
                f"got bc={bc.name}"
            )
        f = self.face(name)
        f.bc = bc
        f.partner = partner
        f.axis_map = axis_map
        f.bc_code = bc_code

    @property
    def periodic(self) -> tuple[bool, bool]:
        return (self.periodic_i, self.periodic_j)
