"""Body: a geometry (STL or CAD) with structured sphere projection."""

from __future__ import annotations

import functools
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import trimesh

from shore.exceptions import MeshInputError
from shore.geometry import CAD_SUFFIXES, AnchoredCaster, Geometry, load_geometry
from shore.geometry.stl import StlGeometry
from shore.surface.anchor import resolve_anchor
from shore.surface.io import load_stl
from shore.surface.projector import Projector
from shore.surface.remesh import (
    project_caps_to_surface,
    sphere_latlon_directions,
)
from shore.volume.hyperbolic import _stl_normals


@dataclass
class Body:
    """Body surface with ray-cast projection capability.

    Parameters
    ----------
    mesh :
        Triangulated body surface. For a CAD body (``geometry`` given) it may
        be omitted and defaults to the geometry's merged tessellation, which
        serves the mesh-based consumers (anchor resolution, C-grid slicing).
    anchor :
        Ray-cast origin, shape (3,).  When *None*, defaults to the mesh
        centroid (set by ``__post_init__``).
    geometry :
        The body as a ``Geometry``. When *None*, an ``StlGeometry`` over
        *mesh*. Projections (k=0 ray casts, marching normals) go through it.
    anchor_error :
        Why no anchor could be resolved (``from_stl`` / ``from_file``): the
        body loads, and the error is raised where the anchor is first used,
        by ``projector`` (ADR 0018). Collars never ray-cast, so an open
        half-model body without an interior point still serves them.
    """

    mesh: trimesh.Trimesh | None = None
    anchor: np.ndarray = field(default=None, compare=False, repr=False)  # type: ignore[assignment]
    geometry: Geometry = field(default=None, compare=False, repr=False)  # type: ignore[assignment]
    anchor_error: MeshInputError | None = field(default=None, compare=False, repr=False)

    def __post_init__(self) -> None:
        if self.geometry is None:
            if self.mesh is None:
                raise MeshInputError("Body needs a mesh or a geometry.")
            self.geometry = StlGeometry(self.mesh)
        if self.mesh is None:
            proxy = getattr(self.geometry, "merged_proxy", None)
            if proxy is None:
                raise MeshInputError(
                    f"{type(self.geometry).__name__} has no tessellation; pass mesh= as well."
                )
            self.mesh = proxy
        if self.anchor is None:
            self.anchor = self.mesh.centroid.copy()

    @classmethod
    def from_stl(
        cls,
        path: str | Path,
        anchor_strategy: str = "auto",
        skip_hygiene: bool = False,
    ) -> Body:
        """Load an STL and resolve the anchor point.

        Parameters
        ----------
        path :
            Path to the STL file.
        anchor_strategy :
            Anchor resolution strategy passed to
            ``shore.surface.anchor.resolve_anchor``.
            ``'auto'`` (default) tries incenter then bbox with containment
            validation; ``'centroid'`` replicates the legacy behaviour.
        skip_hygiene :
            Bypass mesh hygiene (merge_vertices / fix_normals / remove
            degenerate faces).  Preserves legacy behaviour for regression tests.
        """
        mesh = load_stl(path, skip_hygiene=skip_hygiene)
        try:
            return cls(mesh=mesh, anchor=resolve_anchor(mesh, strategy=anchor_strategy))
        except MeshInputError as err:
            return cls(mesh=mesh, anchor_error=err)

    @classmethod
    def from_file(
        cls,
        path: str | Path,
        anchor_strategy: str = "auto",
        skip_hygiene: bool = False,
    ) -> Body:
        """Load an STL or a CAD file (STEP/IGES/BREP) and resolve the anchor.

        STL goes through ``from_stl``. CAD needs ``shore-mesh[cad]``; the
        anchor is resolved on the merged tessellation, and *skip_hygiene*
        does not apply.
        """
        if Path(path).suffix.lower() not in CAD_SUFFIXES:
            return cls.from_stl(path, anchor_strategy=anchor_strategy, skip_hygiene=skip_hygiene)
        geometry = load_geometry(path)
        body = cls(geometry=geometry)
        try:
            body.anchor = resolve_anchor(body.mesh, strategy=anchor_strategy)
        except MeshInputError as err:
            body.anchor_error = err
        return body

    # ------------------------------------------------------------------
    # Cached projector (Phase 3 — single BVH build per Body instance)
    # ------------------------------------------------------------------

    @property
    def is_cad(self) -> bool:
        """True when the body has exact (non-STL) geometry."""
        return not isinstance(self.geometry, StlGeometry)

    @property
    def projector(self) -> Projector | AnchoredCaster:
        """Anchored ray caster for the k=0 builders, one per Body.

        STL: ``geometry.projector(anchor)`` (the legacy ``Projector``, sharing
        the geometry's BVH; results bit-identical to before). CAD: an
        ``AnchoredCaster`` that casts on the tessellation and snaps hits to
        the exact surface.

        Raises the anchor resolution error when the body has no anchor.
        """
        if self.anchor_error is not None:
            raise self.anchor_error
        if isinstance(self.geometry, StlGeometry):
            return self.geometry.projector(self.anchor)
        return self._cad_caster

    @functools.cached_property
    def _cad_caster(self) -> AnchoredCaster:
        return AnchoredCaster(self.geometry, self.anchor, snap=True)

    def normals(self, points: np.ndarray) -> np.ndarray:
        """Unit outward surface normals at the closest body points (..., 3).

        The marcher's normal blending (``blend_normals_k``) reads these. A CAD
        body gives its exact geometry normals. An STL body gives the blend of
        trimesh's vertex normals (``_stl_normals``), not the crease-aware
        normals of ``StlGeometry`` (ADR 0014): those jump by the edge angle
        across a sharp edge, and blended into the first marching steps they
        fold the cells straddling it (ADR 0015).
        """
        if self.is_cad:
            return self.geometry.project(points).normals
        pts = np.asarray(points, dtype=np.float64)
        return _stl_normals(self.mesh, pts.reshape(-1, 3)).reshape(pts.shape)

    # ------------------------------------------------------------------
    # Projection
    # ------------------------------------------------------------------

    def project(
        self,
        ni: int,
        nj: int,
        theta_cap: float | None = None,
        theta_cap_deg: float | None = None,
        i_spacing: str = "uniform",
        i_beta: float = 3.0,
        j_spacing: str = "uniform",
        j_beta: float = 3.0,
        j_cluster_lon: float | None = None,
    ) -> np.ndarray:
        """Ray-cast a structured lat-lon grid onto the body surface.

        Parameters
        ----------
        ni, nj :
            Latitude and longitude resolution.
        theta_cap :
            Polar exclusion angle in radians.
        theta_cap_deg :
            Polar exclusion angle in degrees (mutually exclusive with theta_cap).

        Returns
        -------
        np.ndarray, shape (ni, nj, 3)
            Coordinates of the k=0 structured surface layer.
        """
        from shore.surface.remesh import sphere_latlon_directions

        dirs = sphere_latlon_directions(
            ni,
            nj,
            theta_cap=theta_cap,
            theta_cap_deg=theta_cap_deg,
            i_spacing=i_spacing,
            i_beta=i_beta,
            j_spacing=j_spacing,
            j_beta=j_beta,
            j_cluster_lon=j_cluster_lon,
        )
        result = self.projector.project(dirs.reshape(-1, 3))
        return result.locations.reshape(ni, nj, 3)

    def project_caps(
        self,
        ni: int,
        nj: int,
        theta_cap: float | None = None,
        theta_cap_deg: float | None = None,
        i_spacing: str = "uniform",
        i_beta: float = 3.0,
        j_spacing: str = "uniform",
        j_beta: float = 3.0,
        j_cluster_lon: float | None = None,
        pole_eps_deg: float = 0.5,
    ) -> tuple[np.ndarray, int]:
        """Unified ray-cast covering the full sphere for the cubed-sphere topology.

        Returns the full projected surface plus the split index ``nc_cap``.
        Rows ``0 .. nc_cap-1`` are the south cap, rows ``nc_cap .. ni_total-nc_cap``
        are the equatorial band, rows ``ni_total-nc_cap ..`` are the north cap.

        Parameters
        ----------
        ni :
            Number of equatorial rows.
        nj :
            Longitude resolution (multiple of 4).
        theta_cap :
            Polar split angle in radians.  Default ``0.05 * pi`` (~9 deg).
        theta_cap_deg :
            Convenience alias in degrees.
        pole_eps_deg :
            Angular exclusion at each pole in degrees (default 0.5 deg).

        Returns
        -------
        surface : np.ndarray, shape (ni_total, nj, 3)
        nc_cap : int
        """
        return project_caps_to_surface(
            self.mesh,
            ni,
            nj,
            theta_cap=theta_cap,
            theta_cap_deg=theta_cap_deg,
            i_spacing=i_spacing,
            i_beta=i_beta,
            j_spacing=j_spacing,
            j_beta=j_beta,
            j_cluster_lon=j_cluster_lon,
            pole_eps_deg=pole_eps_deg,
            projector=self.projector,
        )

    def directions(
        self,
        ni: int,
        nj: int,
        theta_cap: float | None = None,
        theta_cap_deg: float | None = None,
        i_spacing: str = "uniform",
        i_beta: float = 3.0,
        j_spacing: str = "uniform",
        j_beta: float = 3.0,
        j_cluster_lon: float | None = None,
    ) -> np.ndarray:
        """Return (ni, nj, 3) unit direction vectors without casting rays."""
        return sphere_latlon_directions(
            ni,
            nj,
            theta_cap=theta_cap,
            theta_cap_deg=theta_cap_deg,
            i_spacing=i_spacing,
            i_beta=i_beta,
            j_spacing=j_spacing,
            j_beta=j_beta,
            j_cluster_lon=j_cluster_lon,
        )

    # ------------------------------------------------------------------
    # Surface I/O
    # ------------------------------------------------------------------

    @staticmethod
    def save_surface(surface: np.ndarray, path: str | Path) -> None:
        """Save a k=0 surface layer to a compressed .npz file."""
        path = Path(path)
        if path.suffix != ".npz":
            path = path.with_suffix(".npz")
        np.savez_compressed(path, k0=surface)

    @staticmethod
    def load_surface(path: str | Path) -> np.ndarray:
        """Load a k=0 surface layer from a .npz file."""
        path = Path(path)
        if not path.exists():
            raise MeshInputError(f"Surface file not found: {path}")
        data = np.load(path)
        if "k0" not in data:
            raise MeshInputError(
                f"'{path}' has no 'k0' array. Was it created by Body.save_surface?"
            )
        return data["k0"]
