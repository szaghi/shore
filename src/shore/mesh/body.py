"""Body: STL surface wrapper with structured sphere projection."""

from __future__ import annotations

import functools
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import trimesh

from shore.exceptions import MeshInputError
from shore.surface.anchor import resolve_anchor
from shore.surface.io import load_stl
from shore.surface.remesh import (
    project_caps_to_surface,
    sphere_latlon_directions,
)


@dataclass
class Body:
    """Triangulated body surface with ray-cast projection capability.

    Parameters
    ----------
    mesh :
        Triangulated body surface.
    anchor :
        Ray-cast origin, shape (3,).  When *None*, defaults to the mesh
        centroid (set by ``__post_init__``).
    """

    mesh: trimesh.Trimesh
    anchor: np.ndarray = field(default=None, compare=False, repr=False)  # type: ignore[assignment]

    def __post_init__(self) -> None:
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
        anchor = resolve_anchor(mesh, strategy=anchor_strategy)
        return cls(mesh=mesh, anchor=anchor)

    # ------------------------------------------------------------------
    # Cached projector (Phase 3 — single BVH build per Body instance)
    # ------------------------------------------------------------------

    @functools.cached_property
    def projector(self) -> object:
        """Lazily constructed Projector — BVH built exactly once per Body."""
        from shore.surface.projector import Projector

        return Projector(self.mesh, self.anchor)

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
        result = self.projector.project(dirs.reshape(-1, 3))  # type: ignore[union-attr]
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
        """Unified ray-cast covering the full sphere for the caps topology.

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
