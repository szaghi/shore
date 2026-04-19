"""Body: STL surface wrapper with structured sphere projection."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import trimesh

from shore.exceptions import MeshInputError
from shore.surface.io import load_stl
from shore.surface.remesh import (
    project_caps_to_surface,
    project_sphere_to_surface,
    sphere_latlon_directions,
)


@dataclass
class Body:
    """Triangulated body surface with ray-cast projection capability."""

    mesh: trimesh.Trimesh

    @classmethod
    def from_stl(cls, path: str | Path) -> Body:
        return cls(mesh=load_stl(path))

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
        numpy.ndarray, shape (ni, nj, 3)
            Coordinates of the k=0 structured surface layer.
        """
        return project_sphere_to_surface(
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
        )

    @staticmethod
    def save_surface(surface: np.ndarray, path: str | Path) -> None:
        """Save a k=0 surface layer to a compressed .npz file.

        Parameters
        ----------
        surface :
            Shape (ni, nj, 3) — the projected k=0 layer.
        path :
            Output path.  The ``.npz`` extension is appended if absent.
        """
        path = Path(path)
        if path.suffix != ".npz":
            path = path.with_suffix(".npz")
        np.savez_compressed(path, k0=surface)

    @staticmethod
    def load_surface(path: str | Path) -> np.ndarray:
        """Load a k=0 surface layer from a .npz file.

        Parameters
        ----------
        path :
            Path to a ``.npz`` file produced by :meth:`save_surface` or
            ``shore project``.

        Returns
        -------
        numpy.ndarray, shape (ni, nj, 3)

        Raises
        ------
        MeshInputError
            If the file does not exist or has no ``'k0'`` key.
        """
        path = Path(path)
        if not path.exists():
            raise MeshInputError(f"Surface file not found: {path}")
        data = np.load(path)
        if "k0" not in data:
            raise MeshInputError(
                f"'{path}' has no 'k0' array. Was it created by Body.save_surface?"
            )
        return data["k0"]

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
            Number of equatorial rows (the value passed to CapsTopology).
        nj :
            Longitude resolution (multiple of 4).
        theta_cap :
            Polar split angle in radians.  Default ``0.05 * pi`` (~9 deg).
        theta_cap_deg :
            Convenience alias in degrees.  Mutually exclusive with theta_cap.
        pole_eps_deg :
            Angular exclusion at each pole in degrees (default 0.5 deg).

        Returns
        -------
        surface : numpy.ndarray, shape (ni_total, nj, 3)
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
