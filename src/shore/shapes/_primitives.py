"""Analytic primitive bodies (sphere, cylinder, torus, box, capsule).

Thin wrappers around ``trimesh.creation`` that return
``trimesh.Trimesh`` instances ready to be inspected, exported as
STL, or fed into SHORE's meshing pipeline.
"""

from __future__ import annotations

import trimesh

from shore.exceptions import ParameterError


def sphere(radius: float = 1.0, subdivisions: int = 3) -> trimesh.Trimesh:
    """Icosahedral sphere centred at the origin.

    Parameters
    ----------
    radius :
        Sphere radius.  Must be positive.
    subdivisions :
        Icosphere subdivision level.  ``0`` = 20 triangles (raw
        icosahedron); ``3`` (default) = 1280 triangles; each level
        quadruples the count.  ``4`` is dense enough for most CFD
        wraps.

    Returns
    -------
    trimesh.Trimesh
        Watertight icosahedral sphere.
    """
    if radius <= 0:
        raise ParameterError(f"sphere: radius must be positive, got {radius}")
    if subdivisions < 0:
        raise ParameterError(f"sphere: subdivisions must be >= 0, got {subdivisions}")
    return trimesh.creation.icosphere(subdivisions=subdivisions, radius=radius)


def cylinder(
    radius: float = 1.0,
    height: float = 1.0,
    sections: int = 32,
) -> trimesh.Trimesh:
    """Right circular cylinder centred at the origin, axis along z.

    Parameters
    ----------
    radius :
        Cylinder radius.  Must be positive.
    height :
        Total length along the z-axis (the cylinder spans
        ``[-height/2, +height/2]``).
    sections :
        Number of segments around the circumference.  Higher =
        smoother.

    Returns
    -------
    trimesh.Trimesh
        Watertight cylinder with circular caps.
    """
    if radius <= 0:
        raise ParameterError(f"cylinder: radius must be positive, got {radius}")
    if height <= 0:
        raise ParameterError(f"cylinder: height must be positive, got {height}")
    if sections < 3:
        raise ParameterError(f"cylinder: sections must be >= 3, got {sections}")
    return trimesh.creation.cylinder(radius=radius, height=height, sections=sections)


def torus(
    r_major: float = 1.0,
    r_minor: float = 0.25,
    major_sections: int = 32,
    minor_sections: int = 32,
) -> trimesh.Trimesh:
    """Torus centred at the origin, axis of revolution = z.

    Parameters
    ----------
    r_major :
        Major radius (centre of the tube to centre of the torus).
    r_minor :
        Minor radius (tube radius).  Must be smaller than ``r_major``
        for a self-non-intersecting torus.
    major_sections :
        Segments around the major circle.
    minor_sections :
        Segments around the tube (minor circle).

    Returns
    -------
    trimesh.Trimesh
        Watertight ring torus (genus-1 surface).
    """
    if r_major <= 0:
        raise ParameterError(f"torus: r_major must be positive, got {r_major}")
    if r_minor <= 0:
        raise ParameterError(f"torus: r_minor must be positive, got {r_minor}")
    if r_minor >= r_major:
        raise ParameterError(
            f"torus: r_minor ({r_minor}) must be < r_major ({r_major}) "
            "for a self-non-intersecting (ring) torus."
        )
    if major_sections < 3 or minor_sections < 3:
        raise ParameterError(
            f"torus: major_sections / minor_sections must be >= 3, "
            f"got ({major_sections}, {minor_sections})"
        )
    return trimesh.creation.torus(
        major_radius=r_major,
        minor_radius=r_minor,
        major_sections=major_sections,
        minor_sections=minor_sections,
    )


def box(extents: tuple[float, float, float] = (1.0, 1.0, 1.0)) -> trimesh.Trimesh:
    """Axis-aligned rectangular box centred at the origin.

    Parameters
    ----------
    extents :
        ``(dx, dy, dz)`` total side lengths.  All must be positive.

    Returns
    -------
    trimesh.Trimesh
        12-triangle watertight box.
    """
    if len(extents) != 3 or any(e <= 0 for e in extents):
        raise ParameterError(f"box: extents must be 3 positive floats, got {extents}")
    return trimesh.creation.box(extents=list(extents))


def capsule(
    radius: float = 0.5,
    height: float = 1.0,
    sections: int = 16,
) -> trimesh.Trimesh:
    """Capsule (cylinder with hemispherical caps), axis along z.

    Parameters
    ----------
    radius :
        Hemispherical cap radius and cylinder radius.
    height :
        Length of the cylindrical middle section (caps add
        ``2 * radius`` to total length).
    sections :
        Segments around the circumference; the same value is used
        for both circumferential and meridional resolution on the
        caps.

    Returns
    -------
    trimesh.Trimesh
        Watertight capsule.
    """
    if radius <= 0:
        raise ParameterError(f"capsule: radius must be positive, got {radius}")
    if height <= 0:
        raise ParameterError(f"capsule: height must be positive, got {height}")
    if sections < 3:
        raise ParameterError(f"capsule: sections must be >= 3, got {sections}")
    return trimesh.creation.capsule(height=height, radius=radius, count=[sections, sections])
