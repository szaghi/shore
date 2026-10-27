"""Section extrusion: 2D Selig profile -> 3D wing STL.

Given a 2D cross-section polyline (e.g. from
:func:`shore.profiles.naca`), build a wing STL by stacking copies of
the section along a span axis.  Per-station scale (taper), in-plane
shift (sweep), and in-plane rotation (twist) are supported, so a
single call produces a tapered swept twisted wing.

Conventions
-----------
- The 2D section lives in the plane of the two non-span axes
  (``span_axis`` defaults to ``"z"`` -> section in xy).  The
  section's first coordinate is along the chord axis (``+x`` by
  Selig convention); the second coordinate is along the thickness
  axis.
- Span stations are equally spaced along ``span_axis`` between
  ``z = 0`` and ``z = span``.
- Taper: each station's section is scaled by a per-station factor
  computed by linear interpolation between ``1.0`` (root) and
  ``taper`` (tip).
- Sweep: the leading edge of each station is shifted along the
  chord axis by ``z * tan(sweep_deg)``.
- Twist: each station's section is rotated about its leading edge
  by a per-station angle interpolated linearly from ``0.0`` (root)
  to ``twist_deg`` (tip).
- End caps: by default the wing has open ends — useful for C-grid
  topology where the slicer expects no end-cap triangles at the
  span extremes.  Set ``end_caps=True`` to add centroid-fan caps
  for a watertight body.
"""

from __future__ import annotations

import numpy as np
import trimesh

from shore.exceptions import ParameterError

_AXIS_INDEX = {"x": 0, "y": 1, "z": 2}


def extrude_section(
    section_2d: np.ndarray,
    *,
    span: float,
    n_stations: int,
    taper: float = 1.0,
    sweep_deg: float = 0.0,
    twist_deg: float = 0.0,
    span_axis: str = "z",
    end_caps: bool = False,
) -> trimesh.Trimesh:
    """Extrude a 2D section along a span axis into a 3D wing.

    Parameters
    ----------
    section_2d :
        ``(N, 2)`` Selig-ordered closed-or-open section polyline.
        For closed sections (sharp TE) ``section_2d[0]`` and
        ``section_2d[-1]`` should coincide; for open / blunt TE
        they differ.  No reordering is applied.
    span :
        Total wing span (length along ``span_axis``).  Must be
        positive.
    n_stations :
        Number of span stations including root and tip.  ``>= 2``.
    taper :
        Tip-to-root chord ratio.  ``1.0`` (default) = constant
        chord; ``0.5`` = tip chord half the root chord.
    sweep_deg :
        Sweep angle of the leading edge, in degrees, about the
        thickness axis.  Positive = sweepback.  ``0.0`` (default).
    twist_deg :
        Geometric twist (washout / washin) at the tip, in degrees,
        about the LE-thickness axis.  Linearly interpolated from
        ``0.0`` at the root.  Positive = nose-up.
    span_axis :
        ``"x"``, ``"y"``, or ``"z"``.  Default ``"z"``.
    end_caps :
        If ``True``, add a centroid-fan cap at each span end so the
        wing is watertight.  Default ``False`` — the SHORE C-grid
        topology relies on open ends to slice cleanly.

    Returns
    -------
    trimesh.Trimesh
        Wing STL.

    Raises
    ------
    ParameterError
        On invalid arguments.
    """
    sec = np.asarray(section_2d, dtype=np.float64)
    if sec.ndim != 2 or sec.shape[1] != 2:
        raise ParameterError(f"extrude_section: section_2d must have shape (N, 2), got {sec.shape}")
    if sec.shape[0] < 4:
        raise ParameterError(
            f"extrude_section: section_2d must have at least 4 points, got {sec.shape[0]}"
        )
    # Drop a duplicate closing vertex if present (Selig sharp-TE
    # convention has section[0] == section[-1]).  Without this, the
    # i -> (i+1) % N wrap below produces a zero-length edge that
    # degenerates two triangles per span station.
    if np.allclose(sec[0], sec[-1], atol=1e-12):
        sec = sec[:-1]
    if span <= 0:
        raise ParameterError(f"extrude_section: span must be positive, got {span}")
    if n_stations < 2:
        raise ParameterError(f"extrude_section: n_stations must be >= 2, got {n_stations}")
    if taper <= 0:
        raise ParameterError(f"extrude_section: taper must be positive, got {taper}")
    if span_axis not in _AXIS_INDEX:
        raise ParameterError(
            f"extrude_section: span_axis must be 'x', 'y', or 'z', got {span_axis!r}"
        )

    n_sec = sec.shape[0]
    span_idx = _AXIS_INDEX[span_axis]
    non_span = [i for i in range(3) if i != span_idx]
    chord_idx, thick_idx = non_span  # chord is the first non-span axis

    # ----- per-station section transforms ------------------------------------
    span_coords = np.linspace(0.0, span, n_stations)
    sweep_per_unit = float(np.tan(np.radians(sweep_deg)))
    le_x = sec[:, 0].min()  # leading-edge x in the local 2D frame

    stations = np.empty((n_stations, n_sec, 3), dtype=np.float64)
    for k, z in enumerate(span_coords):
        frac = z / span if span > 0 else 0.0
        scale = 1.0 + frac * (taper - 1.0)
        twist = np.radians(frac * twist_deg)
        cos_t, sin_t = np.cos(twist), np.sin(twist)

        # scale about the LE
        sec_k = (sec - np.array([le_x, 0.0])) * scale + np.array([le_x, 0.0])
        # rotate about (le_x, 0) by `twist`
        x_rel = sec_k[:, 0] - le_x
        y_rel = sec_k[:, 1]
        sec_k = np.column_stack(
            [le_x + cos_t * x_rel - sin_t * y_rel, sin_t * x_rel + cos_t * y_rel]
        )
        # sweep: shift the chord coord by tan(sweep) * z
        sec_k[:, 0] += sweep_per_unit * z

        # lift to 3D
        v3 = np.zeros((n_sec, 3), dtype=np.float64)
        v3[:, chord_idx] = sec_k[:, 0]
        v3[:, thick_idx] = sec_k[:, 1]
        v3[:, span_idx] = z
        stations[k] = v3

    vertices = stations.reshape(-1, 3)

    # ----- side triangles ----------------------------------------------------
    faces: list[list[int]] = []
    for k in range(n_stations - 1):
        for i in range(n_sec):
            i_next = (i + 1) % n_sec
            v00 = k * n_sec + i
            v01 = k * n_sec + i_next
            v10 = (k + 1) * n_sec + i
            v11 = (k + 1) * n_sec + i_next
            faces.append([v00, v11, v10])
            faces.append([v00, v01, v11])

    # ----- end caps (optional) -----------------------------------------------
    if end_caps:
        for cap_k, outward_negative in ((0, True), (n_stations - 1, False)):
            ring = stations[cap_k]
            centroid = ring.mean(axis=0)
            vertices = np.vstack([vertices, centroid])
            cid = vertices.shape[0] - 1
            for i in range(n_sec):
                i_next = (i + 1) % n_sec
                v_a = cap_k * n_sec + i
                v_b = cap_k * n_sec + i_next
                if outward_negative:
                    faces.append([cid, v_b, v_a])
                else:
                    faces.append([cid, v_a, v_b])

    # outward normals: the faces above are outward for a counter-clockwise section in a
    # right-handed (chord, thickness, span) frame; flip them otherwise
    area = 0.5 * float(
        np.sum(sec[:, 0] * np.roll(sec[:, 1], -1) - np.roll(sec[:, 0], -1) * sec[:, 1])
    )
    frame = np.zeros((3, 3))
    frame[0, chord_idx] = frame[1, thick_idx] = frame[2, span_idx] = 1.0
    face_array = np.asarray(faces, dtype=np.int64)
    if area * np.linalg.det(frame) < 0.0:
        face_array = face_array[:, ::-1]
    mesh = trimesh.Trimesh(vertices=vertices, faces=face_array)
    mesh.merge_vertices()
    return mesh
