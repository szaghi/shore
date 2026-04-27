"""Equator band construction from cap seam boundaries.

Given the two cap k=0 layers (north and south), the equator's k=0 layer is
built so its first/last rows match the cap seams bit-exactly, and its
interior follows great-circle meridians from each south seam node to the
corresponding north seam node.

The cap is built first (gnomonic flat-square seed); the equator inherits the
cap-edge node distribution along its seam (this is how ``j_spacing`` flows
from cap to equator).  The equator's ``i_spacing`` controls only the
distribution along each meridian, so the cap and equator are decoupled in
the i-direction and coupled in the j-direction.

Public API
----------
build_equator(cap_south, cap_north, ni, nj, projector, i_spacing, i_beta)
    Return ``(ni, nj, 3)`` equator k=0 layer.
"""

from __future__ import annotations

import numpy as np

from shore.volume.spacing import build_i_coords, tanh_one_sided, tanh_two_sided


def _seam_ring(cap: np.ndarray, pole: str) -> np.ndarray:
    """Concatenate the 4 cap edges in sub-block order to form a (nj, 3) ring.

    The mapping mirrors ``shore.topology.sphere_cap._overwrite_cap_boundaries``:

    South cap (z>0):
        sub0 = cap[:, 0, :]                (phi 0°→90°)
        sub1 = cap[-1, :, :]               (phi 90°→180°)
        sub2 = cap[:, -1, :][::-1]         (phi 180°→270°)
        sub3 = cap[0, :, :][::-1]          (phi 270°→360°)

    North cap (z<0):
        sub0 = cap[0, :, :]                (phi 0°→90°)
        sub1 = cap[:, -1, :]               (phi 90°→180°)
        sub2 = cap[-1, :, :][::-1]         (phi 180°→270°)
        sub3 = cap[:, 0, :][::-1]          (phi 270°→360°)

    Concatenation drops the duplicated corner at the start of each successive
    sub-row, yielding ``nj = 4 * (nc - 1)`` distinct seam nodes.
    """
    nc = cap.shape[0]
    if pole == "south":
        sub_seams = [
            cap[:, 0, :],
            cap[-1, :, :],
            cap[:, -1, :][::-1],
            cap[0, :, :][::-1],
        ]
    else:  # north
        sub_seams = [
            cap[0, :, :],
            cap[:, -1, :],
            cap[-1, :, :][::-1],
            cap[:, 0, :][::-1],
        ]

    nj = 4 * (nc - 1)
    seam = np.empty((nj, 3), dtype=np.float64)
    for q, s in enumerate(sub_seams):
        seam[q * (nc - 1) : (q + 1) * (nc - 1)] = s[:-1]
    return seam


def _constrained_distribution(
    ni: int,
    start: float,
    stop: float,
    first_step: float,
    last_step: float,
    spacing: str,
    beta: float,
) -> np.ndarray:
    """Return ``(ni,)`` coords with pinned first and last intervals.

    ``out[1] - out[0] == first_step`` and ``out[-1] - out[-2] == last_step``;
    the interior ``ni-2`` coords lie on ``[start + first_step, stop - last_step]``
    distributed by the chosen law.

    Constraints:
        - ``ni >= 2``
        - ``first_step > 0``, ``last_step > 0``
        - ``first_step + last_step <= stop - start``  (interior must fit)

    For ``ni == 2`` the law is unused; ``out = [start, stop]`` and the caller
    must ensure ``first_step == stop - start`` (no constraint on last).  For
    ``ni == 3`` interior collapses to a single coord at ``start + first_step``,
    requiring ``first_step + last_step == stop - start``.
    """
    if ni < 2:
        raise ValueError(f"ni must be >= 2, got {ni}")
    out = np.empty(ni, dtype=np.float64)
    out[0] = start
    out[-1] = stop
    if ni == 2:
        return out
    out[1] = start + first_step
    out[-2] = stop - last_step
    if ni == 3:
        return out
    # Interior coords 1..ni-2 lie on [start+first_step, stop-last_step] with
    # `ni-2` points (already pinned at endpoints).  The law shapes the gap.
    inner_start = start + first_step
    inner_stop = stop - last_step
    if spacing == "uniform":
        out[1:-1] = np.linspace(inner_start, inner_stop, ni - 2)
    elif spacing == "tanh":
        out[1:-1] = tanh_one_sided(ni - 2, inner_start, inner_stop, beta)
    elif spacing == "tanh2":
        out[1:-1] = tanh_two_sided(ni - 2, inner_start, inner_stop, beta)
    else:
        raise ValueError(f"unknown spacing law: {spacing!r}")
    return out


def _colatitude(point: np.ndarray, anchor: np.ndarray, pole_axis: np.ndarray) -> np.ndarray:
    """Return colatitude (rad) of ``point`` from ``pole_axis`` viewed from ``anchor``.

    ``point`` may be ``(3,)`` or ``(..., 3)``.  Returns matching shape minus the
    trailing 3.
    """
    d = point - anchor
    norm = np.linalg.norm(d, axis=-1, keepdims=True)
    cos_t = np.einsum("...i,i->...", d / np.maximum(norm, 1e-30), pole_axis)
    return np.arccos(np.clip(cos_t, -1.0, 1.0))


def _slerp_meridian(
    p_south: np.ndarray,
    p_north: np.ndarray,
    anchor: np.ndarray,
    i_coords: np.ndarray,
) -> np.ndarray:
    """Return ``(ni, 3)`` directions from anchor along a great-circle arc.

    The arc passes through ``p_south`` at ``i_coords[0]`` and ``p_north`` at
    ``i_coords[-1]`` (interpreted as colatitude); intermediate directions are
    spherical-linearly interpolated between the two endpoint directions.

    ``i_coords`` is assumed monotone in ``[theta_south_corner,
    theta_north_corner]`` but only its normalised position along the interval
    matters for the slerp parameter.
    """
    d_s = p_south - anchor
    d_n = p_north - anchor
    u_s = d_s / np.linalg.norm(d_s)
    u_n = d_n / np.linalg.norm(d_n)

    cos_omega = float(np.clip(np.dot(u_s, u_n), -1.0, 1.0))
    if cos_omega > 1.0 - 1e-12:
        # Endpoints essentially coincide — degenerate meridian.
        return np.broadcast_to(u_s, (i_coords.size, 3)).copy()

    omega = float(np.arccos(cos_omega))
    sin_omega = float(np.sin(omega))

    # Normalised slerp parameter t ∈ [0, 1].
    t = (i_coords - i_coords[0]) / (i_coords[-1] - i_coords[0])
    a = np.sin((1.0 - t) * omega) / sin_omega
    b = np.sin(t * omega) / sin_omega
    return a[:, None] * u_s[None, :] + b[:, None] * u_n[None, :]


def build_equator(
    cap_south: np.ndarray,
    cap_north: np.ndarray,
    ni: int,
    theta_cap: float,
    projector: object,
    i_spacing: str = "uniform",
    i_beta: float = 3.0,
) -> np.ndarray:
    """Return the ``(ni, nj, 3)`` equator k=0 layer.

    Row ``i=0`` is the south seam (cap_south boundary, traversed in sub-block
    order); row ``i=ni-1`` is the north seam.  Interior rows are built by
    slerp along great-circle meridians from each south seam node to the
    corresponding north seam node.  Each meridian is then ray-cast through
    the cached Projector so its interior nodes lie on the body surface.

    Parameters
    ----------
    cap_south, cap_north :
        Cap k=0 layers, shape ``(nc, nc, 3)`` each.
    ni :
        Number of equator latitude rows (including both seam rows).
    theta_cap :
        Polar exclusion angle in radians used to build the cap blocks.  Drives the
        meridional slerp interval ``[theta_cap, pi - theta_cap]`` for the
        i-spacing law.
    projector :
        ``Projector`` with cached BVH and anchor.
    i_spacing :
        Latitude spacing law (``'uniform'``, ``'tanh'``, ``'tanh2'``).
    i_beta :
        Clustering strength for tanh/tanh2.

    Returns
    -------
    numpy.ndarray, shape ``(ni, nj, 3)`` with ``nj = 4*(nc-1)``
    """
    seam_south = _seam_ring(cap_south, "south")
    seam_north = _seam_ring(cap_north, "north")
    nj = seam_south.shape[0]
    assert seam_north.shape[0] == nj

    anchor = np.asarray(projector.anchor, dtype=np.float64)  # type: ignore[attr-defined]

    # Per-j cap-seam first step: angular distance (viewed from anchor)
    # between seam node j and seam node (j+1) mod nj on each cap.  This is
    # the first cell length along the cap-edge multi-block edges meeting at
    # seam node j.  The equator's meridional first interval (also a multi-
    # block edge starting at the same seam vertex) must equal this so that
    # all 3 edges meeting at the cap-corner vertex share a single first-step
    # value (and at edge-midpoint vertices the 2 edges share theirs).
    def _angular_step(seam_ring: np.ndarray) -> np.ndarray:
        d = seam_ring - anchor
        u = d / np.linalg.norm(d, axis=-1, keepdims=True)
        u_next = np.roll(u, -1, axis=0)
        cos_t = np.einsum("ij,ij->i", u, u_next)
        return np.arccos(np.clip(cos_t, -1.0, 1.0))  # (nj,)

    cap_seam_step_s = _angular_step(seam_south)  # arc angle to next seam node
    cap_seam_step_n = _angular_step(seam_north)

    # Pole-axis colatitudes (used to set per-j slerp interval [start_j, stop_j]).
    pole_south_axis = np.array([0.0, 0.0, 1.0])
    pole_north_axis = np.array([0.0, 0.0, -1.0])
    theta_seam_s = _colatitude(seam_south, anchor, pole_south_axis)
    theta_seam_n = _colatitude(seam_north, anchor, pole_north_axis)

    # Convert the cap-seam arc step (which is in *anchor-angle* units, i.e.
    # the angle subtended at the anchor) to the colatitude-step units used
    # by _constrained_distribution.  The conversion is per-j:
    #   first_step_colat = cap_seam_step * interval_j / omega_j
    # where interval_j is the colatitude span [theta_seam_s, pi - theta_seam_n]
    # and omega_j is the anchor-angle subtended by the full meridian (south
    # to north).  On a sphere with anchor at origin, interval_j ≈ omega_j and
    # the conversion is essentially 1:1; the formula keeps it exact for any
    # body shape and anchor placement.
    u_s_seam = (seam_south - anchor) / np.linalg.norm(seam_south - anchor, axis=-1, keepdims=True)
    u_n_seam = (seam_north - anchor) / np.linalg.norm(seam_north - anchor, axis=-1, keepdims=True)
    omega = np.arccos(np.clip(np.einsum("ij,ij->i", u_s_seam, u_n_seam), -1.0, 1.0))  # (nj,)
    interval = (np.pi - theta_seam_n) - theta_seam_s  # (nj,)
    # Guard against degenerate omega (≈0): in that case fall through to the
    # unpinned-fallback branch below by setting first/last to a value the
    # compatibility check rejects.
    safe_omega = np.where(omega > 1e-9, omega, 1.0)
    scale = interval / safe_omega
    first_step = cap_seam_step_s * scale  # (nj,) — colatitude units
    last_step = cap_seam_step_n * scale  # (nj,)

    # Build (ni, nj, 3) of directions per meridian.  Each meridian uses its
    # own constrained-endpoint i-coord sequence so the seam steps match the
    # cap on both sides — but only when the cap and equator radial step
    # magnitudes are comparable.  When they aren't (e.g. coarse cap with fine
    # equator), pinning would produce a sharp interior cell-size jump; in
    # that regime we fall back to the unpinned i_spacing law.
    #
    # The interval is per-j: the seam-node colatitude at j is not theta_cap
    # (only at corners; edge midpoints dip poleward), and similarly at the
    # north seam.  The natural first/last interval of the unpinned law on
    # [start_j, stop_j] is used as the comparability scale.
    dirs = np.empty((ni, nj, 3), dtype=np.float64)
    pin_warned = False
    for j in range(nj):
        start_j = float(theta_seam_s[j])
        stop_j = float(np.pi - theta_seam_n[j])
        interval_j = stop_j - start_j
        fs = float(first_step[j])
        ls = float(last_step[j])

        nat = build_i_coords(ni, start_j, stop_j, spacing=i_spacing, beta=i_beta)
        nat_first = float(nat[1] - nat[0])
        nat_last = float(nat[-1] - nat[-2])

        # Pin only if cap-step is comparable to the natural law's seam step
        # (within a factor of 2 each side).  Beyond that, pinning forces a
        # discontinuity in radial cell size that destabilises marching.
        within_first = 0.5 <= fs / max(nat_first, 1e-30) <= 2.0
        within_last = 0.5 <= ls / max(nat_last, 1e-30) <= 2.0
        compatible = within_first and within_last
        sane = (fs > 0) and (ls > 0) and (fs + ls < 0.9 * interval_j)

        if compatible and sane:
            i_coords = _constrained_distribution(
                ni,
                start_j,
                stop_j,
                first_step=fs,
                last_step=ls,
                spacing=i_spacing,
                beta=i_beta,
            )
        else:
            i_coords = nat
            if not pin_warned:
                import warnings

                warnings.warn(
                    f"Cap radial step (~{fs:.3g} rad) is incompatible with "
                    f"equator '{i_spacing}' law (natural first step "
                    f"~{nat_first:.3g} rad); falling back to unpinned "
                    "spacing.  C1 match at the cap seam is approximate.  "
                    "Choose ni/nj/theta_cap so cap and equator radial steps "
                    "are within ~2x of each other for exact C1 match.",
                    stacklevel=3,
                )
                pin_warned = True
        dirs[:, j, :] = _slerp_meridian(
            seam_south[j],
            seam_north[j],
            anchor,
            i_coords,
        )

    # The two seam rows are pinned bit-exactly to the cap edges (avoid a
    # round-trip through the projector for those nodes).  Interior rows go
    # through the projector.
    out = np.empty((ni, nj, 3), dtype=np.float64)
    out[0] = seam_south
    out[-1] = seam_north
    if ni > 2:
        interior_dirs = dirs[1:-1].reshape(-1, 3)
        result = projector.project(interior_dirs, on_failure="raise")  # type: ignore[attr-defined]
        out[1:-1] = result.locations.reshape(ni - 2, nj, 3)
    return out
