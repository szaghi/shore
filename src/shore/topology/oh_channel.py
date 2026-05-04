"""Channel-block frustum builder for the OH (cubed-sphere + H-channel) topology.

A *channel block* sits between one wrap block's outer envelope (k=nk_w-1) and
one face of the user-supplied bounding box.  The two end faces are:

- ``k=0``     ``inner_face`` — bit-copy of the wrap's ``k_hi`` slice.
- ``k=nk-1``  bilinear (a, b) sweep of ``outer_quad`` (the box face's 4 corners).

Interior k-layers are produced by ruled-surface interpolation between the
two end faces under a 1D spacing law (``uniform`` / ``tanh`` / ``tanh2``).
The outer-face parameterisation is index-for-index aligned with the inner
face so that the channel-channel lateral seams (j_hi ↔ j_lo on adjacent
sub-blocks, and i-seams to the cap channels) inherit the wrap's seam
topology with identity AxisMaps.

C1 continuity at the wrap-channel seam is optional via ``c1_first_step``:
when set, the first interior layer is placed at distance ``c1_first_step``
from k=0 along the local k-direction, matching the wrap's last marching
step exactly so that the spacing has no jump across the seam.
"""

from __future__ import annotations

import numpy as np

from shore.exceptions import ParameterError
from shore.volume.spacing import build_i_coords


def build_channel_block(
    inner_face: np.ndarray,
    outer_quad: np.ndarray,
    nk_channel: int,
    k_law: str = "uniform",
    k_beta: float = 3.0,
    c1_first_step: float | np.ndarray | None = None,
) -> np.ndarray:
    """Build a channel frustum block, shape ``(ni, nj, nk_channel, 3)``.

    Parameters
    ----------
    inner_face :
        Wrap block's k_hi outer envelope, shape ``(ni, nj, 3)``.  This
        becomes the channel's k=0 row, bit-exactly.
    outer_quad :
        Four 3D corners of the box face this channel block "owns",
        shape ``(4, 3)``.  Order: ``(a=0, b=0)``, ``(a=1, b=0)``,
        ``(a=1, b=1)``, ``(a=0, b=1)``.  The same parametric (a, b)
        used by ``surface/cap.py``, so that diagonal sub-block
        junctions in the wrap map onto box-face corners.
    nk_channel :
        Total number of layers in the channel block (k=0 + interior + k=nk-1).
        Must be ``>= 2``.
    k_law :
        Spacing law along k: ``'uniform'``, ``'tanh'`` (cluster near
        inner face / wrap seam), or ``'tanh2'`` (cluster at both ends).
    k_beta :
        Clustering strength for tanh / tanh2.
    c1_first_step :
        When supplied as a scalar, the distance from k=0 to k=1 at every
        (i, j) equals this value exactly.  When supplied as an
        ``(ni, nj)`` array, the distance at each (i, j) equals the
        corresponding entry — useful for matching the wrap's per-(i, j)
        last marching step exactly so the seam is C1 at every node, not
        just on average.  Must be positive everywhere it's supplied.

    Returns
    -------
    numpy.ndarray, shape ``(ni, nj, nk_channel, 3)``
        Frustum block in block-local (i, j, k) layout.

    Raises
    ------
    ParameterError
        On invalid shapes, unknown spacing laws, or non-positive
        ``c1_first_step``.
    """
    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------
    if not (
        isinstance(inner_face, np.ndarray) and inner_face.ndim == 3 and inner_face.shape[-1] == 3
    ):
        raise ParameterError(
            f"build_channel_block: inner_face must have shape (ni, nj, 3); "
            f"got {None if not isinstance(inner_face, np.ndarray) else inner_face.shape}"
        )
    if not (isinstance(outer_quad, np.ndarray) and outer_quad.shape == (4, 3)):
        raise ParameterError(
            f"build_channel_block: outer_quad must have shape (4, 3); "
            f"got {None if not isinstance(outer_quad, np.ndarray) else outer_quad.shape}"
        )
    if nk_channel < 2:
        raise ParameterError(f"build_channel_block: nk_channel must be >= 2, got {nk_channel}")
    if k_law not in ("uniform", "tanh", "tanh2"):
        raise ParameterError(
            f"build_channel_block: k_law must be one of 'uniform' / 'tanh' / 'tanh2', got {k_law!r}"
        )
    if c1_first_step is not None:
        c1_arr = np.asarray(c1_first_step, dtype=np.float64)
        if c1_arr.ndim not in (0, 2):
            raise ParameterError(
                f"build_channel_block: c1_first_step must be a scalar or "
                f"(ni, nj) array, got shape {c1_arr.shape}"
            )
        if c1_arr.ndim == 2 and c1_arr.shape != inner_face.shape[:2]:
            raise ParameterError(
                f"build_channel_block: c1_first_step shape {c1_arr.shape} "
                f"does not match inner_face shape "
                f"{inner_face.shape[:2]}"
            )
        if np.any(c1_arr <= 0):
            raise ParameterError(
                f"build_channel_block: c1_first_step must be positive "
                f"everywhere; got min {float(c1_arr.min())}"
            )

    ni, nj = inner_face.shape[:2]

    # ------------------------------------------------------------------
    # Build the outer face: bilinear sweep of outer_quad over (a, b).
    # ------------------------------------------------------------------
    # (a, b) are uniform in [0, 1]^2 — same convention as cap.py's
    # tangent-plane parameterisation, so the (i, j) → (a, b) mapping
    # matches the wrap's outer envelope and seam topology is preserved.
    a = np.linspace(0.0, 1.0, ni)[:, None, None]  # (ni, 1, 1)
    b = np.linspace(0.0, 1.0, nj)[None, :, None]  # (1, nj, 1)

    c00, c10, c11, c01 = outer_quad[0], outer_quad[1], outer_quad[2], outer_quad[3]
    outer_face = (
        (1.0 - a) * (1.0 - b) * c00 + a * (1.0 - b) * c10 + a * b * c11 + (1.0 - a) * b * c01
    )  # (ni, nj, 3)

    # ------------------------------------------------------------------
    # Build the 1D k-parameter t in [0, 1] under the requested law.
    # ------------------------------------------------------------------
    # build_i_coords handles uniform / tanh / tanh2 on a fixed [start, stop].
    # We use [0, 1] so the parameter is a direct interpolation weight.
    t = build_i_coords(nk_channel, 0.0, 1.0, spacing=k_law, beta=k_beta)  # (nk,)

    # ------------------------------------------------------------------
    # Optional C1 first-step pin: per-(i, j) t-grid so that
    # ||x[i, j, k=1] - x[i, j, k=0]|| == c1_first_step exactly at every (i, j).
    # ------------------------------------------------------------------
    # Without the pin, t is a single 1D array shared across all (i, j) and
    # the channel column at any fixed (i, j) is interpolated as
    #     x = (1 - t) * inner + t * outer.
    # The first cell's length is t[1] * ||outer - inner||, which depends on
    # the local gap.  When the wrap envelope is non-spherical or the box
    # face is non-parallel to the envelope, the gap varies across (i, j)
    # and a single t[1] cannot deliver a uniform first step.
    #
    # The C1 pin solves this by using a per-(i, j) t-grid that
    #   (a) anchors t[1] = c1_first_step / gap[i, j]  — bit-exact seam, and
    #   (b) distributes the remaining nk-2 cells under the chosen k_law,
    #       *rescaled* so the tail intervals preserve the law's shape but
    #       sum to (1 - t[1]).
    # That keeps the tanh / tanh2 spacing variation smooth across the seam
    # — no second-cell jump on top of the pinned first cell.
    # ------------------------------------------------------------------
    inner = inner_face[:, :, None, :]  # (ni, nj, 1, 3)
    outer = outer_face[:, :, None, :]  # (ni, nj, 1, 3)

    if c1_first_step is None:
        # Standard 1D-t path under the user-chosen law.
        t = build_i_coords(nk_channel, 0.0, 1.0, spacing=k_law, beta=k_beta)
        t_arr = t[None, None, :, None]  # (1, 1, nk, 1)
    else:
        # With the C1 pin, the first cell size is fixed at c1_first_step.
        # The remaining nk-2 cells grow geometrically with a per-(i, j)
        # ratio r[i, j] solved from the geometric-series constraint
        #     pin + pin*r + pin*r^2 + ... + pin*r^(nk-2) = gap
        # so that the cells fill the gap exactly while preserving smooth
        # growth from the pin out to the box face.  The k_law / k_beta
        # arguments are ignored on this branch — the pin determines the
        # first cell, and a geometric tail is the natural smooth
        # continuation (no clustering at the seam, no mid-channel jump).
        gap = np.linalg.norm(outer_face - inner_face, axis=-1)  # (ni, nj)
        first_step_2d = (
            np.broadcast_to(c1_arr, gap.shape).astype(np.float64, copy=False)
            if c1_arr.ndim == 0
            else c1_arr
        )
        if np.any(first_step_2d >= gap):
            bad = float(first_step_2d.max()) / float(gap.min())
            raise ParameterError(
                f"build_channel_block: c1_first_step exceeds the local "
                f"inner->outer gap somewhere (worst ratio first_step/gap "
                f"= {bad:.3f}); the channel cannot contain a first step "
                f"that large.  Reduce c1_first_step or enlarge the box."
            )
        # Solve gap = pin * (1 - r^(nk-1)) / (1 - r) for r per (i, j).
        # Equivalently: f(r) = r^(nk-1) - (gap/pin) * r + ((gap/pin) - 1) = 0.
        # We need r >= 1 (or r > 0 if the gap is small enough that uniform/
        # contracting works).  Newton iteration converges quickly from r=1.05.
        n_tail = nk_channel - 1  # number of intervals
        ratio_target = gap / first_step_2d  # gap_in_pin_units; > 1 always
        r = _solve_geometric_ratio(ratio_target, n_tail)  # (ni, nj)
        # Build cumulative parametric positions: t[k] = pin/gap * (1-r^k)/(1-r)
        # for k=0..n_tail.  Special-case r ~ 1 to avoid division by ~0.
        k_idx = np.arange(n_tail + 1, dtype=np.float64)  # 0..nk-1
        # Broadcast r to (ni, nj, 1).
        r3 = r[:, :, None]
        with np.errstate(divide="ignore", invalid="ignore"):
            geom = np.where(
                np.abs(r3 - 1.0) < 1e-12,
                k_idx[None, None, :] / n_tail,
                (1.0 - r3 ** k_idx[None, None, :]) / (1.0 - r3**n_tail),
            )  # (ni, nj, nk) in [0, 1]
        # Scale by gap / |outer - inner| (which is 1 since we are in the
        # parametric t-space already; geom is already on [0, 1]).  geom[..., 0] = 0,
        # geom[..., -1] = 1 by construction.
        t_arr = geom[..., None]  # (ni, nj, nk, 1)

    # ------------------------------------------------------------------
    # Ruled-surface interpolation: x[i, j, k] = (1 - t) * inner + t * outer.
    # ------------------------------------------------------------------
    block = (1.0 - t_arr) * inner + t_arr * outer
    return np.ascontiguousarray(block)


def _solve_geometric_ratio(
    ratio_target: np.ndarray,
    n_tail: int,
    tol: float = 1e-12,
    max_iter: int = 60,
) -> np.ndarray:
    """Solve sum_{k=0}^{n-1} r^k = R for r > 0, vectorised over (i, j).

    Equivalently: ``(1 - r^n) / (1 - r) = R`` (with the limit ``r=1 → n``).
    Newton iteration on ``g(r) = (1 - r^n)/(1 - r) - R``, robust because
    g is monotonically increasing in r for r > 0.

    Parameters
    ----------
    ratio_target :
        ``(ni, nj)`` array of gap / first_step values (always >= 1).
    n_tail :
        Number of channel intervals (``nk_channel - 1``).
    tol, max_iter :
        Newton convergence parameters.
    """
    R = np.asarray(ratio_target, dtype=np.float64)
    # Handle the uniform-spacing limit (R == n_tail → r == 1) explicitly so
    # Newton doesn't divide by zero in g'(1) = n*(n-1)/2 (well-defined, but
    # the closed-form g(r) = (1-r^n)/(1-r) hits 0/0 there).  Start Newton
    # away from r=1 to keep the formula numerically stable.
    r = np.where(R > n_tail, 1.5, np.where(R < n_tail, 0.5, 1.0))
    n = float(n_tail)
    for _ in range(max_iter):
        # Skip cells already at the r=1 limit (uniform).
        active = np.abs(r - 1.0) > 1e-10
        if not np.any(active):
            break
        ra = r[active]
        rn = ra**n
        # g(r) = (1 - r^n) / (1 - r) - R   (well-defined for r != 1)
        g = (1.0 - rn) / (1.0 - ra) - R[active]
        # Quotient rule on u/v with u = 1 - r^n, v = 1 - r:
        #   g'(r) = (u' v - u v') / v^2
        #         = (-n*r^{n-1} * (1 - r) - (1 - r^n) * (-1)) / (1 - r)^2
        #         = (-n*r^{n-1} * (1 - r) + (1 - r^n)) / (1 - r)^2
        gp = (-n * ra ** (n - 1.0) * (1.0 - ra) + (1.0 - rn)) / (1.0 - ra) ** 2
        delta = g / gp
        ra_new = ra - delta
        # Guard against negative r — clamp to a small positive value.
        ra_new = np.where(ra_new <= 0.0, 0.5 * ra, ra_new)
        r[active] = ra_new
        if np.max(np.abs(delta)) < tol:
            break
    return r
