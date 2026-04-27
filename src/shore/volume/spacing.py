"""1D grid spacing laws for surface angle redistribution and wall-normal extrusion.

All functions are pure numpy — no scipy dependency.

Public API
----------
build_i_coords         : latitude parameter redistribution
build_j_coords         : longitude parameter redistribution (periodic axis)
build_k_steps          : wall-normal step sizes
pinned_distribution    : 1D coords with first/last interval pinned to caller-given
                         step sizes (used for C1 matching at multi-block seams)

Low-level primitives (reusable):
uniform, tanh_one_sided, tanh_two_sided, geometric
"""

from __future__ import annotations

import numpy as np

from shore.exceptions import ParameterError

# ---------------------------------------------------------------------------
# Low-level 1D laws
# ---------------------------------------------------------------------------


def uniform(n: int, start: float, stop: float) -> np.ndarray:
    """Uniform distribution: n points in [start, stop].

    Parameters
    ----------
    n :
        Number of points.
    start, stop :
        Interval endpoints (inclusive).

    Returns
    -------
    numpy.ndarray, shape (n,)
    """
    return np.linspace(start, stop, n)


def tanh_one_sided(n: int, start: float, stop: float, beta: float) -> np.ndarray:
    """Roberts one-sided tanh clustering toward *start*.

    Dense near ``start``, coarse near ``stop``.

    Uses the reversed Roberts formula so that the fine spacing is at
    ``tau=0`` (the wall / cap boundary):

        t(tau) = 1 - tanh(beta * (1 - tau)) / tanh(beta)

    Parameters
    ----------
    n :
        Number of points.
    start, stop :
        Interval endpoints.
    beta :
        Clustering strength (> 0).  Larger values → tighter clustering.
        beta=1 is mild; beta=5 is aggressive.

    Returns
    -------
    numpy.ndarray, shape (n,)
    """
    tau = np.linspace(0.0, 1.0, n)
    t = 1.0 - np.tanh(beta * (1.0 - tau)) / np.tanh(beta)
    return start + (stop - start) * t


def tanh_two_sided(
    n: int,
    start: float,
    stop: float,
    beta: float,
    cluster_at: float | None = None,
) -> np.ndarray:
    """Roberts two-sided tanh clustering.

    Parameters
    ----------
    n :
        Number of points.
    start, stop :
        Interval endpoints.
    beta :
        Clustering strength (> 0).
    cluster_at :
        If None: symmetric clustering at both ``start`` and ``stop``.
        If a value in ``[start, stop]``: asymmetric — the dense region is
        shifted toward that value.  Used for j-direction longitude clustering
        (e.g. leading-edge longitude on an airfoil).

    Returns
    -------
    numpy.ndarray, shape (n,)
    """
    tau = np.linspace(0.0, 1.0, n)

    if cluster_at is None:
        # Symmetric: cluster at both ends
        t = 0.5 * (1.0 + np.tanh(beta * (tau - 0.5)) / np.tanh(0.5 * beta))
    else:
        total = stop - start
        frac = (cluster_at - start) / total  # map cluster_at → [0, 1]
        t = 0.5 * (1.0 + np.tanh(beta * (tau - frac)) / np.tanh(0.5 * beta))
        t = np.clip(t, 0.0, 1.0)

    return start + (stop - start) * t


# ---------------------------------------------------------------------------
# Domain-specific builders
# ---------------------------------------------------------------------------


def build_i_coords(
    ni: int,
    theta_start: float,
    theta_stop: float,
    spacing: str = "uniform",
    beta: float = 3.0,
) -> np.ndarray:
    """Return *ni* latitude coordinates in ``[theta_start, theta_stop]``.

    Parameters
    ----------
    ni :
        Number of latitude points.
    theta_start, theta_stop :
        Colatitude range in radians (e.g. ``theta_cap`` to ``pi - theta_cap``).
    spacing :
        ``'uniform'`` — equal angular spacing (default, backward-compatible).
        ``'tanh'``    — one-sided Roberts: cluster near ``theta_start``
                        (the cap-boundary row; highest-quality zone for caps).
        ``'tanh2'``   — two-sided Roberts: cluster at both ``theta_start``
                        and ``theta_stop``.
    beta :
        Clustering strength for tanh/tanh2 (ignored for uniform).

    Returns
    -------
    numpy.ndarray, shape (ni,)
    """
    if spacing == "uniform":
        return uniform(ni, theta_start, theta_stop)
    if spacing == "tanh":
        return tanh_one_sided(ni, theta_start, theta_stop, beta)
    if spacing == "tanh2":
        return tanh_two_sided(ni, theta_start, theta_stop, beta)
    raise ParameterError(f"Unknown i-spacing law '{spacing}'. Choose uniform, tanh, or tanh2.")


def build_j_coords(
    nj: int,
    spacing: str = "uniform",
    beta: float = 3.0,
    cluster_lon: float | None = None,
) -> np.ndarray:
    """Return *nj* longitude coordinates in ``[0, 2*pi)``, periodic.

    Parameters
    ----------
    nj :
        Number of longitude points (periodic: endpoint excluded).
    spacing :
        ``'uniform'`` — equal angular spacing (default, backward-compatible).
        ``'tanh2'``   — two-sided Roberts: cluster toward ``cluster_lon``
                        (defaults to ``pi``, the midpoint, if not supplied).
        Note: ``'tanh'`` (one-sided) is not meaningful for a periodic axis.
    beta :
        Clustering strength for tanh2 (ignored for uniform).
    cluster_lon :
        Target longitude in radians for tanh2 clustering.  Only used when
        ``spacing='tanh2'``.  Default ``pi``.

    Returns
    -------
    numpy.ndarray, shape (nj,)
        Values in ``[0, 2*pi)`` with endpoint excluded.

    Notes
    -----
    For tanh2, the algorithm generates nj+1 symmetric tanh-spaced points on
    ``[0, 2*pi]``, then rolls the array so the dense region is centred at
    ``cluster_lon``.  The last (repeated) point is dropped to give a periodic
    grid.
    """
    start, stop = 0.0, 2.0 * np.pi

    if spacing == "uniform":
        return np.linspace(start, stop, nj, endpoint=False)

    if spacing == "tanh2":
        target = cluster_lon if cluster_lon is not None else np.pi
        # 1. Build nj+1 symmetric tanh2 points on [0, 2*pi] (endpoint included)
        coords = tanh_two_sided(nj + 1, start, stop, beta)  # dense at both ends
        # The symmetric law has its dense regions at coords[0] and coords[-1].
        # We want the dense region at `target`.  Roll so that the midpoint of
        # the dense zone aligns with `target`.
        # Dense zone centre is at coords[0] = 0.  Shift by target.
        coords = (coords + target) % stop
        # Re-sort to restore monotone order in [0, 2*pi)
        coords = np.sort(coords[:-1])  # drop duplicate endpoint after wrap
        return coords

    raise ParameterError(f"Unknown j-spacing law '{spacing}'. Choose uniform or tanh2.")


def build_k_steps(
    nk: int,
    ds: float,
    growth: float = 1.15,
    spacing: str = "geometric",
    beta: float = 3.0,
    total_thickness: float | None = None,
) -> np.ndarray:
    """Return ``(nk-1,)`` array of wall-normal step sizes.

    Parameters
    ----------
    nk :
        Total number of layers including k=0.
    ds :
        First-layer thickness (geometric law) or ignored for tanh/tanh2.
    growth :
        Geometric growth ratio (geometric law only).
    spacing :
        ``'geometric'`` — ``ds * growth**k`` (default, backward-compatible).
        ``'tanh'``      — one-sided Roberts from 0 to ``total_thickness``:
                          cluster near the wall.
        ``'tanh2'``     — two-sided Roberts from 0 to ``total_thickness``:
                          cluster at both wall and far-field.
    beta :
        Clustering strength for tanh/tanh2.
    total_thickness :
        Total wall-normal extent.  Required for tanh/tanh2; ignored for
        geometric.

    Returns
    -------
    numpy.ndarray, shape (nk-1,)
        Physical step sizes, all positive.
    """
    if spacing == "geometric":
        k = np.arange(nk - 1)
        return ds * growth**k

    if total_thickness is None:
        raise ParameterError("total_thickness is required for tanh/tanh2 k-spacing.")

    if spacing == "tanh":
        coords = tanh_one_sided(nk, 0.0, total_thickness, beta)
    elif spacing == "tanh2":
        coords = tanh_two_sided(nk, 0.0, total_thickness, beta)
    else:
        raise ParameterError(
            f"Unknown k-spacing law '{spacing}'. Choose geometric, tanh, or tanh2."
        )

    return np.diff(coords)  # (nk-1,) step sizes, all positive


# ---------------------------------------------------------------------------
# Pinned-endpoint distribution
# ---------------------------------------------------------------------------


def pinned_distribution(
    n: int,
    start: float,
    stop: float,
    first_step: float,
    last_step: float,
    spacing: str = "uniform",
    beta: float = 3.0,
) -> np.ndarray:
    """Return ``(n,)`` coords with pinned first and last intervals.

    ``out[1] - out[0] == first_step`` and ``out[-1] - out[-2] == last_step``;
    the interior ``n-2`` coords lie on ``[start + first_step, stop - last_step]``
    distributed by the chosen law.

    Used for C1 spacing matching at multi-block seams: each edge whose two
    endpoints share a step size with neighbouring blocks needs the first and
    last cells pinned to those neighbour steps, while the interior is free
    to follow a clustering law.

    Parameters
    ----------
    n :
        Number of points (>= 2).
    start, stop :
        Interval endpoints (inclusive).
    first_step, last_step :
        Required first and last interval lengths (both > 0).
    spacing :
        Interior law: ``'uniform'``, ``'tanh'``, or ``'tanh2'``.
    beta :
        Clustering strength for tanh/tanh2 (ignored for uniform).

    Returns
    -------
    numpy.ndarray, shape ``(n,)``

    Raises
    ------
    ParameterError
        If ``n < 2``, ``first_step`` or ``last_step`` is non-positive, or the
        two pinned intervals don't fit inside ``[start, stop]``.

    Notes
    -----
    Edge cases:
        - ``n == 2``: law is unused; ``out = [start, stop]`` and the caller
          must ensure ``first_step == stop - start`` (no constraint on last).
        - ``n == 3``: interior collapses to a single coord at
          ``start + first_step``, requiring
          ``first_step + last_step == stop - start``.
    """
    if n < 2:
        raise ParameterError(f"n must be >= 2 for pinned_distribution, got {n}")
    if first_step <= 0.0 or last_step <= 0.0:
        raise ParameterError(
            f"first_step and last_step must be positive (got {first_step}, {last_step})"
        )

    out = np.empty(n, dtype=np.float64)
    out[0] = start
    out[-1] = stop
    if n == 2:
        # Both pins collapse onto the single interval [start, stop].  The
        # caller is responsible for ensuring first_step matches the span;
        # last_step is unconstrained at n==2.
        return out

    span = stop - start
    if first_step + last_step > span + 1e-12:
        raise ParameterError(
            f"first_step + last_step ({first_step + last_step:.6g}) exceeds "
            f"interval span ({span:.6g})"
        )
    out[1] = start + first_step
    out[-2] = stop - last_step
    if n == 3:
        return out
    inner_start = start + first_step
    inner_stop = stop - last_step
    if spacing == "uniform":
        out[1:-1] = np.linspace(inner_start, inner_stop, n - 2)
    elif spacing == "tanh":
        out[1:-1] = tanh_one_sided(n - 2, inner_start, inner_stop, beta)
    elif spacing == "tanh2":
        out[1:-1] = tanh_two_sided(n - 2, inner_start, inner_stop, beta)
    else:
        raise ParameterError(
            f"Unknown interior spacing law '{spacing}' for pinned_distribution."
            " Choose uniform, tanh, or tanh2."
        )
    return out
