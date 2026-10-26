"""Meshing pipeline steps shared by the CLI and the command layer (ADR 0004).

``build_topology`` turns a topology name plus parameters into wired blocks.
It is the single code path behind ``shore mesh`` and
``Project.build_topology``, so a journal replay produces the same blocks as
the CLI run that recorded it.
"""

from __future__ import annotations

import inspect
from collections.abc import Callable
from typing import Any

import numpy as np

from shore.exceptions import ParameterError
from shore.mesh.block import HexBlock
from shore.mesh.body import Body
from shore.mesh.spacing import Spacing1D
from shore.mesh.topology import (
    CGridTopology,
    CHTopology,
    CubedSphereTopology,
    OGridTopology,
    OHTopology,
    Topology,
)

TOPOLOGIES = ("cubed_sphere", "ogrid", "cgrid", "ch", "oh")


def k_spacing(
    law: str = "geometric",
    *,
    ds: float = 1e-3,
    growth: float = 1.15,
    beta: float = 3.0,
    thickness: float | None = None,
) -> Spacing1D:
    """Wall-normal spacing law, as ``shore mesh --volume-k-*`` builds it."""
    if law in {"tanh", "tanh2"} and thickness is None:
        raise ParameterError(f"k-spacing law '{law}' needs a total thickness.")
    return Spacing1D(law=law, ds=ds, growth=growth, beta=beta, total_thickness=thickness)


def make_topology(name: str, options: dict[str, Any] | None = None) -> Topology:
    """Topology object for *name*; *options* are its constructor arguments.

    ``ch`` takes ``chord_min/max``, ``normal_min/max``, ``nk_channel``,
    ``freestream_world``, ``span_axis``, ``k_law_channel``,
    ``k_beta_channel``, ``i_split_body_upper/lower``, ``ni_outflow``,
    ``c1_seam``; ``oh`` takes ``box_vmin``, ``box_vmax``, ``nk_channel``,
    ``k_law_channel``, ``k_beta_channel``, ``c1_seam``; the others take none.
    """
    opts = dict(options or {})
    if name == "cubed_sphere":
        cls: type[Topology] = CubedSphereTopology
    elif name == "ogrid":
        cls = OGridTopology
    elif name == "cgrid":
        cls = CGridTopology
    elif name == "ch":
        cls = CHTopology
        if opts.get("freestream_world") is not None:
            opts["freestream_world"] = np.array(opts["freestream_world"], dtype=float)
    elif name == "oh":
        cls = OHTopology
    else:
        raise ParameterError(f"Unknown topology '{name}' (expected one of {TOPOLOGIES}).")
    _check_kwargs(cls, opts, f"topology '{name}' options")
    return cls(**opts)


def build_topology(
    body: Body,
    name: str,
    *,
    nk: int,
    spacing: Spacing1D,
    options: dict[str, Any] | None = None,
    **params: Any,
) -> tuple[Topology, list[HexBlock]]:
    """Build the wired blocks of topology *name* on *body*.

    *params* go to the topology's ``build``: ``ni``, ``nj``, ``theta_cap``,
    ``i_spacing``, ``i_beta``, ``j_spacing``, ``j_beta`` for the ray-cast
    topologies (``cubed_sphere``, ``ogrid``, ``oh``); ``ni_body``,
    ``ni_wake``, ``n_stations``, ``wake_length``, ``wake_growth``,
    ``span_axis``, ``te_dihedral_deg``, ``blocks`` and the ``cap_*``
    options for ``cgrid``; the C-grid subset for ``ch``.
    """
    topo = make_topology(name, options)
    _check_kwargs(topo.build, params, f"topology '{name}' build")
    if name in ("cgrid", "ch"):
        blocks = topo.build(surface=None, nk=nk, mesh=body.mesh, spacing_k=spacing, **params)
    else:
        blocks = topo.build(
            surface=None, nk=nk, spacing_k=spacing, projector=body.projector, **params
        )
    return topo, blocks


def _check_kwargs(fn: Callable[..., Any], kwargs: dict[str, Any], what: str) -> None:
    """Reject keyword names that are not named parameters of *fn*.

    The topology ``build`` methods end in a ``**kwargs`` catch-all; checking
    against the named parameters keeps a misspelt journal argument from
    being silently ignored.
    """
    named = {
        n
        for n, p in inspect.signature(fn).parameters.items()
        if p.kind not in (inspect.Parameter.VAR_KEYWORD, inspect.Parameter.VAR_POSITIONAL)
    }
    unknown = sorted(set(kwargs) - named)
    if unknown:
        raise ParameterError(f"{what}: unknown parameter(s) {', '.join(unknown)}.")
