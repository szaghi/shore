"""Data charts of the docs, computed with SHORE itself.

* ``spacing-laws``: node positions and cell sizes of the 1-D spacing laws
  (``Spacing1D``) on ``[0, 1]`` with 21 nodes.
* ``quality-histogram``: scaled Jacobian and wall orthogonality of the sphere
  tutorial grid (cubed sphere, ``ni=33, nj=64, nk=25``, geometric k-steps
  ``ds=0.01, growth=1.15``), for both marching schemes.
* ``load-balance``: per-rank cells of that grid on 8 ranks, before and after
  ``plan_balance`` + ``split_blocks`` (assignment by ``write_proc_input``).

Grids and intermediate files go to a temporary directory, never to the repo.
"""

from __future__ import annotations

import atexit
import shutil
import tempfile
from functools import cache
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from theme import save, style

from shore import quality
from shore.mesh.spacing import Spacing1D

N_NODES = 21
NP = 8
BLOCKS_PER_RANK = 4.0


# ---------------------------------------------------------------------------
# spacing-laws
# ---------------------------------------------------------------------------


def _laws() -> list[tuple[str, np.ndarray]]:
    n = N_NODES
    growth = 1.15
    ds = (growth - 1.0) / (growth ** (n - 1) - 1.0)  # geometric steps summing to 1
    geometric = np.r_[0.0, np.cumsum(Spacing1D("geometric", ds=ds, growth=growth).build_steps(n))]
    return [
        ("uniform", Spacing1D("uniform").build(n)),
        ("tanh  $\\beta = 3$", Spacing1D("tanh", beta=3.0).build(n)),
        ("tanh2  $\\beta = 3$", Spacing1D("tanh2", beta=3.0).build(n)),
        (f"geometric  $g = {growth}$", geometric),
        (
            "pinned  0.02 | 0.08",
            Spacing1D("pinned", first_step=0.02, last_step=0.08, interior="uniform").build(n),
        ),
    ]


def spacing_laws(mode: str) -> None:
    with style(mode) as pal:
        laws = _laws()
        colours = [pal.slate, pal.ocean, pal.sea, pal.sunset, pal.plum]
        fig, (strips, sizes) = plt.subplots(
            2, 1, figsize=(7.2, 5.0), sharex=True, gridspec_kw={"height_ratios": [1.0, 1.15]}
        )
        for row, ((name, x), colour) in enumerate(zip(laws, colours, strict=True)):
            y = len(laws) - 1 - row
            strips.plot([0, 1], [y, y], color=pal.faint, lw=1.0, zorder=1)
            strips.vlines(x, y - 0.24, y + 0.24, color=colour, lw=1.3, zorder=2)
            strips.text(-0.02, y, name, ha="right", va="center", fontsize=9.5, color=pal.ink)
            mid = 0.5 * (x[1:] + x[:-1])
            sizes.plot(mid, np.diff(x), "o-", color=colour, lw=1.5, ms=3.2, label=name)
        strips.set_ylim(-0.6, len(laws) - 0.4)
        strips.set_yticks([])
        strips.grid(False)
        for side in ("left", "bottom"):
            strips.spines[side].set_visible(False)
        strips.tick_params(axis="x", length=0)
        sizes.set_xlim(0, 1)
        sizes.set_ylim(0, None)
        sizes.set_xlabel("position $x$ (cell sizes at cell centres)")
        sizes.set_ylabel("cell size $\\Delta x$")
        sizes.grid(axis="x", visible=False)
        fig.subplots_adjust(hspace=0.08, left=0.22)
        save(fig, "spacing-laws", mode, "svg")


# ---------------------------------------------------------------------------
# The sphere tutorial grid
# ---------------------------------------------------------------------------


@cache
def _workdir() -> Path:
    work = Path(tempfile.mkdtemp(prefix="shore-charts-"))
    atexit.register(shutil.rmtree, work, ignore_errors=True)
    return work


@cache
def _sphere_mesh(scheme: str):
    from shore.mesh import CubedSphereTopology, Mesh
    from shore.shapes import sphere

    stl = _workdir() / "sphere.stl"
    if not stl.exists():
        sphere(radius=1.0, subdivisions=4).export(stl)
    mesh = Mesh.from_stl(
        stl,
        topology=CubedSphereTopology(),
        ni=33,
        nj=64,
        nk=25,
        theta_cap=float(np.deg2rad(35.0)),  # as tutorials/sphere
        spacing_k=Spacing1D(law="geometric", ds=0.01, growth=1.15),
    )
    mesh.march(ds=0.01, growth=1.15, march_scheme=scheme)
    return mesh


# ---------------------------------------------------------------------------
# quality-histogram
# ---------------------------------------------------------------------------


def quality_histogram(mode: str) -> None:
    with style(mode) as pal:
        fig, (left, right) = plt.subplots(1, 2, figsize=(7.4, 3.0))
        sj_bins = np.linspace(-1.0, 1.0, 81)
        ko_bins = np.linspace(0.0, 8.0, 65)
        series = (
            ("algebraic", pal.coral, "algebraic (default)"),
            ("hyperbolic", pal.ocean, "hyperbolic"),
        )
        for scheme, colour, label in series:
            blocks = _sphere_mesh(scheme).blocks
            sj = np.concatenate([quality.scaled_jacobian(b.nodes).ravel() for b in blocks])
            ko = np.concatenate([quality.k_orthogonality_deg(b.nodes).ravel() for b in blocks])
            folded = int((sj <= 0.0).sum())
            for ax, data, bins in ((left, sj, sj_bins), (right, ko, ko_bins)):
                ax.hist(
                    data,
                    bins=bins,
                    histtype="stepfilled",
                    color=pal.tint(colour, 0.8),
                    alpha=0.6,
                    lw=0,
                )
                ax.hist(
                    data,
                    bins=bins,
                    histtype="step",
                    color=colour,
                    lw=1.3,
                    label=f"{label}: {folded} folded" if ax is left else label,
                )
        left.axvspan(-1.0, 0.0, color=pal.tint(pal.coral, 0.9), lw=0, zorder=0)
        left.axvline(0.0, color=pal.muted, lw=0.8)
        left.text(
            -0.5,
            0.97,
            "folded",
            transform=left.get_xaxis_transform(),
            ha="center",
            va="top",
            fontsize=8.8,
            color=pal.muted,
        )
        left.set_xlabel("scaled Jacobian (min corner per cell)")
        left.set_ylabel("cells")
        left.set_yscale("log")
        left.set_xlim(-1.0, 1.0)
        right.set_xlabel("wall orthogonality (deg)")
        right.set_ylabel("wall nodes")
        right.set_xlim(0.0, 8.0)
        for ax in (left, right):
            ax.grid(axis="x", visible=False)
        left.legend(loc="upper left", fontsize=8.6, handlelength=1.4, bbox_to_anchor=(0.0, 0.9))
        fig.subplots_adjust(wspace=0.3)
        save(fig, "quality-histogram", mode, "svg")


# ---------------------------------------------------------------------------
# load-balance
# ---------------------------------------------------------------------------


@cache
def _balance() -> tuple:
    from shore.balance import plan_balance
    from shore.io.adjacency import write_adjacency_json
    from shore.io.geo import write_geo
    from shore.io.grd import write_grd
    from shore.io.proc_input import write_proc_input
    from shore.split import BlockSplit, split_blocks

    work = _workdir()
    blocks = _sphere_mesh("algebraic").blocks

    def pack(blocks, stem):
        geos = []
        for b in blocks:
            geo = work / f"{stem}_{b.label}.geo"
            write_geo(geo, b.nodes)
            geos.append(geo)
        return write_grd(work / f"{stem}.grd", geos)

    grd = pack(blocks, "whole")
    adjacency = write_adjacency_json(blocks, work / "whole.adjacency.json")
    _, before = write_proc_input(grd, work / "whole.proc.input", np=NP)
    plan = plan_balance(
        grd, np=NP, blocks_per_rank=BLOCKS_PER_RANK, adjacency_path=adjacency, blocks=blocks
    )
    cuts: dict[tuple[str, str], list[int]] = {}
    for label, axis, at in plan.splits:
        cuts.setdefault((label, axis), []).append(at)
    split = split_blocks(blocks, [BlockSplit(lab, ax, at) for (lab, ax), at in cuts.items()])
    _, after = write_proc_input(pack(split, "split"), work / "split.proc.input", np=NP)
    return before, after


def load_balance(mode: str) -> None:
    before, after = _balance()
    with style(mode) as pal:
        fig, ax = plt.subplots(figsize=(7.2, 3.0))
        ranks = np.arange(NP)
        width = 0.38
        mean = sum(before.per_rank_cells) / NP / 1e3
        for offset, summary, colour, label in (
            (-width / 2, before, pal.slate, "whole blocks"),
            (width / 2, after, pal.ocean, f"after balance, {BLOCKS_PER_RANK:g} blocks/rank"),
        ):
            cells = np.array(summary.per_rank_cells) / 1e3
            ax.bar(
                ranks + offset,
                cells,
                width=width * 0.92,
                color=pal.tint(colour, 0.55),
                ec=colour,
                lw=1.0,
                label=f"{label}: {summary.n_blocks} blocks, "
                f"imbalance {summary.imbalance * 100:.0f} %",
            )
        ax.axhline(mean, color=pal.ink, lw=1.0, ls=(0, (4, 3)))
        ax.text(NP - 0.45, mean, " mean", va="center", ha="left", fontsize=8.8, color=pal.ink)
        ax.set_xticks(ranks)
        ax.set_xlim(-0.6, NP - 0.4)
        ax.set_xlabel("MPI rank")
        ax.set_ylabel("cells per rank (thousands)")
        ax.grid(axis="x", visible=False)
        ax.legend(loc="upper right", fontsize=8.8, handlelength=1.2)
        ax.set_ylim(0, max(before.per_rank_cells) / 1e3 * 1.35)
        save(fig, "load-balance", mode, "svg")


FIGURES = {
    "spacing-laws": spacing_laws,
    "quality-histogram": quality_histogram,
    "load-balance": load_balance,
}
