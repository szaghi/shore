"""Generate all figures for the SHORE docs.

Outputs land under ``docs/public/figures/`` as PNGs.  Run from the project
root with the venv active::

    .venv/bin/python docs/scripts/make_figures.py

The script is idempotent — re-running overwrites the figures.  Heavy 3D
renders (cap surfaces, exploded views, marching layers) are produced with
PyVista in off-screen mode; schematic 2D figures are produced with
matplotlib.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

# Force off-screen rendering for headless / WSL environments before any
# pyvista import.
os.environ.setdefault("DISPLAY", "")

import matplotlib.pyplot as plt
import numpy as np
import trimesh

# Project the docs/public/figures path so we can dump there.
ROOT = Path(__file__).resolve().parents[2]
FIG = ROOT / "docs" / "public" / "figures"
FIG.mkdir(parents=True, exist_ok=True)


# =========================================================================
# Schematic figures (matplotlib)
# =========================================================================


def _save(fig: plt.Figure, name: str) -> Path:
    out = FIG / f"{name}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  wrote {out.relative_to(ROOT)}")
    return out


def fig_pipeline() -> None:
    """High-level pipeline flowchart: STL → cap blocks + equator → march → 6 .geo."""
    fig, ax = plt.subplots(figsize=(10, 3.4))
    ax.axis("off")
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 3.4)

    boxes = [
        (0.2, 1.2, 1.6, 1.0, "STL\n(triangulated\nbody surface)", "#e8f0fc"),
        (2.4, 2.0, 1.6, 1.0, "Cap k=0\n(gnomonic\nflat-square)", "#f9d6c5"),
        (2.4, 0.4, 1.6, 1.0, "Equator k=0\n(great-circle\nmeridians)", "#cfe9c9"),
        (4.6, 1.2, 1.6, 1.0, "6-block O-grid\nk=0 surface", "#dfe8d2"),
        (6.8, 1.2, 1.6, 1.0, "Hyperbolic\nmarching\n(nk layers)", "#f5e1a4"),
        (8.6, 1.2, 1.2, 1.0, "6 .geo\nfiles", "#cfd8e0"),
    ]
    for (x, y, w, h, label, c) in boxes:
        ax.add_patch(plt.Rectangle((x, y), w, h, facecolor=c, edgecolor="#3a3a3a", linewidth=1.2))
        ax.text(x + w / 2, y + h / 2, label, ha="center", va="center", fontsize=9)

    arrows = [
        (1.8, 1.7, 2.4, 2.5),     # STL → cap
        (1.8, 1.7, 2.4, 0.9),     # STL → equator
        (4.0, 2.5, 4.6, 1.85),    # cap → 6-block (top)
        (4.0, 0.9, 4.6, 1.55),    # equator → 6-block (bottom)
        (6.2, 1.7, 6.8, 1.7),     # 6-block → march
        (8.4, 1.7, 8.6, 1.7),     # march → .geo
    ]
    for (x0, y0, x1, y1) in arrows:
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle="->", color="#3a3a3a", lw=1.4))

    ax.text(2.6, 3.2, "k=0 construction (cubed-sphere topology)", fontsize=8.5,
            color="#3a3a3a", style="italic")
    ax.text(7.0, 3.2, "volume extrusion", fontsize=8.5,
            color="#3a3a3a", style="italic")

    _save(fig, "pipeline")


def fig_anchor_raycast() -> None:
    """Anchor-based ray-cast onto a star-shaped body (2D cross-section)."""
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.set_aspect("equal"); ax.axis("off")

    # Body silhouette: a slightly distorted ellipse-like blob.
    t = np.linspace(0, 2 * np.pi, 200)
    r = 2.0 + 0.25 * np.cos(3 * t) - 0.15 * np.sin(2 * t)
    bx, by = r * np.cos(t), r * np.sin(t)
    ax.fill(bx, by, color="#e8e8e8", edgecolor="#555", linewidth=1.4, zorder=1)

    anchor = np.array([0.1, -0.05])
    ax.plot(*anchor, "o", color="#c33", markersize=8, zorder=4)
    ax.text(anchor[0] + 0.15, anchor[1] - 0.15, "anchor", color="#c33", fontsize=10)

    # Fan of rays.
    n_rays = 18
    for i in range(n_rays):
        ang = 2 * np.pi * i / n_rays + 0.1
        # Find intersection with body silhouette analytically (sort-of):
        # interpolate r(ang) by lookup.
        r_hit = 2.0 + 0.25 * np.cos(3 * ang) - 0.15 * np.sin(2 * ang)
        end = anchor + r_hit * np.array([np.cos(ang), np.sin(ang)])
        ax.plot([anchor[0], end[0]], [anchor[1], end[1]],
                color="#666", lw=0.7, alpha=0.6, zorder=2)
        ax.plot(end[0], end[1], "o", color="#1a73e8", markersize=4, zorder=3)

    ax.text(0, 2.7, "k=0 surface nodes",
            ha="center", color="#1a73e8", fontsize=10)
    ax.text(0, -2.9, "Body STL (star-shaped from anchor)",
            ha="center", color="#444", fontsize=9, style="italic")

    ax.set_xlim(-3.0, 3.0); ax.set_ylim(-3.2, 3.0)
    _save(fig, "projection-raycast")


def fig_gnomonic_cap() -> None:
    """Tangent-plane diamond + radial projection schematic."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))

    # Left panel: tangent-plane diamond (top-down).
    ax = axes[0]; ax.set_aspect("equal"); ax.axis("off")
    h = 1.0  # = tan(theta_cap), schematic
    corners = np.array([(h, 0), (0, h), (-h, 0), (0, -h), (h, 0)])
    ax.plot(corners[:, 0], corners[:, 1], color="#2a4d8f", lw=2)
    ax.fill(corners[:, 0], corners[:, 1], color="#dfe8f5", alpha=0.5)
    # Bilinear grid of dots, nc=7.
    nc = 7
    for i in range(nc):
        for j in range(nc):
            a = i / (nc - 1); b = j / (nc - 1)
            x = h * (1 - a - b); y = h * (a - b)
            ax.plot(x, y, "o", color="#2a4d8f", markersize=3.5)
    # Annotate corners.
    for (cx, cy, label, off) in [
        (h, 0, "C0\n(φ=0°)",  (0.1, -0.18)),
        (0, h, "C1\n(φ=90°)", (-0.05, 0.12)),
        (-h, 0, "C2\n(φ=180°)", (-0.55, -0.18)),
        (0, -h, "C3\n(φ=270°)", (-0.05, -0.4)),
    ]:
        ax.text(cx + off[0], cy + off[1], label, fontsize=9, color="#2a4d8f")
    ax.set_xlim(-1.4, 1.4); ax.set_ylim(-1.5, 1.4)
    ax.set_title("Tangent-plane diamond at the cap pole",
                 fontsize=10, color="#222222")
    ax.text(0, -1.45, "$h = \\tan(\\theta_{cap})$, half-side along ±x and ±y",
            ha="center", fontsize=9, color="#444", style="italic")

    # Right panel: side view with radial projection onto the sphere.
    ax = axes[1]; ax.set_aspect("equal"); ax.axis("off")
    th = np.linspace(0, np.pi, 200)
    R = 1.0
    sx = R * np.sin(th); sy = R * np.cos(th)
    ax.plot(sx, sy, color="#444", lw=1.4)
    ax.plot(-sx, sy, color="#444", lw=1.4)
    ax.fill(np.concatenate([sx, -sx[::-1]]),
            np.concatenate([sy, sy[::-1]]),
            color="#f3f3f3", zorder=0)

    # Anchor
    ax.plot(0, 0, "o", color="#c33", markersize=7)
    ax.text(0.05, -0.12, "anchor", color="#c33", fontsize=9)

    # Tangent square at top.
    th_cap = np.deg2rad(30)
    h_tan = np.tan(th_cap)  # half-side
    # Square at z = 1 (in pole-local, which in our 2D side view sits at y=cos(th_cap))
    z_plane = 1.0
    ax.plot([-h_tan, h_tan], [z_plane, z_plane], color="#2a4d8f", lw=2)
    ax.plot([h_tan, h_tan], [z_plane, z_plane], "o", color="#2a4d8f", markersize=5)
    ax.plot([-h_tan, -h_tan], [z_plane, z_plane], "o", color="#2a4d8f", markersize=5)
    # Sample interior nodes on tangent line.
    for u in np.linspace(-h_tan, h_tan, 7):
        ax.plot(u, z_plane, "o", color="#2a4d8f", markersize=4)
    # Rays from anchor through tangent nodes; intersect sphere of radius 1.
    for u in np.linspace(-h_tan, h_tan, 7):
        d = np.array([u, z_plane])
        d_unit = d / np.linalg.norm(d)
        hit = R * d_unit
        ax.plot([0, u], [0, z_plane], color="#888", lw=0.6, alpha=0.6)
        ax.plot(hit[0], hit[1], "o", color="#1a73e8", markersize=5)

    # Annotations
    ax.annotate("flat tangent square\n(z = 1, half-side $h = \\tan \\theta_{cap}$)",
                xy=(0, z_plane + 0.02), xytext=(-1.6, 1.35),
                fontsize=9, color="#2a4d8f",
                arrowprops=dict(arrowstyle="->", color="#2a4d8f", lw=1.0))
    ax.annotate("cap k=0 patch\non body surface",
                xy=(0.45, 0.85), xytext=(1.15, 0.55),
                fontsize=9, color="#1a73e8",
                arrowprops=dict(arrowstyle="->", color="#1a73e8", lw=1.0))
    ax.text(0, -1.35, "Right: radial projection of the diamond onto the body. "
                       "Corners at $\\theta_{cap}$; edges dip poleward at midpoints.",
            ha="center", fontsize=9, color="#444", style="italic")

    ax.set_xlim(-1.7, 1.7); ax.set_ylim(-1.4, 1.5)
    ax.set_title("Radial gnomonic projection (side view)",
                 fontsize=10, color="#222222")

    _save(fig, "gnomonic-cap")


def fig_3edge_junction() -> None:
    """3-edge cap-corner junction with first-cell lengths annotated."""
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    ax.set_aspect("equal"); ax.axis("off")

    # Draw the cap region as a small upper-left square with great-circle-ish curves.
    corner = np.array([0.0, 0.0])
    # Cap edge 1 (along +x of the cap, going east along the seam).
    e1_end = np.array([2.5, 0.3])
    e1_first = np.array([0.7, 0.07])
    ax.plot([corner[0], e1_end[0]], [corner[1], e1_end[1]],
            color="#c33", lw=1.5)
    ax.plot([corner[0], e1_first[0]], [corner[1], e1_first[1]],
            color="#c33", lw=4)
    # Cap edge 2 (along +y of the cap, perpendicular along the seam).
    e2_end = np.array([-0.3, 2.5])
    e2_first = np.array([-0.07, 0.7])
    ax.plot([corner[0], e2_end[0]], [corner[1], e2_end[1]],
            color="#c33", lw=1.5)
    ax.plot([corner[0], e2_first[0]], [corner[1], e2_first[1]],
            color="#c33", lw=4)
    # Equator meridional edge (going outward / inward into the equator band).
    e3_end = np.array([2.0, -2.5])
    e3_first = np.array([0.55, -0.7])
    ax.plot([corner[0], e3_end[0]], [corner[1], e3_end[1]],
            color="#1a73e8", lw=1.5)
    ax.plot([corner[0], e3_first[0]], [corner[1], e3_first[1]],
            color="#1a73e8", lw=4)

    # The 3-edge vertex.
    ax.plot(*corner, "o", color="#222222", markersize=12, zorder=5)
    # Mark first cells.
    ax.plot(*e1_first, "o", color="#c33", markersize=7, zorder=4)
    ax.plot(*e2_first, "o", color="#c33", markersize=7, zorder=4)
    ax.plot(*e3_first, "o", color="#1a73e8", markersize=7, zorder=4)

    # Length annotations.
    ax.text(0.4, 0.22, "$\\Delta s$", fontsize=12, color="#c33")
    ax.text(0.16, 0.4, "$\\Delta s$", fontsize=12, color="#c33")
    ax.text(0.4, -0.5, "$\\Delta s$", fontsize=12, color="#1a73e8")

    # Edge labels.
    ax.text(2.0, 0.55, "cap edge\n(seam, +x)", color="#c33", fontsize=10)
    ax.text(-1.7, 2.4, "cap edge\n(seam, +y)", color="#c33", fontsize=10)
    ax.text(2.0, -2.7, "equator meridian\n(into equator)", color="#1a73e8", fontsize=10)

    ax.text(0, 3.1, "3-edge cap-corner junction\n(all three first-cells share one length)",
            ha="center", fontsize=11, color="#222222")

    ax.set_xlim(-2.5, 3.0); ax.set_ylim(-3.2, 3.2)
    _save(fig, "junction-3edge")


def fig_per_j_pinning() -> None:
    """Meridional profile: pinned first/last interval vs. fallback."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.6))
    ni = 10
    L = 1.0  # total interval length

    # ---- Compatible: cap first step ~ natural step; pinning active ----
    # Mimics _constrained_distribution: pin first and last interval;
    # interior ni-3 intervals share the residual evenly.
    ax = axes[0]
    nat_step = L / (ni - 1)             # natural uniform step
    cap_step = 0.85 * nat_step          # within the 2x compatibility band
    interior_total = L - 2 * cap_step
    interior_step = interior_total / (ni - 3)
    diffs = np.empty(ni - 1)
    diffs[0] = cap_step
    diffs[-1] = cap_step
    diffs[1:-1] = interior_step
    ax.plot(diffs, "o-", color="#1a73e8", ms=6, lw=1.5,
            label="equator meridional step")
    ax.axhline(cap_step, color="#c33", lw=1.2, ls="--",
               label="cap first-step (pin target)")
    ax.axhline(nat_step, color="#888", lw=1.0, ls=":",
               label="natural uniform step")
    ax.set_title("Compatible sizes — pinning active",
                 fontsize=10)
    ax.set_xlabel("interval index")
    ax.set_ylabel("$\\Delta\\theta$ (rad)")
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, max(diffs) * 1.6)

    # ---- Incompatible (fallback): uniform natural everywhere ----
    # Cap step is much smaller than the natural step -> pinning would
    # cause a huge interior jump, so the equator falls back to the
    # unpinned i_spacing law.  The cap step is shown for reference.
    ax = axes[1]
    cap_step_bad = 0.20 * nat_step      # ~5x smaller than natural
    diffs = np.full(ni - 1, nat_step)
    ax.plot(diffs, "o-", color="#1a73e8", ms=6, lw=1.5,
            label="equator step (unpinned fallback)")
    ax.axhline(cap_step_bad, color="#c33", lw=1.2, ls="--",
               label="cap first-step (incompatible)")
    ax.axhline(nat_step, color="#888", lw=1.0, ls=":",
               label="natural uniform step")
    ax.set_title("Incompatible sizes — fallback to natural law",
                 fontsize=10)
    ax.set_xlabel("interval index")
    ax.set_ylabel("$\\Delta\\theta$ (rad)")
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, nat_step * 1.6)

    fig.suptitle("Per-j C1 seam pinning: compatible vs. fallback regime",
                 fontsize=11)
    _save(fig, "per-j-pinning")


def fig_seam_aware_normals() -> None:
    """Sub-block seam: one-sided vs centred difference for the j-tangent."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))

    for ax, title, mode in zip(axes, ["one-sided (problem)", "centred via ghost row (fix)"],
                                ["onesided", "centred"]):
        ax.axis("off"); ax.set_aspect("equal")
        ax.set_xlim(-1, 5); ax.set_ylim(-1.0, 4.0)

        # Sub_q grid (left)
        for j in range(0, 4):
            for i in range(0, 4):
                ax.plot(j * 0.6, i * 0.7, "o", color="#1a73e8", ms=4)
        # Seam column (shared)
        for i in range(0, 4):
            ax.plot(2.0, i * 0.7, "s", color="#c33", ms=8)
        # Sub_{q+1} grid (right)
        for j in range(0, 4):
            for i in range(0, 4):
                ax.plot(2.0 + (j + 1) * 0.6, i * 0.7, "o", color="#0a8a3a", ms=4)

        ax.text(0.4, 3.4, "sub_q", color="#1a73e8", fontsize=11, ha="center")
        ax.text(3.5, 3.4, "sub_{q+1}", color="#0a8a3a", fontsize=11, ha="center")
        ax.text(2.0, 3.4, "seam", color="#c33", fontsize=11, ha="center")

        # Tangent vectors at the seam node (i=1).
        seam_y = 0.7
        if mode == "onesided":
            # Sub_q: backward only (uses j-1 inside its block).
            ax.annotate("", xy=(2.0, seam_y), xytext=(1.4, seam_y),
                        arrowprops=dict(arrowstyle="->", color="#1a73e8", lw=2))
            ax.text(1.7, seam_y + 0.15, "$t_q$ (one-sided)", color="#1a73e8", fontsize=9)
            ax.annotate("", xy=(2.6, seam_y), xytext=(2.0, seam_y),
                        arrowprops=dict(arrowstyle="->", color="#0a8a3a", lw=2))
            ax.text(2.2, seam_y + 0.15, "$t_{q+1}$ (one-sided)",
                    color="#0a8a3a", fontsize=9)
            ax.text(2.0, -0.5, "Different tangents → different normals → spacing jump",
                    ha="center", color="#c33", fontsize=10)
        else:
            # Both blocks compute the same centred difference.
            ax.annotate("", xy=(2.6, seam_y), xytext=(1.4, seam_y),
                        arrowprops=dict(arrowstyle="->", color="#222222", lw=2.5))
            ax.text(1.7, seam_y + 0.18, "$t_{centred}$",
                    color="#222222", fontsize=10)
            # Ghost row indicator on each side.
            ax.plot(2.6, seam_y, "o", mec="#222", mfc="none", ms=14, lw=2)
            ax.plot(1.4, seam_y, "o", mec="#222", mfc="none", ms=14, lw=2)
            ax.text(2.0, -0.5, "Ghost rows from partners → both blocks agree",
                    ha="center", color="#0a8a3a", fontsize=10)

        ax.set_title(title, fontsize=11)

    fig.suptitle("Seam-aware normals at the sub-sub j-seam", fontsize=11)
    _save(fig, "seam-aware-normals")


# =========================================================================
# 3D renders (PyVista off-screen)
# =========================================================================


def _make_mesh(theta_cap_deg: float = 30.0, ni: int = 40, nj: int = 60, nk: int = 20):
    """Build a unit-sphere mesh with sensible defaults for figures."""
    from shore.mesh import CubedSphereTopology, Mesh, Spacing1D
    stl = ROOT / "docs" / "scripts" / "_sphere_for_figs.stl"
    if not stl.exists():
        stl.parent.mkdir(parents=True, exist_ok=True)
        trimesh.creation.icosphere(subdivisions=4, radius=1.0).export(str(stl))
    m = Mesh.from_stl(
        str(stl),
        topology=CubedSphereTopology(),
        ni=ni, nj=nj, nk=nk,
        theta_cap=float(np.deg2rad(theta_cap_deg)),
        spacing_k=Spacing1D(law="geometric", ds=0.02, growth=1.1),
    )
    m.march(ds=0.02, growth=1.1)
    return m


def fig_caps_topology_blocks() -> None:
    """k=0 surface coloured by block: 4 sub-blocks + 2 caps."""
    import pyvista as pv
    pv.OFF_SCREEN = True

    m = _make_mesh()
    colours = ["#1a73e8", "#0a8a3a", "#f5a623", "#9c27b0", "#e74c3c", "#16a085"]
    pl = pv.Plotter(off_screen=True, window_size=(1100, 900))
    for blk, c in zip(m.blocks, colours):
        surf = blk.nodes[:, :, 0, :]
        sg = pv.StructuredGrid(surf[..., 0], surf[..., 1], surf[..., 2])
        pl.add_mesh(sg, color=c, show_edges=True, edge_color="#222222",
                    line_width=0.6)
    pl.camera_position = [(2.4, -2.4, 1.6), (0, 0, 0.1), (0, 0, 1)]
    pl.add_axes()
    out = FIG / "cubed-sphere-topology-blocks.png"
    pl.screenshot(str(out))
    pl.close()
    print(f"  wrote {out.relative_to(ROOT)}")


def fig_caps_topology_topdown() -> None:
    """Top-down view of the cap region — square cells, 3-edge junctions."""
    import pyvista as pv
    pv.OFF_SCREEN = True
    m = _make_mesh()
    colours = ["#cfe1f7", "#cfeec3", "#fce6c4", "#e6cef0", "#f7c4be", "#bfe5db"]
    pl = pv.Plotter(off_screen=True, window_size=(900, 900))
    for blk, c in zip(m.blocks, colours):
        surf = blk.nodes[:, :, 0, :]
        sg = pv.StructuredGrid(surf[..., 0], surf[..., 1], surf[..., 2])
        pl.add_mesh(sg, color=c, show_edges=True, edge_color="#1a3e8b",
                    line_width=0.7)
    pl.camera_position = [(0, 0, 4), (0, 0, 0), (0, 1, 0)]
    pl.add_axes()
    out = FIG / "cubed-sphere-topology-topdown.png"
    pl.screenshot(str(out))
    pl.close()
    print(f"  wrote {out.relative_to(ROOT)}")


def fig_marching_layers() -> None:
    """Cross-section showing wall-normal layers expanding outward."""
    import pyvista as pv
    pv.OFF_SCREEN = True
    m = _make_mesh()

    pl = pv.Plotter(off_screen=True, window_size=(1100, 800))
    # Slice each block by a thin XY band (z near 0) and render each layer surface.
    nk = m.blocks[0].nodes.shape[2]
    cmap = plt.cm.viridis
    for k in range(0, nk, 2):
        for blk in m.blocks:
            surf = blk.nodes[:, :, k, :]
            sg = pv.StructuredGrid(surf[..., 0], surf[..., 1], surf[..., 2])
            colour = cmap(k / max(nk - 1, 1))
            rgb = (colour[0], colour[1], colour[2])
            pl.add_mesh(sg, color=rgb, show_edges=True,
                        edge_color="#222222", line_width=0.4, opacity=0.55)
    pl.camera_position = [(3.5, -3.5, 2.4), (0, 0, 0), (0, 0, 1)]
    pl.add_axes()
    out = FIG / "marching-layers.png"
    pl.screenshot(str(out))
    pl.close()
    print(f"  wrote {out.relative_to(ROOT)}")


def fig_volume_seam_continuity() -> None:
    """Outer-layer surface — smooth across block seams (the slice.png fix)."""
    import pyvista as pv
    pv.OFF_SCREEN = True
    m = _make_mesh()

    pl = pv.Plotter(off_screen=True, window_size=(1100, 800))
    for blk in m.blocks:
        surf_in = blk.nodes[:, :, 0, :]
        sg_in = pv.StructuredGrid(surf_in[..., 0], surf_in[..., 1], surf_in[..., 2])
        pl.add_mesh(sg_in, color="#1a73e8", show_edges=True,
                    edge_color="#cc3333", line_width=0.5)
        surf_out = blk.nodes[:, :, -1, :]
        sg_out = pv.StructuredGrid(surf_out[..., 0], surf_out[..., 1], surf_out[..., 2])
        pl.add_mesh(sg_out, color="#f5a623", show_edges=True,
                    edge_color="#16a085", line_width=0.5, opacity=0.6)
    pl.camera_position = [(3.8, -3.8, 2.4), (0, 0, 0), (0, 0, 1)]
    pl.add_axes()
    out = FIG / "volume-seam-continuity.png"
    pl.screenshot(str(out))
    pl.close()
    print(f"  wrote {out.relative_to(ROOT)}")


# =========================================================================
# Main
# =========================================================================


def main() -> int:
    print("Generating SHORE documentation figures...\n")
    print("Schematic figures (matplotlib):")
    fig_pipeline()
    fig_anchor_raycast()
    fig_gnomonic_cap()
    fig_3edge_junction()
    fig_per_j_pinning()
    fig_seam_aware_normals()

    print("\n3D renders (pyvista off-screen):")
    fig_caps_topology_blocks()
    fig_caps_topology_topdown()
    fig_marching_layers()
    fig_volume_seam_continuity()

    print(f"\nAll figures written to {FIG.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
