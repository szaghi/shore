"""2-D schematics of the docs: ray casting, the cubed-sphere cap, seams, the block model.

Each function draws one figure for one site theme (``"light"`` or ``"dark"``)
and writes it as SVG through :func:`theme.save`. Colours come from the theme
palette only; the drawings follow the concept pages that show them
(``docs/concepts/projection.md``, ``k0-construction.md``,
``seam-aware-normals.md``, ``block-model.md``).
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Arc, FancyArrowPatch, Polygon
from theme import Palette, save, style

from shore.mesh.face import AxisMap
from shore.volume.spacing import build_i_coords, pinned_distribution

LW = 1.5  # outlines and main strokes
THIN = 0.8  # grid lines, rays, guides
SMALL = 8.8  # annotation font size


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _canvas(width: float, height: float, ncols: int = 1) -> tuple[plt.Figure, list[plt.Axes]]:
    fig, axes = plt.subplots(1, ncols, figsize=(width, height))
    axes = list(np.atleast_1d(axes))
    for ax in axes:
        ax.set_aspect("equal")
        ax.axis("off")
    return fig, axes


def _arrow(ax, p0, p1, colour, lw=LW, head=9, style_="-|>", z=4, ls="-"):
    ax.add_patch(
        FancyArrowPatch(
            tuple(p0),
            tuple(p1),
            arrowstyle=f"{style_},head_length=0.7,head_width=0.32",
            mutation_scale=head * 1.25,
            color=colour,
            lw=lw,
            linestyle=ls,
            shrinkA=0,
            shrinkB=0,
            zorder=z,
            joinstyle="miter",
            capstyle="butt",
        )
    )


def _dot(ax, p, colour, size=4.5, z=5, face=None, pal: Palette | None = None):
    ax.plot(
        *p,
        "o",
        ms=size,
        mfc=colour if face is None else face,
        mec=colour,
        mew=1.0,
        zorder=z,
    )


def _label(ax, p, text, colour, size=SMALL, **kw):
    kw.setdefault("ha", "center")
    kw.setdefault("va", "center")
    ax.text(p[0], p[1], text, color=colour, fontsize=size, **kw)


def _v(*xs: float) -> np.ndarray:
    """A coordinate tuple as an array (offsets added to points)."""
    return np.array(xs, dtype=float)


def _oblique(p, depth=0.5, angle_deg=32.0) -> np.ndarray:
    """Cabinet-style oblique projection of 3-D points (x right, y depth, z up)."""
    angle = np.deg2rad(angle_deg)
    p = np.asarray(p, dtype=float)
    x = p[..., 0] + depth * p[..., 1] * np.cos(angle)
    y = p[..., 2] + depth * p[..., 1] * np.sin(angle)
    return np.stack([x, y], axis=-1)


# ---------------------------------------------------------------------------
# projection-raycast
# ---------------------------------------------------------------------------


def _first_hits(curve: np.ndarray, anchor: np.ndarray, angles: np.ndarray):
    """First intersection of rays from *anchor* with a closed polyline, and all crossings."""
    a, b = curve, np.roll(curve, -1, axis=0)
    first, every = [], []
    for ang in angles:
        d = np.array([np.cos(ang), np.sin(ang)])
        e = b - a
        w = a - anchor
        den = d[0] * e[:, 1] - d[1] * e[:, 0]
        with np.errstate(divide="ignore", invalid="ignore"):
            t = (w[:, 0] * e[:, 1] - w[:, 1] * e[:, 0]) / den
            s = (w[:, 0] * d[1] - w[:, 1] * d[0]) / den
        ok = (t > 1e-9) & (s >= 0) & (s < 1)
        ts = np.sort(t[ok])
        first.append(anchor + ts[0] * d)
        every.append([anchor + tt * d for tt in ts])
    return np.array(first), every


def projection_raycast(mode: str) -> None:
    with style(mode) as pal:
        fig, (left, right) = _canvas(7.4, 3.5, 2)
        body_fill = pal.tint(pal.slate, 0.84)
        n_rays = 32

        # Left: a star-shaped section, every ray hits once.
        t = np.linspace(0, 2 * np.pi, 400, endpoint=False)
        r = 1.0 + 0.16 * np.cos(3 * t) - 0.09 * np.sin(2 * t)
        curve = np.c_[r * np.cos(t), r * np.sin(t)]
        anchor = np.zeros(2)
        angles = np.linspace(0, 2 * np.pi, n_rays, endpoint=False) + 0.05
        hits, _ = _first_hits(curve, anchor, angles)
        left.add_patch(
            Polygon(
                curve, closed=True, fc=body_fill, ec=pal.ink, lw=LW, zorder=1, joinstyle="round"
            )
        )
        for h in hits:
            left.plot(*np.c_[anchor, h], color=pal.muted, lw=THIN, zorder=2)
        left.plot(*hits.T, "o", ms=4.2, color=pal.ocean, zorder=4)
        _dot(left, anchor, pal.coral, size=6.5)
        _label(left, (-1.3, -1.05), "anchor", pal.coral, ha="left")
        left.plot(*np.c_[(-0.98, -1.0), anchor - (0.05, 0.05)], color=pal.coral, lw=THIN)
        _label(left, (0, -1.42), "star-shaped: one hit per ray", pal.ink, size=9.5)
        _label(left, (1.38, 0.82), "k = 0 nodes", pal.ocean, ha="left")
        left.annotate(
            "",
            xy=hits[3] + 0.02,
            xytext=(1.36, 0.82),
            arrowprops={"arrowstyle": "-", "color": pal.ocean, "lw": THIN},
        )
        left.set_xlim(-1.35, 1.95)
        left.set_ylim(-1.55, 1.3)

        # Right: a kidney section, the anchor cannot see the whole surface.
        t = np.linspace(0, 2 * np.pi, 600, endpoint=False)
        r = 1.0 + 0.10 * np.cos(2 * t) - 0.75 * np.exp(-(((t - np.pi / 2) / 0.3) ** 2))
        curve = np.c_[r * np.cos(t), r * np.sin(t)]
        anchor = np.array([-0.55, 0.0])
        hits, every = _first_hits(curve, anchor, angles)
        # Surface points that are a first hit of some ray: dense sampling.
        dense = np.linspace(0, 2 * np.pi, 2000, endpoint=False)
        seen, _ = _first_hits(curve, anchor, dense)
        dist = np.min(np.linalg.norm(curve[:, None, :] - seen[None, :, :], axis=-1), axis=1)
        hidden = dist > 0.01
        right.add_patch(
            Polygon(
                curve, closed=True, fc=body_fill, ec=pal.ink, lw=LW, zorder=1, joinstyle="round"
            )
        )
        hc = np.where(hidden[:, None], curve, np.nan)
        right.plot(*hc.T, color=pal.coral, lw=3.2, zorder=3, solid_capstyle="round")
        for h, ev in zip(hits, every, strict=True):
            right.plot(*np.c_[anchor, h], color=pal.muted, lw=THIN, zorder=2)
            if len(ev) > 1:
                far = ev[-1]
                right.plot(*np.c_[h, far], color=pal.muted, lw=THIN, ls=(0, (2, 2)), zorder=2)
                for p in ev[1:]:
                    _dot(right, p, pal.muted, size=4.2, face="none", z=4)
        right.plot(*hits.T, "o", ms=4.2, color=pal.ocean, zorder=4)
        _dot(right, anchor, pal.coral, size=6.5)
        _label(right, (0, -1.42), "not star-shaped: first hit kept", pal.ink, size=9.5)
        hid = curve[hidden]
        tip = hid[np.argmax(hid[:, 1])]
        _label(
            right, (1.02, 1.12), "hidden from the anchor:\nno node lands here", pal.coral, ha="left"
        )
        right.annotate(
            "",
            xy=tip + _v(0.04, 0.03),
            xytext=(1.0, 1.06),
            arrowprops={"arrowstyle": "-", "color": pal.coral, "lw": THIN},
        )
        right.set_xlim(-1.35, 1.95)
        right.set_ylim(-1.55, 1.3)
        fig.subplots_adjust(wspace=0.02)
        save(fig, "projection-raycast", mode, "svg")


# ---------------------------------------------------------------------------
# gnomonic-cap
# ---------------------------------------------------------------------------


def gnomonic_cap(mode: str) -> None:
    with style(mode) as pal:
        fig, (top, side) = _canvas(7.4, 3.6, 2)
        theta = np.deg2rad(30.0)  # schematic theta_cap; the default is 0.05 pi
        h = np.tan(theta)
        nc = 9

        # Left: the diamond in the tangent plane, seen from above the pole.
        a = np.linspace(0, 1, nc)
        ag, bg = np.meshgrid(a, a, indexing="ij")
        x = h * (1 - ag - bg)
        y = h * (ag - bg)
        top.add_patch(
            Polygon(
                [(h, 0), (0, h), (-h, 0), (0, -h)],
                closed=True,
                fc=pal.tint(pal.sky, 0.82),
                ec="none",
                zorder=0,
            )
        )
        for i in range(nc):
            w = LW if i in (0, nc - 1) else THIN
            c = pal.ink if i in (0, nc - 1) else pal.ocean
            top.plot(x[i, :], y[i, :], color=c, lw=w, zorder=2)
            top.plot(x[:, i], y[:, i], color=c, lw=w, zorder=2)
        # The seam parallel (colatitude theta_cap) is the circle of radius h.
        ring = np.linspace(0, 2 * np.pi, 300)
        top.plot(h * np.cos(ring), h * np.sin(ring), color=pal.muted, lw=THIN, ls=(0, (4, 3)))
        top.plot([-1.35 * h, 1.35 * h], [0, 0], color=pal.faint, lw=THIN, zorder=-1)
        top.plot([0, 0], [-1.35 * h, 1.35 * h], color=pal.faint, lw=THIN, zorder=-1)
        corners = [
            ((h, 0), "cap[0, 0]\n$\\varphi = 0°$", (0.1, 0.0), "left"),
            ((0, h), "cap[-1, 0]\n$\\varphi = 90°$", (0.0, 0.17), "center"),
            ((-h, 0), "cap[-1, -1]\n$\\varphi = 180°$", (-0.1, 0.0), "right"),
            ((0, -h), "cap[0, -1]\n$\\varphi = 270°$", (0.0, -0.19), "center"),
        ]
        for p, text, off, ha in corners:
            _dot(top, p, pal.ink, size=5)
            _label(top, (p[0] + off[0], p[1] + off[1]), text, pal.ink, ha=ha, linespacing=1.3)
        # Half-diagonal h, dimensioned below the diagonal.
        yd = -0.62 * h
        for xe in (0.0, h):
            top.plot(
                [xe, xe], [yd - 0.05, -0.02 if xe else -h], color=pal.faint, lw=THIN, zorder=-1
            )
        top.annotate(
            "",
            xy=(h, yd),
            xytext=(0, yd),
            arrowprops={
                "arrowstyle": "<|-|>,head_length=0.5,head_width=0.18",
                "color": pal.coral,
                "lw": 1.0,
                "shrinkA": 0,
                "shrinkB": 0,
                "mutation_scale": 10,
            },
            zorder=6,
        )
        _label(top, (h * 0.62, yd - 0.12), "$h = \\tan\\theta_{cap}$", pal.coral, ha="left")
        _label(
            top, (-0.62, -0.74), "seam parallel\n$\\theta = \\theta_{cap}$", pal.muted, ha="right"
        )
        top.set_xlim(-1.05, 1.05)
        top.set_ylim(-0.95, 0.95)

        # Right: side view through two opposite corners, radial projection.
        th = np.linspace(-0.62 * np.pi, 0.62 * np.pi, 300)
        side.plot(np.sin(th), np.cos(th), color=pal.ink, lw=LW, zorder=1)
        side.fill(np.sin(th), np.cos(th), color=pal.tint(pal.slate, 0.86), zorder=0, lw=0)
        side.plot([-1.25 * h, 1.25 * h], [1, 1], color=pal.muted, lw=THIN, zorder=1)
        u = np.linspace(-h, h, nc)
        for ui in u:
            d = np.array([ui, 1.0])
            hit = d / np.linalg.norm(d)
            side.plot([0, ui], [0, 1], color=pal.muted, lw=THIN, zorder=1)
            _dot(side, (ui, 1.0), pal.ink, size=3.8, face=pal.background)
            _dot(side, hit, pal.ocean, size=4.2)
        side.plot([-h, h], [1, 1], color=pal.ink, lw=LW, zorder=2)
        _dot(side, (0, 0), pal.coral, size=6.5)
        _label(side, (0.06, -0.12), "anchor", pal.coral, ha="left")
        # theta_cap angle at the anchor.
        side.add_patch(
            Arc(
                (0, 0),
                0.62,
                0.62,
                theta1=90 - np.rad2deg(theta),
                theta2=90,
                color=pal.coral,
                lw=1.1,
            )
        )
        side.plot([0, 0], [0, 1], color=pal.faint, lw=THIN, zorder=0)
        _label(
            side,
            (0.13, 0.43),
            "$\\theta_{cap}$",
            pal.coral,
            bbox={"fc": pal.tint(pal.slate, 0.86), "ec": "none", "pad": 0.6},
        )
        _label(side, (-h - 0.07, 1.0), "tangent plane", pal.ink, ha="right")
        _label(side, (-0.98, 0.42), "body", pal.muted, ha="right")
        _label(side, (0.86, 0.80), "cap k = 0", pal.ocean, ha="left")
        side.set_xlim(-1.35, 1.35)
        side.set_ylim(-0.75, 1.25)
        fig.subplots_adjust(wspace=0.08)
        save(fig, "gnomonic-cap", mode, "svg")


# ---------------------------------------------------------------------------
# junction-3edge
# ---------------------------------------------------------------------------


def _graded(n: int, first: float, growth: float) -> np.ndarray:
    steps = first * growth ** np.arange(n)
    return np.r_[0.0, np.cumsum(steps)]


def junction_3edge(mode: str) -> None:
    with style(mode) as pal:
        fig, (ax,) = _canvas(5.6, 4.6)
        s = _graded(5, 0.3, 1.22)  # node positions along each edge from the junction
        dirs = {
            "seam_q": np.deg2rad(210.0),
            "seam_q1": np.deg2rad(330.0),
            "meridian": np.deg2rad(90.0),
        }
        unit = {k: np.array([np.cos(v), np.sin(v)]) for k, v in dirs.items()}
        sectors = [
            ("cap", "seam_q", "seam_q1", pal.sand),
            ("sub$_q$", "meridian", "seam_q", pal.ocean),
            ("sub$_{q+1}$", "seam_q1", "meridian", pal.sea),
        ]
        for name, ea, eb, colour in sectors:
            u, v = unit[ea], unit[eb]
            far = s[-1]
            quad = np.array([(0, 0), far * u, far * (u + v), far * v])
            ax.add_patch(Polygon(quad, closed=True, fc=pal.tint(colour, 0.8), ec="none", zorder=0))
            for si in s[1:-1]:
                ax.plot(*np.c_[si * u, si * u + far * v], color=pal.tint(pal.ink, 0.55), lw=THIN)
                ax.plot(*np.c_[si * v, si * v + far * u], color=pal.tint(pal.ink, 0.55), lw=THIN)
            ax.plot(*np.c_[far * u, far * (u + v), far * v], color=pal.tint(pal.ink, 0.4), lw=THIN)
            centre = 0.55 * far * (u + v)
            _label(ax, centre, name, pal.ink, size=10.5)
        for key, colour in (("seam_q", pal.ink), ("seam_q1", pal.ink), ("meridian", pal.ink)):
            u = unit[key]
            ax.plot(*np.c_[(0, 0), s[-1] * u], color=colour, lw=LW, zorder=3)
            ax.plot(
                *np.c_[(0, 0), s[1] * u], color=pal.coral, lw=4.0, zorder=4, solid_capstyle="butt"
            )
            for si in s[1:]:
                _dot(ax, si * u, pal.ink, size=3.4, z=5)
            n = np.array([-u[1], u[0]])
            _label(ax, 0.5 * s[1] * u + 0.17 * n, "$\\Delta s$", pal.coral, size=10)
        _dot(ax, (0, 0), pal.ink, size=7.5, z=6)
        # Edge names along each edge, outside the patches.
        _label(
            ax, s[-1] * unit["meridian"] + (0, 0.14), "meridian (sub$_q$ | sub$_{q+1}$)", pal.ink
        )
        _label(
            ax,
            s[-1] * unit["seam_q"] + (-0.08, -0.14),
            "cap edge (cap | sub$_q$)",
            pal.ink,
            ha="center",
        )
        _label(
            ax,
            s[-1] * unit["seam_q1"] + (0.08, -0.14),
            "cap edge (cap | sub$_{q+1}$)",
            pal.ink,
            ha="center",
        )
        far = s[-1]
        ax.set_xlim(-0.87 * far - 0.45, 0.87 * far + 0.45)
        ax.set_ylim(-0.5 * far - 0.3, far + 0.3)
        save(fig, "junction-3edge", mode, "svg")


# ---------------------------------------------------------------------------
# per-j-pinning
# ---------------------------------------------------------------------------


def per_j_pinning(mode: str) -> None:
    with style(mode) as pal:
        fig, axes = plt.subplots(1, 2, figsize=(7.4, 2.9), sharey=True)
        ni = 17
        theta_cap = 0.05 * np.pi
        start, stop = theta_cap, np.pi - theta_cap
        natural = np.diff(build_i_coords(ni, start, stop, spacing="uniform"))
        nat = natural[0]
        idx = np.arange(1, ni)
        cases = [
            (axes[0], 0.7 * nat, True, "cap step $0.7\\times$ natural: pinned"),
            (axes[1], 0.15 * nat, False, "cap step $0.15\\times$ natural: fallback"),
        ]
        for ax, cap, compatible, title in cases:
            ax.axhspan(0.5 * nat, 2.0 * nat, color=pal.tint(pal.sea, 0.85), lw=0, zorder=0)
            ax.axhline(nat, color=pal.muted, lw=THIN, ls=(0, (1.5, 2)), zorder=1)
            ax.axhline(cap, color=pal.coral, lw=1.2, ls=(0, (5, 3)), zorder=1)
            pinned = np.diff(pinned_distribution(ni, start, stop, cap, cap, spacing="uniform"))
            if compatible:
                ax.plot(idx, pinned, "o-", color=pal.ocean, ms=4, lw=LW, zorder=3)
            else:
                ax.plot(
                    idx,
                    pinned,
                    marker="o",
                    color=pal.muted,
                    ms=3,
                    lw=THIN,
                    zorder=2,
                    mfc="none",
                    ls=(0, (3, 2)),
                )
                ax.plot(idx, natural, "o-", color=pal.ocean, ms=4, lw=LW, zorder=3)
                _label(
                    ax,
                    (8.5, pinned[1] + 0.035),
                    "pinned would jump here: rejected",
                    pal.muted,
                    size=SMALL,
                )
            ax.set_title(title, fontsize=9.8, fontweight="normal", loc="left", pad=6)
            ax.set_xlabel("meridional interval $i$")
            ax.set_xlim(0.3, ni - 0.3)
            ax.grid(False)
            ax.set_xticks([1, 4, 8, 12, 16])
        axes[0].set_ylabel("$\\Delta\\theta$ (rad)")
        axes[0].set_ylim(0, 2.15 * nat)
        _label(
            axes[0],
            (8.5, 1.75 * nat),
            "compatibility window $[0.5, 2]\\times$ natural",
            pal.sea,
            size=SMALL,
        )
        _label(axes[0], (8.5, nat - 0.022), "natural step", pal.muted, size=SMALL)
        _label(axes[0], (2.6, 0.7 * nat - 0.022), "cap step", pal.coral, size=SMALL, ha="left")
        fig.subplots_adjust(wspace=0.08)
        save(fig, "per-j-pinning", mode, "svg")


# ---------------------------------------------------------------------------
# seam-aware-normals
# ---------------------------------------------------------------------------


def seam_aware_normals(mode: str) -> None:
    with style(mode) as pal:
        fig, axes = _canvas(7.4, 3.0, 2)
        # Wall row along j on a curved surface: sub_q fine, sub_{q+1} coarse.
        ang_q = np.deg2rad(np.array([-40, -32, -24, -16, -8, 0]))
        ang_q1 = np.deg2rad(np.array([0, 15, 30, 45]))
        r = 2.4

        def pt(a):
            return np.c_[r * np.sin(a), r * np.cos(a)] - (0, r)

        pq, pq1 = pt(ang_q), pt(ang_q1)
        seam = pq[-1]
        for ax, centred in zip(axes, (False, True), strict=True):
            arc = pt(np.linspace(np.deg2rad(-44), np.deg2rad(49), 200))
            ax.plot(*arc.T, color=pal.faint, lw=5.0, zorder=0, solid_capstyle="round")
            ax.plot(*pq.T, "-", color=pal.ocean, lw=LW, zorder=2)
            ax.plot(*pq1.T, "-", color=pal.sea, lw=LW, zorder=2)
            ax.plot(*pq[:-1].T, "o", ms=4.2, color=pal.ocean, zorder=3)
            ax.plot(*pq1[1:].T, "o", ms=4.2, color=pal.sea, zorder=3)
            ax.plot(*seam, "s", ms=6.5, color=pal.ink, zorder=5)
            _label(ax, pq[2] + (-0.1, -0.32), "sub$_q$", pal.ocean, size=10)
            _label(ax, pq1[2] + (0.1, -0.32), "sub$_{q+1}$", pal.sea, size=10)
            length = 1.15
            if not centred:
                for a, b, colour in ((pq[-2], seam, pal.ocean), (seam, pq1[1], pal.sea)):
                    t = (b - a) / np.linalg.norm(b - a)
                    n = np.array([-t[1], t[0]])
                    _arrow(ax, seam, seam + length * n, colour, lw=LW)
                ax.text(
                    0,
                    -1.25,
                    "one-sided tangents: two normals at one node",
                    color=pal.ink,
                    fontsize=9.5,
                    ha="center",
                )
            else:
                for g, colour in ((pq1[1], pal.sea), (pq[-2], pal.ocean)):
                    _dot(ax, g, colour, size=11, face="none", z=4)
                ax.plot(*np.c_[pq[-2], pq1[1]], color=pal.ink, lw=THIN, ls=(0, (3, 2)), zorder=1)
                t = (pq1[1] - pq[-2]) / np.linalg.norm(pq1[1] - pq[-2])
                n = np.array([-t[1], t[0]])
                _arrow(ax, seam, seam + length * n, pal.ink, lw=LW)
                _label(ax, pq1[1] + (0.08, 0.3), "ghost", pal.muted, ha="left")
                _label(ax, pq[-2] + (-0.08, 0.3), "ghost", pal.muted, ha="right")
                ax.text(
                    0,
                    -1.25,
                    "centred against ghost rows: one normal",
                    color=pal.ink,
                    fontsize=9.5,
                    ha="center",
                )
            _label(ax, seam + _v(0.0, -0.3), "seam", pal.ink, size=SMALL)
            if not centred:
                _label(ax, seam + _v(-0.35, 0.95), "$n_q$", pal.ocean, size=10)
                _label(ax, seam + _v(0.4, 0.95), "$n_{q+1}$", pal.sea, size=10)
            ax.set_xlim(-2.0, 2.05)
            ax.set_ylim(-1.4, 1.3)
        fig.subplots_adjust(wspace=0.04)
        save(fig, "seam-aware-normals", mode, "svg")


# ---------------------------------------------------------------------------
# block-model
# ---------------------------------------------------------------------------

_BOX = np.array([2.2, 1.4, 1.1])  # drawn extents along i, j, k


def _vertex_xyz(index: int, origin=(0.0, 0.0, 0.0), size=_BOX) -> np.ndarray:
    """Corner *index* = 4 i + 2 j + k of a box (HexBlock.vertices order)."""
    i, j, k = (index >> 2) & 1, (index >> 1) & 1, index & 1
    return np.asarray(origin) + np.array([i, j, k]) * size


def block_model(mode: str) -> None:
    with style(mode) as pal:
        fig, (ax,) = _canvas(6.4, 4.4)
        axis_colour = {"i": pal.ocean, "j": pal.plum, "k": pal.sea}
        v3 = np.array([_vertex_xyz(n) for n in range(8)])
        v = _oblique(v3)
        hidden_vertex = 2  # (0, 1, 0): behind the front faces
        faces = {
            "i_lo": [0, 1, 3, 2],
            "i_hi": [4, 5, 7, 6],
            "j_lo": [0, 1, 5, 4],
            "j_hi": [2, 3, 7, 6],
            "k_lo": [0, 2, 6, 4],
            "k_hi": [1, 3, 7, 5],
        }
        # WALL face first, under everything.
        ax.add_patch(
            Polygon(
                v[faces["k_lo"]], closed=True, fc=pal.tint(pal.coral, 0.75), ec="none", zorder=0
            )
        )
        for name in ("j_lo", "i_hi", "k_hi"):
            ax.add_patch(
                Polygon(
                    v[faces[name]],
                    closed=True,
                    fc=pal.tint(pal.slate, 0.9),
                    ec="none",
                    alpha=0.55,
                    zorder=1,
                )
            )
        # 12 edges, coloured by axis, hidden ones dashed.
        for a in range(8):
            for bit, axis in ((4, "i"), (2, "j"), (1, "k")):
                if a & bit:
                    continue
                b = a | bit
                hidden = hidden_vertex in (a, b)
                ax.plot(
                    *v[[a, b]].T,
                    color=axis_colour[axis],
                    lw=LW,
                    ls=(0, (3.5, 2.5)) if hidden else "-",
                    zorder=2,
                    solid_capstyle="round",
                )
        # Vertices: index and (i, j, k) corner.
        mid = v.mean(axis=0)
        for n in range(8):
            ax.plot(*v[n], "o", ms=13, mfc=pal.background, mec=pal.ink, mew=1.0, zorder=6)
            _label(ax, v[n], str(n), pal.ink, size=8.6, zorder=7)
            ijk = ",".join(str((n >> s) & 1) for s in (2, 1, 0))
            out = (v[n] - mid) / np.linalg.norm(v[n] - mid)
            at = v[n] + {2: (0.0, -0.27), 5: (0.4, 0.06)}.get(n, 0.34 * out)
            _label(ax, at, f"({ijk})", pal.muted, size=7.8)

        # Face names: visible faces inside, hidden ones with a leader.
        def centre(name):
            return v[faces[name]].mean(axis=0)

        _label(ax, centre("j_lo"), "j_lo", pal.ink, size=9.5)
        _label(ax, centre("i_hi"), "i_hi", pal.ink, size=9.5)
        _label(ax, centre("k_hi"), "k_hi", pal.ink, size=9.5)
        for name, at, ha in (
            ("i_lo", (-0.75, 1.15), "right"),
            ("j_hi", (3.55, 1.95), "left"),
        ):
            c = centre(name)
            ax.plot(*c, "o", ms=2.6, color=pal.muted, zorder=3)
            ax.plot(*np.c_[c, at], color=pal.muted, lw=THIN, zorder=3)
            _label(
                ax,
                (at[0] + (-0.05 if ha == "right" else 0.05), at[1]),
                name,
                pal.muted,
                size=9.5,
                ha=ha,
            )
        c = centre("k_lo")
        ax.plot(*c, "o", ms=2.6, color=pal.coral, zorder=3)
        at = (3.55, -0.32)
        ax.plot(*np.c_[c, (c[0] + 0.5, at[1]), at], color=pal.coral, lw=THIN, zorder=3)
        _label(ax, (at[0] + 0.05, at[1]), "k_lo  WALL", pal.coral, size=9.5, ha="left")
        # Axis triad, offset from vertex 0.
        o = np.array([-1.45, -0.55, 0.0])
        for e, axis in (
            (np.array([0.55, 0, 0]), "i"),
            (np.array([0, 0.8, 0]), "j"),
            (np.array([0, 0, 0.55]), "k"),
        ):
            p0, p1 = _oblique(o), _oblique(o + e)
            _arrow(ax, p0, p1, axis_colour[axis], lw=1.3, head=8)
            d = (p1 - p0) / np.linalg.norm(p1 - p0)
            _label(ax, p1 + 0.13 * d, axis, axis_colour[axis], size=10.5, style="italic")
        ax.set_xlim(-1.8, 4.4)
        ax.set_ylim(-0.75, 2.15)
        save(fig, "block-model", mode, "svg")


# ---------------------------------------------------------------------------
# axis-map
# ---------------------------------------------------------------------------


def axis_map(mode: str) -> None:
    """Block A's i_hi shared with block B's j_lo: A's j runs against B's i."""
    # The AxisMaps drawn below, checked against the code.
    amap = AxisMap.from_pairs(("i", True), ("k", False))
    amap.validate_against_partner("i_hi", "j_lo")
    inv = amap.inverse("i_hi", "j_lo")
    assert inv.along == (("j", True), ("k", False))

    with style(mode) as pal:
        fig, (ax,) = _canvas(7.2, 3.6)
        axis_colour = {"i": pal.ocean, "j": pal.plum, "k": pal.sea}
        size = np.array([1.3, 1.2, 1.0])
        gap = 1.15
        # Block A: i -> +x, j -> +y, k -> +z.  Block B: i -> -y, j -> +x, k -> +z.
        a_corner = np.array(
            [[x, y, z] for x in (0, size[0]) for y in (0, size[1]) for z in (0, size[2])]
        )
        b0 = np.array([size[0] + gap, 0, 0])
        b_size = np.array([1.2, size[1], size[2]])
        b_corner = np.array(
            [
                b0 + _v(x, y, z)
                for x in (0, b_size[0])
                for y in (0, b_size[1])
                for z in (0, b_size[2])
            ]
        )
        shared = pal.tint(pal.sand, 0.72)
        for corners in (a_corner, b_corner):
            p = _oblique(corners)
            # corner order here: index = 4 x + 2 y + z (x, y, z in {0, 1})
            for a in range(8):
                for bit in (4, 2, 1):
                    if a & bit:
                        continue
                    b = a | bit
                    hidden = 2 in (a, b)
                    ax.plot(
                        *p[[a, b]].T,
                        color=pal.ink if not hidden else pal.muted,
                        lw=1.1 if not hidden else THIN,
                        ls="-" if not hidden else (0, (3, 2.5)),
                        zorder=2,
                    )
        pa, pb = _oblique(a_corner), _oblique(b_corner)
        face_a = [4, 5, 7, 6]  # x = max of A: its i_hi
        face_b = [0, 1, 3, 2]  # x = min of B: its j_lo
        ax.add_patch(Polygon(pa[face_a], closed=True, fc=shared, ec="none", zorder=1))
        ax.add_patch(Polygon(pb[face_b], closed=True, fc=shared, ec="none", zorder=1, alpha=0.8))
        # Along-axes on each face, drawn from the index origin of the face.
        inset = 0.16
        # A.i_hi: origin at (j=0, k=0) = (x=1.3, y=0, z=0); j -> +y, k -> +z.
        oa = np.array([size[0], inset, inset])
        _arrow(ax, _oblique(oa), _oblique(oa + _v(0, size[1] - 2 * inset, 0)), axis_colour["j"])
        ax.plot(*_oblique(oa), "o", ms=3.5, color=pal.ink, zorder=5)
        _arrow(ax, _oblique(oa), _oblique(oa + _v(0, 0, size[2] - 2 * inset)), axis_colour["k"])
        # B.j_lo: origin at (i=0, k=0) = (x=b0, y=size_y, z=0); i -> -y, k -> +z.
        ob = b0 + np.array([0, size[1] - inset, inset])
        _arrow(ax, _oblique(ob), _oblique(ob - _v(0, size[1] - 2 * inset, 0)), axis_colour["i"])
        ax.plot(*_oblique(ob), "o", ms=3.5, color=pal.ink, zorder=5)
        _arrow(ax, _oblique(ob), _oblique(ob + _v(0, 0, size[2] - 2 * inset)), axis_colour["k"])

        # Axis names next to the arrows.
        def name(p, text, axis, off):
            _label(ax, _oblique(p) + off, text, axis_colour[axis], size=10.5, style="italic")

        name(oa + _v(0, size[1] - 2 * inset, 0), "j", "j", (0.12, 0.02))
        name(oa + _v(0, 0, size[2] - 2 * inset), "k", "k", (0.0, 0.13))
        name(ob - _v(0, size[1] - 2 * inset, 0), "i", "i", (-0.1, -0.08))
        name(ob + _v(0, 0, size[2] - 2 * inset), "k", "k", (0.0, 0.13))
        # Block and face names.
        _label(ax, _oblique([size[0] / 2, 0, size[2] / 2]), "A", pal.ink, size=13)
        _label(ax, _oblique(b0 + _v(b_size[0] / 2, 0, size[2] / 2)), "B", pal.ink, size=13)
        _label(
            ax, _oblique([size[0], 0, 0]) + _v(0.05, -0.18), "A.i_hi", pal.ink, size=9.5, ha="left"
        )
        _label(ax, _oblique(b0) + _v(-0.05, -0.18), "B.j_lo", pal.ink, size=9.5, ha="right")
        # The mapping between the two faces.
        mid = (_oblique([size[0], size[1], size[2]]) + _oblique(b0 + _v(0, 0, size[2]))) / 2
        _label(
            ax,
            mid + _v(0.0, 0.62),
            "A.j = B.i, reversed\nA.k = B.k",
            pal.ink,
            size=SMALL,
            linespacing=1.4,
        )
        # Block-local triads under each block.
        for origin, axes_dirs in (
            (np.array([0.0, -0.05, -0.55]), {"i": [1, 0, 0], "j": [0, 1, 0], "k": [0, 0, 1]}),
            (b0 + np.array([0.0, 0.85, -0.55]), {"i": [0, -1, 0], "j": [1, 0, 0], "k": [0, 0, 1]}),
        ):
            for name, d in axes_dirs.items():
                d = np.array(d, dtype=float) * (0.6 if name != "k" else 0.0)
                if name == "k":
                    continue
                p0, p1 = _oblique(origin), _oblique(origin + d)
                _arrow(ax, p0, p1, axis_colour[name], lw=1.1, head=7)
                u = (p1 - p0) / np.linalg.norm(p1 - p0)
                _label(ax, p1 + 0.12 * u, name, axis_colour[name], size=9.5, style="italic")
        # The records, as AxisMap prints them.
        mono = {"family": "monospace", "size": 8.6}
        ax.text(
            -0.2,
            -1.05,
            "A.i_hi  along (j, k)  ->  (('i', True), ('k', False))",
            color=pal.ink,
            fontdict=mono,
            ha="left",
        )
        ax.text(
            -0.2,
            -1.32,
            "B.j_lo  along (i, k)  ->  (('j', True), ('k', False))",
            color=pal.ink,
            fontdict=mono,
            ha="left",
        )
        ax.set_xlim(-0.35, 4.6)
        ax.set_ylim(-1.45, 2.05)
        save(fig, "axis-map", mode, "svg")


FIGURES = {
    "projection-raycast": projection_raycast,
    "gnomonic-cap": gnomonic_cap,
    "junction-3edge": junction_3edge,
    "per-j-pinning": per_j_pinning,
    "seam-aware-normals": seam_aware_normals,
    "block-model": block_model,
    "axis-map": axis_map,
}
