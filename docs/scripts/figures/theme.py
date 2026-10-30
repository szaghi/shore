"""The SHORE figure theme: one palette for matplotlib and PyVista, light and dark.

Every figure is drawn twice, once per site theme, and written to
``docs/public/figures/{light,dark}/<name>.<ext>`` (SVG charts and
schematics, WebP renders); the ``<Figure>`` component
of the docs theme shows the one matching the reader's theme. The colours
mirror ``docs/.vitepress/theme/custom.css``.
"""

from __future__ import annotations

from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap, to_rgb

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "docs" / "public" / "figures"
MODES = ("light", "dark")


@dataclass(frozen=True)
class Palette:
    """Colours of one site theme."""

    mode: str
    ink: str  # text, outlines
    muted: str  # secondary text, guides
    faint: str  # grid lines, light fills
    panel: str  # box fills in schematics
    background: str  # the page (figures themselves are transparent)
    edge: str  # grid lines drawn over 3-D surfaces
    blocks: tuple[str, ...]  # categorical: one colour per block or role

    @property
    def ocean(self) -> str:
        return self.blocks[0]

    @property
    def sea(self) -> str:
        return self.blocks[1]

    @property
    def sand(self) -> str:
        return self.blocks[2]

    @property
    def sunset(self) -> str:
        return self.blocks[3]

    @property
    def coral(self) -> str:
        return self.blocks[4]

    @property
    def plum(self) -> str:
        return self.blocks[5]

    @property
    def sky(self) -> str:
        return self.blocks[6]

    @property
    def slate(self) -> str:
        return self.blocks[7]

    def block(self, k: int) -> str:
        return self.blocks[k % len(self.blocks)]

    def tint(self, colour: str, amount: float) -> tuple[float, float, float]:
        """*colour* blended towards the background by *amount* (0..1)."""
        c, b = to_rgb(colour), to_rgb(self.background)
        return tuple(x + (y - x) * amount for x, y in zip(c, b, strict=True))

    @property
    def quality(self) -> LinearSegmentedColormap:
        """Sequential map for quality metrics: bad (coral) to good (ocean)."""
        return LinearSegmentedColormap.from_list(
            f"shore-quality-{self.mode}", [self.coral, self.sunset, self.sand, self.sea, self.ocean]
        )

    @property
    def depth(self) -> LinearSegmentedColormap:
        """Sequential map for layer index or distance: near wall to far field."""
        return LinearSegmentedColormap.from_list(
            f"shore-depth-{self.mode}", [self.ocean, self.sea, self.sand, self.sunset]
        )


LIGHT = Palette(
    mode="light",
    ink="#1d2733",
    muted="#5b6b7b",
    faint="#d5dde5",
    panel="#f3f6f9",
    background="#ffffff",
    edge="#24313d",
    blocks=(
        "#1f6f8b",  # ocean
        "#2a9d8f",  # sea
        "#e9c46a",  # sand
        "#f4a261",  # sunset
        "#e76f51",  # coral
        "#9d6b9b",  # plum
        "#6fb3d6",  # sky
        "#56697a",  # slate
    ),
)

DARK = Palette(
    mode="dark",
    ink="#e3e8ee",
    muted="#9aa8b6",
    faint="#33404d",
    panel="#232a31",
    background="#1b1b1f",
    edge="#0f1418",
    blocks=(
        "#4aa3c2",
        "#43c2b3",
        "#f0d283",
        "#f7b77f",
        "#ef8b72",
        "#c494c2",
        "#93cbe8",
        "#8fa3b5",
    ),
)

PALETTES = {"light": LIGHT, "dark": DARK}


def _font() -> str:
    names = {f.name for f in font_manager.fontManager.ttflist}
    return next((f for f in ("Inter", "Source Sans 3", "DejaVu Sans") if f in names), "sans-serif")


@contextmanager
def style(mode: str) -> Iterator[Palette]:
    """matplotlib rcParams for *mode*; yields its palette."""
    pal = PALETTES[mode]
    rc = {
        "font.family": _font(),
        "font.size": 10.5,
        "mathtext.fontset": "dejavusans",
        "text.color": pal.ink,
        "axes.labelcolor": pal.ink,
        "axes.edgecolor": pal.muted,
        "axes.facecolor": "none",
        "axes.titleweight": "semibold",
        "axes.titlesize": 11.5,
        "axes.titlecolor": pal.ink,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "axes.prop_cycle": mpl.cycler(color=list(pal.blocks)),
        "grid.color": pal.faint,
        "grid.linewidth": 0.6,
        "xtick.color": pal.muted,
        "ytick.color": pal.muted,
        "xtick.labelcolor": pal.ink,
        "ytick.labelcolor": pal.ink,
        "legend.frameon": False,
        "figure.facecolor": "none",
        "savefig.facecolor": "none",
        "savefig.transparent": True,
        "lines.linewidth": 1.8,
        "patch.edgecolor": pal.ink,
    }
    with mpl.rc_context(rc):
        yield pal


def save(fig: plt.Figure, name: str, mode: str, ext: str = "svg") -> Path:
    """Write *fig* to ``figures/<mode>/<name>.<ext>`` (transparent background)."""
    out = OUT / mode / f"{name}.{ext}"
    out.parent.mkdir(parents=True, exist_ok=True)
    kw = {"metadata": {"Date": None}} if ext == "svg" else {"dpi": 200}
    fig.savefig(out, bbox_inches="tight", pad_inches=0.08, transparent=True, **kw)
    plt.close(fig)
    if ext == "svg":  # matplotlib leaves trailing spaces in path data
        out.write_text("\n".join(line.rstrip() for line in out.read_text().splitlines()) + "\n")
    return out


def plotter(mode: str, size: tuple[int, int] = (1200, 900)):
    """An off-screen PyVista plotter lit for *mode*."""
    import pyvista as pv

    pv.OFF_SCREEN = True
    pl = pv.Plotter(off_screen=True, window_size=size, lighting="none")
    pl.set_background(PALETTES[mode].background)
    key = pv.Light(position=(3, -4, 6), focal_point=(0, 0, 0), intensity=0.75)
    fill = pv.Light(position=(-5, 2, 1), focal_point=(0, 0, 0), intensity=0.35)
    head = pv.Light(light_type="headlight", intensity=0.25)
    for light in (key, fill, head):
        pl.add_light(light)
    pl.enable_anti_aliasing("ssaa")
    return pl


def shot(pl, name: str, mode: str) -> Path:
    """Screenshot *pl* to ``figures/<mode>/<name>.webp`` with a transparent background."""
    out = OUT / mode / f"{name}.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    pl.screenshot(str(out), transparent_background=True)
    pl.close()
    _trim(out)
    return webp(out)


def webp(png: Path) -> Path:
    """Replace *png* by a WebP (quality 90, alpha kept): renders 4x smaller."""
    from PIL import Image

    out = png.with_suffix(".webp")
    Image.open(png).save(out, "WEBP", quality=90, method=6)
    png.unlink()
    return out


def _trim(path: Path, pad: int = 12) -> None:
    """Crop the transparent margin of a PNG."""
    import numpy as np
    from PIL import Image

    img = Image.open(path)
    alpha = np.asarray(img)[..., 3]
    ys, xs = np.nonzero(alpha)
    if len(xs):
        box = (
            max(xs.min() - pad, 0),
            max(ys.min() - pad, 0),
            min(xs.max() + pad, img.width),
            min(ys.max() + pad, img.height),
        )
        img.crop(box).save(path)
