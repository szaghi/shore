"""Tutorial 2: C-grids around a NACA 0012 wing, with tip caps and a CH far field.

Writes naca0012.dat, wing.stl, the C-grid blocks wing_<label>.geo, the
conformal capped C-grid capped_<label>.geo, the chimera tip cap
cap_jlo_cap_<label>.geo and the CH hybrid ch_<label>.geo to the current
directory.
"""

from __future__ import annotations

import numpy as np

from shore import quality
from shore.project import Project


def report(blocks) -> None:
    """One line per block: label, node counts, scaled-Jacobian range."""
    for block in blocks:
        sj = quality.scaled_jacobian(block.nodes)
        print(f"{block.label:24s} {block.nodes.shape[:3]}  sJ [{sj.min():+.3f}, {sj.max():+.3f}]")


# #region body
from shore.profiles import naca, read_section, write_section  # noqa: E402
from shore.shapes import extrude_section  # noqa: E402

write_section("naca0012.dat", naca("0012", n_pts=40), name="NACA 0012")
section = read_section("naca0012.dat")  # (79, 2): upper TE -> LE -> lower TE, sharp TE
wing = extrude_section(section, span=0.5, n_stations=8, span_axis="z", end_caps=True)
wing.export("wing.stl")
print(f"section {section.shape}, wing {len(wing.faces)} triangles")
# #endregion body

# #region cgrid
p = Project()
p.load_geometry("wing.stl", skip_hygiene=True)
p.build_topology(
    "cgrid",
    nk=20,
    ds=5e-4,
    growth=1.15,
    ni_body=120,  # nodes round the section, TE -> LE -> TE
    ni_wake=30,  # nodes along each wake branch, the TE node included
    n_stations=5,  # span stations (j)
    wake_length=15.0,  # chords
    span_axis="z",
)
p.march(smooth_strength=0.1, smooth_iters=1)
p.write_geo("wing")  # wing_upper.geo, wing_body.geo, wing_lower.geo
p.write_adjacency("wing.adjacency.json")
report(p.mesh.blocks)
# #endregion cgrid

# #region capped
p = Project()
p.load_geometry("wing.stl", skip_hygiene=True)
p.build_topology(
    "cgrid",
    nk=20,
    ds=5e-4,
    growth=1.15,
    ni_body=120,
    ni_wake=30,
    n_stations=5,
    wake_length=15.0,
    span_axis="z",
    cap_tips="both",  # 5 cap blocks + 1 wake-base cap per tip; the C-grid becomes 4-block
    cap_tip_inset=0.02,  # spanwise thickness of the cap
)
p.march(smooth_strength=0.1, smooth_iters=1)
p.write_geo("capped")
report(p.mesh.blocks)
# #endregion capped

# #region tipcap
p = Project()
p.build_tip_cap("wing_body.geo", tip="jlo")  # from the sharp C-grid of the cgrid region
p.write_geo("cap_jlo", per_block=True)
p.write_adjacency("cap_jlo.adjacency.json")
report(p.mesh.blocks)
# #endregion tipcap

# #region ch
p = Project()
p.load_geometry("wing.stl", skip_hygiene=True)
p.build_topology(
    "ch",
    nk=12,
    ds=5e-4,
    growth=1.08,
    options={  # the far-field frame and the channel blocks
        "chord_min": -2.0,
        "chord_max": 20.0,
        "normal_min": -2.5,
        "normal_max": 2.5,
        "nk_channel": 8,
    },
    ni_body=80,
    ni_wake=20,
    n_stations=5,
    wake_length=15.0,
)
p.march(smooth_strength=0.0, smooth_iters=0)
p.write_geo("ch")
report(p.mesh.blocks)
lo = np.min([b.nodes.reshape(-1, 3).min(axis=0) for b in p.mesh.blocks], axis=0)
hi = np.max([b.nodes.reshape(-1, 3).max(axis=0) for b in p.mesh.blocks], axis=0)
print(f"far-field box {lo.round(3)} .. {hi.round(3)}")
# #endregion ch
