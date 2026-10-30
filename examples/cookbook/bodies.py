"""Make analytic test bodies, NACA sections and extruded wings.

Writes sphere.stl, capsule.stl, torus.stl, naca0012.dat, naca2412_blunt.dat,
wing_straight.stl and wing_swept.stl to the current directory.
"""

from __future__ import annotations

import numpy as np

from shore import profiles, shapes

# #region analytic
# Analytic bodies are trimesh.Trimesh objects; export() writes the STL.
shapes.sphere(radius=1.0, subdivisions=3).export("sphere.stl")
shapes.capsule(radius=0.5, height=2.0, sections=24).export("capsule.stl")
torus = shapes.torus(r_major=1.0, r_minor=0.25)
torus.export("torus.stl")
print(f"torus: {len(torus.faces)} triangles, watertight={torus.is_watertight}")
# #endregion analytic

# #region naca
sharp = profiles.naca("0012", n_pts=80)  # (159, 2), upper TE -> LE -> lower TE
blunt = profiles.naca("2412", n_pts=80, closed_te=False)
profiles.write_section("naca0012.dat", sharp, name="NACA 0012")
profiles.write_section("naca2412_blunt.dat", blunt, name="NACA 2412 blunt TE")
gap = np.linalg.norm(blunt[0] - blunt[-1])
print(f"NACA 0012: {sharp.shape[0]} points; NACA 2412 blunt TE gap {gap:.2e} chord")
# #endregion naca

# #region wing
section = profiles.read_section("naca0012.dat")
# Straight wing with open ends: the C-grid topology slices it per span station.
shapes.extrude_section(section, span=1.0, n_stations=8).export("wing_straight.stl")
# Tapered, swept, twisted wing with watertight end caps.
swept = shapes.extrude_section(
    section,
    span=2.0,
    n_stations=12,
    taper=0.6,
    sweep_deg=25.0,
    twist_deg=-3.0,
    span_axis="y",
    end_caps=True,
)
swept.export("wing_swept.stl")
print(f"swept wing: extents {np.round(swept.extents, 3)}, watertight={swept.is_watertight}")
# #endregion wing
