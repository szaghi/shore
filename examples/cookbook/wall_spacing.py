"""Choose the wall-normal spacing: geometric, tanh/tanh2, and a first cell from a target y+.

Writes sphere.stl and the six sphere_*.geo blocks of a marched cubed sphere to the
current directory.
"""

from __future__ import annotations

import math

import numpy as np

from shore import shapes
from shore.project import Project
from shore.volume.spacing import build_k_steps

# #region geometric
ds, growth, thickness = 1e-3, 1.15, 0.5
# Smallest nk whose geometric layers reach the thickness: ds * (g^(nk-1) - 1) / (g - 1) >= T.
nk = 1 + math.ceil(math.log(1.0 + thickness * (growth - 1.0) / ds) / math.log(growth))
steps = build_k_steps(nk, ds, growth)  # (nk - 1,) step sizes ds * growth**k
print(f"geometric: nk = {nk}, reach {steps.sum():.3f}, last step {steps[-1]:.3e}")
# #endregion geometric


# #region tanh
def first_step(beta: float, law: str, nk: int, thickness: float) -> float:
    return float(build_k_steps(nk, 0.0, spacing=law, beta=beta, total_thickness=thickness)[0])


def beta_for(target: float, law: str, nk: int, thickness: float) -> float:
    """Clustering strength whose first step is *target* (bisection; it falls as beta grows)."""
    lo, hi = 1e-3, 20.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if first_step(mid, law, nk, thickness) > target else (lo, mid)
    return 0.5 * (lo + hi)


nk, thickness, target = 33, 0.5, 1e-3
for law in ("tanh", "tanh2"):
    beta = beta_for(target, law, nk, thickness)
    s = build_k_steps(nk, 0.0, spacing=law, beta=beta, total_thickness=thickness)
    ratio = (s[1:] / s[:-1]).max()
    print(f"{law:5s}: beta = {beta:.3f}, first {s[0]:.2e}, last {s[-1]:.2e}, max ratio {ratio:.3f}")
# #endregion tanh


# #region yplus
def first_cell(y_plus: float, reynolds: float, length: float = 1.0) -> float:
    """Flat-plate ESTIMATE of the wall distance of y+ = *y_plus* at x = *length*.

    Turbulent skin friction Cf = 0.026 Re^(-1/7), u_tau = U sqrt(Cf / 2),
    y = y+ nu / u_tau, i.e. y / L = y+ sqrt(2 / Cf) / Re (Re = U L / nu).
    """
    cf = 0.026 * reynolds ** (-1.0 / 7.0)
    return y_plus * length * math.sqrt(2.0 / cf) / reynolds


for re in (1e5, 1e6, 1e7):
    print(f"Re = {re:.0e}: first cell for y+ = 1 is {first_cell(1.0, re):.2e} L")
# #endregion yplus

# #region march
ds = first_cell(1.0, 1e6, length=2.0)  # sphere diameter 2 as the reference length
shapes.sphere(radius=1.0, subdivisions=4).export("sphere.stl")
p = Project()
p.load_geometry("sphere.stl")
p.build_topology("cubed_sphere", nk=41, ds=ds, growth=1.2, ni=33, nj=64, theta_cap=np.deg2rad(35.0))
p.march()
wall = p.mesh.blocks[0].nodes
first = np.linalg.norm(wall[:, :, 1] - wall[:, :, 0], axis=-1)
print(f"ds = {ds:.2e}: first cell {first.min():.2e}..{first.max():.2e}")
p.write_geo("sphere")
# #endregion march

# #region components
# One k-schedule per component: every block of a component shares it, but a second
# component (here an O-grid on the same sphere) can use another law.
p.use_component("ogrid")
p.build_topology(
    "ogrid", nk=33, k_law="tanh", k_beta=2.582, k_thickness=0.5, ni=33, nj=49,
    theta_cap=np.deg2rad(10.0),
)  # fmt: skip
p.march()
for name, mesh in p.components.items():
    k_line = mesh.blocks[0].nodes[16, 0, :]
    s = np.linalg.norm(np.diff(k_line, axis=0), axis=-1)
    print(f"{name}: {len(mesh.blocks)} blocks, nk = {len(k_line)}, steps {s[1]:.2e}..{s[-1]:.2e}")
# #endregion components
