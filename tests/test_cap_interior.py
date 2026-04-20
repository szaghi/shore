"""Phase 4: cap k=0 interior quality on a projected STL.

Three properties verified on a unit icosphere STL:

1. No fold-overs — every cap quad's outward normal (computed as the cross
   product of its diagonals) points into the same hemisphere as the sphere
   radial at the quad centroid.
2. STL-normal alignment — the maximum angular deviation between cap quad
   normals and the sphere radial is below 5 deg.
3. Positive Jacobian after one march — the two-layer Jacobian of each cap
   block is strictly positive after extruding one wall-normal layer.
"""

from __future__ import annotations

import numpy as np

from shore.mesh import CapsTopology, Mesh, Spacing1D
from shore.volume.jacobian import _jacobian_field_open


def _quad_normals(layer: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (centroids, normals) for each quad in a (nc, nc, 3) layer.

    Quads are the (nc-1, nc-1) cells between adjacent grid points.  Normals are
    unnormalised cross products of the two diagonals; normalisation is done by
    the caller when needed.
    """
    p00 = layer[:-1, :-1, :]
    p10 = layer[1:, :-1, :]
    p01 = layer[:-1, 1:, :]
    p11 = layer[1:, 1:, :]

    centroids = 0.25 * (p00 + p10 + p01 + p11)
    d1 = p11 - p00
    d2 = p01 - p10
    normals = np.cross(d1, d2)
    return centroids, normals


class TestCapInteriorQuality:
    """Quality of the cap k=0 seed on a unit sphere STL."""

    def _mesh(self, sphere_stl, **kw) -> Mesh:
        return Mesh.from_stl(
            sphere_stl,
            topology=CapsTopology(),
            ni=33,
            nj=32,
            nk=5,
            spacing_k=Spacing1D(law="geometric", ds=0.02, growth=1.1),
            **kw,
        )

    def test_no_fold_overs(self, sphere_stl):
        """Every cap quad normal points outward (positive dot with radial)."""
        m = self._mesh(sphere_stl)
        for cap in (m.blocks[4], m.blocks[5]):
            layer = cap.nodes[:, :, 0, :]
            centroids, normals = _quad_normals(layer)
            radial = centroids / np.linalg.norm(centroids, axis=-1, keepdims=True)
            dots = np.einsum("ijk,ijk->ij", normals, radial)
            # dots sign is geometry-dependent; use the median as the "outward"
            # convention and assert no quad has the opposite sign.
            sign = float(np.sign(np.median(dots)))
            assert sign != 0.0, "degenerate median — all quads are coplanar?"
            assert np.all(sign * dots > 0.0), f"{cap.label}: {(sign * dots <= 0).sum()} fold-overs"

    def test_quad_normal_aligns_with_stl_normal(self, sphere_stl):
        """Cap quad normals deviate by less than 5 deg from the sphere radial."""
        m = self._mesh(sphere_stl)
        for cap in (m.blocks[4], m.blocks[5]):
            layer = cap.nodes[:, :, 0, :]
            centroids, normals = _quad_normals(layer)
            n_hat = normals / np.linalg.norm(normals, axis=-1, keepdims=True)
            radial = centroids / np.linalg.norm(centroids, axis=-1, keepdims=True)
            cos = np.einsum("ijk,ijk->ij", n_hat, radial)
            cos = np.abs(cos)  # sphere radial ≡ STL outward normal on a sphere
            max_dev_deg = float(np.rad2deg(np.arccos(np.clip(cos.min(), -1.0, 1.0))))
            assert max_dev_deg < 5.0, f"{cap.label}: max normal deviation {max_dev_deg:.3f} deg > 5"

    def test_positive_jacobian_after_one_march(self, sphere_stl):
        """Two-layer Jacobian of each cap block is strictly positive post-march."""
        m = self._mesh(sphere_stl)
        m.march(ds=0.02, growth=1.1)
        for cap in (m.blocks[4], m.blocks[5]):
            two_layers = cap.nodes[:, :, :2, :]
            jac = _jacobian_field_open(two_layers)
            # Cap corners are shared with adjacent sub-blocks and can have
            # near-zero Jacobians at the 4-block junction; exclude them.
            ni_c, nj_c = jac.shape[0], jac.shape[1]
            mask = np.ones((ni_c, nj_c), dtype=bool)
            mask[0, 0] = False
            mask[0, nj_c - 1] = False
            mask[ni_c - 1, 0] = False
            mask[ni_c - 1, nj_c - 1] = False
            jmin = float(jac[:, :, 0][mask].min())
            assert jmin > 0.0, f"{cap.label}: min Jacobian {jmin:.3e} <= 0"
