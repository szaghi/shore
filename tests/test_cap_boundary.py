"""Phase 4: cap/equator boundary conformance after the bit-exact override.

``build_cap_k0`` finishes with ``_override_boundary``, which copies the 4 sub-
block pole rows into the cap's 4 edges.  After override:

- The 4 cap corners coincide (bit-exact) with the 4 sub-block seam corners.
- Each of the 4 cap edges coincides (bit-exact) with the corresponding sub-
  block pole row.

These properties are the prerequisite for SHARED/DIRICHLET face BCs in the
mesh layer: without them the cap/equator seam would crack under marching.
"""

from __future__ import annotations

import numpy as np

from shore.mesh import CapsTopology
from shore.surface.cap import build_cap_k0
from tests.conftest import sphere_layer

# ---------------------------------------------------------------------------
# Standalone build_cap_k0 override tests
# ---------------------------------------------------------------------------


def _pole_rows(surface: np.ndarray, nj: int, pole: str) -> list[np.ndarray]:
    sector = nj // 4
    ni = surface.shape[0]
    i_pick = ni - 1 if pole == "north" else 0
    rows = []
    for q in range(4):
        j0 = q * sector
        j_idx = [*range(j0, j0 + sector), (j0 + sector) % nj]
        rows.append(surface[i_pick, j_idx, :])
    return rows


class TestBuildCapK0Boundary:
    """Boundary conformance direct from ``build_cap_k0``."""

    def _build(self, pole: str, ni: int = 33, nj: int = 32, theta_cap_deg: float = 15.0):
        theta_cap = np.deg2rad(theta_cap_deg)
        surface = sphere_layer(ni, nj, theta_cap=theta_cap)
        rows = _pole_rows(surface, nj, pole)
        cap = build_cap_k0(rows, theta_cap=theta_cap, ni=ni, pole=pole)
        return cap, rows

    def test_north_edges_match_pole_rows(self):
        cap, rows = self._build("north")
        # Mapping mirrors shore.topology.sphere_cap._overwrite_cap_boundaries
        np.testing.assert_array_equal(cap[0, :, :], rows[0])
        np.testing.assert_array_equal(cap[:, -1, :], rows[1])
        np.testing.assert_array_equal(cap[-1, :, :], rows[2][::-1])
        np.testing.assert_array_equal(cap[:, 0, :], rows[3][::-1])

    def test_south_edges_match_pole_rows(self):
        cap, rows = self._build("south")
        np.testing.assert_array_equal(cap[:, 0, :], rows[0])
        np.testing.assert_array_equal(cap[-1, :, :], rows[1])
        np.testing.assert_array_equal(cap[:, -1, :], rows[2][::-1])
        np.testing.assert_array_equal(cap[0, :, :], rows[3][::-1])

    def test_north_corners_are_shared_junctions(self):
        """The 4 cap corners must be the 4 sub-block pole-row endpoints."""
        cap, rows = self._build("north")
        # cap[0, 0] = sub0[0] ; cap[0, -1] = sub0[-1] = sub1[0]
        np.testing.assert_array_equal(cap[0, 0, :], rows[0][0])
        np.testing.assert_array_equal(cap[0, -1, :], rows[0][-1])
        np.testing.assert_array_equal(cap[0, -1, :], rows[1][0])
        np.testing.assert_array_equal(cap[-1, -1, :], rows[1][-1])
        np.testing.assert_array_equal(cap[-1, -1, :], rows[2][0])
        np.testing.assert_array_equal(cap[-1, 0, :], rows[2][-1])
        np.testing.assert_array_equal(cap[-1, 0, :], rows[3][0])
        np.testing.assert_array_equal(cap[0, 0, :], rows[3][-1])


# ---------------------------------------------------------------------------
# CapsTopology end-to-end: k=0 seam conformity
# ---------------------------------------------------------------------------


class TestCapsTopologyBoundary:
    """The full topology builder must preserve the bit-exact seam."""

    def test_k0_cap_edges_match_sub_pole_rows(self):
        ni, nj = 33, 32
        topo = CapsTopology()
        blocks = topo.build(sphere_layer(ni, nj), nk=3)
        subs, cap_north, cap_south = blocks[:4], blocks[4], blocks[5]

        # North cap (shore: z<0, i=ni-1 of each sub-block)
        for q, sub in enumerate(subs):
            pole_row = sub.nodes[ni - 1, :, 0, :]
            if q == 0:
                np.testing.assert_array_equal(cap_north.nodes[0, :, 0, :], pole_row)
            elif q == 1:
                np.testing.assert_array_equal(cap_north.nodes[:, -1, 0, :], pole_row)
            elif q == 2:
                np.testing.assert_array_equal(cap_north.nodes[-1, :, 0, :], pole_row[::-1])
            else:
                np.testing.assert_array_equal(cap_north.nodes[:, 0, 0, :], pole_row[::-1])

        # South cap (z>0, i=0)
        for q, sub in enumerate(subs):
            pole_row = sub.nodes[0, :, 0, :]
            if q == 0:
                np.testing.assert_array_equal(cap_south.nodes[:, 0, 0, :], pole_row)
            elif q == 1:
                np.testing.assert_array_equal(cap_south.nodes[-1, :, 0, :], pole_row)
            elif q == 2:
                np.testing.assert_array_equal(cap_south.nodes[:, -1, 0, :], pole_row[::-1])
            else:
                np.testing.assert_array_equal(cap_south.nodes[0, :, 0, :], pole_row[::-1])

    def test_sub_block_periodic_seam_k0(self):
        """sub3.j_hi ↔ sub0.j_lo (the periodic wrap) must be node-equal at k=0."""
        topo = CapsTopology()
        blocks = topo.build(sphere_layer(33, 32), nk=3)
        np.testing.assert_array_equal(
            blocks[3].nodes[:, -1, 0, :],
            blocks[0].nodes[:, 0, 0, :],
        )
