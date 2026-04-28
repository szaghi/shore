"""Tests for the Xall cc.par writer.

Covers:
- header field formatting
- block section
- BC code mapping (WALL/FREE/SHARED/DIRICHLET/PERIODIC)
- patch index assignment + symmetric family pairing
- cubed-sphere cap-equator orientation derivation from axis_map
  (8 cross-axis seams)
- BlockMeta per-block free_face_bc override
- end-to-end via shore cc-par CLI
"""

from __future__ import annotations

import json

import pytest
from typer.testing import CliRunner

from shore.cli import app
from shore.exceptions import ParameterError
from shore.io.adjacency import ADJACENCY_VERSION, write_adjacency_json
from shore.io.cc_par import (
    BlockMeta,
    _orientation_code_from_axis_map,
    write_cc_par,
)
from shore.mesh.face import AxisMap

# All payloads in this file use the current adjacency version.
_VERSION = ADJACENCY_VERSION


def _face(
    bc: str,
    partner: tuple[str, str] | None = None,
    axis_map: list[list] | None = None,
) -> dict:
    """Build a face entry for a hand-crafted v2 adjacency payload."""
    return {
        "bc": bc,
        "partner": ({"label": partner[0], "face": partner[1]} if partner else None),
        "axis_map": axis_map,
    }


def _identity_map(face_name: str) -> list[list]:
    """Same-axis no-flip axis_map JSON for *face_name* (for tests only)."""
    a, b = {
        "i_lo": ("j", "k"),
        "i_hi": ("j", "k"),
        "j_lo": ("i", "k"),
        "j_hi": ("i", "k"),
        "k_lo": ("i", "j"),
        "k_hi": ("i", "j"),
    }[face_name]
    return [[a, False], [b, False]]


# ---------------------------------------------------------------------------
# Tiny payload fixtures (no real meshing required)
# ---------------------------------------------------------------------------


def _two_block_shared_payload() -> dict:
    """Two adjacent blocks sharing a SHARED j_hi/j_lo seam (same-axis, no flip)."""
    return {
        "version": _VERSION,
        "blocks": [
            {
                "label": "a",
                "shape": [4, 4, 4],
                "faces": {
                    "i_lo": _face("WALL"),
                    "i_hi": _face("FREE"),
                    "j_lo": _face("FREE"),
                    "j_hi": _face("SHARED", ("b", "j_lo"), _identity_map("j_hi")),
                    "k_lo": _face("FREE"),
                    "k_hi": _face("FREE"),
                },
            },
            {
                "label": "b",
                "shape": [4, 4, 4],
                "faces": {
                    "i_lo": _face("WALL"),
                    "i_hi": _face("FREE"),
                    "j_lo": _face("SHARED", ("a", "j_hi"), _identity_map("j_lo")),
                    "j_hi": _face("FREE"),
                    "k_lo": _face("FREE"),
                    "k_hi": _face("FREE"),
                },
            },
        ],
    }


def _parse_patches(text: str) -> list[dict]:
    """Pull out every patch line into a {block, face, bc, family, ext, idx} dict."""
    out: list[dict] = []
    in_patches = False
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if "patches number" in stripped:
            in_patches = True
            continue
        if in_patches and "edges number" in stripped:
            break
        if not in_patches:
            continue
        # Patch line: "block face bc family I0 I1 J0 J1 K0 K1 ! patch_idx"
        if "!" not in stripped:
            continue
        head, _tail = stripped.split("!", 1)
        parts = head.split()
        if len(parts) != 10:
            continue
        idx = int(_tail.strip())
        out.append(
            {
                "block": int(parts[0]),
                "face": int(parts[1]),
                "bc": int(parts[2]),
                "family": int(parts[3]),
                "ext": tuple(int(x) for x in parts[4:10]),
                "idx": idx,
            }
        )
    return out


# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------


class TestHeader:
    def test_header_fields(self, tmp_path):
        payload = _two_block_shared_payload()
        out = write_cc_par(
            payload,
            tmp_path / "cc.par",
            grd_filename="my.grd",
            output_basename="bench.",
            save_ghost_cells=True,
            increase_overlap=False,
            extend_internal_wall=True,
            multigrid_levels=(5, 2),
            boundary_layer_thickness=-2.5e-3,
            numerical_beach=-0.5,
        )
        text = out.read_text()
        lines = text.splitlines()
        assert lines[0] == "'my.grd'"
        assert lines[1] == "'bench.'"
        assert ".true." in lines[2]
        assert ".false." in lines[3]
        assert ".true." in lines[4]
        assert lines[5] == ""
        assert lines[6].split() == ["5", "2"]
        assert "-2.5000E-03" in lines[8]
        assert "-5.0000E-01" in lines[9]


# ---------------------------------------------------------------------------
# Blocks section
# ---------------------------------------------------------------------------


class TestBlocksSection:
    def test_default_metadata(self, tmp_path):
        payload = _two_block_shared_payload()
        out = write_cc_par(payload, tmp_path / "cc.par", grd_filename="x.grd")
        text = out.read_text()
        # Default level=1 group=0 priority=1, comment defaults to label.
        assert "        1         0         1 ! 1 : a" in text
        assert "        1         0         1 ! 2 : b" in text

    def test_blocks_count_line(self, tmp_path):
        payload = _two_block_shared_payload()
        out = write_cc_par(payload, tmp_path / "cc.par", grd_filename="x.grd")
        assert "2 ! blocks number" in out.read_text()

    def test_custom_metadata(self, tmp_path):
        payload = _two_block_shared_payload()
        out = write_cc_par(
            payload,
            tmp_path / "cc.par",
            grd_filename="x.grd",
            block_metadata={
                "a": BlockMeta(level=2, group=1, priority=7, comment="aero block"),
                # b uses defaults
            },
        )
        text = out.read_text()
        assert "        2         1         7 ! 1 : aero block" in text
        assert "        1         0         1 ! 2 : b" in text


# ---------------------------------------------------------------------------
# Patch BC mapping and symmetry
# ---------------------------------------------------------------------------


class TestPatchBCMapping:
    def test_wall_emits_bc_1(self, tmp_path):
        payload = _two_block_shared_payload()
        out = write_cc_par(payload, tmp_path / "cc.par", grd_filename="x.grd")
        patches = _parse_patches(out.read_text())
        # Block 1 (a) face 1 (i_lo) is WALL.
        wall = next(p for p in patches if p["block"] == 1 and p["face"] == 1)
        assert wall["bc"] == 1
        assert wall["family"] == 0

    def test_free_emits_default_offgen(self, tmp_path):
        payload = _two_block_shared_payload()
        out = write_cc_par(payload, tmp_path / "cc.par", grd_filename="x.grd")
        patches = _parse_patches(out.read_text())
        # Block 1 face 2 (i_hi) is FREE; default free_face_bc = 40.
        free = next(p for p in patches if p["block"] == 1 and p["face"] == 2)
        assert free["bc"] == 40
        assert free["family"] == 0

    def test_free_per_block_override(self, tmp_path):
        payload = _two_block_shared_payload()
        out = write_cc_par(
            payload,
            tmp_path / "cc.par",
            grd_filename="x.grd",
            block_metadata={
                "a": BlockMeta(comment="a", free_face_bc=-9),
                "b": BlockMeta(comment="b", free_face_bc=40),
            },
        )
        patches = _parse_patches(out.read_text())
        # a's FREE faces -> -9, b's FREE faces -> 40.
        a_free = next(p for p in patches if p["block"] == 1 and p["face"] == 2)
        b_free = next(p for p in patches if p["block"] == 2 and p["face"] == 4)
        assert a_free["bc"] == -9
        assert b_free["bc"] == 40

    def test_shared_emits_orientation_codes_with_paired_family(self, tmp_path):
        payload = _two_block_shared_payload()
        out = write_cc_par(payload, tmp_path / "cc.par", grd_filename="x.grd")
        patches = _parse_patches(out.read_text())
        a_seam = next(p for p in patches if p["block"] == 1 and p["face"] == 4)  # a.j_hi
        b_seam = next(p for p in patches if p["block"] == 2 and p["face"] == 3)  # b.j_lo
        # a.j_hi → b.j_lo: this perp at max (j_hi), partner perp at min (j_lo,
        # increasing inward) → lj=3 → 135.  The reverse side b.j_lo → a.j_hi:
        # partner perp at max (j_hi, decreasing inward) → lj=4 → 145.
        # Both sides describe the same seam but from their own frame.
        assert a_seam["bc"] == 135
        assert b_seam["bc"] == 145
        # Cross-referenced patch indices (symmetric).
        assert a_seam["family"] == b_seam["idx"]
        assert b_seam["family"] == a_seam["idx"]

    def test_periodic_warns_and_emits_offgen(self, tmp_path):
        payload = {
            "version": _VERSION,
            "blocks": [
                {
                    "label": "a",
                    "shape": [4, 4, 4],
                    "faces": {
                        "i_lo": _face("WALL"),
                        "i_hi": _face("FREE"),
                        "j_lo": _face("PERIODIC"),
                        "j_hi": _face("PERIODIC"),
                        "k_lo": _face("FREE"),
                        "k_hi": _face("FREE"),
                    },
                }
            ],
        }
        with pytest.warns(UserWarning, match="PERIODIC"):
            out = write_cc_par(payload, tmp_path / "cc.par", grd_filename="x.grd")
        patches = _parse_patches(out.read_text())
        per = next(p for p in patches if p["face"] == 3)
        assert per["bc"] == 40
        assert per["family"] == 0

    def test_unknown_bc_raises(self, tmp_path):
        payload = {
            "version": _VERSION,
            "blocks": [
                {
                    "label": "a",
                    "shape": [4, 4, 4],
                    "faces": {
                        "i_lo": _face("BOGUS"),
                        "i_hi": _face("FREE"),
                        "j_lo": _face("FREE"),
                        "j_hi": _face("FREE"),
                        "k_lo": _face("FREE"),
                        "k_hi": _face("FREE"),
                    },
                }
            ],
        }
        with pytest.raises(ParameterError, match="unknown FaceBC"):
            write_cc_par(payload, tmp_path / "cc.par", grd_filename="x.grd")


# ---------------------------------------------------------------------------
# Symmetry validation
# ---------------------------------------------------------------------------


class TestSymmetry:
    def test_every_connection_is_symmetric(self, tmp_path):
        """Run on the cubed sphere; every BC>99 patch must be paired."""
        from shore.io.cc_par import write_cc_par

        # Build a synthetic 6-block cubed-sphere adjacency by hand.
        payload = _build_cubed_sphere_payload()
        out = write_cc_par(payload, tmp_path / "cc.par", grd_filename="x.grd")
        patches = _parse_patches(out.read_text())
        by_idx = {p["idx"]: p for p in patches}
        for p in patches:
            if p["bc"] <= 99:
                continue
            assert p["family"] in by_idx, f"patch {p['idx']} family={p['family']} out of range"
            partner = by_idx[p["family"]]
            assert partner["family"] == p["idx"], (
                f"patch {p['idx']} -> {p['family']}, but partner -> {partner['family']}"
            )

    def test_asymmetric_partner_raises(self, tmp_path):
        """A hand-constructed broken pairing must be rejected."""
        # Block a says j_hi -> b.j_lo, but block b's j_lo doesn't reference a.
        payload = {
            "version": _VERSION,
            "blocks": [
                {
                    "label": "a",
                    "shape": [4, 4, 4],
                    "faces": {
                        "i_lo": _face("WALL"),
                        "i_hi": _face("FREE"),
                        "j_lo": _face("FREE"),
                        "j_hi": _face("SHARED", ("b", "j_lo"), _identity_map("j_hi")),
                        "k_lo": _face("FREE"),
                        "k_hi": _face("FREE"),
                    },
                },
                {
                    "label": "b",
                    "shape": [4, 4, 4],
                    "faces": {
                        "i_lo": _face("WALL"),
                        "i_hi": _face("FREE"),
                        # b.j_lo points at b.j_hi instead of a.j_hi — broken.
                        "j_lo": _face("SHARED", ("b", "j_hi"), _identity_map("j_lo")),
                        "j_hi": _face("SHARED", ("b", "j_lo"), _identity_map("j_hi")),
                        "k_lo": _face("FREE"),
                        "k_hi": _face("FREE"),
                    },
                },
            ],
        }
        with pytest.raises(ParameterError, match=r"symmetric|asymmetric"):
            write_cc_par(payload, tmp_path / "cc.par", grd_filename="x.grd")


# ---------------------------------------------------------------------------
# Cubed-sphere cap-equator cross-axis pairs
# ---------------------------------------------------------------------------


def _build_cubed_sphere_payload() -> dict:
    """Build a minimal 6-block cubed-sphere adjacency JSON by hand,
    matching what shore.mesh.CubedSphereTopology produces.

    Sub-block ↔ sub-block: same-axis SHARED ring (j_hi ↔ j_lo).
    Sub-block ↔ cap: face-to-face DIRICHLET, with cap-side partner being
    one of the 4 cap lateral faces (i_lo/i_hi/j_lo/j_hi).  The cap-side
    AxisMap is the inverse of the sub-side AxisMap.  Sourced from
    ``shore.mesh.topology._CAP_EQUATOR_NORTH/_SOUTH``.
    """
    # (sub_q, sub_face) -> (cap_face, sub_side_axis_map_pairs)
    # Mirrors shore.mesh.topology._CAP_EQUATOR_{NORTH,SOUTH}.
    cap_pairs: dict[tuple[int, str], tuple[str, str, list[list]]] = {
        (0, "i_hi"): ("i_lo", "cap_north", [["j", False], ["k", False]]),
        (1, "i_hi"): ("j_hi", "cap_north", [["i", False], ["k", False]]),
        (2, "i_hi"): ("i_hi", "cap_north", [["j", True], ["k", False]]),
        (3, "i_hi"): ("j_lo", "cap_north", [["i", True], ["k", False]]),
        (0, "i_lo"): ("j_lo", "cap_south", [["i", False], ["k", False]]),
        (1, "i_lo"): ("i_hi", "cap_south", [["j", False], ["k", False]]),
        (2, "i_lo"): ("j_hi", "cap_south", [["i", True], ["k", False]]),
        (3, "i_lo"): ("i_lo", "cap_south", [["j", True], ["k", False]]),
    }

    def _inverse_map(sub_face: str, cap_face: str, sub_pairs: list[list]) -> list[list]:
        """Return the cap-side axis_map (inverse of sub_pairs)."""
        sub_along = {"i_lo": ("j", "k"), "i_hi": ("j", "k")}[sub_face]
        cap_along = {
            "i_lo": ("j", "k"),
            "i_hi": ("j", "k"),
            "j_lo": ("i", "k"),
            "j_hi": ("i", "k"),
        }[cap_face]
        # rev: partner_axis -> (my_axis_label, flipped)
        rev: dict[str, tuple[str, bool]] = {}
        for my_ax, (p_ax, flip) in zip(sub_along, sub_pairs, strict=True):
            rev[p_ax] = (my_ax, flip)
        return [
            [rev[cap_along[0]][0], rev[cap_along[0]][1]],
            [rev[cap_along[1]][0], rev[cap_along[1]][1]],
        ]

    blocks: list[dict] = []
    # 4 equatorial sub-blocks.
    for q in range(4):
        q_prev = (q - 1) % 4
        q_next = (q + 1) % 4
        cap_n_face, cap_n_label, sub_n_pairs = cap_pairs[(q, "i_hi")]
        cap_s_face, cap_s_label, sub_s_pairs = cap_pairs[(q, "i_lo")]
        blocks.append(
            {
                "label": f"sub{q}",
                "shape": [9, 5, 5],
                "faces": {
                    "i_lo": _face("DIRICHLET", (cap_s_label, cap_s_face), sub_s_pairs),
                    "i_hi": _face("DIRICHLET", (cap_n_label, cap_n_face), sub_n_pairs),
                    "j_lo": _face("SHARED", (f"sub{q_prev}", "j_hi"), _identity_map("j_lo")),
                    "j_hi": _face("SHARED", (f"sub{q_next}", "j_lo"), _identity_map("j_hi")),
                    "k_lo": _face("WALL"),
                    "k_hi": _face("FREE"),
                },
            }
        )

    # Cap blocks: lateral faces partner the matching sub block with the
    # inverse axis_map.  k_lo is the cap's wall, k_hi is its outer.
    for pole in ("north", "south"):
        cap_label = f"cap_{pole}"
        # Build per-cap-face partner lookup from cap_pairs.
        partner_of: dict[str, tuple[int, str, list[list]]] = {}
        for (q, sub_face), (
            cap_face,
            cap_label_t,
            sub_pairs,
        ) in cap_pairs.items():
            if cap_label_t != cap_label:
                continue
            partner_of[cap_face] = (q, sub_face, sub_pairs)

        face_dict: dict = {}
        for cap_face in ("i_lo", "i_hi", "j_lo", "j_hi"):
            q, sub_face, sub_pairs = partner_of[cap_face]
            cap_pairs_inv = _inverse_map(sub_face, cap_face, sub_pairs)
            face_dict[cap_face] = _face("DIRICHLET", (f"sub{q}", sub_face), cap_pairs_inv)
        face_dict["k_lo"] = _face("WALL")
        face_dict["k_hi"] = _face("FREE")
        blocks.append(
            {
                "label": cap_label,
                "shape": [5, 5, 5],
                "faces": face_dict,
            }
        )
    return {"version": _VERSION, "blocks": blocks}


class TestOrientationCodeDerivation:
    """Pure-function tests for :func:`_orientation_code_from_axis_map`.

    Replaces the old `_CAP_EQUATOR_PARTNER` table tests: the writer no
    longer has a per-topology table — orientation codes are derived
    from ``(face_name, partner_face_name, axis_map)``.
    """

    def test_same_axis_no_flip_yields_135(self):
        # Generic same-axis seam: i_hi <-> i_lo, no flips on j or k.
        # Perpendicular: partner is i_lo (Imin) -> li=1.
        # Along-j: partner j, no flip -> lj=3.
        # Along-k: partner k, no flip -> lk=5.
        # Code = 135.
        amap = AxisMap.from_pairs(("j", False), ("k", False))
        assert _orientation_code_from_axis_map("i_hi", "i_lo", amap) == 135

    def test_partner_at_max_yields_2(self):
        # i_hi <-> i_hi: perpendicular partner at Imax -> li=2.
        amap = AxisMap.from_pairs(("j", False), ("k", False))
        assert _orientation_code_from_axis_map("i_hi", "i_hi", amap) == 235

    def test_along_axis_flipped_yields_4(self):
        # i_hi <-> i_hi with j-flip -> lj=4.
        amap = AxisMap.from_pairs(("j", True), ("k", False))
        assert _orientation_code_from_axis_map("i_hi", "i_hi", amap) == 245

    @pytest.mark.parametrize(
        ("face", "partner", "pairs", "expected"),
        [
            # Cubed-sphere cap-equator seams (sub side codes from the
            # historical table).
            ("i_hi", "i_lo", (("j", False), ("k", False)), 135),  # sub0/N
            ("i_hi", "j_hi", (("i", False), ("k", False)), 415),  # sub1/N
            ("i_hi", "i_hi", (("j", True), ("k", False)), 245),  # sub2/N
            ("i_hi", "j_lo", (("i", True), ("k", False)), 325),  # sub3/N
            ("i_lo", "j_lo", (("i", False), ("k", False)), 315),  # sub0/S
            ("i_lo", "i_hi", (("j", False), ("k", False)), 235),  # sub1/S
            ("i_lo", "j_hi", (("i", True), ("k", False)), 425),  # sub2/S
            ("i_lo", "i_lo", (("j", True), ("k", False)), 145),  # sub3/S
        ],
    )
    def test_cubed_sphere_sub_side_codes(self, face, partner, pairs, expected):
        amap = AxisMap.from_pairs(pairs[0], pairs[1])
        assert _orientation_code_from_axis_map(face, partner, amap) == expected

    @pytest.mark.parametrize(
        ("face", "partner", "pairs", "expected"),
        [
            # Cubed-sphere cap-side codes.  Cap faces have along-axes
            # (j, k) for i_lo/i_hi and (i, k) for j_lo/j_hi; the AxisMap
            # below uses the canonical cap-side along-axis order, then
            # maps to the sub partner's along-axes (j, k).
            ("i_lo", "i_hi", (("j", False), ("k", False)), 235),  # cap_n.i_lo ↔ sub0
            ("j_hi", "i_hi", (("j", False), ("k", False)), 325),  # cap_n.j_hi ↔ sub1
            ("i_hi", "i_hi", (("j", True), ("k", False)), 245),  # cap_n.i_hi ↔ sub2
            ("j_lo", "i_hi", (("j", True), ("k", False)), 425),  # cap_n.j_lo ↔ sub3
            ("j_lo", "i_lo", (("j", False), ("k", False)), 315),  # cap_s.j_lo ↔ sub0
            ("i_hi", "i_lo", (("j", False), ("k", False)), 135),  # cap_s.i_hi ↔ sub1
            ("j_hi", "i_lo", (("j", True), ("k", False)), 415),  # cap_s.j_hi ↔ sub2
            ("i_lo", "i_lo", (("j", True), ("k", False)), 145),  # cap_s.i_lo ↔ sub3
        ],
    )
    def test_cubed_sphere_cap_side_codes(self, face, partner, pairs, expected):
        amap = AxisMap.from_pairs(pairs[0], pairs[1])
        assert _orientation_code_from_axis_map(face, partner, amap) == expected


class TestCubedSphereCapEquator:
    def test_cubed_sphere_writes_36_patches(self, tmp_path):
        payload = _build_cubed_sphere_payload()
        out = write_cc_par(payload, tmp_path / "cc.par", grd_filename="x.grd")
        text = out.read_text()
        assert "36 ! patches number" in text
        patches = _parse_patches(text)
        assert len(patches) == 36

    def test_cubed_sphere_sub0_north_cap_pair(self, tmp_path):
        """sub0.i_hi (block 1, face 2) <-> cap_north.i_lo (block 5, face 1).
        Per the table: sub-side BC=135, cap-side BC=235."""
        payload = _build_cubed_sphere_payload()
        out = write_cc_par(payload, tmp_path / "cc.par", grd_filename="x.grd")
        patches = _parse_patches(out.read_text())
        sub0_ihi = next(p for p in patches if p["block"] == 1 and p["face"] == 2)
        cap_n_ilo = next(p for p in patches if p["block"] == 5 and p["face"] == 1)
        assert sub0_ihi["bc"] == 135
        assert cap_n_ilo["bc"] == 235
        assert sub0_ihi["family"] == cap_n_ilo["idx"]
        assert cap_n_ilo["family"] == sub0_ihi["idx"]

    def test_cubed_sphere_sub2_north_cap_pair(self, tmp_path):
        """sub2.i_hi (block 3, face 2) <-> cap_north.i_hi (block 5, face 2).
        Per the table: both sides BC=245 (j-axis flipped on both sides)."""
        payload = _build_cubed_sphere_payload()
        out = write_cc_par(payload, tmp_path / "cc.par", grd_filename="x.grd")
        patches = _parse_patches(out.read_text())
        sub2_ihi = next(p for p in patches if p["block"] == 3 and p["face"] == 2)
        cap_n_ihi = next(p for p in patches if p["block"] == 5 and p["face"] == 2)
        assert sub2_ihi["bc"] == 245
        assert cap_n_ihi["bc"] == 245
        assert sub2_ihi["family"] == cap_n_ihi["idx"]


# ---------------------------------------------------------------------------
# CLI end-to-end
# ---------------------------------------------------------------------------


class TestCLI:
    def test_cc_par_cli_round_trip(self, sphere_stl, tmp_path):
        """Build a real cubed-sphere mesh, write its adjacency JSON, run
        ``shore cc-par``, and parse the result."""
        from shore.mesh import CubedSphereTopology, Mesh

        m = Mesh.from_stl(sphere_stl, topology=CubedSphereTopology(), ni=9, nj=16, nk=5)
        adj = tmp_path / "sphere.adjacency.json"
        write_adjacency_json(m.blocks, adj)
        ccpar = tmp_path / "sphere.cc.par"

        runner = CliRunner()
        r = runner.invoke(app, ["cc-par", str(adj), str(ccpar), "--grd", "sphere.grd"])
        assert r.exit_code == 0, r.output
        assert ccpar.exists()
        text = ccpar.read_text()
        assert "'sphere.grd'" in text
        assert "6 ! blocks number" in text
        assert "36 ! patches number" in text

    def test_cc_par_cli_free_bc_override(self, tmp_path):
        """--free-bc -9 should propagate to every block's FREE faces."""
        adj = tmp_path / "two.adjacency.json"
        adj.write_text(json.dumps(_two_block_shared_payload()))
        ccpar = tmp_path / "two.cc.par"

        runner = CliRunner()
        r = runner.invoke(
            app,
            [
                "cc-par",
                str(adj),
                str(ccpar),
                "--grd",
                "x.grd",
                "--free-bc",
                "-9",
            ],
        )
        assert r.exit_code == 0, r.output
        patches = _parse_patches(ccpar.read_text())
        free_patches = [p for p in patches if p["bc"] == -9]
        # 4 FREE faces per block * 2 blocks = 8 patches
        assert len(free_patches) == 8

    def test_cc_par_cli_missing_input(self, tmp_path):
        runner = CliRunner()
        r = runner.invoke(
            app,
            [
                "cc-par",
                str(tmp_path / "nope.json"),
                str(tmp_path / "cc.par"),
                "--grd",
                "x.grd",
            ],
        )
        assert r.exit_code == 1
        # Click may wrap the error line at the terminal width; collapse
        # whitespace before matching so the test is invariant to wrap pos.
        assert "not found" in " ".join(r.output.split())

    def test_cc_par_cli_meta_inline_json(self, tmp_path):
        """--meta as an inline JSON string overrides per-block free_face_bc."""
        adj = tmp_path / "two.adjacency.json"
        adj.write_text(json.dumps(_two_block_shared_payload()))
        ccpar = tmp_path / "two.cc.par"

        meta = json.dumps({"a": {"free_face_bc": -9}, "b": {"free_face_bc": 40}})
        runner = CliRunner()
        r = runner.invoke(
            app,
            ["cc-par", str(adj), str(ccpar), "--grd", "x.grd", "--meta", meta],
        )
        assert r.exit_code == 0, r.output
        patches = _parse_patches(ccpar.read_text())
        # Block 1 (a) FREE faces -> -9; block 2 (b) FREE faces -> 40.
        a_free = [p for p in patches if p["block"] == 1 and p["bc"] == -9]
        b_free = [p for p in patches if p["block"] == 2 and p["bc"] == 40]
        assert len(a_free) == 4
        assert len(b_free) == 4

    def test_cc_par_cli_meta_file(self, tmp_path):
        """--meta as a path to a JSON file produces the same result as inline."""
        adj = tmp_path / "two.adjacency.json"
        adj.write_text(json.dumps(_two_block_shared_payload()))
        meta_file = tmp_path / "meta.json"
        meta_file.write_text(json.dumps({"a": {"free_face_bc": -9}, "b": {"free_face_bc": 40}}))
        ccpar = tmp_path / "two.cc.par"

        runner = CliRunner()
        r = runner.invoke(
            app,
            ["cc-par", str(adj), str(ccpar), "--grd", "x.grd", "--meta", str(meta_file)],
        )
        assert r.exit_code == 0, r.output
        patches = _parse_patches(ccpar.read_text())
        a_free = [p for p in patches if p["block"] == 1 and p["bc"] == -9]
        b_free = [p for p in patches if p["block"] == 2 and p["bc"] == 40]
        assert len(a_free) == 4
        assert len(b_free) == 4

    def test_cc_par_cli_meta_level_and_comment(self, tmp_path):
        """--meta level and comment fields are written to the blocks section."""
        adj = tmp_path / "two.adjacency.json"
        adj.write_text(json.dumps(_two_block_shared_payload()))
        ccpar = tmp_path / "two.cc.par"

        meta = json.dumps({"a": {"level": 3, "group": 1, "priority": 7, "comment": "body"}})
        runner = CliRunner()
        r = runner.invoke(
            app,
            ["cc-par", str(adj), str(ccpar), "--grd", "x.grd", "--meta", meta],
        )
        assert r.exit_code == 0, r.output
        text = ccpar.read_text()
        assert "        3         1         7 ! 1 : body" in text

    def test_cc_par_cli_meta_not_json_not_file(self, tmp_path):
        """--meta with a value that is neither valid JSON nor an existing file exits 1."""
        adj = tmp_path / "two.adjacency.json"
        adj.write_text(json.dumps(_two_block_shared_payload()))
        runner = CliRunner()
        r = runner.invoke(
            app,
            [
                "cc-par",
                str(adj),
                str(tmp_path / "two.cc.par"),
                "--grd",
                "x.grd",
                "--meta",
                "not-json-not-a-file",
            ],
        )
        assert r.exit_code == 1
        assert "not valid JSON" in " ".join(r.output.split())

    def test_cc_par_cli_meta_invalid_json_file(self, tmp_path):
        """--meta pointing at a file whose content is not valid JSON exits 1."""
        adj = tmp_path / "two.adjacency.json"
        adj.write_text(json.dumps(_two_block_shared_payload()))
        bad_file = tmp_path / "bad.json"
        bad_file.write_text("this is not json {{{")
        runner = CliRunner()
        r = runner.invoke(
            app,
            [
                "cc-par",
                str(adj),
                str(tmp_path / "two.cc.par"),
                "--grd",
                "x.grd",
                "--meta",
                str(bad_file),
            ],
        )
        assert r.exit_code == 1
        assert "not valid JSON" in " ".join(r.output.split())

    def test_cc_par_cli_meta_non_object_json(self, tmp_path):
        """--meta whose JSON value is not an object (e.g. a list) exits 1."""
        adj = tmp_path / "two.adjacency.json"
        adj.write_text(json.dumps(_two_block_shared_payload()))
        runner = CliRunner()
        r = runner.invoke(
            app,
            [
                "cc-par",
                str(adj),
                str(tmp_path / "two.cc.par"),
                "--grd",
                "x.grd",
                "--meta",
                "[1, 2, 3]",
            ],
        )
        assert r.exit_code == 1
        assert "object" in " ".join(r.output.split()).lower()
