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

import numpy as np
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
    bc_code: int | None = None,
) -> dict:
    """Build a face entry for a hand-crafted v2 adjacency payload."""
    entry = {
        "bc": bc,
        "partner": ({"label": partner[0], "face": partner[1]} if partner else None),
        "axis_map": axis_map,
    }
    if bc_code is not None:
        entry["bc_code"] = bc_code
    return entry


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
                "a": BlockMeta(comment="a", free_face_bc=9),
                "b": BlockMeta(comment="b", free_face_bc=40),
            },
        )
        patches = _parse_patches(out.read_text())
        # a's FREE faces -> 9 (natural far-field), b's FREE faces -> 40.
        a_free = next(p for p in patches if p["block"] == 1 and p["face"] == 2)
        b_free = next(p for p in patches if p["block"] == 2 and p["face"] == 4)
        assert a_free["bc"] == 9
        assert b_free["bc"] == 40

    def test_per_face_bc_code_overrides_free_default(self, tmp_path):
        """A per-face ``bc_code`` on a FREE face wins over free_face_bc, so a
        single block can carry inflow (3) on one outer face and extrapolation
        (9) on the rest."""
        payload = _two_block_shared_payload()
        # block "a": i_lo is WALL; put inflow=3 on i_hi, leave the other FREE
        # faces to the per-block default of 9.
        payload["blocks"][0]["faces"]["i_hi"] = _face("FREE", bc_code=3)
        out = write_cc_par(
            payload,
            tmp_path / "cc.par",
            grd_filename="x.grd",
            block_metadata={"a": BlockMeta(comment="a", free_face_bc=9)},
        )
        patches = _parse_patches(out.read_text())
        a_inflow = next(p for p in patches if p["block"] == 1 and p["face"] == 2)  # i_hi
        a_other = next(p for p in patches if p["block"] == 1 and p["face"] == 5)  # k_lo
        assert a_inflow["bc"] == 3  # per-face override
        assert a_other["bc"] == 9  # per-block default still applies

    def test_per_face_bc_code_survives_sidecar_round_trip(self, tmp_path):
        """A FREE face's bc_code set on a HexBlock survives write -> read."""
        from shore.io.adjacency import read_adjacency_json
        from shore.mesh.block import HexBlock
        from shore.mesh.face import FaceBC

        blk = HexBlock(label="bg", nodes=np.zeros((4, 4, 4, 3)), periodic_j=False)
        blk._set_face_bc("i_lo", FaceBC.FREE, bc_code=3)  # inflow
        path = write_adjacency_json([blk], tmp_path / "bg.adjacency.json")
        data = read_adjacency_json(path)
        face = data["blocks"][0]["faces"]["i_lo"]
        assert face["bc"] == "FREE"
        assert face["bc_code"] == 3
        # a FREE face with no override emits no bc_code key
        assert "bc_code" not in data["blocks"][0]["faces"]["i_hi"]

    def test_bc_code_rejected_on_non_free_face(self):
        from shore.mesh.block import HexBlock
        from shore.mesh.face import FaceBC

        blk = HexBlock(label="b", nodes=np.zeros((4, 4, 4, 3)), periodic_j=False)
        with pytest.raises(ValueError, match="only valid on FREE"):
            blk._set_face_bc("k_lo", FaceBC.WALL, bc_code=3)

    def test_bc_code_out_of_range_rejected_by_reader(self, tmp_path):
        payload = _two_block_shared_payload()
        payload["blocks"][0]["faces"]["i_hi"] = _face("FREE", bc_code=99)
        adj = tmp_path / "bad.adjacency.json"
        adj.write_text(json.dumps(payload))
        from shore.exceptions import MeshInputError
        from shore.io.adjacency import read_adjacency_json

        with pytest.raises(MeshInputError, match=r"natural-BC code in \[1, 19\]"):
            read_adjacency_json(adj)

    def test_shared_emits_orientation_codes_with_paired_family(self, tmp_path):
        payload = _two_block_shared_payload()
        out = write_cc_par(payload, tmp_path / "cc.par", grd_filename="x.grd")
        patches = _parse_patches(out.read_text())
        a_seam = next(p for p in patches if p["block"] == 1 and p["face"] == 4)  # a.j_hi
        b_seam = next(p for p in patches if p["block"] == 2 and p["face"] == 3)  # b.j_lo
        # A j_hi <-> j_lo seam: the faces are on opposite sides of their blocks,
        # so Xall's pre-processor reads an odd normal digit (normind, ADR 0036):
        # 135 from both sides, as in a solver-run INSEAN cc.par.
        assert a_seam["bc"] == 135
        assert b_seam["bc"] == 135
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
            # Cubed-sphere cap-equator seams, sub side; the normal digit is odd
            # for opposite sides (lo <-> hi), even for the same side (ADR 0036).
            ("i_hi", "i_lo", (("j", False), ("k", False)), 135),  # sub0/N
            ("i_hi", "j_hi", (("i", False), ("k", False)), 415),  # sub1/N
            ("i_hi", "i_hi", (("j", True), ("k", False)), 245),  # sub2/N
            ("i_hi", "j_lo", (("i", True), ("k", False)), 325),  # sub3/N
            ("i_lo", "j_lo", (("i", False), ("k", False)), 415),  # sub0/S
            ("i_lo", "i_hi", (("j", False), ("k", False)), 135),  # sub1/S
            ("i_lo", "j_hi", (("i", True), ("k", False)), 325),  # sub2/S
            ("i_lo", "i_lo", (("j", True), ("k", False)), 245),  # sub3/S
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
            ("i_lo", "i_hi", (("j", False), ("k", False)), 135),  # cap_n.i_lo ↔ sub0
            ("j_hi", "i_hi", (("j", False), ("k", False)), 325),  # cap_n.j_hi ↔ sub1
            ("i_hi", "i_hi", (("j", True), ("k", False)), 245),  # cap_n.i_hi ↔ sub2
            ("j_lo", "i_hi", (("j", True), ("k", False)), 415),  # cap_n.j_lo ↔ sub3
            ("j_lo", "i_lo", (("j", False), ("k", False)), 325),  # cap_s.j_lo ↔ sub0
            ("i_hi", "i_lo", (("j", False), ("k", False)), 135),  # cap_s.i_hi ↔ sub1
            ("j_hi", "i_lo", (("j", True), ("k", False)), 415),  # cap_s.j_hi ↔ sub2
            ("i_lo", "i_lo", (("j", True), ("k", False)), 245),  # cap_s.i_lo ↔ sub3
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
        Opposite sides: BC=135 on both (ADR 0036)."""
        payload = _build_cubed_sphere_payload()
        out = write_cc_par(payload, tmp_path / "cc.par", grd_filename="x.grd")
        patches = _parse_patches(out.read_text())
        sub0_ihi = next(p for p in patches if p["block"] == 1 and p["face"] == 2)
        cap_n_ilo = next(p for p in patches if p["block"] == 5 and p["face"] == 1)
        assert sub0_ihi["bc"] == 135
        assert cap_n_ilo["bc"] == 135
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

        meta = json.dumps({"a": {"free_face_bc": 9}, "b": {"free_face_bc": 40}})
        runner = CliRunner()
        r = runner.invoke(
            app,
            ["cc-par", str(adj), str(ccpar), "--grd", "x.grd", "--meta", meta],
        )
        assert r.exit_code == 0, r.output
        patches = _parse_patches(ccpar.read_text())
        # Block 1 (a) FREE faces -> 9; block 2 (b) FREE faces -> 40.
        a_free = [p for p in patches if p["block"] == 1 and p["bc"] == 9]
        b_free = [p for p in patches if p["block"] == 2 and p["bc"] == 40]
        assert len(a_free) == 4
        assert len(b_free) == 4

    def test_cc_par_cli_meta_file(self, tmp_path):
        """--meta as a path to a JSON file produces the same result as inline."""
        adj = tmp_path / "two.adjacency.json"
        adj.write_text(json.dumps(_two_block_shared_payload()))
        meta_file = tmp_path / "meta.json"
        meta_file.write_text(json.dumps({"a": {"free_face_bc": 9}, "b": {"free_face_bc": 40}}))
        ccpar = tmp_path / "two.cc.par"

        runner = CliRunner()
        r = runner.invoke(
            app,
            ["cc-par", str(adj), str(ccpar), "--grd", "x.grd", "--meta", str(meta_file)],
        )
        assert r.exit_code == 0, r.output
        patches = _parse_patches(ccpar.read_text())
        a_free = [p for p in patches if p["block"] == 1 and p["bc"] == 9]
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


# ---------------------------------------------------------------------------
# merge_cc_par
# ---------------------------------------------------------------------------


def _write_two_cc_pars(tmp_path) -> tuple:
    """Generate two cc.par files from two independent two-block payloads.

    Returns (a_path, b_path, n_patches_each).  Each input has 2 blocks
    and 12 patches (one SHARED seam pair on j; everything else
    WALL/FREE).  The two files are independent — no cross-file
    adjacency, exactly the expected use case for merge_cc_par.
    """
    a_payload = _two_block_shared_payload()  # blocks "a", "b"
    a_path = tmp_path / "a.cc.par"
    write_cc_par(a_payload, a_path, grd_filename="a.grd")

    # Build a second payload with different labels so we can verify
    # the auto-prefix behaviour distinguishes them.
    b_payload = {
        "version": _VERSION,
        "blocks": [
            {
                "label": "x",
                "shape": [4, 4, 4],
                "faces": {
                    "i_lo": _face("WALL"),
                    "i_hi": _face("FREE"),
                    "j_lo": _face("FREE"),
                    "j_hi": _face("SHARED", ("y", "j_lo"), _identity_map("j_hi")),
                    "k_lo": _face("FREE"),
                    "k_hi": _face("FREE"),
                },
            },
            {
                "label": "y",
                "shape": [4, 4, 4],
                "faces": {
                    "i_lo": _face("WALL"),
                    "i_hi": _face("FREE"),
                    "j_lo": _face("SHARED", ("x", "j_hi"), _identity_map("j_lo")),
                    "j_hi": _face("FREE"),
                    "k_lo": _face("FREE"),
                    "k_hi": _face("FREE"),
                },
            },
        ],
    }
    b_path = tmp_path / "b.cc.par"
    write_cc_par(b_payload, b_path, grd_filename="b.grd")
    return a_path, b_path


def _parse_block_lines(text: str) -> list[tuple]:
    """Pull out (level, group, priority, comment) tuples from block section."""
    out = []
    in_blocks = False
    n_expected = None
    for line in text.splitlines():
        stripped = line.strip()
        if "blocks number" in stripped:
            n_expected = int(stripped.split("!", 1)[0].strip())
            in_blocks = True
            continue
        if in_blocks and "patches number" in stripped:
            break
        if not in_blocks or not stripped or "!" not in stripped:
            continue
        head, _, tail = stripped.partition("!")
        parts = head.split()
        if len(parts) != 3:
            continue
        level, group, priority = (int(x) for x in parts)
        comment = tail.split(":", 1)[-1].strip() if ":" in tail else ""
        out.append((level, group, priority, comment))
        if n_expected is not None and len(out) == n_expected:
            break
    return out


class TestMergeCCPar:
    def test_merged_block_count_is_sum(self, tmp_path):
        from shore.io.cc_par import merge_cc_par

        a, b = _write_two_cc_pars(tmp_path)
        out = merge_cc_par([a, b], tmp_path / "merged.cc.par", grd_filename="m.grd")
        text = out.read_text()
        assert "4 ! blocks number" in text

    def test_merged_patch_count_is_sum(self, tmp_path):
        from shore.io.cc_par import merge_cc_par

        a, b = _write_two_cc_pars(tmp_path)
        out = merge_cc_par([a, b], tmp_path / "merged.cc.par", grd_filename="m.grd")
        text = out.read_text()
        # 12 patches per file (6 per block * 2 blocks).
        assert "24 ! patches number" in text

    def test_grd_filename_overrides_inputs0(self, tmp_path):
        from shore.io.cc_par import merge_cc_par

        a, b = _write_two_cc_pars(tmp_path)
        out = merge_cc_par([a, b], tmp_path / "merged.cc.par", grd_filename="merged.grd")
        first_line = out.read_text().splitlines()[0]
        assert first_line == "'merged.grd'"

    def test_block_comments_auto_prefixed(self, tmp_path):
        from shore.io.cc_par import merge_cc_par

        a, b = _write_two_cc_pars(tmp_path)
        out = merge_cc_par([a, b], tmp_path / "merged.cc.par", grd_filename="m.grd")
        block_lines = _parse_block_lines(out.read_text())
        comments = [bl[3] for bl in block_lines]
        # File stems are 'a' and 'b' (the trailing .cc is stripped).
        assert comments == ["a.a", "a.b", "b.x", "b.y"]

    def test_block_indices_are_sequential(self, tmp_path):
        from shore.io.cc_par import merge_cc_par

        a, b = _write_two_cc_pars(tmp_path)
        out = merge_cc_par([a, b], tmp_path / "merged.cc.par", grd_filename="m.grd")
        text = out.read_text()
        # Block-section trailing comments include "1 :", "2 :", ... in order.
        for n in (1, 2, 3, 4):
            assert f"! {n} :" in text

    def test_patch_block_index_shifted(self, tmp_path):
        from shore.io.cc_par import merge_cc_par

        a, b = _write_two_cc_pars(tmp_path)
        out = merge_cc_par([a, b], tmp_path / "merged.cc.par", grd_filename="m.grd")
        patches = _parse_patches(out.read_text())
        # b's blocks (originally 1, 2) are now blocks 3, 4.
        b_blocks = {p["block"] for p in patches if p["idx"] > 12}
        assert b_blocks == {3, 4}

    def test_patch_family_pointers_remain_symmetric(self, tmp_path):
        """Every connection (BC≥100) patch must still pair symmetrically
        after the merge offset is applied."""
        from shore.io.cc_par import merge_cc_par

        a, b = _write_two_cc_pars(tmp_path)
        out = merge_cc_par([a, b], tmp_path / "merged.cc.par", grd_filename="m.grd")
        patches = _parse_patches(out.read_text())
        by_idx = {p["idx"]: p for p in patches}
        for p in patches:
            if p["bc"] < 100 or p["family"] == 0:
                continue
            partner = by_idx[p["family"]]
            assert partner["family"] == p["idx"], (
                f"asymmetric pairing: patch {p['idx']} -> {p['family']}, "
                f"but {p['family']} -> {partner['family']}"
            )

    def test_a_side_family_unchanged(self, tmp_path):
        """Patches inside inputs[0] reference partners inside inputs[0],
        so their family pointers must NOT be shifted."""
        from shore.io.cc_par import merge_cc_par

        a, b = _write_two_cc_pars(tmp_path)
        # Read the original a-only patches.
        a_only = _parse_patches(a.read_text())
        a_family_by_idx = {p["idx"]: p["family"] for p in a_only}

        out = merge_cc_par([a, b], tmp_path / "merged.cc.par", grd_filename="m.grd")
        merged = _parse_patches(out.read_text())
        for p in merged:
            if p["idx"] > 12:
                break
            assert p["family"] == a_family_by_idx[p["idx"]], (
                f"a-side patch {p['idx']} family changed unexpectedly"
            )

    def test_b_side_family_offset_by_a_count(self, tmp_path):
        """Patches inside inputs[1] reference partners inside inputs[1];
        their family pointers must shift by len(a's patches) = 12."""
        from shore.io.cc_par import merge_cc_par

        a, b = _write_two_cc_pars(tmp_path)
        b_only = _parse_patches(b.read_text())
        b_family_by_idx = {p["idx"]: p["family"] for p in b_only}

        out = merge_cc_par([a, b], tmp_path / "merged.cc.par", grd_filename="m.grd")
        merged = _parse_patches(out.read_text())
        for p in merged:
            if p["idx"] <= 12:
                continue
            original_idx = p["idx"] - 12
            original_family = b_family_by_idx[original_idx]
            expected = original_family + 12 if original_family > 0 else 0
            assert p["family"] == expected

    def test_warns_on_header_drift(self, tmp_path):
        """Disagreement between input headers (other than grd filename)
        emits a UserWarning but the merge proceeds."""
        from shore.io.cc_par import merge_cc_par

        a_payload = _two_block_shared_payload()
        b_payload = _two_block_shared_payload()
        a = tmp_path / "a.cc.par"
        b = tmp_path / "b.cc.par"
        write_cc_par(a, a, grd_filename="a.grd") if False else None
        write_cc_par(a_payload, a, grd_filename="a.grd", multigrid_levels=(4, 1))
        write_cc_par(b_payload, b, grd_filename="b.grd", multigrid_levels=(3, 1))
        with pytest.warns(UserWarning, match="header line"):
            merge_cc_par([a, b], tmp_path / "merged.cc.par", grd_filename="m.grd")

    def test_rejects_single_input(self, tmp_path):
        from shore.io.cc_par import merge_cc_par

        a, _ = _write_two_cc_pars(tmp_path)
        with pytest.raises(ValueError, match="at least 2"):
            merge_cc_par([a], tmp_path / "merged.cc.par", grd_filename="m.grd")

    def test_rejects_missing_file(self, tmp_path):
        from shore.io.cc_par import merge_cc_par

        a, _ = _write_two_cc_pars(tmp_path)
        with pytest.raises(FileNotFoundError, match="not found"):
            merge_cc_par(
                [a, tmp_path / "nope.cc.par"], tmp_path / "merged.cc.par", grd_filename="m.grd"
            )


class TestCLICCParMerge:
    def test_basic_invocation(self, tmp_path):
        a, b = _write_two_cc_pars(tmp_path)
        runner = CliRunner()
        r = runner.invoke(
            app,
            [
                "cc-par-merge",
                str(a),
                str(b),
                "-o",
                str(tmp_path / "merged.cc.par"),
                "--grd",
                "merged.grd",
            ],
        )
        assert r.exit_code == 0, r.output
        text = (tmp_path / "merged.cc.par").read_text()
        assert "4 ! blocks number" in text
        assert "'merged.grd'" in text.splitlines()[0]

    def test_single_input_rejected(self, tmp_path):
        a, _ = _write_two_cc_pars(tmp_path)
        runner = CliRunner()
        r = runner.invoke(
            app,
            [
                "cc-par-merge",
                str(a),
                "-o",
                str(tmp_path / "merged.cc.par"),
                "--grd",
                "m.grd",
            ],
        )
        assert r.exit_code == 1
        assert "at least 2" in r.output

    def test_missing_input_exits_1(self, tmp_path):
        a, _ = _write_two_cc_pars(tmp_path)
        runner = CliRunner()
        r = runner.invoke(
            app,
            [
                "cc-par-merge",
                str(a),
                str(tmp_path / "nope.cc.par"),
                "-o",
                str(tmp_path / "merged.cc.par"),
                "--grd",
                "m.grd",
            ],
        )
        assert r.exit_code == 1
        assert "not found" in r.output.lower()


# ---------------------------------------------------------------------------
# write_cc_par(boxes=...) — boxes section emission
# ---------------------------------------------------------------------------


def _parse_box_records(text: str) -> list[dict]:
    """Pull box records out of the boxes section.

    Each record: ``{"type": int, "block": int, "vertices": (8, 3) array}``.
    Skips non-type-9 boxes (returns them as ``{"type": int, "block": int}``).
    Returns an empty list when ``boxes number`` = 0.
    """
    import numpy as np

    out: list[dict] = []
    lines = text.splitlines()
    in_boxes = False
    nbox = None
    cursor = 0
    for idx, line in enumerate(lines):
        if "boxes number" in line:
            nbox = int(line.split("!", 1)[0].strip())
            in_boxes = True
            cursor = idx + 1
            break
    if not in_boxes or nbox == 0:
        return out

    while len(out) < nbox and cursor < len(lines):
        line = lines[cursor].strip()
        cursor += 1
        if not line:
            continue
        head = line.split("!", 1)[0]
        parts = head.split()
        if len(parts) != 2:
            continue
        box_type, block_idx = int(parts[0]), int(parts[1])
        if box_type == 9:
            verts = np.empty((8, 3), dtype=np.float64)
            v = 0
            while v < 8 and cursor < len(lines):
                vline = lines[cursor].strip()
                cursor += 1
                if not vline:
                    continue
                vparts = vline.split("!", 1)[0].split()
                verts[v] = [float(x) for x in vparts]
                v += 1
            out.append({"type": 9, "block": block_idx, "vertices": verts})
        else:
            out.append({"type": box_type, "block": block_idx})
    return out


class TestWriteCCParBoxes:
    def test_no_boxes_emits_zero_count(self, tmp_path):
        """Default behaviour (boxes=None) emits '0 ! boxes number'."""
        payload = _two_block_shared_payload()
        out = write_cc_par(payload, tmp_path / "cc.par", grd_filename="x.grd")
        text = out.read_text()
        assert "0 ! boxes number" in text
        assert _parse_box_records(text) == []

    def test_single_box_emits_correct_section(self, tmp_path):
        from shore.io.boxes import Box

        payload = _two_block_shared_payload()
        boxes = [Box.from_aabb((-0.5, -0.5, -0.5), (0.5, 0.5, 0.5), "a")]
        out = write_cc_par(
            payload,
            tmp_path / "cc.par",
            grd_filename="x.grd",
            boxes=boxes,
        )
        text = out.read_text()
        assert "1 ! boxes number" in text
        records = _parse_box_records(text)
        assert len(records) == 1
        # Block 'a' is index 1 in the payload.
        assert records[0]["type"] == 9
        assert records[0]["block"] == 1
        # Canonical AABB vertex ordering.
        assert records[0]["vertices"][0].tolist() == [-0.5, -0.5, -0.5]
        assert records[0]["vertices"][7].tolist() == [0.5, 0.5, 0.5]

    def test_multiple_boxes_block_index_resolution(self, tmp_path):
        """Block labels resolve to their 1-based positions in the payload."""
        from shore.io.boxes import Box

        payload = _two_block_shared_payload()  # blocks "a", "b"
        boxes = [
            Box.from_aabb((0, 0, 0), (1, 1, 1), "a"),  # block index 1
            Box.from_aabb((10, 0, 0), (11, 1, 1), "b"),  # block index 2
            Box.from_aabb((-1, -1, -1), (0, 0, 0), "a"),  # block index 1 again
        ]
        out = write_cc_par(
            payload,
            tmp_path / "cc.par",
            grd_filename="x.grd",
            boxes=boxes,
        )
        records = _parse_box_records(out.read_text())
        assert [r["block"] for r in records] == [1, 2, 1]

    def test_unknown_block_label_rejected(self, tmp_path):
        """Box pointing at a label not in the adjacency must raise."""
        from shore.io.boxes import Box

        payload = _two_block_shared_payload()
        boxes = [Box.from_aabb((0, 0, 0), (1, 1, 1), "ghost")]
        with pytest.raises(ParameterError, match="unknown block"):
            write_cc_par(
                payload,
                tmp_path / "cc.par",
                grd_filename="x.grd",
                boxes=boxes,
            )

    def test_boxes_dont_affect_other_sections(self, tmp_path):
        """Adding boxes must not change the blocks / patches sections."""
        from shore.io.boxes import Box

        payload = _two_block_shared_payload()
        out_no = write_cc_par(payload, tmp_path / "no.cc.par", grd_filename="x.grd")
        out_yes = write_cc_par(
            payload,
            tmp_path / "yes.cc.par",
            grd_filename="x.grd",
            boxes=[Box.from_aabb((0, 0, 0), (1, 1, 1), "a")],
        )
        # Same blocks + patches; differ only after the patches section.
        text_no = out_no.read_text()
        text_yes = out_yes.read_text()
        edges_no = text_no.index("0 ! edges number")
        edges_yes = text_yes.index("0 ! edges number")
        assert text_no[:edges_no] == text_yes[:edges_yes]


class TestMergeCCParBoxes:
    """The merger must concatenate boxes too, with each box's
    associated-block index shifted by the cumulative block offset of
    the preceding inputs.
    """

    def _write_pair(self, tmp_path, *, boxes_in_a: int, boxes_in_b: int):
        from shore.io.boxes import Box

        a_payload = _two_block_shared_payload()  # blocks a, b
        a_boxes = [
            Box.from_aabb((i, 0, 0), (i + 0.5, 0.5, 0.5), "a" if i % 2 == 0 else "b")
            for i in range(boxes_in_a)
        ]
        a_path = tmp_path / "A.cc.par"
        write_cc_par(a_payload, a_path, grd_filename="A.grd", boxes=a_boxes or None)

        b_payload = {
            "version": _VERSION,
            "blocks": [
                {
                    "label": "x",
                    "shape": [4, 4, 4],
                    "faces": {
                        "i_lo": _face("WALL"),
                        "i_hi": _face("FREE"),
                        "j_lo": _face("FREE"),
                        "j_hi": _face("SHARED", ("y", "j_lo"), _identity_map("j_hi")),
                        "k_lo": _face("FREE"),
                        "k_hi": _face("FREE"),
                    },
                },
                {
                    "label": "y",
                    "shape": [4, 4, 4],
                    "faces": {
                        "i_lo": _face("WALL"),
                        "i_hi": _face("FREE"),
                        "j_lo": _face("SHARED", ("x", "j_hi"), _identity_map("j_lo")),
                        "j_hi": _face("FREE"),
                        "k_lo": _face("FREE"),
                        "k_hi": _face("FREE"),
                    },
                },
            ],
        }
        b_boxes = [
            Box.from_aabb(
                (100 + i, 0, 0),
                (100 + i + 0.5, 0.5, 0.5),
                "x" if i % 2 == 0 else "y",
            )
            for i in range(boxes_in_b)
        ]
        b_path = tmp_path / "B.cc.par"
        write_cc_par(b_payload, b_path, grd_filename="B.grd", boxes=b_boxes or None)
        return a_path, b_path

    def test_zero_boxes_each_side_passes_through(self, tmp_path):
        from shore.io.cc_par import merge_cc_par

        a, b = self._write_pair(tmp_path, boxes_in_a=0, boxes_in_b=0)
        out = merge_cc_par([a, b], tmp_path / "merged.cc.par", grd_filename="m.grd")
        text = out.read_text()
        assert "0 ! boxes number" in text

    def test_boxes_only_in_a(self, tmp_path):
        """When only inputs[0] has boxes, they appear unchanged in the merge."""
        from shore.io.cc_par import merge_cc_par

        a, b = self._write_pair(tmp_path, boxes_in_a=2, boxes_in_b=0)
        out = merge_cc_par([a, b], tmp_path / "merged.cc.par", grd_filename="m.grd")
        records = _parse_box_records(out.read_text())
        assert len(records) == 2
        # A's block indices (1, 2) should be unchanged in the merge.
        assert {r["block"] for r in records} == {1, 2}

    def test_boxes_only_in_b_get_offset(self, tmp_path):
        """B's box block indices shift by len(A.blocks) = 2."""
        from shore.io.cc_par import merge_cc_par

        a, b = self._write_pair(tmp_path, boxes_in_a=0, boxes_in_b=2)
        out = merge_cc_par([a, b], tmp_path / "merged.cc.par", grd_filename="m.grd")
        records = _parse_box_records(out.read_text())
        assert len(records) == 2
        # B's blocks (1, 2) become (3, 4) in the merged file.
        assert {r["block"] for r in records} == {3, 4}

    def test_boxes_in_both_concatenated_and_offset(self, tmp_path):
        from shore.io.cc_par import merge_cc_par

        a, b = self._write_pair(tmp_path, boxes_in_a=1, boxes_in_b=1)
        out = merge_cc_par([a, b], tmp_path / "merged.cc.par", grd_filename="m.grd")
        records = _parse_box_records(out.read_text())
        assert len(records) == 2
        # First box from A: block 1 (i=0 → "a").
        assert records[0]["block"] == 1
        # Second box from B: B's "x" was 1, shifted by 2 → 3.
        assert records[1]["block"] == 3

    def test_box_vertices_byte_identical(self, tmp_path):
        """Vertex coordinates pass through untouched."""
        from shore.io.boxes import Box
        from shore.io.cc_par import merge_cc_par

        a, b = self._write_pair(tmp_path, boxes_in_a=0, boxes_in_b=0)
        # Re-write A with a known box so we can compare.
        payload = _two_block_shared_payload()
        verts = np.array(
            [
                [-0.5, -0.5, -0.5],
                [0.5, -0.5, -0.5],
                [-0.5, 0.5, -0.5],
                [0.5, 0.5, -0.5],
                [-0.5, -0.5, 0.5],
                [0.5, -0.5, 0.5],
                [-0.5, 0.5, 0.5],
                [0.5, 0.5, 0.5],
            ]
        )
        write_cc_par(
            payload, a, grd_filename="A.grd", boxes=[Box(vertices=verts, associated_block="a")]
        )

        out = merge_cc_par([a, b], tmp_path / "merged.cc.par", grd_filename="m.grd")
        records = _parse_box_records(out.read_text())
        assert len(records) == 1
        np.testing.assert_array_equal(records[0]["vertices"], verts)


# ---------------------------------------------------------------------------
# CLI: --boxes flag on shore cc-par
# ---------------------------------------------------------------------------


class TestCLICCParBoxes:
    def _payload_path(self, tmp_path):
        adj = tmp_path / "two.adjacency.json"
        adj.write_text(json.dumps(_two_block_shared_payload()))
        return adj

    def test_inline_boxes_json(self, tmp_path):
        adj = self._payload_path(tmp_path)
        boxes_json = json.dumps(
            [
                {
                    "associated_block": "a",
                    "vertices": [
                        [-0.5, -0.5, -0.5],
                        [0.5, -0.5, -0.5],
                        [-0.5, 0.5, -0.5],
                        [0.5, 0.5, -0.5],
                        [-0.5, -0.5, 0.5],
                        [0.5, -0.5, 0.5],
                        [-0.5, 0.5, 0.5],
                        [0.5, 0.5, 0.5],
                    ],
                }
            ]
        )
        runner = CliRunner()
        r = runner.invoke(
            app,
            [
                "cc-par",
                str(adj),
                str(tmp_path / "out.cc.par"),
                "--grd",
                "x.grd",
                "--boxes",
                boxes_json,
            ],
        )
        assert r.exit_code == 0, r.output
        text = (tmp_path / "out.cc.par").read_text()
        assert "1 ! boxes number" in text

    def test_boxes_from_file(self, tmp_path):
        adj = self._payload_path(tmp_path)
        boxes_file = tmp_path / "boxes.json"
        boxes_file.write_text(
            json.dumps(
                [
                    {
                        "associated_block": "b",
                        "vertices": [
                            [0, 0, 0],
                            [1, 0, 0],
                            [0, 1, 0],
                            [1, 1, 0],
                            [0, 0, 1],
                            [1, 0, 1],
                            [0, 1, 1],
                            [1, 1, 1],
                        ],
                    }
                ]
            )
        )
        runner = CliRunner()
        r = runner.invoke(
            app,
            [
                "cc-par",
                str(adj),
                str(tmp_path / "out.cc.par"),
                "--grd",
                "x.grd",
                "--boxes",
                str(boxes_file),
            ],
        )
        assert r.exit_code == 0, r.output
        records = _parse_box_records((tmp_path / "out.cc.par").read_text())
        assert records[0]["block"] == 2  # "b" is the second block

    def test_unknown_block_label_exits_1(self, tmp_path):
        adj = self._payload_path(tmp_path)
        runner = CliRunner()
        r = runner.invoke(
            app,
            [
                "cc-par",
                str(adj),
                str(tmp_path / "out.cc.par"),
                "--grd",
                "x.grd",
                "--boxes",
                json.dumps(
                    [
                        {
                            "associated_block": "ghost",
                            "vertices": [
                                [0, 0, 0],
                                [1, 0, 0],
                                [0, 1, 0],
                                [1, 1, 0],
                                [0, 0, 1],
                                [1, 0, 1],
                                [0, 1, 1],
                                [1, 1, 1],
                            ],
                        }
                    ]
                ),
            ],
        )
        assert r.exit_code == 1
        assert "unknown block" in r.output

    def test_invalid_box_geometry_exits_1(self, tmp_path):
        """Coplanar vertices are rejected at Box construction; the CLI
        surfaces the ParameterError as a non-zero exit."""
        adj = self._payload_path(tmp_path)
        # Eight points all on z=0 (degenerate).
        flat = [
            [0, 0, 0],
            [1, 0, 0],
            [0, 1, 0],
            [1, 1, 0],
            [0.1, 0.1, 0],
            [0.9, 0.1, 0],
            [0.1, 0.9, 0],
            [0.9, 0.9, 0],
        ]
        runner = CliRunner()
        r = runner.invoke(
            app,
            [
                "cc-par",
                str(adj),
                str(tmp_path / "out.cc.par"),
                "--grd",
                "x.grd",
                "--boxes",
                json.dumps([{"associated_block": "a", "vertices": flat}]),
            ],
        )
        assert r.exit_code == 1
        assert "Box" in r.output or "degenerate" in r.output.lower()

    def test_non_list_json_exits_1(self, tmp_path):
        adj = self._payload_path(tmp_path)
        runner = CliRunner()
        r = runner.invoke(
            app,
            [
                "cc-par",
                str(adj),
                str(tmp_path / "out.cc.par"),
                "--grd",
                "x.grd",
                "--boxes",
                json.dumps({"not": "a list"}),
            ],
        )
        assert r.exit_code == 1
        assert "list" in r.output.lower()

    def test_invalid_json_not_a_file_exits_1(self, tmp_path):
        adj = self._payload_path(tmp_path)
        runner = CliRunner()
        r = runner.invoke(
            app,
            [
                "cc-par",
                str(adj),
                str(tmp_path / "out.cc.par"),
                "--grd",
                "x.grd",
                "--boxes",
                "not-json-not-a-file",
            ],
        )
        assert r.exit_code == 1
        assert "not valid JSON" in r.output
