"""Tests for the Xall cc.par writer.

Covers:
- header field formatting
- block section
- BC code mapping (WALL/FREE/SHARED/DIRICHLET/PERIODIC)
- patch index assignment + symmetric family pairing
- cubed-sphere cap-equator orientation table (8 cross-axis pairs)
- BlockMeta per-block free_face_bc override
- end-to-end via shore cc-par CLI
"""

from __future__ import annotations

import json

import pytest
from typer.testing import CliRunner

from shore.cli import app
from shore.exceptions import ParameterError
from shore.io.adjacency import write_adjacency_json
from shore.io.cc_par import (
    _CAP_EQUATOR_PARTNER,
    BlockMeta,
    write_cc_par,
)

# ---------------------------------------------------------------------------
# Tiny payload fixtures (no real meshing required)
# ---------------------------------------------------------------------------


def _two_block_shared_payload() -> dict:
    """Two adjacent blocks sharing a SHARED j_hi/j_lo seam."""
    return {
        "version": 1,
        "blocks": [
            {
                "label": "a",
                "shape": [4, 4, 4],
                "faces": {
                    "i_lo": {"bc": "WALL", "partner": None},
                    "i_hi": {"bc": "FREE", "partner": None},
                    "j_lo": {"bc": "FREE", "partner": None},
                    "j_hi": {"bc": "SHARED", "partner": {"label": "b", "face": "j_lo"}},
                    "k_lo": {"bc": "FREE", "partner": None},
                    "k_hi": {"bc": "FREE", "partner": None},
                },
            },
            {
                "label": "b",
                "shape": [4, 4, 4],
                "faces": {
                    "i_lo": {"bc": "WALL", "partner": None},
                    "i_hi": {"bc": "FREE", "partner": None},
                    "j_lo": {"bc": "SHARED", "partner": {"label": "a", "face": "j_hi"}},
                    "j_hi": {"bc": "FREE", "partner": None},
                    "k_lo": {"bc": "FREE", "partner": None},
                    "k_hi": {"bc": "FREE", "partner": None},
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
        assert "        1         0         1 ! a" in text
        assert "        1         0         1 ! b" in text

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
        assert "        2         1         7 ! aero block" in text
        assert "        1         0         1 ! b" in text


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

    def test_shared_emits_135_with_paired_family(self, tmp_path):
        payload = _two_block_shared_payload()
        out = write_cc_par(payload, tmp_path / "cc.par", grd_filename="x.grd")
        patches = _parse_patches(out.read_text())
        a_seam = next(p for p in patches if p["block"] == 1 and p["face"] == 4)  # a.j_hi
        b_seam = next(p for p in patches if p["block"] == 2 and p["face"] == 3)  # b.j_lo
        assert a_seam["bc"] == 135
        assert b_seam["bc"] == 135
        # Cross-referenced patch indices (symmetric).
        assert a_seam["family"] == b_seam["idx"]
        assert b_seam["family"] == a_seam["idx"]

    def test_periodic_warns_and_emits_offgen(self, tmp_path):
        payload = {
            "version": 1,
            "blocks": [
                {
                    "label": "a",
                    "shape": [4, 4, 4],
                    "faces": {
                        "i_lo": {"bc": "WALL", "partner": None},
                        "i_hi": {"bc": "FREE", "partner": None},
                        "j_lo": {"bc": "PERIODIC", "partner": None},
                        "j_hi": {"bc": "PERIODIC", "partner": None},
                        "k_lo": {"bc": "FREE", "partner": None},
                        "k_hi": {"bc": "FREE", "partner": None},
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
            "version": 1,
            "blocks": [
                {
                    "label": "a",
                    "shape": [4, 4, 4],
                    "faces": {
                        "i_lo": {"bc": "BOGUS", "partner": None},
                        "i_hi": {"bc": "FREE", "partner": None},
                        "j_lo": {"bc": "FREE", "partner": None},
                        "j_hi": {"bc": "FREE", "partner": None},
                        "k_lo": {"bc": "FREE", "partner": None},
                        "k_hi": {"bc": "FREE", "partner": None},
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
            "version": 1,
            "blocks": [
                {
                    "label": "a",
                    "shape": [4, 4, 4],
                    "faces": {
                        "i_lo": {"bc": "WALL", "partner": None},
                        "i_hi": {"bc": "FREE", "partner": None},
                        "j_lo": {"bc": "FREE", "partner": None},
                        "j_hi": {"bc": "SHARED", "partner": {"label": "b", "face": "j_lo"}},
                        "k_lo": {"bc": "FREE", "partner": None},
                        "k_hi": {"bc": "FREE", "partner": None},
                    },
                },
                {
                    "label": "b",
                    "shape": [4, 4, 4],
                    "faces": {
                        "i_lo": {"bc": "WALL", "partner": None},
                        "i_hi": {"bc": "FREE", "partner": None},
                        # b.j_lo points at b.j_hi instead of a.j_hi — broken.
                        "j_lo": {"bc": "SHARED", "partner": {"label": "b", "face": "j_hi"}},
                        "j_hi": {"bc": "SHARED", "partner": {"label": "b", "face": "j_lo"}},
                        "k_lo": {"bc": "FREE", "partner": None},
                        "k_hi": {"bc": "FREE", "partner": None},
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
    matching what shore.mesh.CubedSphereTopology produces (sub_q.j_lo/j_hi
    SHARED ring + sub_q.i_lo/i_hi DIRICHLET to caps)."""

    def face(bc: str, partner: tuple[str, str] | None) -> dict:
        return {
            "bc": bc,
            "partner": ({"label": partner[0], "face": partner[1]} if partner else None),
        }

    blocks = []
    # Sub-blocks
    for q in range(4):
        q_prev = (q - 1) % 4
        q_next = (q + 1) % 4
        blocks.append(
            {
                "label": f"sub{q}",
                "shape": [9, 5, 5],
                "faces": {
                    "i_lo": face("DIRICHLET", ("cap_south", "k_lo")),
                    "i_hi": face("DIRICHLET", ("cap_north", "k_lo")),
                    "j_lo": face("SHARED", (f"sub{q_prev}", "j_hi")),
                    "j_hi": face("SHARED", (f"sub{q_next}", "j_lo")),
                    "k_lo": face("WALL", None),
                    "k_hi": face("FREE", None),
                },
            }
        )
    # Cap blocks: SHORE leaves their lateral faces DIRICHLET-without-partner;
    # the writer reconstructs the symmetric link from _CAP_EQUATOR_PARTNER.
    for pole in ("north", "south"):
        blocks.append(
            {
                "label": f"cap_{pole}",
                "shape": [5, 5, 5],
                "faces": {
                    "i_lo": face("DIRICHLET", None),
                    "i_hi": face("DIRICHLET", None),
                    "j_lo": face("DIRICHLET", None),
                    "j_hi": face("DIRICHLET", None),
                    "k_lo": face("WALL", None),
                    "k_hi": face("FREE", None),
                },
            }
        )
    return {"version": 1, "blocks": blocks}


class TestCubedSphereCapEquator:
    def test_table_covers_all_8_pairs(self):
        """The hard-coded table has exactly 4 sub-blocks * 2 poles = 8 entries."""
        assert len(_CAP_EQUATOR_PARTNER) == 8
        assert all(0 <= q <= 3 for q, _ in _CAP_EQUATOR_PARTNER.keys())
        assert all(p in ("north", "south") for _, p in _CAP_EQUATOR_PARTNER.keys())

    def test_orientation_codes_decode_canonically(self):
        """Every code's three digits must each be in {1..6}."""
        for entry in _CAP_EQUATOR_PARTNER.values():
            _cap_face, _sub_face, sub_code, cap_code = entry
            for code in (sub_code, cap_code):
                li = code // 100
                lj = (code // 10) % 10
                lk = code % 10
                assert 1 <= li <= 6, f"{code}: li={li} out of range"
                assert 1 <= lj <= 6, f"{code}: lj={lj} out of range"
                assert 1 <= lk <= 6, f"{code}: lk={lk} out of range"

    def test_each_cap_face_used_exactly_once_per_pole(self):
        """Cap_north has 4 lateral faces, each must bind to exactly one sub."""
        for pole in ("north", "south"):
            cap_faces = [entry[0] for (q, p), entry in _CAP_EQUATOR_PARTNER.items() if p == pole]
            assert sorted(cap_faces) == ["i_hi", "i_lo", "j_hi", "j_lo"], (
                f"{pole} cap-face usage: {cap_faces}"
            )

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
        assert "        3         1         7 ! body" in text

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
