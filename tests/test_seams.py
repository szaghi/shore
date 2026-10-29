"""Unit tests of the seam alias groups and their sync (ADR 0040)."""

from __future__ import annotations

import numpy as np

from shore.mesh.seams import Aliases


def _flat_caps_copies():
    """The flat-caps primitive's 5 blocks, each nodes array an independent copy."""
    from shore.project import Project

    p = Project()
    p.build_primitive("flat_caps", r_out=3.0, z0=-1.0, z1=4.0, ni=5, nc=6, nk=5, inner_frac=0.4)
    blocks = p.mesh.blocks
    for b in blocks:
        b.nodes = b.nodes.copy()
    return blocks


class TestAliases:
    def test_groups_join_coincident_nodes(self):
        blocks = _flat_caps_copies()
        aliases = Aliases(blocks, 1e-9)
        assert aliases.groups
        sizes = {len(g) for g in aliases.groups}
        assert min(sizes) == 2 and max(sizes) >= 3  # block corners join three blocks
        for group in aliases.groups:  # every member is the same point
            points = [blocks[b].nodes[idx][0] for b, idx in aliases._locate(group) for _ in [0]]
            assert np.ptp(np.array(points), axis=0).max() <= 1e-9

    def test_sync_last_writer_wins(self):
        blocks = _flat_caps_copies()
        aliases = Aliases(blocks, 1e-9)
        group = next(g for g in aliases.groups if len(g) >= 3)
        located = aliases._locate(group)
        changed = [np.zeros(b.nodes.shape[:3], bool) for b in blocks]
        for value, (b, idx) in zip((10.0, 20.0), located[:2], strict=False):
            blocks[b].nodes[idx] = value
            changed[b][idx] = True
        aliases.sync(changed)
        winner = 20.0 if len(located) >= 2 else 10.0
        for b, idx in located:
            assert np.all(blocks[b].nodes[idx] == winner)

    def test_unchanged_groups_untouched(self):
        blocks = _flat_caps_copies()
        before = [b.nodes.copy() for b in blocks]
        aliases = Aliases(blocks, 1e-9)
        aliases.sync([np.zeros(b.nodes.shape[:3], bool) for b in blocks])
        for b, old in zip(blocks, before, strict=True):
            assert np.array_equal(b.nodes, old)
