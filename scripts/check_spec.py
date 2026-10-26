"""Check that the frozen acceptance specs are unmodified (ADR 0005).

Every path listed in ``tests/spec/SHA256SUMS`` must match one of its listed
SHA-256 hashes. Spec modules (``tests/spec/test_*.py``) are hashed with their
``@PENDING`` lines removed: deleting such a line once a case passes is the
only permitted edit.

Usage: python scripts/check_spec.py   (exit 1 on any mismatch)
"""

from __future__ import annotations

import hashlib
import sys
from collections import defaultdict
from pathlib import Path


def _normalised(path: Path) -> bytes:
    """File bytes; spec modules without their ``@PENDING`` marker lines."""
    data = path.read_bytes()
    if path.parent.name == "spec" and path.name.startswith("test_"):
        lines = data.decode().splitlines(keepends=True)
        data = "".join(ln for ln in lines if ln.strip() != "@PENDING").encode()
    return data


def main() -> int:
    root = Path(__file__).resolve().parent.parent
    allowed: dict[str, set[str]] = defaultdict(set)
    for line in (root / "tests/spec/SHA256SUMS").read_text().splitlines():
        if line.strip() and not line.startswith("#"):
            digest, path = line.split()[:2]
            allowed[path].add(digest)
    bad = 0
    for path, digests in sorted(allowed.items()):
        actual = hashlib.sha256(_normalised(root / path)).hexdigest()
        ok = actual in digests
        bad += not ok
        print(f"{'ok      ' if ok else 'CHANGED '} {path}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
