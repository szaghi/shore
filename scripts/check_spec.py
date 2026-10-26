"""Check that the frozen acceptance specs are unmodified (ADR 0005).

Every path listed in ``tests/spec/SHA256SUMS`` must match one of its listed
SHA-256 hashes (a spec file may list a ``frozen`` and an ``accepted`` hash:
the only permitted edit is removing its ``pytestmark`` xfail line).

Usage: python scripts/check_spec.py   (exit 1 on any mismatch)
"""

from __future__ import annotations

import hashlib
import sys
from collections import defaultdict
from pathlib import Path


def main() -> int:
    root = Path(__file__).resolve().parent.parent
    allowed: dict[str, set[str]] = defaultdict(set)
    for line in (root / "tests/spec/SHA256SUMS").read_text().splitlines():
        if line.strip() and not line.startswith("#"):
            digest, path = line.split()[:2]
            allowed[path].add(digest)
    bad = 0
    for path, digests in sorted(allowed.items()):
        actual = hashlib.sha256((root / path).read_bytes()).hexdigest()
        ok = actual in digests
        bad += not ok
        print(f"{'ok      ' if ok else 'CHANGED '} {path}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
