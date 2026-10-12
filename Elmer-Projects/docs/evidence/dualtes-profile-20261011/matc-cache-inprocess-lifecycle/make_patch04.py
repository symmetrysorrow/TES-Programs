#!/usr/bin/env python3
"""Regenerate patch 04 from a patched-03 tree and a diagnostic source tree."""
from __future__ import annotations

import argparse
import difflib
import hashlib
import json
from pathlib import Path

FILES = ("fem/src/ElmerSolver.F90", "fem/src/modules/HeatSolve.F90")


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pre04", required=True, type=Path, help="Elmer 26.1 source with patches 01-03 applied")
    ap.add_argument("--diagnostic", required=True, type=Path, help="test-only diagnostic source tree")
    ap.add_argument("--out", type=Path, default=Path(__file__).resolve().parent / "04_matc_cache_meshmode_diagnostics.patch")
    args = ap.parse_args()
    chunks = []
    changed = []
    for rel in FILES:
        before, after = args.pre04 / rel, args.diagnostic / rel
        if before.read_bytes() == after.read_bytes():
            continue
        body = "".join(difflib.unified_diff(
            before.read_text(encoding="utf-8").splitlines(keepends=True),
            after.read_text(encoding="utf-8").splitlines(keepends=True),
            fromfile="a/" + rel, tofile="b/" + rel, n=3))
        chunks.append("diff --git a/" + rel + " b/" + rel + "\n" + body)
        changed.append({"path": rel, "pre04_sha256": sha(before), "diagnostic_sha256": sha(after)})
    if len(changed) != len(FILES):
        raise SystemExit(f"expected edits in {FILES}; found {[x['path'] for x in changed]}")
    args.out.write_text("".join(chunks), encoding="utf-8", newline="")
    manifest = {"schema": "phase24-matc-cache-meshmode-diagnostics-v1",
                "patch": args.out.name, "patch_sha256": sha(args.out), "files": changed}
    (args.out.parent / "patch04_source_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
