"""Verify the native Stage 11 source patch and recorded runtime identity."""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import tempfile
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run(*args: str, cwd: Path | None = None) -> str:
    return subprocess.check_output(args, cwd=cwd, text=True).strip()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=Path("artifacts/phase24_stage11_native_provenance_manifest.json"))
    parser.add_argument("--source", type=Path)
    parser.add_argument("--skip-clean-patch-check", action="store_true")
    args = parser.parse_args()

    manifest_path = args.manifest.resolve()
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    root = manifest_path.parents[1]
    source = (args.source or root / manifest["native"]["source_relative_to_repo"]).resolve()
    native = manifest["native"]
    if not (source / ".git").exists():
        raise SystemExit(f"native Git checkout not found: {source}")

    base = native["base_commit_sha"]
    actual_head = run("git", "rev-parse", base, cwd=source)
    if actual_head != base:
        raise SystemExit(f"base commit mismatch: expected {base}, resolved {actual_head}")

    patch = (root / native["patch_relative_to_repo"]).resolve()
    if sha256(patch) != native["patch_sha256"]:
        raise SystemExit("native patch SHA256 mismatch")

    actual_patch = subprocess.check_output(("git", "diff", "--binary", base), cwd=source)
    actual_patch_sha = hashlib.sha256(actual_patch).hexdigest()
    if actual_patch_sha != native["worktree_diff_sha256"]:
        raise SystemExit("native worktree diff SHA256 mismatch")

    for relative, expected in native["patched_source_sha256"].items():
        actual = sha256(source / relative)
        if actual != expected:
            raise SystemExit(f"source SHA256 mismatch: {relative}")

    if not args.skip_clean_patch_check:
        with tempfile.TemporaryDirectory(prefix="phase24-native-verify-") as temporary:
            clone = Path(temporary) / "source"
            run("git", "clone", "--no-hardlinks", "--quiet", str(source), str(clone))
            run("git", "checkout", "--quiet", base, cwd=clone)
            subprocess.run(("git", "apply", "--check", str(patch)), cwd=clone, check=True)

    for relative, expected in manifest["runtime_artifacts"].items():
        path = (root / relative).resolve()
        if not path.is_file():
            raise SystemExit(f"runtime artifact not found: {path}")
        if sha256(path) != expected:
            raise SystemExit(f"runtime artifact SHA256 mismatch: {relative}")

    print("PASS: native base, complete patch, source hashes, and runtime artifacts verified")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
