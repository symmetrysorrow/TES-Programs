#!/usr/bin/env python3
from __future__ import annotations

import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parent
MANIFEST = ROOT / "SHA256SUMS"
rows = []
for path in sorted(p for p in ROOT.rglob("*") if p.is_file() and p != MANIFEST):
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    rows.append(f"{digest}  {path.relative_to(ROOT).as_posix()}")
MANIFEST.write_text("\n".join(rows) + "\n", encoding="ascii", newline="\n")
print(f"wrote {MANIFEST} with {len(rows)} entries")
