#!/usr/bin/env python3
from pathlib import Path

ROOT = Path(__file__).resolve().parent
TEXT_SUFFIXES = {".md", ".json", ".py", ".ps1", ".sif", ".patch", ".txt", ".log",
                 ".msh", ".header", ".nodes", ".elements", ".boundary", ".gitattributes"}
changed = 0
for path in ROOT.rglob("*"):
    if not path.is_file() or (path.suffix.lower() not in TEXT_SUFFIXES and path.name != "SHA256SUMS"):
        continue
    raw = path.read_bytes()
    text = raw.decode("utf-8-sig") if raw else ""
    text = text.replace("\r\n", "\n").replace("\r", "\n")
    if path.suffix.lower() == ".sif":
        text = text.rstrip("\n") + ("\n" if text else "")
    path.write_bytes(text.encode("utf-8"))
    changed += 1
print(f"normalized {changed} text files to UTF-8/LF")
