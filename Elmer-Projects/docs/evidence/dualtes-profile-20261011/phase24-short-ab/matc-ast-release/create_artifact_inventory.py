from pathlib import Path
import hashlib, json

HERE = Path(__file__).resolve().parent
excluded = {"artifact_inventory.json", "evidence_manifest.json"}
files = []
for path in sorted(p for p in HERE.rglob("*") if p.is_file() and p.name not in excluded and "__pycache__" not in p.parts):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    files.append({"path": path.relative_to(HERE).as_posix(), "bytes": path.stat().st_size,
                  "sha256": h.hexdigest()})
(HERE / "artifact_inventory.json").write_text(json.dumps({"schema":"phase24-ast-artifacts-v1", "files":files}, indent=2)+"\n",encoding="utf-8")
print(f"INVENTORY_PASS files={len(files)}")
