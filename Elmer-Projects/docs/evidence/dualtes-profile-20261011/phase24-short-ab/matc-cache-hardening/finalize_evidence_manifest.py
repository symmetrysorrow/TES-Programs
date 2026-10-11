from __future__ import annotations
import hashlib
import json
from pathlib import Path

root = Path(__file__).resolve().parent
files = []
for path in sorted(root.rglob('*')):
    if path.is_file() and path.name != 'evidence_manifest.json':
        files.append({
            'path': path.relative_to(root).as_posix(),
            'bytes': path.stat().st_size,
            'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
        })
manifest = {
    'repository': 'symmetrysorrow/TES-Programs',
    'branch': 'codex/dualtes-profile-evidence-20261011',
    'input_commit': 'e519ad74ce3f9d942075d3b0728d34b4c9309327',
    'decision': 'BLOCKED_NOT_SAFE_TO_ENABLE',
    'cache_default': 'OFF',
    'vector_assembly': 'OFF',
    'production_dll_or_sif_modified': False,
    'full_75ms_simulation_run': False,
    'runtime_negative_test_count': 1,
    'files': files,
}
(root / 'evidence_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
print(json.dumps({'files': len(files), 'decision': manifest['decision']}))
