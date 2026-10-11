from __future__ import annotations
import argparse
import hashlib
import json
import re
import subprocess
from pathlib import Path

BASE = '5a8de867068be0568f09af40fb90fee300dfbede'

def git(root: Path, *args: str) -> str:
    return subprocess.check_output(['git', '-C', str(root), *args], text=True).strip()

def filtered_blob(root: Path, relative: str, file: Path) -> str:
    return git(root, 'hash-object', '--path=' + relative, str(file))

parser = argparse.ArgumentParser()
parser.add_argument('--reference-source', type=Path, required=True)
parser.add_argument('--target-source', type=Path, required=True)
parser.add_argument('--replay-source', type=Path)
parser.add_argument('--output', type=Path)
args = parser.parse_args()

head = git(args.reference_source, 'rev-parse', 'HEAD')
if head != BASE:
    raise SystemExit(f'Expected reference {BASE}, got {head}')
patch = Path(__file__).with_name('01_phase24_source_from_release_26_1.patch')
text = patch.read_text(encoding='utf-8', errors='replace')
paths = re.findall(r'^diff --git a/(.*?) b/', text, re.M)
rows = []
for relative in paths:
    target = args.target_source / relative
    before = git(args.reference_source, 'rev-parse', f'{BASE}:{relative}')
    after = filtered_blob(args.reference_source, relative, target)
    row = {
        'path': relative,
        'before_git_blob': before,
        'target_git_blob': after,
        'target_raw_sha256': hashlib.sha256(target.read_bytes()).hexdigest(),
        'target_bytes': target.stat().st_size,
    }
    if args.replay_source:
        replayed = args.replay_source / relative
        actual = filtered_blob(args.reference_source, relative, replayed)
        row['replayed_git_blob'] = actual
        row['replay_matches_target'] = actual == after
    rows.append(row)

result = {
    'schema': 'phase24-source-reconstruction-hashes-v1',
    'reference': {'repository': 'ElmerCSC/elmerfem', 'tag': 'release-26.1', 'commit': BASE},
    'patch': {
        'path': patch.name,
        'sha256': hashlib.sha256(patch.read_bytes()).hexdigest(),
        'bytes': patch.stat().st_size,
        'apply_order': 1,
    },
    'files': rows,
    'replay_all_match': all(x.get('replay_matches_target', True) for x in rows),
}
encoded = json.dumps(result, indent=2) + '\n'
if args.output:
    args.output.write_text(encoded, encoding='utf-8')
print(encoded, end='')
