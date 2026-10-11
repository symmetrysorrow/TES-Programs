#!/usr/bin/env python3
from __future__ import annotations
import hashlib, json
from pathlib import Path
here=Path(__file__).resolve().parent
root=here.parents[5]
source=Path(r'C:\Users\Nao\AppData\Local\Temp\phase24-debug-src\fem\src\modules\HeatSolve.F90')
build=Path(r'C:\Users\Nao\AppData\Local\Temp\phase24-release-build')
def h(p): return hashlib.sha256(p.read_bytes()).hexdigest()
files=[]
for p in sorted(here.rglob('*')):
    if p.is_file() and p.name != 'safety_manifest.json':
        files.append({'path':p.relative_to(here).as_posix(),'bytes':p.stat().st_size,'sha256':h(p)})
manifest={
 'repository':'symmetrysorrow/TES-Programs',
 'branch':'codex/dualtes-profile-evidence-20261011',
 'base_commit':'329b27b9a196175efe48a91bc6e3083d91f3233a',
 'decision':'BLOCKED_NOT_SAFE_TO_ENABLE',
 'production_files_modified':False,
 'runtime_negative_tests_executed':False,
 'isolated_source':{'path':str(source),'git_managed':False,'sha256':h(source) if source.exists() else None},
 'previous_release_build':{
  'solver':str(build/'fem/src/ElmerSolver_mpi.exe'),
  'solver_sha256':h(build/'fem/src/ElmerSolver_mpi.exe') if (build/'fem/src/ElmerSolver_mpi.exe').exists() else None,
  'HeatSolve':str(build/'fem/src/modules/HeatSolve.dll'),
  'HeatSolve_sha256':h(build/'fem/src/modules/HeatSolve.dll') if (build/'fem/src/modules/HeatSolve.dll').exists() else None,
  'purpose':'hash-only identification of pre-existing Release validation artifacts; no binary copied or modified'},
 'files':files}
(here/'safety_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'file_entries':len(files),'source_hash':manifest['isolated_source']['sha256'],'decision':manifest['decision']}))
