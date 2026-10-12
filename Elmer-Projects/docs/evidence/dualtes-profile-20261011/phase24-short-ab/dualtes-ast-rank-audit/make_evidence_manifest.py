#!/usr/bin/env python3
"""Hash tracked evidence files and record the isolated Release runtime."""
from __future__ import annotations
import hashlib, json, subprocess
from pathlib import Path

HERE=Path(__file__).resolve().parent
BUILD=Path(r'C:\Users\Nao\AppData\Local\Temp\phase24-repro-release-build')
RUNTIME=Path(r'C:\Users\Nao\AppData\Local\Temp\phase24-cache-runtime')
def sha(p: Path) -> str:
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1<<20),b''): h.update(block)
    return h.hexdigest()
def git(*args: str) -> str:
    return subprocess.check_output(['git',*args],cwd=HERE.parents[5],text=True).strip()
files=[]
for p in sorted(HERE.rglob('*')):
    if p.is_file() and p.name!='evidence_manifest.json':
        files.append({'path':p.relative_to(HERE).as_posix(),'bytes':p.stat().st_size,'sha256':sha(p)})
manifest={
 'schema':'phase24-dualtes-ast-rank-audit-v1',
 'repository':'symmetrysorrow/TES-Programs',
 'branch':'codex/dualtes-profile-evidence-20261011',
 'requested_base_commit':'5f341c76bb886c87401e29acc4247d17bdc1f3c5',
 'recorded_from_commit':git('rev-parse','HEAD'),
 'runtime':{'configuration':'Release','mpi':'Microsoft MPI','ranks':6,'steps':4,'vector_assembly':False,'verify_run_cache':True,'performance_run_cache_verify':False,'default_cache_setting':'OFF','full_75ms_run':False},
 'binaries':{str(p.relative_to(BUILD)).replace('\\','/'):sha(p) for p in [BUILD/'fem/src/ElmerSolver_mpi.exe',BUILD/'fem/src/modules/HeatSolve.dll',BUILD/'matc/src/libmatc.dll']},
 'ast_test':{'log':'matc_ast_purity_test.log','passed':'RESULT failures=0' in (HERE/'matc_ast_purity_test.log').read_text(encoding='utf-8',errors='replace')},
 'diagnostic_returncode':0 if 'MAIN: Time: 4/4:' in (HERE/'dualtes_rankaudit_verify.solver.log').read_text(encoding='utf-8',errors='replace') else None,
 'performance_run_returncodes':[x['returncode'] for x in json.loads((HERE/'performance_repeats.json').read_text(encoding='utf-8'))['records']],
 'files':files}
(HERE/'evidence_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
print(f"hashed_files={len(files)} base={manifest['requested_base_commit']} current={manifest['recorded_from_commit']}")
