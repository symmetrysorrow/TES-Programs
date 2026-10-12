#!/usr/bin/env python3
from __future__ import annotations
import difflib, hashlib, json, subprocess
from pathlib import Path
HERE=Path(__file__).resolve().parent
PRE=Path(r'C:\Users\Nao\AppData\Local\Temp\phase24-source-replay-worktree\fem\src\modules\HeatSolve.F90')
POST=Path(r'C:\Users\Nao\AppData\Local\Temp\phase24-03-source\fem\src\modules\HeatSolve.F90')
REL='fem/src/modules/HeatSolve.F90'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def blob(p): return subprocess.check_output(['git','-C',r'C:\Users\Nao\AppData\Local\Temp\phase24-03-source','hash-object',f'--path={REL}',str(p)],text=True).strip()
pre=PRE.read_text(encoding='utf-8').splitlines(keepends=True)
post=POST.read_text(encoding='utf-8').splitlines(keepends=True)
body=''.join(difflib.unified_diff(pre,post,fromfile='a/'+REL,tofile='b/'+REL,n=3))
if not body.startswith('--- a/'+REL): raise SystemExit('expected unified diff was not generated')
before=blob(PRE); after=blob(POST)
patch='diff --git a/'+REL+' b/'+REL+'\nindex '+before+'..'+after+' 100644\n'+body
out=HERE/'03_matc_cache_lifecycle.patch'; out.write_text(patch,encoding='utf-8',newline='')
meta={'schema':'phase24-lifecycle-patch-v1','patch':{'path':out.name,'sha256':sha(out),'order':3},'files':[{'path':REL,'before_git_blob':before,'after_git_blob':after,'before_raw_sha256':sha(PRE),'after_raw_sha256':sha(POST)}]}
(HERE/'03_source_hash_manifest.json').write_text(json.dumps(meta,indent=2)+'\n',encoding='utf-8')
print(json.dumps(meta,indent=2))
