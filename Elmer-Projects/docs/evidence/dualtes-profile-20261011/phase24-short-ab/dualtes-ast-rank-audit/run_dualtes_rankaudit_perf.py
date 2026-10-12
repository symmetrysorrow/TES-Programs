#!/usr/bin/env python3
"""Repeat 4-step DualTES cache OFF/ON using the AST Release build, verify disabled."""
from __future__ import annotations
import hashlib, json, os, re, shutil, subprocess, time
from pathlib import Path
HERE = Path(__file__).resolve().parent
AST = HERE.parent / "matc-ast-release"
RUNTIME = Path(r"C:\Users\Nao\AppData\Local\Temp\phase24-cache-runtime")
BUILD = Path(r"C:\Users\Nao\AppData\Local\Temp\phase24-repro-release-build")
SOLVER = BUILD / "fem" / "src" / "ElmerSolver_mpi.exe"
MPI = Path(r"C:\Program Files\Microsoft MPI\Bin\mpiexec.exe")
MESH = RUNTIME / "workspace" / "work" / "meshes" / "mesh_toml_95333801e1a3"
RESULT_DIR = RUNTIME / "mesh_toml_95333801e1a3"
STATE_DIR = AST.parent / "bulk-matc-cache"
PATHS = [RUNTIME / "build" / "udf", BUILD / "fem" / "src" / "modules", BUILD / "fem" / "src",
 BUILD / "matc" / "src", BUILD / "fhutiter" / "src", BUILD / "fem" / "src" / "binio",
 BUILD / "mathlibs" / "src" / "arpack", BUILD / "mathlibs" / "src" / "parpack", RUNTIME / "bin",
 RUNTIME / "share" / "elmersolver" / "lib", Path(r"G:\Elmer-Migration\msys64\ucrt64\bin"),
 Path(r"G:\Elmer-Migration\msys64\usr\bin"), Path(r"G:\Elmer-Migration\elmer\bin"), MPI.parent]
def sha(p: Path) -> str:
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''): h.update(b)
 return h.hexdigest()
def main():
 env=os.environ.copy(); env.update(ELMER_HOME=str(RUNTIME),OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE',PHASE24_HYPRE_REPORT='1')
 env['PATH']=os.pathsep.join(map(str,PATHS))+os.pathsep+env.get('PATH','')
 for i in (1,2): shutil.copy2(STATE_DIR/f'phase24_short_ab_dualtes_vector_{i}.state', MESH)
 order=[('off',1),('on',1),('on',2),('off',2),('off',3),('on',3)]
 records=[]
 for mode,rep in order:
  raw=(AST/f'ast_dualtes_{mode}.sif').read_text(encoding='utf-8-sig')
  prefix=f'rankaudit_perf_{mode}_{rep}'
  raw=raw.replace(f'ast_dualtes_{mode}',prefix)
  if 'Phase24 Dynamic K Node Cache Verify = Logical True' in raw:
   raise SystemExit('verification must be disabled during timing runs')
  sif=HERE/f'{prefix}.sif'; sif.write_text(raw,encoding='utf-8',newline='\n'); shutil.copy2(sif,RUNTIME/sif.name)
  log=HERE/f'{prefix}.solver.log'; err=HERE/f'{prefix}.stderr.log'
  t=time.perf_counter()
  with log.open('w',encoding='utf-8',newline='\n') as out, err.open('w',encoding='utf-8',newline='\n') as eout:
   p=subprocess.run([str(MPI),'-n','6',str(SOLVER),sif.name],cwd=RUNTIME,env=env,stdout=out,stderr=eout,check=False)
  wall=time.perf_counter()-t
  txt=log.read_text(encoding='utf-8',errors='replace')
  norms=[float(x) for x in re.findall(r'HeatSolve:\s+Result Norm\s*:\s*([0-9.Ee+-]+)',txt)]
  hypre=[[int(n),float(r)] for n,r in re.findall(r'Phase24 PCG iterations=(\d+) final_rel_res=([0-9.Ee+-]+)',txt)]
  outputs=[]
  for tes in (1,2):
   for kind in ('series','iterations'):
    src=RUNTIME/f'{prefix}_{tes}_{kind}.csv'
    if src.is_file():
     dst=HERE/src.name; shutil.copy2(src,dst); outputs.append({'file':dst.name,'sha256':sha(dst)})
  shards=[RESULT_DIR/f'{prefix}.result']+sorted(RESULT_DIR.glob(f'{prefix}.result.*'))
  outputs += [{'file':x.name,'bytes':x.stat().st_size,'sha256':sha(x)} for x in shards if x.is_file()]
  rec={'case':'DualTES','mode':mode,'repeat':rep,'returncode':p.returncode,'wall_seconds':wall,'ranks':6,'steps':4,'vector_assembly':False,'cache_verify':False,'cache_enabled':mode=='on','solver_sha256':sha(SOLVER),'heat_solve_sha256':sha(BUILD/'fem'/'src'/'modules'/'HeatSolve.dll'),'matc_sha256':sha(BUILD/'matc'/'src'/'libmatc.dll'),'sif_sha256':sha(sif),'stdout_log':log.name,'stderr_log':err.name,'result_norms':norms,'hypre_pcg':hypre,'outputs':outputs}
  records.append(rec); (HERE/'performance_repeats.json').write_text(json.dumps({'build':'AST-purity reconstructed Release; same binary for all six runs','order':order,'records':records},indent=2)+'\n',encoding='utf-8')
  print(json.dumps(rec,separators=(',',':')),flush=True)
  if p.returncode: raise SystemExit(f'failed {prefix}: {p.returncode}')
if __name__=='__main__': main()
