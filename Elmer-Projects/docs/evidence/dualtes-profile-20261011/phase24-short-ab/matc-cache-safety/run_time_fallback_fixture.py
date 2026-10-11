#!/usr/bin/env python3
"""Run one isolated explicit-Time fallback fixture using the prior Release DLL."""
from __future__ import annotations
import datetime, hashlib, json, os, re, shutil, subprocess, time
from pathlib import Path
HERE=Path(__file__).resolve().parent
EVIDENCE=HERE.parent
INPUT=EVIDENCE/'matc-cache-release/dualtes_on_1.sif'
STATE_DIR=EVIDENCE/'bulk-matc-cache'
RUNTIME=Path(r'C:\Users\Nao\AppData\Local\Temp\phase24-cache-runtime')
BUILD=Path(r'C:\Users\Nao\AppData\Local\Temp\phase24-release-build')
MESH_TAG='mesh_toml_95333801e1a3'
MESH=RUNTIME/'workspace/work/meshes'/MESH_TAG
MPI=Path(r'C:\Program Files\Microsoft MPI\Bin\mpiexec.exe')
SOLVER=BUILD/'fem/src/ElmerSolver_mpi.exe'
LABEL='matc_safety_time_dependency'

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    for p in (INPUT, MPI, SOLVER, MESH/'partitioning.6'):
        if not p.exists(): raise SystemExit(f'missing isolated input: {p}')
    text=INPUT.read_text(encoding='utf-8-sig')
    text=text.replace('Heat Conductivity = Variable Temperature\n    Real MATC', 'Heat Conductivity = Variable Time\n    Real MATC', 1)
    if 'Heat Conductivity = Variable Time' not in text: raise SystemExit('time fixture edit failed')
    text=text.replace('dualtes_release_on_1', LABEL)
    text=text.replace('dualtes_cache_verify_1_series.csv', LABEL+'_series.csv')
    text=text.replace('dualtes_cache_verify_1_iterations.csv', LABEL+'_iterations.csv')
    sif=RUNTIME/(LABEL+'.sif')
    sif.write_text(text,encoding='utf-8')
    for i in (1,2):
        src=STATE_DIR/f'phase24_short_ab_dualtes_vector_{i}.state'
        if not src.exists(): raise SystemExit(f'missing same-state restart input: {src}')
        shutil.copy2(src,MESH/src.name)
    paths=[RUNTIME/'build/udf',BUILD/'fem/src/modules',BUILD/'fem/src',BUILD/'matc/src',BUILD/'fhutiter/src',BUILD/'fem/src/binio',BUILD/'mathlibs/src/arpack',BUILD/'mathlibs/src/parpack',RUNTIME/'bin',RUNTIME/'share/elmersolver/lib',Path(r'G:\Elmer-Migration\msys64\ucrt64\bin'),Path(r'G:\Elmer-Migration\msys64\usr\bin'),Path(r'G:\Elmer-Migration\elmer\bin'),MPI.parent]
    env=os.environ.copy(); env.update(ELMER_HOME=str(RUNTIME),OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE',PHASE24_HYPRE_REPORT='1'); env['PATH']=os.pathsep.join([str(p) for p in paths]+[env.get('PATH','')])
    log=HERE/(LABEL+'.solver.log')
    start=datetime.datetime.now().astimezone().isoformat(timespec='seconds'); t0=time.perf_counter()
    with log.open('w',encoding='utf-8',newline='\n') as f:
        proc=subprocess.run([str(MPI),'-n','6',str(SOLVER),sif.name],cwd=RUNTIME,env=env,stdout=f,stderr=subprocess.STDOUT,check=False)
    wall=time.perf_counter()-t0; raw=log.read_text(encoding='utf-8',errors='replace')
    matches=[line.strip() for line in raw.splitlines() if 'dynamic-k generic cache: no safe temperature-only MATC prisms' in line]
    result={'case':LABEL,'fixture':'Material 7 test-copy changes declared conductivity dependency Temperature -> Time; cache option remains True','runtime_isolated':True,'ranks':6,'steps':4,'returncode':proc.returncode,'started':start,'wall_seconds':wall,'sif_sha256':sha(sif),'solver_sha256':sha(SOLVER),'HeatSolve_sha256':sha(BUILD/'fem/src/modules/HeatSolve.dll'),'log':log.name,'log_sha256':sha(log),'fallback_messages':matches,'fallback_observed':len(matches)>0,'production_files_modified':False,'interpretation':'fallback behavior only; altered-material outputs are not numerical or performance evidence'}
    (HERE/(LABEL+'_result.json')).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    shutil.copy2(sif,HERE/(LABEL+'.sif'))
    print(json.dumps(result,indent=2))
    return 0 if proc.returncode==0 and matches else 1
if __name__=='__main__': raise SystemExit(main())
