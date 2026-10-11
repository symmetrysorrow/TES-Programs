from pathlib import Path
import os, shutil, subprocess, time, hashlib, json, re
HERE=Path(__file__).resolve().parent
PARENT=HERE.parent/'bulk-matc-cache'
RUNTIME=Path(r'C:\Users\Nao\AppData\Local\Temp\phase24-cache-runtime')
BUILD=Path(r'C:\Users\Nao\AppData\Local\Temp\phase24-release-build')
MESH=RUNTIME/'workspace'/'work'/'meshes'/'mesh_toml_95333801e1a3'
MPI=Path(r'C:\Program Files\Microsoft MPI\Bin\mpiexec.exe')
SOLVER=BUILD/'fem'/'src'/'ElmerSolver_mpi.exe'
STATES=[f'phase24_short_ab_dualtes_vector_{i}.state' for i in (1,2)]
PATHS=[RUNTIME/'build'/'udf',BUILD/'fem'/'src'/'modules',BUILD/'fem'/'src',BUILD/'matc'/'src',BUILD/'fhutiter'/'src',BUILD/'fem'/'src'/'binio',BUILD/'mathlibs'/'src'/'arpack',BUILD/'mathlibs'/'src'/'parpack',RUNTIME/'bin',RUNTIME/'share'/'elmersolver'/'lib',Path(r'G:\Elmer-Migration\msys64\ucrt64\bin'),Path(r'G:\Elmer-Migration\msys64\usr\bin'),Path(r'G:\Elmer-Migration\elmer\bin'),MPI.parent]
env=os.environ.copy(); env.update(ELMER_HOME=str(RUNTIME),OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE',PHASE24_HYPRE_REPORT='1'); env['PATH']=os.pathsep.join(map(str,PATHS))+os.pathsep+env.get('PATH','')
records=[]
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''): h.update(b)
 return h.hexdigest()
for mode in ('off','on'):
 for step in range(1,5):
  label=f'{mode}_ts{step}'
  prefix=RUNTIME/f'release_capture_{label}'
  for nlit in range(1,121): (prefix/f'ts{step:04d}_nl{nlit:04d}').mkdir(parents=True,exist_ok=True)
  for state in STATES: shutil.copy2(PARENT/state,MESH/state)
  src=PARENT/f'dualtes_field_{mode}_6rank.sif'; s=src.read_text(encoding='utf-8-sig')
  s=s.replace('Phase24 Outer HeatSolve Capture = Logical False','Phase24 Outer HeatSolve Capture = Logical True')
  s=re.sub(r'(?m)^\s*Phase24 Outer HeatSolve Capture Timestep\s*=.*$',f'  Phase24 Outer HeatSolve Capture Timestep = Integer {step}',s)
  s=re.sub(r'(?m)^\s*Phase24 Outer HeatSolve Capture Max Iterations\s*=.*$','  Phase24 Outer HeatSolve Capture Max Iterations = Integer 120',s)
  s=re.sub(r'(?m)^\s*Phase24 Outer HeatSolve Capture Prefix\s*=.*$',f'  Phase24 Outer HeatSolve Capture Prefix = String "release_capture_{label}"',s)
  s=s.replace(f'dualtes_field_{mode}_',f'dualtes_release_capture_{label}_')
  s=s.replace(f'  Output File = dualtes_field_{mode}.result',f'  Output File = dualtes_release_capture_{label}.result')
  sif=HERE/f'{label}.sif'; sif.write_text(s,encoding='utf-8',newline=''); shutil.copy2(sif,RUNTIME/sif.name)
  log=HERE/f'{label}.solver.log'; t=time.perf_counter()
  with log.open('w',encoding='utf-8',newline='\n') as f: p=subprocess.run([str(MPI),'-n','6',str(SOLVER),sif.name],cwd=RUNTIME,env=env,stdout=f,stderr=subprocess.STDOUT)
  prefix=RUNTIME/f'release_capture_{label}'; files=sorted(x for x in prefix.rglob('*') if x.is_file()) if prefix.exists() else []
  records.append({'label':label,'mode':mode,'timestep':step,'returncode':p.returncode,'wall_seconds':time.perf_counter()-t,'log':log.name,'capture_files':[{'relative':x.relative_to(prefix).as_posix(),'bytes':x.stat().st_size,'sha256':sha(x)} for x in files]})
  (HERE/'release_all_assembly_capture.json').write_text(json.dumps({'compiler':'GNU Fortran/C/C++ 16.1 Release -O3','ranks':6,'vector_assembly':False,'capture':'before/after DefaultSolve; max nonlinear iterations 120','runs':records},indent=2)+'\n',encoding='utf-8')
  if p.returncode: raise SystemExit(f'{label} failed rc={p.returncode}')
