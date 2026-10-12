from pathlib import Path
import os,re,subprocess,time,json,hashlib,shutil
here=Path(__file__).resolve().parent
rt=Path(r'C:\Users\Nao\AppData\Local\Temp\phase24-cache-runtime')
build=Path(r'C:\Users\Nao\AppData\Local\Temp\phase24-03-release-build')
mpi=Path(r'C:\Program Files\Microsoft MPI\Bin\mpiexec.exe')
mesh=rt/'workspace'/'work'/'meshes'/'mesh_toml_95333801e1a3'
for i in (1,2): shutil.copy2(here.parent/'bulk-matc-cache'/f'phase24_short_ab_dualtes_vector_{i}.state',mesh/f'phase24_short_ab_dualtes_vector_{i}.state')
s=(here/'final_dualtes_on.sif').read_text(encoding='utf-8-sig')
s=s.replace('final_dualtes_on','lifecycle_dualtes_verify')
s=re.sub(r'(?m)^\s*Phase24 Dynamic K Node Cache Verify\s*=.*$', '  Phase24 Dynamic K Node Cache Verify = Logical True',s)
p=here/'lifecycle_dualtes_verify.sif'; p.write_text(s,encoding='utf-8',newline='\n'); shutil.copy2(p,rt/p.name)
env=os.environ.copy(); env.update(ELMER_HOME=str(rt),OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE',PHASE24_HYPRE_REPORT='1')
paths=[rt/'build'/'udf',build/'fem'/'src'/'modules',build/'fem'/'src',build/'matc'/'src',build/'fhutiter'/'src',build/'fem'/'src'/'binio',build/'mathlibs'/'src'/'arpack',build/'mathlibs'/'src'/'parpack',rt/'bin',rt/'share'/'elmersolver'/'lib',Path(r'G:\Elmer-Migration\msys64\ucrt64\bin'),Path(r'G:\Elmer-Migration\msys64\usr\bin'),Path(r'G:\Elmer-Migration\elmer\bin'),mpi.parent]
env['PATH']=os.pathsep.join(map(str,paths))+os.pathsep+env.get('PATH','')
log=here/'lifecycle_dualtes_verify.solver.log'; t=time.perf_counter()
with log.open('w',encoding='utf-8',newline='\n') as out: proc=subprocess.run([str(mpi),'-n','6',str(build/'fem'/'src'/'ElmerSolver_mpi.exe'),p.name],cwd=rt,env=env,stdout=out,stderr=subprocess.STDOUT)
wall=time.perf_counter()-t
sha=lambda x:hashlib.sha256(x.read_bytes()).hexdigest()
rec={'returncode':proc.returncode,'wall_seconds':wall,'ranks':6,'steps':4,'cache':True,'verify':True,'solver_sha256':sha(build/'fem'/'src'/'ElmerSolver_mpi.exe'),'heat_solve_sha256':sha(build/'fem'/'src'/'modules'/'HeatSolve.dll'),'sif_sha256':sha(p),'log':log.name}
(here/'lifecycle_verify_run.json').write_text(json.dumps(rec,indent=2)+'\n',encoding='utf-8');print(json.dumps(rec));
if proc.returncode: raise SystemExit(proc.returncode)
