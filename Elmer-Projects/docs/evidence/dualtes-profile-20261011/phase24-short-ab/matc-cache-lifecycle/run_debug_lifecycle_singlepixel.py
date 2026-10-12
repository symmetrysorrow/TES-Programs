from pathlib import Path
import os,re,subprocess,time,shutil,hashlib,json
h=Path(__file__).resolve().parent;rt=Path(r'C:\Users\Nao\AppData\Local\Temp\phase24-cache-runtime');b=Path(r'C:\Users\Nao\AppData\Local\Temp\phase24-final-lifecycle-debug');mpi=Path(r'C:\Program Files\Microsoft MPI\Bin\mpiexec.exe')
mesh=rt/'workspace'/'work'/'meshes'/'mesh_hybrid_fullconf_h8';state=h.parent/'phase24_short_ab_singlepixel_vector_1.state';shutil.copy2(state,mesh/state.name)
raw=(h/'final_singlepixel_on.sif').read_text(encoding='utf-8-sig').replace('final_singlepixel_on','lifecycle_debug_singlepixel')
raw=re.sub(r'(?m)^\s*Phase24 Dynamic K Node Cache Verify\s*=.*$','  Phase24 Dynamic K Node Cache Verify = Logical True',raw)
sif=h/'lifecycle_debug_singlepixel.sif';sif.write_text(raw,encoding='utf-8',newline='\n');shutil.copy2(sif,rt/sif.name)
paths=[rt/'build'/'udf',b/'fem'/'src'/'modules',b/'fem'/'src',b/'matc'/'src',b/'fhutiter'/'src',b/'fem'/'src'/'binio',b/'mathlibs'/'src'/'arpack',b/'mathlibs'/'src'/'parpack',rt/'bin',rt/'share'/'elmersolver'/'lib',Path(r'G:\Elmer-Migration\msys64\ucrt64\bin'),Path(r'G:\Elmer-Migration\msys64\usr\bin'),Path(r'G:\Elmer-Migration\elmer\bin'),mpi.parent]
env=os.environ.copy();env.update(ELMER_HOME=str(rt),OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE',GFORTRAN_ERROR_BACKTRACE='1',GFORTRAN_SHOW_LOCUS='1',GFORTRAN_UNBUFFERED_ALL='y');env['PATH']=os.pathsep.join(map(str,paths))+os.pathsep+env.get('PATH','')
log=h/'debug_lifecycle_singlepixel.solver.log';t=time.perf_counter()
with log.open('w',encoding='utf-8',newline='\n') as o:p=subprocess.run([str(mpi),'-n','6',str(b/'fem'/'src'/'ElmerSolver_mpi.exe'),sif.name],cwd=rt,env=env,stdout=o,stderr=subprocess.STDOUT)
rec={'returncode':p.returncode,'wall_seconds':time.perf_counter()-t,'ranks':6,'steps':4,'cache':True,'verify':True,'fortran_checks':'-fcheck=all -finit-real=snan -finit-integer=-999 -ffpe-trap=invalid,zero,overflow','solver_sha256':hashlib.sha256((b/'fem'/'src'/'ElmerSolver_mpi.exe').read_bytes()).hexdigest(),'heat_solve_sha256':hashlib.sha256((b/'fem'/'src'/'modules'/'HeatSolve.dll').read_bytes()).hexdigest(),'log':log.name}
(h/'debug_lifecycle_run.json').write_text(json.dumps(rec,indent=2)+'\n',encoding='utf-8');print(json.dumps(rec));
if p.returncode:raise SystemExit(p.returncode)
