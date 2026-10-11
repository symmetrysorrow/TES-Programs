from pathlib import Path
import os,shutil,subprocess,time,json,hashlib
HERE=Path(__file__).resolve().parent; R=Path(r'C:\Users\Nao\AppData\Local\Temp\phase24-cache-runtime'); B=Path(r'C:\Users\Nao\AppData\Local\Temp\phase24-release-build'); M=R/'workspace'/'work'/'meshes'/'mesh_hybrid_fullconf_h8'; initial=HERE.parent/'phase24_short_ab_singlepixel_vector_1.state'; mpi=Path(r'C:\Program Files\Microsoft MPI\Bin\mpiexec.exe')
paths=[R/'build'/'udf',B/'fem'/'src'/'modules',B/'fem'/'src',B/'matc'/'src',B/'fhutiter'/'src',B/'fem'/'src'/'binio',B/'mathlibs'/'src'/'arpack',B/'mathlibs'/'src'/'parpack',R/'bin',R/'share'/'elmersolver'/'lib',Path(r'G:\Elmer-Migration\msys64\ucrt64\bin'),Path(r'G:\Elmer-Migration\msys64\usr\bin'),Path(r'G:\Elmer-Migration\elmer\bin'),mpi.parent]
env=os.environ.copy();env.update(ELMER_HOME=str(R),OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE',PHASE24_HYPRE_REPORT='1');env['PATH']=os.pathsep.join(map(str,paths))+os.pathsep+env.get('PATH','')
s=(HERE/'singlepixel_release_on_6rank.sif').read_text().replace('Dynamic K Node Cache Verify = Logical True','Dynamic K Node Cache Verify = Logical False').replace('singlepixel_release_on','singlepixel_release_on_perf')
sif=HERE/'singlepixel_on_perf_6rank.sif';sif.write_text(s,encoding='utf-8',newline='');shutil.copy2(sif,R/sif.name);shutil.copy2(initial,M/'phase24_short_ab_singlepixel_vector_1.state')
log=HERE/'singlepixel_on_perf.solver.log';t=time.perf_counter()
with log.open('w',encoding='utf-8',newline='\n') as f:p=subprocess.run([str(mpi),'-n','6',str(B/'fem'/'src'/'ElmerSolver_mpi.exe'),sif.name],cwd=R,env=env,stdout=f,stderr=subprocess.STDOUT)
rec={'mode':'on','verification':False,'returncode':p.returncode,'wall_seconds':time.perf_counter()-t,'initial_state_sha256':hashlib.sha256(initial.read_bytes()).hexdigest(),'sif_sha256':hashlib.sha256(sif.read_bytes()).hexdigest(),'log':log.name}
(HERE/'singlepixel_on_perf_run.json').write_text(json.dumps(rec,indent=2)+'\n',encoding='utf-8')
if p.returncode:raise SystemExit(1)
