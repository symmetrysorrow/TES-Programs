from pathlib import Path
import json, os, shutil, subprocess, time, hashlib
HERE=Path(__file__).resolve().parent
RUNTIME=Path(r'C:\Users\Nao\AppData\Local\Temp\phase24-cache-runtime')
BUILD=Path(r'C:\Users\Nao\AppData\Local\Temp\phase24-release-build')
MESH=RUNTIME/'workspace'/'work'/'meshes'/'mesh_hybrid_fullconf_h8'
INITIAL=HERE.parent/'phase24_short_ab_singlepixel_vector_1.state'
MPI=Path(r'C:\Program Files\Microsoft MPI\Bin\mpiexec.exe')
SOLVER=BUILD/'fem'/'src'/'ElmerSolver_mpi.exe'
PATHS=[RUNTIME/'build'/'udf',BUILD/'fem'/'src'/'modules',BUILD/'fem'/'src',BUILD/'matc'/'src',BUILD/'fhutiter'/'src',BUILD/'fem'/'src'/'binio',BUILD/'mathlibs'/'src'/'arpack',BUILD/'mathlibs'/'src'/'parpack',RUNTIME/'bin',RUNTIME/'share'/'elmersolver'/'lib',Path(r'G:\\Elmer-Migration\\msys64\\ucrt64\\bin'),Path(r'G:\\Elmer-Migration\\msys64\\usr\\bin'),Path(r'G:\\Elmer-Migration\\elmer\\bin'),MPI.parent]
env=os.environ.copy(); env.update(ELMER_HOME=str(RUNTIME),OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE',PHASE24_HYPRE_REPORT='1'); env['PATH']=os.pathsep.join(map(str,PATHS))+os.pathsep+env.get('PATH','')
records=[]
for mode in ('off','on'):
    shutil.copy2(INITIAL,MESH/'phase24_short_ab_singlepixel_vector_1.state')
    sif=HERE/f'singlepixel_release_{mode}_6rank.sif'; shutil.copy2(sif,RUNTIME/sif.name)
    log=HERE/f'singlepixel_{mode}.solver.log'; t=time.perf_counter()
    with log.open('w',encoding='utf-8',newline='\n') as f: p=subprocess.run([str(MPI),'-n','6',str(SOLVER),sif.name],cwd=RUNTIME,env=env,stdout=f,stderr=subprocess.STDOUT)
    records.append({'mode':mode,'returncode':p.returncode,'wall_seconds':time.perf_counter()-t,'initial_state_sha256':hashlib.sha256(INITIAL.read_bytes()).hexdigest(),'log':log.name,'sif_sha256':hashlib.sha256(sif.read_bytes()).hexdigest()})
    (HERE/'singlepixel_release_runs.json').write_text(json.dumps({'ranks':6,'steps':4,'vector_assembly':False,'cache_default':False,'build':'Release -O3','state_file':'phase24_short_ab_singlepixel_vector_1.state','runs':records},indent=2)+'\n',encoding='utf-8')
    if p.returncode: raise SystemExit(1)
