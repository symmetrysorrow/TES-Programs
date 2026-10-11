from pathlib import Path
import os,shutil,subprocess,time,hashlib,json,re
HERE=Path(__file__).resolve().parent; P=HERE.parent/'bulk-matc-cache'; R=Path(r'C:\\Users\\Nao\\AppData\\Local\\Temp\\phase24-cache-runtime'); B=Path(r'C:\\Users\\Nao\\AppData\\Local\\Temp\\phase24-release-build'); M=R/'workspace'/'work'/'meshes'/'mesh_toml_95333801e1a3'; MPI=Path(r'C:\\Program Files\\Microsoft MPI\\Bin\\mpiexec.exe')
paths=[R/'build'/'udf',B/'fem'/'src'/'modules',B/'fem'/'src',B/'matc'/'src',B/'fhutiter'/'src',B/'fem'/'src'/'binio',B/'mathlibs'/'src'/'arpack',B/'mathlibs'/'src'/'parpack',R/'bin',R/'share'/'elmersolver'/'lib',Path(r'G:\\Elmer-Migration\\msys64\\ucrt64\\bin'),Path(r'G:\\Elmer-Migration\\msys64\\usr\\bin'),Path(r'G:\\Elmer-Migration\\elmer\\bin'),MPI.parent]
env=os.environ.copy();env.update(ELMER_HOME=str(R),OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE',PHASE24_HYPRE_REPORT='1');env['PATH']=os.pathsep.join(map(str,paths))+os.pathsep+env.get('PATH','')
s=(P/'dualtes_cache_verify_final_6rank.sif').read_text(encoding='utf-8-sig').replace('dualtes_cache_verify','dualtes_release_verify')
sif=HERE/'dualtes_release_verify_6rank.sif';sif.write_text(s,encoding='utf-8',newline='');shutil.copy2(sif,R/sif.name)
for n in (1,2):shutil.copy2(P/f'phase24_short_ab_dualtes_vector_{n}.state',M/f'phase24_short_ab_dualtes_vector_{n}.state')
log=HERE/'dualtes_release_verify.solver.log';t=time.perf_counter()
with log.open('w',encoding='utf-8',newline='\n') as f:proc=subprocess.run([str(MPI),'-n','6',str(B/'fem'/'src'/'ElmerSolver_mpi.exe'),sif.name],cwd=R,env=env,stdout=f,stderr=subprocess.STDOUT)
rec={'returncode':proc.returncode,'wall_seconds':time.perf_counter()-t,'ranks':6,'steps':4,'vector_assembly':False,'cache':True,'cache_verify':True,'log':log.name,'sif_sha256':hashlib.sha256(sif.read_bytes()).hexdigest()}
(HERE/'dualtes_release_verify_run.json').write_text(json.dumps(rec,indent=2)+'\n',encoding='utf-8')
if proc.returncode:raise SystemExit(1)
