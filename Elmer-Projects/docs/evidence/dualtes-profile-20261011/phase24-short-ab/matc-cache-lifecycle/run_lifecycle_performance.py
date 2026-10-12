from pathlib import Path
import os,re,subprocess,time,json,hashlib,shutil,statistics
here=Path(__file__).resolve().parent
rt=Path(r'C:\Users\Nao\AppData\Local\Temp\phase24-cache-runtime')
build=Path(r'C:\Users\Nao\AppData\Local\Temp\phase24-03-release-build')
mpi=Path(r'C:\Program Files\Microsoft MPI\Bin\mpiexec.exe')
mesh=rt/'workspace'/'work'/'meshes'/'mesh_toml_95333801e1a3'
base=here.parent/'matc-cache-lifecycle'
source=here
states=[here.parent/'bulk-matc-cache'/f'phase24_short_ab_dualtes_vector_{i}.state' for i in (1,2)]
paths=[rt/'build'/'udf',build/'fem'/'src'/'modules',build/'fem'/'src',build/'matc'/'src',build/'fhutiter'/'src',build/'fem'/'src'/'binio',build/'mathlibs'/'src'/'arpack',build/'mathlibs'/'src'/'parpack',rt/'bin',rt/'share'/'elmersolver'/'lib',Path(r'G:\Elmer-Migration\msys64\ucrt64\bin'),Path(r'G:\Elmer-Migration\msys64\usr\bin'),Path(r'G:\Elmer-Migration\elmer\bin'),mpi.parent]
env=os.environ.copy();env.update(ELMER_HOME=str(rt),OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE',PHASE24_HYPRE_REPORT='1');env['PATH']=os.pathsep.join(map(str,paths))+os.pathsep+env.get('PATH','')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
records=[]
for i in range(1,4):
 for mode in ('off','on'):
  for st in states: shutil.copy2(st,mesh/st.name)
  old=source/f'final_dualtes_{mode}.sif'; prefix=f'lifecycle_perf_dualtes_{mode}_r{i}'
  raw=old.read_text(encoding='utf-8-sig').replace(f'final_dualtes_{mode}',prefix)
  raw=re.sub(r'(?m)^\s*Phase24 Dynamic K Node Cache Verify\s*=.*$','  Phase24 Dynamic K Node Cache Verify = Logical False',raw)
  sif=here/f'{prefix}.sif';sif.write_text(raw,encoding='utf-8',newline='\n');shutil.copy2(sif,rt/sif.name)
  log=here/f'{prefix}.solver.log';t=time.perf_counter()
  with log.open('w',encoding='utf-8',newline='\n') as out: proc=subprocess.run([str(mpi),'-n','6',str(build/'fem'/'src'/'ElmerSolver_mpi.exe'),sif.name],cwd=rt,env=env,stdout=out,stderr=subprocess.STDOUT)
  wall=time.perf_counter()-t;text=log.read_text(encoding='utf-8',errors='replace')
  rec={'case':'dualtes','mode':mode,'repeat':i,'returncode':proc.returncode,'wall_seconds':wall,'solver_sha256':sha(build/'fem'/'src'/'ElmerSolver_mpi.exe'),'heat_solve_sha256':sha(build/'fem'/'src'/'modules'/'HeatSolve.dll'),'result_norms':[float(x) for x in re.findall(r'HeatSolve:\s+Result Norm\s*:\s*([0-9.Ee+-]+)',text)],'hypre_pcg':[[int(a),float(b)] for a,b in re.findall(r'Phase24 PCG iterations=(\d+) final_rel_res=([0-9.Ee+-]+)',text)],'exit_ok':proc.returncode==0,'log':log.name}
  records.append(rec);(here/'performance_runs.json').write_text(json.dumps({'build':'Phase24 official 26.1 Release + 01 + 02 + 03','ranks':6,'steps':4,'verify':False,'records':records},indent=2)+'\n',encoding='utf-8');print(json.dumps(rec),flush=True)
  if proc.returncode: raise SystemExit(proc.returncode)
summary={}
for mode in ('off','on'):
 v=[x['wall_seconds'] for x in records if x['mode']==mode];summary[mode]={'runs':v,'mean':statistics.mean(v),'stdev_sample':statistics.stdev(v),'median':statistics.median(v)}
summary['speedup_off_over_on']=summary['off']['mean']/summary['on']['mean']
(here/'performance_summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8');print(json.dumps(summary,indent=2))
