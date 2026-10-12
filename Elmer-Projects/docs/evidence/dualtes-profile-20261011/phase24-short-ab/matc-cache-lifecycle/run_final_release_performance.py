from pathlib import Path
import os,re,subprocess,time,json,hashlib,shutil,statistics

here=Path(__file__).resolve().parent
rt=Path(r'C:\Users\Nao\AppData\Local\Temp\phase24-cache-runtime')
build=Path(r'C:\Users\Nao\AppData\Local\Temp\phase24-final-lifecycle-release')
mpi=Path(r'C:\Program Files\Microsoft MPI\Bin\mpiexec.exe')
meshroot=rt/'workspace'/'work'/'meshes'
cases=[
 ('dualtes','mesh_toml_95333801e1a3',[here.parent/'bulk-matc-cache'/f'phase24_short_ab_dualtes_vector_{i}.state' for i in (1,2)]),
 ('singlepixel','mesh_hybrid_fullconf_h8',[here.parent/f'phase24_short_ab_singlepixel_vector_1.state']),
]
paths=[rt/'build'/'udf',build/'fem'/'src'/'modules',build/'fem'/'src',build/'matc'/'src',build/'fhutiter'/'src',build/'fem'/'src'/'binio',build/'mathlibs'/'src'/'arpack',build/'mathlibs'/'src'/'parpack',rt/'bin',rt/'share'/'elmersolver'/'lib',Path(r'G:\Elmer-Migration\msys64\ucrt64\bin'),Path(r'G:\Elmer-Migration\msys64\usr\bin'),Path(r'G:\Elmer-Migration\elmer\bin'),mpi.parent]
env=os.environ.copy();env.update(ELMER_HOME=str(rt),OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE',PHASE24_HYPRE_REPORT='1');env['PATH']=os.pathsep.join(map(str,paths))+os.pathsep+env.get('PATH','')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
records=[]
for case,meshname,states in cases:
 mesh=meshroot/meshname
 for repeat in range(1,4):
  for mode in ('off','on'):
   for st in states: shutil.copy2(st,mesh/st.name)
   src=here/f'finalreplay_{case}_{mode}.sif'
   prefix=f'finalperf_{case}_{mode}_r{repeat}'
   raw=src.read_text(encoding='utf-8-sig').replace(f'finalreplay_{case}_{mode}',prefix)
   raw=re.sub(r'(?m)^\s*Phase24 Dynamic K Node Cache Verify\s*=.*$','  Phase24 Dynamic K Node Cache Verify = Logical False',raw)
   sif=here/f'{prefix}.sif';sif.write_text(raw,encoding='utf-8',newline='\n');shutil.copy2(sif,rt/sif.name)
   log=here/f'{prefix}.solver.log';t=time.perf_counter()
   with log.open('w',encoding='utf-8',newline='\n') as out:
    proc=subprocess.run([str(mpi),'-n','6',str(build/'fem'/'src'/'ElmerSolver_mpi.exe'),sif.name],cwd=rt,env=env,stdout=out,stderr=subprocess.STDOUT)
   wall=time.perf_counter()-t;text=log.read_text(encoding='utf-8',errors='replace')
   rec={'case':case,'mode':mode,'repeat':repeat,'returncode':proc.returncode,'wall_seconds':wall,'ranks':6,'steps':4,'verify':False,'cache':mode=='on','vector_assembly':False,'solver_sha256':sha(build/'fem'/'src'/'ElmerSolver_mpi.exe'),'heat_solve_sha256':sha(build/'fem'/'src'/'modules'/'HeatSolve.dll'),'result_norms':[float(x) for x in re.findall(r'HeatSolve:\s+Result Norm\s*:\s*([0-9.Ee+-]+)',text)],'hypre_pcg':[[int(a),float(b)] for a,b in re.findall(r'Phase24 PCG iterations=(\d+) final_rel_res=([0-9.Ee+-]+)',text)],'log':log.name}
   records.append(rec);(here/'final_release_performance_runs.json').write_text(json.dumps({'build':'official Elmer 26.1 Release + 01/02/03; same machine/run sequence','records':records},indent=2)+'\n',encoding='utf-8');print(json.dumps(rec),flush=True)
   if proc.returncode: raise SystemExit(proc.returncode)
summary={}
for case,_,_ in cases:
 summary[case]={}
 for mode in ('off','on'):
  vals=[r['wall_seconds'] for r in records if r['case']==case and r['mode']==mode]
  summary[case][mode]={'runs':vals,'mean':statistics.mean(vals),'stdev_sample':statistics.stdev(vals),'median':statistics.median(vals)}
 summary[case]['speedup']=summary[case]['off']['mean']/summary[case]['on']['mean']
(here/'final_release_performance_summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8');print(json.dumps(summary,indent=2))
