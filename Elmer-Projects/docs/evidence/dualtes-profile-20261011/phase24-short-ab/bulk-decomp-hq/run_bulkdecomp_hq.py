from pathlib import Path
import re,shutil,hashlib,subprocess,os,json,datetime,sys
base=Path(__file__).resolve().parent.parent
out=Path(__file__).resolve().parent
runtime_root=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(os.environ.get('PHASE24_RUNTIME_ROOT',r'D:\\github\\TES-Programs\\Elmer-Projects')).resolve()
build=Path(os.environ['TEMP'])/'phase24-debug-build'; runtime=Path(os.environ['TEMP'])/'phase24-debug-runtime'; mpi=Path(r'C:\Program Files\Microsoft MPI\Bin\mpiexec.exe')
parts=[runtime_root/'build/udf',build/'fem/src/modules',build/'fem/src',build/'matc/src',build/'fhutiter/src',build/'fem/src/binio',build/'mathlibs/src/arpack',build/'mathlibs/src/parpack',runtime/'bin',runtime/'share/elmersolver/lib',Path(r'G:\Elmer-Migration\msys64\ucrt64\bin'),Path(r'G:\Elmer-Migration\msys64\usr\bin'),Path(r'G:\Elmer-Migration\elmer\bin'),mpi.parent]
env=os.environ.copy();env.update(ELMER_HOME=str(runtime),PATH=os.pathsep.join([str(x) for x in parts]+[os.environ.get('PATH','')]),OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE')
runs=[]
for label,src,ntes,statebase in [('dualtes',base/'bulk-decomp/dualtes_6rank.sif',2,'phase24_short_ab_dualtes_vector'),('singlepixel',base/'bulk-decomp/singlepixel_6rank.sif',1,'phase24_short_ab_singlepixel_vector')]:
 s=src.read_text(encoding='utf-8'); mesh=re.search(r'Mesh DB\s+"[^"]+"\s+"([^"]+)"',s).group(1)
 for tes in range(1,ntes+1):
  state=runtime_root/'workspace/work/meshes'/mesh/f'{label}_bulkdecomp_hq_{tes}.state';shutil.copy2(base/f'{statebase}_{tes}.state',state)
  for kind in ['State File','Series File','Iteration Series File']:
   pat=rf'(TES(?: {tes})? {re.escape(kind)} = String ")[^"]+(" )';pat=rf'(TES(?: {tes})? {re.escape(kind)} = String ")[^"]+(\")'
   value=state.relative_to(runtime_root).as_posix() if kind=='State File' else f'{label}_bulkdecomp_hq_{tes}_{"series" if kind=="Series File" else "iterations"}.csv'
   s,n=re.subn(pat,lambda m:m.group(1)+value+m.group(2),s,count=1)
   if n!=1: raise RuntimeError((label,tes,kind,n))
 sif=out/f'{label}_6rank.sif';sif.write_text(s,encoding='utf-8');log=out/f'{label}_6rank.solver.log';cmd=[str(mpi),'-n','6',str(runtime/'bin/ElmerSolver_mpi.exe'),str(sif)]
 start=datetime.datetime.now().astimezone().isoformat(timespec='seconds')
 with log.open('w',encoding='utf-8') as f: proc=subprocess.run(cmd,cwd=runtime_root,env=env,stdout=f,stderr=subprocess.STDOUT,check=False)
 outputs={}
 for t in range(1,ntes+1):
  for kind in ['series','iterations']:
   q=runtime_root/f'{label}_bulkdecomp_hq_{t}_{kind}.csv'
   if q.exists(): shutil.copy2(q,out/q.name);outputs[q.name]=hashlib.sha256((out/q.name).read_bytes()).hexdigest()
 runs.append({'case':label,'ranks':6,'steps':4,'returncode':proc.returncode,'started':start,'finished':datetime.datetime.now().astimezone().isoformat(timespec='seconds'),'sif_sha256':hashlib.sha256(sif.read_bytes()).hexdigest(),'log':log.name,'outputs':outputs})
 print(json.dumps(runs[-1],ensure_ascii=False))
(out/'manifest.json').write_text(json.dumps({'base_commit':'8bc3e403988450ebb86d9d60654f8c42adf00f21','compiler':'GNU Fortran 16.1.0','openmp':'GNU libgomp 5.2','mpi':'MS-MPI','build':'Debug -O0 -g3 -fbacktrace -fno-omit-frame-pointer','omp_num_threads':1,'category_timing':'RealTime; group interval and existing per-element stage counters','csr_rhs_timing':'omp_get_wtime high-resolution, one sample per 32 body/class equation calls, scaled by call/sample','timer_probe':'1,000,000 omp_get_wtime calls = 0.044 s (44 ns/call)','runs':runs},indent=2)+'\n',encoding='utf-8')