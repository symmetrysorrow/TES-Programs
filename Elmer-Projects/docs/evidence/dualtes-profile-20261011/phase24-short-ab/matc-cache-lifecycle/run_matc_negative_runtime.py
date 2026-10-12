"""One-step, isolated SIF probes: cache must reject non-pure Material 7 MATC."""
from pathlib import Path
import hashlib,json,os,re,shutil,subprocess,time,sys

HERE=Path(__file__).resolve().parent
RT=Path(r'C:\Users\Nao\AppData\Local\Temp\phase24-cache-runtime')
BUILD=Path(r'C:\Users\Nao\AppData\Local\Temp\phase24-final-lifecycle-release')
MPI=Path(r'C:\Program Files\Microsoft MPI\Bin\mpiexec.exe')
MESH=RT/'workspace'/'work'/'meshes'/'mesh_hybrid_fullconf_h8'
STATE=HERE.parent/'phase24_short_ab_singlepixel_vector_1.state'
BASE=HERE/'finalreplay_singlepixel_on.sif'
PATHS=[RT/'build'/'udf',BUILD/'fem'/'src'/'modules',BUILD/'fem'/'src',BUILD/'matc'/'src',BUILD/'fhutiter'/'src',BUILD/'fem'/'src'/'binio',BUILD/'mathlibs'/'src'/'arpack',BUILD/'mathlibs'/'src'/'parpack',RT/'bin',RT/'share'/'elmersolver'/'lib',Path(r'G:\Elmer-Migration\msys64\ucrt64\bin'),Path(r'G:\Elmer-Migration\msys64\usr\bin'),Path(r'G:\Elmer-Migration\elmer\bin'),MPI.parent]
CASES=[
 ('time_or_external','tx+ty','non-tx symbol ty'),
 ('coordinate','tx+x','coordinate symbol x'),
 ('unregistered','tx+unlisted_time_state','unregistered symbol'),
 ('unknown_function','unlisted_function(tx)','unknown function'),
 ('assignment','tx=tx+1','assignment side effect'),
 ('parse_error','tx+(','MATC syntax error'),
 ('time_dependency','TIME_DEP','property declared Variable Time'),
 ('coordinate_dependency','COORD_DEP','property declared Variable Coordinate 1'),
 ('procedure','PROCEDURE','procedure instead of MATC'),
 ('anisotropic','ANISOTROPIC','additional directional conductivity entry'),
 ('multiple_conductivity','MULTIPLE','multiple conductivity entries'),
]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
env=os.environ.copy();env.update(ELMER_HOME=str(RT),OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE');env['PATH']=os.pathsep.join(map(str,PATHS))+os.pathsep+env.get('PATH','')
records=[]
base=BASE.read_text(encoding='utf-8-sig')
for name,expr,reason in CASES:
 if len(sys.argv)>1 and name not in sys.argv[1:]: continue
 prefix='negative_'+name
 if expr=='PROCEDURE':
  raw=re.sub(r'(?m)^\s*Real MATC "[^"]*"','    Real Procedure "tes_transient_heat_source_t0" "AbsorberWindowPulseHeatSource"',base,count=1)
 elif expr=='TIME_DEP':
  raw=base.replace('Heat Conductivity = Variable Temperature','Heat Conductivity = Variable Time',1)
 elif expr=='COORD_DEP':
  raw=base.replace('Heat Conductivity = Variable Temperature','Heat Conductivity = Variable Coordinate 1',1)
 elif expr=='ANISOTROPIC':
  raw=base
  marker='    Real MATC "(7.854e-06)*(0.5*((tx+1.0e-12)+abs(tx-1.0e-12)))^(4.252-1)*0.4*((0.0007)-(0.0005))/(8*((1e-06)+(1.5e-05)+(1e-06))*(0.0005))"'
  raw=raw.replace(marker,marker+'\n  Heat Conductivity(3) = 0.0001 0 0',1)
 elif expr=='MULTIPLE':
  marker='    Real MATC "(7.854e-06)*(0.5*((tx+1.0e-12)+abs(tx-1.0e-12)))^(4.252-1)*0.4*((0.0007)-(0.0005))/(8*((1e-06)+(1.5e-05)+(1e-06))*(0.0005))"'
  raw=base.replace(marker,marker+'\n  Heat Conductivity(3) = 0.0001 0 0\n  Heat Conductivity(3) = 0 0.0001 0',1)
 else:
  raw=re.sub(r'(?m)^\s*Real MATC "[^"]*"',f'    Real MATC "{expr}"',base,count=1)
 if raw==base: raise RuntimeError('Material 7 MATC pattern was not replaced')
 raw=re.sub(r'(?m)^\s*Output File\s*=.*$',f'  Output File = {prefix}.result',raw,count=1)
 raw=raw.replace('finalreplay_singlepixel_on',prefix)
 raw=re.sub(r'(?m)^\s*Timestep Intervals\(3\)\s*=.*$','  Timestep Intervals(3) = 1 0 0',raw,count=1)
 raw=re.sub(r'(?m)^\s*Phase24 Dynamic K Node Cache Verify\s*=.*$','  Phase24 Dynamic K Node Cache Verify = Logical False',raw)
 sif=HERE/f'{prefix}.sif';sif.write_text(raw,encoding='utf-8',newline='\n');shutil.copy2(sif,RT/sif.name);shutil.copy2(STATE,MESH/STATE.name)
 log=HERE/f'{prefix}.solver.log';t=time.perf_counter()
 with log.open('w',encoding='utf-8',newline='\n') as out:
  p=subprocess.run([str(MPI),'-n','6',str(BUILD/'fem'/'src'/'ElmerSolver_mpi.exe'),sif.name],cwd=RT,env=env,stdout=out,stderr=subprocess.STDOUT)
 wall=time.perf_counter()-t;txt=log.read_text(encoding='utf-8',errors='replace')
 rejected=('no safe temperature-only MATC prisms' in txt or 'no safe temperature-only MATC' in txt)
 builds=re.findall(r'PHASE24_MATC_CACHE_BUILD rank=(\d+), materials=(\d+)',txt)
 rec={'test':name,'expression':expr,'reason':reason,'returncode':p.returncode,'wall_seconds':wall,'safe_cache_rejected':rejected,'cache_build_records':builds,'solver_sha256':sha(BUILD/'fem'/'src'/'ElmerSolver_mpi.exe'),'heat_solve_sha256':sha(BUILD/'fem'/'src'/'modules'/'HeatSolve.dll'),'sif_sha256':sha(sif),'log':log.name}
 records.append(rec);(HERE/'negative_runtime_results.json').write_text(json.dumps({'build':'official 26.1 + patches 01/02/03 isolated Release','vector_assembly':False,'steps':1,'records':records},indent=2)+'\n',encoding='utf-8');print(json.dumps(rec),flush=True)
