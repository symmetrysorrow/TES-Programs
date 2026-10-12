from __future__ import annotations
import json, os, re, shutil, signal, subprocess, sys, time
from datetime import datetime, timezone
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import run_experiment as exp
exp.RUNTIME_ROOT = Path(r'C:\Users\Nao\AppData\Local\Temp\dualtes-win-singlepixel-75ms-priorflags-20261013')
input_sif = exp.INPUTS / 'singlepixel_75ms_cache_on.sif'
run_name = 'singlepixel75_full_nl1e6_precondreuseon_cacheon'
runtime = exp.stage_case('singlepixel')
exp.restore_initial_state('singlepixel', runtime)
text = input_sif.read_text(encoding='utf-8')
text = text.replace('singlepixel75_full_cache_on', run_name)
text = re.sub(r'(Nonlinear System Convergence Tolerance\s*=\s*)\S+', r'\g<1>1e-6', text)
text = re.sub(r'(Phase24 Static Matrix Reuse\s*=\s*Logical\s+False\s*)', r'\1  Phase24 Precond Matrix Reuse = Logical True\n', text)
if 'Phase24 Precond Matrix Reuse = Logical True' not in text:
    raise SystemExit('Failed to enable Phase24 preconditioner reuse')
if 'Nonlinear System Convergence Tolerance = 1e-6' not in text:
    raise SystemExit('Failed to set nonlinear tolerance')
if 'Phase24 Dynamic K Node Cache = Logical True' not in text or 'Phase24 Vector Assembly = Logical False' not in text:
    raise SystemExit('Expected cache-on/vector-off settings missing')
local_sif = runtime / 'case.sif'
local_sif.write_text(text, encoding='utf-8')
run_dir = HERE / 'runs' / 'windows-prior-flags-full' / 'singlepixel'
run_dir.mkdir(parents=True, exist_ok=True)
shutil.copy2(local_sif, run_dir / 'case.sif')
log = run_dir / 'solver.log'
env = os.environ.copy()
env.update({'ELMER_HOME': str(runtime), 'OMP_NUM_THREADS': '1', 'OMP_DYNAMIC': 'FALSE', 'PHASE24_HYPRE_REPORT': '1'})
ordered_paths = [runtime/'build/udf', exp.BUILD/'fem/src/modules', exp.BUILD/'fem/src', exp.BUILD/'matc/src', exp.BUILD/'fhutiter/src', exp.BUILD/'fem/src/binio', exp.BUILD/'mathlibs/src/arpack', exp.BUILD/'mathlibs/src/parpack', exp.BUILD/'share/elmersolver/lib', runtime/'share/elmersolver/lib', exp.MSYS_BIN, exp.MSYS_USR_BIN, exp.ELMER_BIN, exp.MPIEXEC.parent]
env['PATH'] = os.pathsep.join(str(p) for p in ordered_paths) + os.pathsep + env.get('PATH', '')
cmd = [str(exp.MPIEXEC), '-n', '6', str(exp.BUILD/'fem/src/ElmerSolver_mpi.exe'), local_sif.name]
start = datetime.now(timezone.utc).isoformat(); t0 = time.monotonic(); timed_out = False
with log.open('w', encoding='utf-8', newline='\n') as f:
    p = subprocess.Popen(cmd, cwd=runtime, env=env, stdout=f, stderr=subprocess.STDOUT, creationflags=subprocess.CREATE_NEW_PROCESS_GROUP)
    try:
        p.wait(timeout=600)
    except subprocess.TimeoutExpired:
        timed_out = True
        try:
            p.send_signal(signal.CTRL_BREAK_EVENT)
            p.wait(timeout=20)
        except Exception:
            subprocess.run(['taskkill', '/PID', str(p.pid), '/T', '/F'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            try: p.wait(timeout=10)
            except subprocess.TimeoutExpired: pass
elapsed = time.monotonic()-t0; end = datetime.now(timezone.utc).isoformat()
logtext = log.read_text(encoding='utf-8', errors='replace')
total = re.findall(r'SOLVER TOTAL TIME\(CPU,REAL\):\s*([0-9.eE+-]+)\s+([0-9.eE+-]+)', logtext)
steps = re.findall(r'MAIN:\s*Time:\s*(\d+)/225:\s*([0-9.Ee+-]+)', logtext)
record = {'case':'singlepixel','target_end_time_s':0.075,'run_mode':'full_75ms','start_utc':start,'end_utc':end,'wall_seconds':elapsed,'timeout_seconds':600,'timed_out':timed_out,'exit_code':p.returncode,'normal_completion':(p.returncode==0 and 'ALL DONE' in logtext),'last_completed_timestep':int(steps[-1][0]) if steps else None,'last_completed_solver_time_s':float(steps[-1][1]) if steps else None,'logged_timestep_markers':len(steps),'ranks':6,'omp_threads_per_rank':1,'matc_cache_on':True,'vector_assembly':False,'nonlinear_tolerance':1e-6,'preconditioner_matrix_reuse':True,'linear_max_iterations':2000,'linear_tolerance':1e-11,'input_sif':'singlepixel_75ms_cache_on.sif (isolated copy; only nonlinear tol/preconditioner reuse/output IDs changed)','input_sif_sha256':exp.sha256(input_sif),'run_sif_sha256':exp.sha256(run_dir/'case.sif'),'solver_sha256':exp.sha256(exp.BUILD/'fem/src/ElmerSolver_mpi.exe'),'heatsolve_sha256':exp.sha256(exp.BUILD/'fem/src/modules/HeatSolve.dll'),'run_record_log':'runs/windows-prior-flags-full/singlepixel/solver.log','runtime_root':str(runtime),'solver_total_cpu_real_s':list(map(float,total[-1])) if total else None}
(run_dir/'run_record.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
for pattern in (f'{run_name}*.csv',f'{run_name}*.result*',f'{run_name}*.state'):
    for src in runtime.glob(pattern):
        if src.is_file() and src.suffix.lower()=='.csv': shutil.copy2(src,run_dir/src.name)
print(json.dumps(record,indent=2))
