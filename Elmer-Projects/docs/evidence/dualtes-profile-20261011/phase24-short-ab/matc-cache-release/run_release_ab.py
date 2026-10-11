#!/usr/bin/env python3
"""Run isolated Release DualTES MATC cache performance trials."""
from __future__ import annotations
import datetime, hashlib, json, os, re, shutil, statistics, subprocess, time
from pathlib import Path

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent / 'bulk-matc-cache'
RUNTIME = Path(r'C:\Users\Nao\AppData\Local\Temp\phase24-cache-runtime')
BUILD = Path(r'C:\Users\Nao\AppData\Local\Temp\phase24-release-build')
SOLVER = BUILD / 'fem' / 'src' / 'ElmerSolver_mpi.exe'
MPI = Path(r'C:\Program Files\Microsoft MPI\Bin\mpiexec.exe')
MESH_TAG = 'mesh_toml_95333801e1a3'
MESH = RUNTIME / 'workspace' / 'work' / 'meshes' / MESH_TAG
RESULTS = RUNTIME / MESH_TAG
STATES = [f'phase24_short_ab_dualtes_vector_{i}.state' for i in (1, 2)]
PATHS = [RUNTIME / 'build' / 'udf', BUILD / 'fem' / 'src' / 'modules',
         BUILD / 'fem' / 'src', BUILD / 'matc' / 'src', BUILD / 'fhutiter' / 'src',
         BUILD / 'fem' / 'src' / 'binio', BUILD / 'mathlibs' / 'src' / 'arpack',
         BUILD / 'mathlibs' / 'src' / 'parpack', RUNTIME / 'bin',
         RUNTIME / 'share' / 'elmersolver' / 'lib',
         Path(r'G:\Elmer-Migration\msys64\ucrt64\bin'),
         Path(r'G:\Elmer-Migration\msys64\usr\bin'),
         Path(r'G:\Elmer-Migration\elmer\bin'), MPI.parent]

def sha(p: Path) -> str:
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()

def main() -> int:
    for p in (SOLVER, MPI, MESH / 'partitioning.6', MESH / 'case_toml_2cc534365acf_steady.result.0'):
        if not p.exists(): raise SystemExit(f'missing required isolated input: {p}')
    env = os.environ.copy()
    env.update(ELMER_HOME=str(RUNTIME), OMP_NUM_THREADS='1', OMP_DYNAMIC='FALSE', PHASE24_HYPRE_REPORT='1')
    env['PATH'] = os.pathsep.join([str(p) for p in PATHS] + [env.get('PATH', '')])
    order = [('warmup_off','off'), ('warmup_on','on'),
             ('off_1','off'), ('on_1','on'), ('on_2','on'), ('off_2','off'),
             ('off_3','off'), ('on_3','on')]
    records = []
    for label, mode in order:
        for state in STATES:
            shutil.copy2(PARENT / state, MESH / state)
        src = PARENT / f'dualtes_field_{mode}_6rank.sif'
        sif_text = src.read_text(encoding='utf-8-sig')
        sif_text = re.sub(r'(?m)^\s*Output File\s*=.*$', f'  Output File = dualtes_release_{label}.result', sif_text, count=1)
        sif_text = sif_text.replace(f'dualtes_field_{mode}_', f'dualtes_release_{label}_')
        sif_text = re.sub(r'dualtes_cache_verify_(\d)_(series|iterations)\.csv', rf'dualtes_release_{label}_\1_\2.csv', sif_text)
        sif_runtime = RUNTIME / f'dualtes_release_{label}.sif'
        sif_runtime.write_text(sif_text, encoding='utf-8')
        shutil.copy2(sif_runtime, HERE / f'{label}.sif')
        log = HERE / f'{label}.solver.log'
        started = datetime.datetime.now().astimezone().isoformat(timespec='seconds')
        t0 = time.perf_counter()
        with log.open('w', encoding='utf-8', newline='\n') as f:
            p = subprocess.run([str(MPI), '-n', '6', str(SOLVER), sif_runtime.name],
                               cwd=RUNTIME, env=env, stdout=f, stderr=subprocess.STDOUT, check=False)
        wall = time.perf_counter() - t0
        log_text = log.read_text(encoding='utf-8', errors='replace')
        outputs = []
        for tes in (1,2):
            for kind in ('series','iterations'):
                old = RUNTIME / f'dualtes_release_{label}_{tes}_{kind}.csv'
                if old.exists():
                    new = HERE / f'{label}_{tes}_{kind}.csv'
                    shutil.copy2(old, new)
                    outputs.append({'file':new.name,'sha256':sha(new)})
        result_base = RESULTS / f'dualtes_release_{label}.result'
        for result in ([result_base] if result_base.exists() else []) + sorted(RESULTS.glob(result_base.name + '.*')):
            if result.is_file(): outputs.append({'file':result.name,'bytes':result.stat().st_size,'sha256':sha(result)})
        bulk = {}
        pat = re.compile(r'PHASE24_BULK_LOCAL: rank=(\d+), iter=(\d+), body=(\d+), material=(\d+), class=(\d+), elements=(\d+), body_class_wall=\s*([0-9.Ee+-]+), material_eval=\s*([0-9.Ee+-]+), local_calc=\s*([0-9.Ee+-]+), global_update=\s*([0-9.Ee+-]+), other=\s*([0-9.Ee+-]+)')
        for m in pat.finditer(log_text):
            rank=m.group(1); bulk.setdefault(rank,[]).append({'iter':int(m.group(2)),'body':int(m.group(3)),'material':int(m.group(4)),'class':int(m.group(5)),'elements':int(m.group(6)),'body_class_wall_s':float(m.group(7)),'material_eval_s':float(m.group(8)),'local_calc_s':float(m.group(9)),'global_update_s':float(m.group(10)),'other_s':float(m.group(11))})
        hypre = re.findall(r'Phase24 PCG iterations=(\d+) final_rel_res=([0-9.Ee+-]+)',log_text)
        norms = [float(x) for x in re.findall(r'HeatSolve:\s+Result Norm\s*:\s*([0-9.Ee+-]+)',log_text)]
        status = {'label':label,'mode':mode,'warmup':label.startswith('warmup'),'returncode':p.returncode,'started':started,'wall_seconds':wall,'sif_sha256':sha(sif_runtime),'log':log.name,'hypre_pcg':[[int(i),float(r)] for i,r in hypre],'result_norms':norms,'bulk_local_by_rank':bulk,'outputs':outputs}
        records.append(status)
        print(json.dumps(status, separators=(',',':')))
    timed = [r for r in records if not r['warmup']]
    summary={}
    for mode in ('off','on'):
        vals=[r['wall_seconds'] for r in timed if r['mode']==mode]
        summary[mode]={'runs_s':vals,'mean_s':statistics.mean(vals),'median_s':statistics.median(vals),'stdev_s':statistics.stdev(vals)}
    (HERE/'release_runs.json').write_text(json.dumps({'compiler_flags':{'Fortran':'-O3','C':'-O3 -DNDEBUG','CXX':'-O3 -DNDEBUG'},'mpi':'MS-MPI','hypre':'same Release build','ranks':6,'steps':4,'omp_threads':1,'vector_assembly':False,'cache_default':False,'order':[x[0] for x in order],'summary':summary,'runs':records},indent=2)+'\n',encoding='utf-8')
    return 0 if all(r['returncode']==0 for r in records) else 1
if __name__=='__main__': raise SystemExit(main())
