from pathlib import Path
import hashlib,json,glob,datetime

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[5]
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
fixed=[
 'report_final_ja.md','final_lifecycle_test_matrix.json','final_release_runtime_runs.json',
 'final_release_full_run_summary.json','final_release_temperature_comparison.json',
 'final_release_performance_runs.json','final_release_performance_summary.json',
 'old03_dualtes_performance_runs.json','old03_dualtes_performance_summary.json',
 'negative_runtime_results.json','matc_ast_unit_finalrelease.log','debug_lifecycle_run.json',
 'final_replay_build.log','debug_build.log','asan_capability_probe.log',
 'replay_01_02_03_and_build.ps1','run_final_release_ab.py','run_final_release_performance.py',
 'run_previous03_dualtes_performance.py','run_matc_negative_runtime.py',
 'run_debug_lifecycle_singlepixel.py','build_debug_diagnostics.ps1','check_asan_capability.ps1',
 'compare_final_release_temperatures.py','make_final_evidence_manifest.py',
 'matc_ast_purity_test.c','asan_capability_probe.c',
 'finalreplay_dualtes_off.sif','finalreplay_dualtes_on.sif',
 'finalreplay_singlepixel_off.sif','finalreplay_singlepixel_on.sif','lifecycle_debug_singlepixel.sif',
 'final_lifecycle_test_matrix.json','final_release_runtime_runs.json',
 'final_release_performance_runs.json','old03_dualtes_performance_runs.json',
 'negative_runtime_results.json'
]
patterns=['finalreplay_*.solver.log','finalperf_*.solver.log','old03_*.solver.log','negative_*.solver.log','negative_*.sif','finalperf_*.sif','old03_*.sif']
files=set(fixed)
for pat in patterns: files.update(p.name for p in HERE.glob(pat))
files.update(p.name for p in HERE.glob('*.log'))
files.add('debug_lifecycle_singlepixel.solver.log')
file_records=[]
for name in sorted(files):
 p=HERE/name
 if p.is_file() and p.suffix.lower()!='.exe': file_records.append({'path':name,'bytes':p.stat().st_size,'sha256':sha(p)})
patches={
 '01':REPO/'Elmer-Projects/docs/evidence/dualtes-profile-20261011/phase24-short-ab/phase24-source-repro/01_phase24_source_from_release_26_1.patch',
 '02':REPO/'Elmer-Projects/docs/evidence/dualtes-profile-20261011/phase24-short-ab/matc-ast-release/02_matc_ast_purity_validator.patch',
 '03':HERE/'03_matc_cache_lifecycle.patch'
}
exe_paths={
 'final_release_solver':Path(r'C:\Users\Nao\AppData\Local\Temp\phase24-final-lifecycle-release\fem\src\ElmerSolver_mpi.exe'),
 'final_release_heat_solve':Path(r'C:\Users\Nao\AppData\Local\Temp\phase24-final-lifecycle-release\fem\src\modules\HeatSolve.dll'),
 'final_release_libmatc':Path(r'C:\Users\Nao\AppData\Local\Temp\phase24-final-lifecycle-release\matc\src\libmatc.dll'),
 'previous03_solver':Path(r'C:\Users\Nao\AppData\Local\Temp\phase24-03-release-build\fem\src\ElmerSolver_mpi.exe'),
 'previous03_heat_solve':Path(r'C:\Users\Nao\AppData\Local\Temp\phase24-03-release-build\fem\src\modules\HeatSolve.dll'),
 'debug_solver':Path(r'C:\Users\Nao\AppData\Local\Temp\phase24-final-lifecycle-debug\fem\src\ElmerSolver_mpi.exe'),
 'debug_heat_solve':Path(r'C:\Users\Nao\AppData\Local\Temp\phase24-final-lifecycle-debug\fem\src\modules\HeatSolve.dll'),
 'hypre_runtime':Path(r'G:\Elmer-Migration\msys64\ucrt64\bin\libHYPRE.dll'),
 'msmpi':Path(r'C:\Program Files\Microsoft MPI\Bin\mpiexec.exe')
}
external={k:({'path':str(p),'sha256':sha(p),'bytes':p.stat().st_size} if p.is_file() else {'path':str(p),'status':'MISSING'}) for k,p in exe_paths.items()}
restart_inputs={
 'dualtes_step_1':HERE.parent/'bulk-matc-cache'/'phase24_short_ab_dualtes_vector_1.state',
 'dualtes_step_2':HERE.parent/'bulk-matc-cache'/'phase24_short_ab_dualtes_vector_2.state',
 'singlepixel_restart':HERE.parent/'phase24_short_ab_singlepixel_vector_1.state'
}
manifest={
 'schema':'phase24-final-lifecycle-evidence-v1',
 'created_local':'2026-10-12',
 'repository_base_commit':'5ef7888636a2d238177eb5234ecb131e4d97651b',
 'target_branch':'codex/dualtes-profile-evidence-20261011',
 'upstream_elmer_commit':'5a8de867068be0568f09af40fb90fee300dfbede',
 'patch_order':[{'patch':k,'path':str(v.relative_to(REPO)).replace('\\','/'),'sha256':sha(v)} for k,v in patches.items()],
 'toolchain':{'GNU Fortran/C/C++':'16.1.0','generator':'CMake + Ninja','release':'optimized Release','debug':'-O0 -g -fcheck=all -finit-real=snan -finit-integer=-999 -fbacktrace -ffpe-trap=invalid,zero,overflow','MPI':'Microsoft MPI 2.0','HYPRE':'libHYPRE.dll; enabled','OpenMP':'enabled','vector_assembly':'OFF','cache_default':'OFF'},
 'external_binaries':external,
 'restart_input_states':{k:({'path':str(v.relative_to(REPO)).replace('\\','/'),'sha256':sha(v),'bytes':v.stat().st_size} if v.is_file() else {'path':str(v),'status':'MISSING'}) for k,v in restart_inputs.items()},
 'runtime_constraints':{'full_75ms_run':'NOT_RUN','production_DLL_modified':False,'production_SIF_modified':False},
 'artifacts':file_records,
 'omitted':'Large runtime binaries are outside the repository; SHA-256 and exact local paths are recorded above. The previous evidence directory and patches 01/02/03 are immutable.'
}
(HERE/'final_evidence_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'artifact_count':len(file_records),'manifest':'final_evidence_manifest.json','missing_binaries':[k for k,v in external.items() if 'status' in v]},indent=2))
