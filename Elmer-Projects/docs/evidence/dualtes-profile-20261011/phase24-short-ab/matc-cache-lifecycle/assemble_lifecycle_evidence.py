from pathlib import Path
import json,re,hashlib,statistics,subprocess
E=Path(__file__).resolve().parent
R=Path(r'C:\Users\Nao\AppData\Local\Temp\phase24-cache-runtime')
B=Path(r'C:\Users\Nao\AppData\Local\Temp\phase24-03-release-build')
S=Path(r'C:\Users\Nao\AppData\Local\Temp\phase24-03-build-source')
repo=E.parents[5]
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
runs=json.loads((E/'runtime_runs.json').read_text()); partial=json.loads((E/'runtime_partial_runs.json').read_text())
rec={(x['case'],x['mode']):x for x in partial['records']}
comp=runs['temperature_comparisons']
# Check output hashes for TES CSVs across each cache pair.
tes={}
for case in ('dualtes','singlepixel'):
 a=rec[(case,'off')]; b=rec[(case,'on')]; amap={x['file'].replace('_off_','_MODE_'):x['sha256'] for x in a['outputs'] if x['file'].endswith('.csv')}; bmap={x['file'].replace('_on_','_MODE_'):x['sha256'] for x in b['outputs'] if x['file'].endswith('.csv')}; tes[case]={'files':len(amap),'same_names':amap.keys()==bmap.keys(),'byte_identical':amap==bmap}
# HYPRE/norm equality and all-temperature coverage.
perf=json.loads((E/'performance_runs.json').read_text()); ps=json.loads((E/'performance_summary.json').read_text())
# Lifecycle per-rank counts from verification log.
log=(E/'lifecycle_dualtes_verify.solver.log').read_text(errors='replace')
build_lines=re.findall(r'PHASE24_MATC_CACHE_BUILD rank=(\d+)(?:, materials=(\d+), unique_nodes=(\d+))?, wall_seconds=\s*([0-9.Ee+-]+)',log)
eval_lines=re.findall(r'PHASE24_MATC_CACHE_EVAL rank=(\d+), materials=(\d+), unique_nodes=(\d+), listgetreal_calls=(\d+), previous_assembly_hits=(\d+), wall_seconds=\s*([0-9.Ee+-]+)',log)
byrank={}
for rank in range(6):
 bitems=[x for x in build_lines if int(x[0])==rank]; eitems=[x for x in eval_lines if int(x[0])==rank]
 byrank[str(rank)]={'map_build_activations':len(bitems),'materials':int(next((x[1] for x in bitems if x[1]),0)),'unique_nodes':int(next((x[2] for x in bitems if x[2]),0)),'matc_eval_assemblies':len(eitems),'listgetreal_calls_per_assembly':sorted(set(int(x[3]) for x in eitems)),'reported_previous_assembly_hits':sorted(set(int(x[4]) for x in eitems if int(x[4])>0)),'mean_build_wall_seconds':statistics.mean(float(x[3]) for x in bitems) if bitems else None}
# Result Norm and PCG sequences equal exact string parsed arrays.
solver_equal={case:{'result_norms_equal':rec[(case,'off')]['result_norms']==rec[(case,'on')]['result_norms'],'hypre_pcg_equal':rec[(case,'off')]['hypre_pcg']==rec[(case,'on')]['hypre_pcg']} for case in ('dualtes','singlepixel')}
mat=json.loads((E/'matrix_rhs_comparison.json').read_text())
num={'build':'Phase24 official Elmer 26.1 Release + patches 01/02/03','solver_sha256':sha(B/'fem/src/ElmerSolver_mpi.exe'),'heat_solve_sha256':sha(B/'fem/src/modules/HeatSolve.dll'),'libmatc_sha256':sha(B/'matc/src/libmatc.dll'),'runtime_runs':[{k:x.get(k) for k in ('case','mode','returncode','wall_seconds','solver_sha256','heat_solve_sha256','result_norms','hypre_pcg')} for x in partial['records']],'temperature_comparisons':comp,'tes_csv_comparisons':tes,'norm_and_hypre_exact_comparisons':solver_equal,'matrix_rhs':mat}
(E/'numeric_comparison_lifecycle.json').write_text(json.dumps(num,indent=2)+'\n',encoding='utf-8')
life={'verify_run':json.loads((E/'lifecycle_verify_run.json').read_text()),'no_safe_temperature_only_message_count':len(re.findall('no safe temperature-only MATC prisms',log)),'map_build_lines':len(build_lines),'matc_eval_lines':len(eval_lines),'rank_summary':byrank,'semantics':'eligible prism element hit counts are observed in previous_assembly_hits; node maps are rebuilt once per HeatSolver activation; nodal conductivity is re-evaluated once per eligible material at each nonlinear assembly.'}
(E/'lifecycle_runtime_summary.json').write_text(json.dumps(life,indent=2)+'\n',encoding='utf-8')
# patches and build provenance
p1=E.parent/'phase24-source-repro/01_phase24_source_from_release_26_1.patch';p2=E.parent/'matc-ast-release/02_matc_ast_purity_validator.patch';p3=E/'03_matc_cache_lifecycle.patch'
manifest={'repository_base_commit':'cca48a8107f6313dbbbaa57267c26a87f1652cfa','elmer_base_commit':'5a8de867068be0568f09af40fb90fee300dfbede','patches':[{'path':str(p.relative_to(repo)).replace('\\','/'),'sha256':sha(p)} for p in (p1,p2,p3)],'patch_01_02_immutable':True,'patch_03_source_manifest':json.loads((E/'03_source_hash_manifest.json').read_text()),'toolchain':{'compiler':'GNU Fortran/C/C++ 16.1.0','generator':'CMake + Ninja','build_type':'Release','MPI':'Microsoft MPI 2.0','HYPRE':'enabled','OpenMP':'enabled','vector_assembly':'OFF','cache_default':'OFF'},'binaries':{'ElmerSolver_mpi.exe':sha(B/'fem/src/ElmerSolver_mpi.exe'),'HeatSolve.dll':sha(B/'fem/src/modules/HeatSolve.dll'),'libmatc.dll':sha(B/'matc/src/libmatc.dll')},'result_summary':{'ast_unit_tests':'PASS','lifecycle_4step_dualtes':'PASS','lifecycle_4step_singlepixel':'PASS','rank_local_cache_usage':'PASS','matrix_rhs_timestep1':'PASS','performance_dualtes_repeated':'PASS','same_process_model_mesh_mutations':'SKIP','same_process_off_on_off':'SKIP','runtime_memory_sanitizer':'SKIP','75ms_full_run':'NOT_RUN'},'artifacts':['report_ja.md','03_matc_cache_lifecycle.patch','03_source_hash_manifest.json','replay_01_02_03_and_build.ps1','run_lifecycle_release_ab.py','run_lifecycle_verify.py','run_lifecycle_matrix_rhs_capture.py','matc_ast_purity_test.c','run_lifecycle_ast_unit_tests.ps1','compare_matrix_rhs_lifecycle.py','run_lifecycle_performance.py','runtime_partial_runs.json','runtime_runs.json','numeric_comparison_lifecycle.json','lifecycle_runtime_summary.json','matrix_rhs_capture_runs.json','matrix_rhs_comparison.json','performance_runs.json','performance_summary.json','lifecycle_dualtes_verify.sif','lifecycle_dualtes_verify.solver.log','matc_ast_unit.log','release_build.log']}
manifest['evidence_files']=[{'path':p.name,'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(E.iterdir()) if p.is_file() and p.name!='evidence_manifest.json']
(E/'evidence_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'temperature':{k:{'max_abs_K':v['max_abs_K'],'max_rel':v['max_rel'],'different':v['total_different_values']} for k,v in comp.items()},'tes':tes,'norm_hypre':solver_equal,'memory_rank':byrank,'perf':ps,'matrix_rhs':mat},indent=2))

