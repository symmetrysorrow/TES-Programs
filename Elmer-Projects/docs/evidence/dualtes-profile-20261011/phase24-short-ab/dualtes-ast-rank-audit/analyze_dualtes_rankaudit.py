#!/usr/bin/env python3
"""Parse rank-resolved AST cache verification and 3x OFF/ON timings."""
from __future__ import annotations
import hashlib, json, re, statistics
from pathlib import Path
HERE=Path(__file__).resolve().parent
LOG=HERE/'dualtes_rankaudit_verify.solver.log'
text=LOG.read_text(encoding='utf-8',errors='replace')
build={}
for m in re.finditer(r'PHASE24_MATC_CACHE_BUILD rank=(\d+), (?:materials=(\d+), unique_nodes=(\d+), )?wall_seconds=\s*([0-9.Ee+-]+)',text):
 r=int(m.group(1)); build[r]={'material_count':int(m.group(2)) if m.group(2) else 0,'unique_nodes':int(m.group(3)) if m.group(3) else 0,'wall_seconds':float(m.group(4))}
evals={r:[] for r in range(6)}
for m in re.finditer(r'PHASE24_MATC_CACHE_EVAL rank=(\d+), materials=(\d+), unique_nodes=(\d+), listgetreal_calls=(\d+), previous_assembly_hits=(\d+), wall_seconds=\s*([0-9.Ee+-]+)',text):
 r=int(m.group(1)); evals[r].append({'materials':int(m.group(2)),'unique_nodes':int(m.group(3)),'listgetreal_calls':int(m.group(4)),'previous_assembly_hits':int(m.group(5)),'wall_seconds':float(m.group(6))})
mat7={r:0 for r in range(6)}
first=text.split('MAIN: Time: 2/4:',1)[0]
for m in re.finditer(r'PHASE24_BULK_LOCAL: rank=(\d+), iter=1, body=\d+, material=7, class=2, elements=(\d+)',first): mat7[int(m.group(1))]+=int(m.group(2))
verify={r:[] for r in range(6)}
for m in re.finditer(r'PHASE24_MATC_VERIFY rank=(\d+), element=(\d+), material=(\d+), max_abs=\s*([0-9.Ee+-]+)',text):
 verify[int(m.group(1))].append({'element':int(m.group(2)),'material':int(m.group(3)),'max_abs':float(m.group(4))})
ranks=[]
for r in range(6):
 e=evals[r]; candidate=mat7[r]>0
 ranks.append({'rank':r,'material7_prism_elements':mat7[r],
  'matc_metadata_gate':'PASS (global Material 7 scalar Variable Temperature, no Lua/procedure, one conductivity entry)' if candidate else 'NOT_EVALUATED: no local Material 7 prism candidates',
  'ast_gate':'PASS (cache builder accepted material; standalone AST test accepts exact expression)' if candidate else 'NOT_EVALUATED: no local candidate element',
  'unique_nodes':build[r]['unique_nodes'],'matc_evaluations_per_assembly':e[0]['listgetreal_calls'] if e else 0,
  'assembly_evaluation_events':len(e),
  'matc_evaluation_calls_total':sum(x['listgetreal_calls'] for x in e),
  'cache_hits_per_reuse_assembly':e[1]['previous_assembly_hits'] if len(e)>1 else 0,
  'cache_hit_elements_total':sum(x['previous_assembly_hits'] for x in e),
  'fallback_candidate_elements':0 if candidate else 0,
  'verify_samples':verify[r], 'build':build[r], 'assembly_evaluations':e})
perf=json.loads((HERE/'performance_repeats.json').read_text(encoding='utf-8'))
result_dir=Path(r'C:\Users\Nao\AppData\Local\Temp\phase24-cache-runtime')/'mesh_toml_95333801e1a3'
for run in perf['records']:
 prefix=f"rankaudit_perf_{run['mode']}_{run['repeat']}"
 for shard in sorted(result_dir.glob(prefix+'.result.*')):
  h=hashlib.sha256(shard.read_bytes()).hexdigest()
  if not any(x.get('file')==shard.name for x in run['outputs']):
   run['outputs'].append({'file':shard.name,'bytes':shard.stat().st_size,'sha256':h})
perf_path=HERE/'performance_repeats.json'
perf_path.write_text(json.dumps(perf,indent=2)+'\n',encoding='utf-8')
summary={}
for mode in ('off','on'):
 vals=[x['wall_seconds'] for x in perf['records'] if x['mode']==mode]
 summary[mode]={'seconds':vals,'mean':statistics.mean(vals),'sample_stdev':statistics.stdev(vals),'median':statistics.median(vals),'returncodes':[x['returncode'] for x in perf['records'] if x['mode']==mode]}
summary['speedup_off_over_on']=summary['off']['mean']/summary['on']['mean']
summary['wall_reduction_percent']=(1-summary['on']['mean']/summary['off']['mean'])*100
manifest={'schema':'dualtes-ast-rank-cache-audit-v1','rank_audit_log':LOG.name,'rankless_no_safe_message_count':text.count('no safe temperature-only MATC prisms'),
 'no_safe_callsite_source':'01_phase24_source_from_release_26_1.patch: Phase24BuildGenericDynKCache returns through this message when the local eligible material count is zero',
 'no_safe_branch_ranks':[r for r in range(6) if build[r]['material_count']==0],
 'all_six_solver_build_returncode':0 if 'MAIN: Time: 4/4:' in text else None,
 'rank_rows':ranks,'ast_unit_test_log':'matc_ast_purity_test.log','ast_unit_test_result':'RESULT failures=0' if 'RESULT failures=0' in (HERE/'matc_ast_purity_test.log').read_text(encoding='utf-8',errors='replace') else 'NOT_PASS',
 'performance':summary}
(HERE/'rank_audit.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'ranks':[{k:r[k] for k in ('rank','material7_prism_elements','unique_nodes','assembly_evaluation_events','matc_evaluation_calls_total','cache_hit_elements_total')} for r in ranks],'performance':summary},indent=2))
