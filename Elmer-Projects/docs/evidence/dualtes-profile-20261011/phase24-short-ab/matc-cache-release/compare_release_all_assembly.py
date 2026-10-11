from pathlib import Path
import json, math
HERE=Path(__file__).resolve().parent
d=json.loads((HERE/'release_all_assembly_capture.json').read_text()); runs={x['label']:x for x in d['runs']}; steps=[]
def sha(x): return x['sha256']
for step in range(1,5):
 a=runs[f'off_ts{step}']; b=runs[f'on_ts{step}']; aa={x['relative']:x for x in a['capture_files']}; bb={x['relative']:x for x in b['capture_files']}; common=sorted(aa.keys()&bb.keys()); missing=sorted(aa.keys()^bb.keys()); nonrhs=[]; rhsmax=0.0; rhsrel=0.0; rhsdiff=0; rhsfiles=0; rhsvalues=0; rowmap=True; rhsdetails=[]
 for rel in common:
  if '/outer_' not in rel: continue
  fname=rel.rsplit('/',1)[-1]
  if '_b.dat.' in fname:
   rhsfiles+=1; pa=Path(r'C:\Users\Nao\AppData\Local\Temp\phase24-cache-runtime')/f'release_capture_off_ts{step}'/rel; pb=Path(r'C:\Users\Nao\AppData\Local\Temp\phase24-cache-runtime')/f'release_capture_on_ts{step}'/rel
   va={int(z.split()[0]):float(z.split()[1]) for z in pa.read_text().splitlines()}; vb={int(z.split()[0]):float(z.split()[1]) for z in pb.read_text().splitlines()}
   if va.keys()!=vb.keys(): rowmap=False; continue
   rhsvalues+=len(va)
   for k,x in va.items():
    y=vb[k]; delta=abs(x-y); relerr=delta/max(abs(x),abs(y),1e-300); rhsmax=max(rhsmax,delta); rhsrel=max(rhsrel,relerr)
    if delta: rhsdiff+=1
   if sha(aa[rel])!=sha(bb[rel]): rhsdetails.append(rel)
  else:
   if aa[rel]['sha256']!=bb[rel]['sha256'] or aa[rel]['bytes']!=bb[rel]['bytes']: nonrhs.append(rel)
 steps.append({'timestep':step,'off_rc':a['returncode'],'on_rc':b['returncode'],'off_files':len(aa),'on_files':len(bb),'matched_files':len(common),'missing_files':missing,'matrix_and_structure_files_byte_identical':not nonrhs,'nonrhs_different_files':nonrhs,'rhs_files_compared':rhsfiles,'rhs_values_compared':rhsvalues,'rhs_values_changed_at_bit_level':rhsdiff,'rhs_shards_with_hash_change':rhsdetails,'rhs_row_maps_identical':rowmap,'rhs_max_abs_difference':rhsmax,'rhs_max_relative_difference':rhsrel})
result={'comparison':'outer DefaultSolve before/after; all nonlinear assemblies at selected timestep in isolated Release runs','matrix_coefficients':'byte-identical','sparsity_and_size_metadata':'byte-identical','rhs':'node-row keyed numeric comparison','steps':steps,'all_pattern_and_matrix_exact':all(x['matrix_and_structure_files_byte_identical'] and not x['missing_files'] for x in steps),'rhs_within_roundoff':all(x['rhs_row_maps_identical'] and x['rhs_max_relative_difference'] <= 1e-15 for x in steps)}
(HERE/'release_all_assembly_comparison.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8'); print(json.dumps(result,indent=2))
