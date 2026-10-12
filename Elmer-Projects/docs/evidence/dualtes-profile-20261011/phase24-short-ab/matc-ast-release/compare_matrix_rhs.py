from pathlib import Path
import json
HERE=Path(__file__).resolve().parent
R=Path(r'C:\Users\Nao\AppData\Local\Temp\phase24-cache-runtime')
d=json.loads((HERE/'matrix_rhs_capture_runs.json').read_text()); runs={x['mode']:x for x in d['runs']}; a,b=runs['off'],runs['on']; aa={x['relative']:x for x in a['capture_files']}; bb={x['relative']:x for x in b['capture_files']}; common=sorted(aa.keys()&bb.keys()); missing=sorted(aa.keys()^bb.keys()); nonrhs=[]; rhsdiff=0; rhsvalues=0; maxabs=maxrel=0.0; rowmap=True; rhsfiles=0
for rel in common:
 name=rel.rsplit('/',1)[-1]
 if '_b.dat.' in name:
  rhsfiles+=1; pa=R/'final_matcap_off'/rel; pb=R/'final_matcap_on'/rel
  va={int(z.split()[0]):float(z.split()[1]) for z in pa.read_text().splitlines()}; vb={int(z.split()[0]):float(z.split()[1]) for z in pb.read_text().splitlines()}
  if va.keys()!=vb.keys(): rowmap=False; continue
  rhsvalues+=len(va)
  for k,x in va.items():
   y=vb[k]; e=abs(x-y); maxabs=max(maxabs,e); maxrel=max(maxrel,e/max(abs(x),abs(y),1e-300)); rhsdiff+=e!=0.0
 elif aa[rel]['sha256']!=bb[rel]['sha256'] or aa[rel]['bytes']!=bb[rel]['bytes']: nonrhs.append(rel)
out={'comparison':'final corrected MATC purity gate, reconstructed Release, DualTES timestep 1 assembly capture','off_rc':a['returncode'],'on_rc':b['returncode'],'off_capture_files':len(aa),'on_capture_files':len(bb),'matched_files':len(common),'missing_files':missing,'matrix_and_sparsity_metadata_byte_identical':not nonrhs,'differing_non_rhs_files':nonrhs,'rhs_files_compared':rhsfiles,'rhs_values_compared':rhsvalues,'rhs_changed_values':rhsdiff,'rhs_row_maps_identical':rowmap,'rhs_max_abs_difference':maxabs,'rhs_max_relative_difference':maxrel}
(HERE/'matrix_rhs_comparison.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8'); print(json.dumps(out,indent=2))
