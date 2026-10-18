import csv
from pathlib import Path
base=Path('dump2event_physics_validation')
rows=list(csv.DictReader((base/'header_identity_matches.csv').open(encoding='utf-8')))
with (base/'observed_nocas_rcasc_mismatches.csv').open('w',newline='',encoding='utf-8') as f:
 w=csv.writer(f);w.writerow(['NCOL4_marker_line','header_data_line','NOCAS','NOBCH','RCASC','RSOUIN','maxcas_from_phits_out','formula_expected_RCASC','observed_minus_expected','scan_scope'])
 for r in rows:
  nocas=int(r['NOCAS']); nobch=int(r['NOBCH']); rcasc=float(r['RCASC'])
  if nobch==1 and nocas!=int(rcasc): w.writerow([r['source_line'],int(r['source_line'])+1,nocas,nobch,rcasc,r['RSOUIN'],1000000,nocas+1000000*(nobch-1),rcasc-(nocas+1000000*(nobch-1)),'headers observed through physical line 5515858; not a whole-file census'])
# Add all observed RCASC-only rows to identity table.
p=base/'history_id_mapping.csv'; old=list(csv.DictReader(p.open(encoding='utf-8')))
known={(x['raw_NOCAS'],x['header_data_line']) for x in old}
for r in rows:
 nocas=r['NOCAS']; data=str(int(r['source_line'])+1)
 if float(r['RCASC'])==318209 and nocas!='318209' and (nocas,data) not in known:
  old.append({'audit_event_id':'318209','converter_event_id_source':'NCOL=4 column 1 (NOCAS)','raw_NOCAS':nocas,'NOBCH':r['NOBCH'],'raw_RCASC':r['RCASC'],'RSOUIN':r['RSOUIN'],'NCOL4_marker_line':r['source_line'],'header_data_line':data,'mapping_status':'RCASC-only false match','notes':'Observed in scanned prefix; NOCAS differs, cannot map to audit event_id by converter rule'})
with p.open('w',newline='',encoding='utf-8') as f:
 w=csv.DictWriter(f,fieldnames=['audit_event_id','converter_event_id_source','raw_NOCAS','NOBCH','raw_RCASC','RSOUIN','NCOL4_marker_line','header_data_line','mapping_status','notes']);w.writeheader();w.writerows(old)
