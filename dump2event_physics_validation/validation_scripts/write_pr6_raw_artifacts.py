import csv, json
from pathlib import Path
base=Path('dump2event_physics_validation')
# Original line numbers and exact reaction blocks from the independently extracted history files.
def raw_map(event):
 d={}
 for s in (base/f'raw_history_NOCAS_{event}.txt').read_text(encoding='ascii').splitlines()[1:]:
  try:n,t=s.split(': ',1); d[int(n)]=t
  except:pass
 return d
r318=raw_map(318209); r119=raw_map(119953)
early={}
with Path(r'H:\\hata\\1332_142_136\\1332keV_6\\dumpall.dat').open('r',encoding='ascii') as source:
 for i,line in enumerate(source,1):
  if i in (2,3,63,64): early[i]=line.rstrip('\\r\\n')
  if i>=64: break
with (base/'raw_dumpall_pr6_electron_annihilation_excerpts.txt').open('w',encoding='ascii',newline='\n') as f:
 f.write('PHITS Dumpall raw excerpts with original 1-based line numbers\n')
 f.write('Source: H:\\hata\\1332_142_136\\1332keV_6\\dumpall.dat (116,820,667,876 bytes)\n')
 f.write('\n## RCASC-only candidate headers around start of file\n')
 for n in (2,3,63,64): f.write(f'{n}: {early.get(n,"")}\n')
 f.write('\n## Audit event_id 318209: exact NOCAS header and NCOL=14 reaction 7\n')
 for n in range(612,614): f.write(f'{n}: {r318[n]}\n')
 for n in range(849,885): f.write(f'{n}: {r318[n]}\n')
 f.write('\n## Audit event_id 119953: exact NOCAS header and NCOL=13 reaction 404\n')
 for n in range(5464052,5464054): f.write(f'{n}: {r119.get(n,'')}\n')
 for n in range(5475195,5475241): f.write(f'{n}: {r119.get(n,'')}\n')

# Identity map, explicitly distinguishing NOCAS match from RCASC-only aliases.
with (base/'history_id_mapping.csv').open('w',newline='',encoding='utf-8') as f:
 w=csv.writer(f); w.writerow(['audit_event_id','converter_event_id_source','raw_NOCAS','NOBCH','raw_RCASC','RSOUIN','NCOL4_marker_line','header_data_line','mapping_status','notes'])
 w.writerow([318209,'NCOL=4 column 1 (NOCAS)',318209,1,318209,12,612,613,'exact NOCAS match','7th NCOL=13/14 reaction in this history is particle 2, line 849'])
 w.writerow([318209,'NCOL=4 column 1 (NOCAS)',191679,1,318209,2,2,3,'RCASC-only false match','Prior extraction candidate; NOCAS differs'])
 w.writerow([318209,'NCOL=4 column 1 (NOCAS)',156535,1,318209,1,63,64,'RCASC-only false match','Prior extraction candidate; NOCAS differs'])
 w.writerow([119953,'NCOL=4 column 1 (NOCAS)',119953,1,119953,1102,5464052,5464053,'exact NOCAS match','404th NCOL=13/14 reaction is particle 3, line 5475195'])

# Parent/child comparison. Energies are QCLUSTS(7) kinetic energy in MeV; weights are QCLUSTS(8).
parents={
 '318209:2:7':dict(event_id=318209,particle_id=2,reaction_index=7,ncol=14,jcoll=14,parent_ityp=12,parent_E=1.00643799692798,parent_EC=1.00643799692798,parent_WT=1.0,parent_OLDWT=1.0,parent_marker_line=849,parent_weight_line=857,parent_E_line=859,parent_EC_line=861,audit_parent_E=1.00643799692798,audit_reported_deposit=-0.042497415724038756,formula_note='E and EC are kinetic energy for this electron track; child kinetic energies are compared at the same NCOL=14 reaction point'),
 '119953:3:404':dict(event_id=119953,particle_id=3,reaction_index=404,ncol=13,jcoll=2,parent_ityp=13,parent_E=0.001000000000001,parent_EC=0.001,parent_WT=1.0,parent_OLDWT=1.0,parent_marker_line=5475195,parent_weight_line=5475203,parent_E_line=5475205,parent_EC_line=5475207,audit_parent_E=0.001000000000001,audit_reported_deposit=-1.020997803999999,formula_note='Parent and daughter values are kinetic energy; two photons carry annihilation rest-mass energy; dead positron has QCLUSTS(8)=1 but is excluded from the transported-only sum'),
}
children={
 '318209:2:7':[(12,0,1.0475541999824101,1.0,870,872),(12,0,0.0013812126696086671,1.0,876,881)],
 '119953:3:404':[(14,0,0.510998902,1.0,5213,5218),(14,0,0.510998902,1.0,5222,5227),(13,-1,0.001,1.0,5231,5236)],
}
audit_children={
 '318209:2:7':[1.0475541999824101,0.0013812126696086671],
 '119953:3:404':[0.510998902,0.510998902,0.001],
}
fields=['audit_event_id','particle_id','history_reaction_index','NCOL','JCOLL','parent_ITYP','parent_E_MeV','parent_EC_MeV','parent_WT','parent_OLDWT','child_index','child_ITYP','child_JCLUSTS4_status','child_QCLUSTS7_kinetic_MeV','child_QCLUSTS8_weight','child_weighted_kinetic_MeV','parent_E_times_OLDWT_MeV_weighted','parent_EC_times_WT_MeV_weighted','all_child_energy_unweighted_sum_MeV','all_child_energy_weighted_sum_MeV','incoming_E_minus_all_children_unweighted_MeV','incoming_E_times_OLDWT_minus_all_children_weighted_MeV','collision_EC_minus_all_children_unweighted_MeV','collision_EC_times_WT_minus_all_children_weighted_MeV','transported_child_energy_unweighted_sum_MeV','transported_child_energy_weighted_sum_MeV','audit_parent_E_MeV','audit_child_energy_MeV','audit_reported_deposit_MeV','raw_vs_audit_parent_delta_MeV','raw_vs_audit_children_delta_MeV','raw_vs_audit_deposit_delta_MeV','parent_weight_raw_line','parent_E_raw_line','parent_EC_raw_line','child_status_raw_line','child_energy_weight_raw_line','energy_definition']
with (base/'parent_child_energy_weight_comparison.csv').open('w',newline='',encoding='utf-8') as f:
 w=csv.DictWriter(f,fieldnames=fields); w.writeheader()
 for key,p in parents.items():
  ch=children[key]; all_u=sum(x[2] for x in ch); all_w=sum(x[2]*x[3] for x in ch)
  trans=[x for x in ch if x[1]>=0]
  tu=sum(x[2] for x in trans); tw=sum(x[2]*x[3] for x in trans)
  aa=audit_children[key]; ae=sum(aa)
  for idx,(ityp,status,e,wt,status_line,energy_line) in enumerate(ch):
   row={k:p[k] for k in ['event_id','particle_id','reaction_index','ncol','jcoll','parent_ityp','parent_E','parent_EC','parent_WT','parent_OLDWT']}
   row.update({'child_index':idx,'child_ITYP':ityp,'child_JCLUSTS4_status':status,'child_QCLUSTS7_kinetic_MeV':e,'child_QCLUSTS8_weight':wt,'child_weighted_kinetic_MeV':e*wt,'parent_E_times_OLDWT_MeV_weighted':p['parent_E']*p['parent_OLDWT'],'parent_EC_times_WT_MeV_weighted':p['parent_EC']*p['parent_WT'],'all_child_energy_unweighted_sum_MeV':all_u,'all_child_energy_weighted_sum_MeV':all_w,'incoming_E_minus_all_children_unweighted_MeV':p['parent_E']-all_u,'incoming_E_times_OLDWT_minus_all_children_weighted_MeV':p['parent_E']*p['parent_OLDWT']-all_w,'collision_EC_minus_all_children_unweighted_MeV':p['parent_EC']-all_u,'collision_EC_times_WT_minus_all_children_weighted_MeV':p['parent_EC']*p['parent_WT']-all_w,'transported_child_energy_unweighted_sum_MeV':tu,'transported_child_energy_weighted_sum_MeV':tw,'audit_parent_E_MeV':p['audit_parent_E'],'audit_child_energy_MeV':ae,'audit_reported_deposit_MeV':p['audit_reported_deposit'],'raw_vs_audit_parent_delta_MeV':p['parent_E']-p['audit_parent_E'],'raw_vs_audit_children_delta_MeV':all_u-ae,'raw_vs_audit_deposit_delta_MeV':(p['parent_E']-tw)-p['audit_reported_deposit'],'parent_weight_raw_line':p['parent_weight_line'],'parent_E_raw_line':p['parent_E_line'],'parent_EC_raw_line':p['parent_EC_line'],'child_status_raw_line':status_line,'child_energy_weight_raw_line':energy_line,'energy_definition':p['formula_note']})
   row['event_id']=p['event_id']; row['reaction_index']=p['reaction_index']; row['ncol']=p['ncol']; row['jcoll']=p['jcoll']; row['parent_ityp']=p['parent_ityp']; row['parent_E']=p['parent_E']; row['parent_EC']=p['parent_EC']; row['parent_WT']=p['parent_WT']; row['parent_OLDWT']=p['parent_OLDWT']
   for key2,val in list(row.items()):
    aliases={'event_id':'audit_event_id','reaction_index':'history_reaction_index','ncol':'NCOL','jcoll':'JCOLL','parent_ityp':'parent_ITYP','parent_E':'parent_E_MeV','parent_EC':'parent_EC_MeV','parent_WT':'parent_WT','parent_OLDWT':'parent_OLDWT'}
    if key2 in aliases: row[aliases[key2]]=row.pop(key2)
   w.writerow(row)



