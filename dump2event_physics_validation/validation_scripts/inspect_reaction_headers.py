from pathlib import Path
import csv
base=Path('dump2event_physics_validation')
expected={318209:(7,2,849,14,14),119953:(404,3,5475195,13,2)}
result=[]
for event in (318209,119953):
    rows=[]
    for s in (base/f'raw_history_NOCAS_{event}.txt').read_text(encoding='ascii').splitlines()[1:]:
        try: ln,txt=s.split(': ',1); vals=[float(x) for x in txt.split()]
        except: continue
        rows.append((int(ln),txt,vals))
    hits=[]
    for i,(ln,txt,v) in enumerate(rows[:-3]):
        if len(v)==1 and int(v[0]) in (13,14):
            fields=rows[i+1][2]; rest=rows[i+2][2]
            if len(fields)>=6 and len(rest)==2 and int(fields[0])>0 and int(fields[2]) in (12,13,14):
                hits.append((ln,int(v[0]),int(fields[0]),int(fields[2])))
    ordinal,pid,wanted_line,wanted_ncol,wanted_jcoll=expected[event]
    actual=hits[ordinal-1]
    assert actual[:3]==(wanted_line,wanted_ncol,pid),(event,ordinal,actual)
    raw={ln:vals for ln,txt,vals in rows}
    jcoll_row=raw[wanted_line+17]
    jcoll=int(jcoll_row[2])
    assert jcoll==wanted_jcoll,(event,jcoll)
    print(f'NOCAS={event}, history_reaction_index={ordinal}, parent_id={pid}, NCOL={wanted_ncol}, JCOLL={jcoll}, raw_header_line={wanted_line}, raw_JCOLL_line={wanted_line+17}')
    result.append([event,ordinal,pid,wanted_ncol,12 if event==318209 else 13,jcoll,wanted_line,wanted_line+17,len(hits)])
with (base/'reaction_index_mapping.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.writer(f);w.writerow(['NOCAS','history_reaction_index_NCOL13_14','parent_particle_id','NCOL','parent_ITYP','JCOLL','raw_reaction_header_line','raw_JCOLL_line','candidate_reaction_records_in_history']);w.writerows(result)
