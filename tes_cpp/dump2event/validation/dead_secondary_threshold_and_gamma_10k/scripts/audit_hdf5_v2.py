#!/usr/bin/env python3
"""Bounded-memory completeness audit for an experimental HDF5 v2 file."""
import argparse,json
from collections import Counter
from pathlib import Path
import h5py

def s(x): return x.decode('utf-8') if isinstance(x,bytes) else str(x)
def main():
 p=argparse.ArgumentParser();p.add_argument('--hdf5',type=Path,required=True);p.add_argument('--output-json',type=Path,required=True);a=p.parse_args()
 c=Counter(); expected=1; max_id_gap=0; min_id=None; max_id=None
 with h5py.File(a.hdf5,'r') as h:
  c['events']=len(h['events/event_id']);c['particles']=len(h['particles/particle_id']);c['tracks']=len(h['tracks/particle_id'])
  for base in range(0,c['events'],65536):
   end=min(base+65536,c['events'])
   for complete,unresolved,absorption in zip(h['events/complete'][base:end],h['events/unresolved_count'][base:end],h['events/all_absorption'][base:end]):
    c['complete_events' if bool(complete) else 'incomplete_events']+=1
    c['events_with_unresolved']+=bool(unresolved)
    c['all_absorption_'+s(absorption)]+=1
  d=h['deposits'];n=len(d['addition_id']);c['deposits']=n
  for base in range(0,n,65536):
   end=min(base+65536,n)
   ids=d['addition_id'][base:end]; energies=d['energy_MeV'][base:end];valid=d['energy_valid'][base:end]; kinds=d['kind'][base:end]
   for x,e,v,k in zip(ids,energies,valid,kinds):
    value=int(s(x));c['addition_ids_checked']+=1
    min_id=value if min_id is None else min(min_id,value)
    max_id=value if max_id is None else max(max_id,value)
    if value!=expected:c['addition_id_sequence_mismatch']+=1;max_id_gap=max(max_id_gap,abs(value-expected))
    expected+=1
    if not bool(v):c['invalid_energy']+=1
    elif float(e)<0:c['negative_energy']+=1
    elif float(e)==0:c['zero_energy']+=1
    else:c['positive_energy']+=1
    c['kind_'+s(k)]+=1
 result=dict(c);result['max_addition_id_index_difference']=max_id_gap
 result['addition_id_min']=min_id;result['addition_id_max']=max_id
 result['addition_id_sequence_expected_count']=expected-1
 a.output_json.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps(result,indent=2))
 return 0 if c['deposits']==c['addition_ids_checked'] and c['addition_id_sequence_mismatch']==0 and c['negative_energy']==0 and c['invalid_energy']==0 else 1
if __name__=='__main__':raise SystemExit(main())
