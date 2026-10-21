#!/usr/bin/env python3
"""Sequentially compare candidate raw-only additions with PHITS PRE trace."""
import argparse,csv,math
from pathlib import Path

FIELDS=("label nocas nobch rcasc particle_id ityp jcoll ncol material_id tally_region_index "
        "tally_region_id region_before region_after oldwt E EC ein ecc edep tlv tally_unit "
        "species_slot lattice_slot time_slot x_start y_start z_start x_end y_end z_end tally_before").split()

def internal(path):
    with path.open(encoding='utf-8',errors='replace') as f:
        for ln,line in enumerate(f,1):
            if not line.startswith(' PRE'): continue
            vals=line.split()
            if len(vals)!=len(FIELDS): raise ValueError(f'{path}:{ln}: expected {len(FIELDS)} fields, got {len(vals)}')
            r=dict(zip(FIELDS,vals)); yield ln,r

def main():
    p=argparse.ArgumentParser(); p.add_argument('--ledger',type=Path,required=True); p.add_argument('--internal',type=Path,required=True); p.add_argument('--out',type=Path,required=True); p.add_argument('--energy-tol',type=float,default=1e-12); p.add_argument('--position-tol',type=float,default=1e-12); a=p.parse_args()
    cnt=unresolved=energy_fail=region_fail=position_fail=positive=0; maxe=maxp=0.; sums={}; byproc={}; examples=[]
    it=iter(internal(a.internal))
    with a.ledger.open(newline='',encoding='utf-8') as lf:
      rd=csv.DictReader(lf)
      for row in rd:
        try: line,ir=next(it)
        except StopIteration: raise ValueError(f'internal trace ended before addition {row["addition_id"]}')
        cnt+=1; proc=row['process']; group=byproc.setdefault(proc,[0,0,0.,0.,0])
        group[0]+=1
        if row['status']!='reconstructed' or not row['E_deposit_MeV']:
            unresolved+=1; continue
        val=float(row['E_deposit_MeV']); tlv=float(ir['tlv']); err=abs(val-tlv); maxe=max(maxe,err); group[2]+=val; group[3]=max(group[3],err)
        if err>a.energy_tol:
            energy_fail+=1
            if len(examples)<30: examples.append((row['addition_id'],row['event_id'],proc,val,tlv,err,line))
        else: group[1]+=1
        if int(row['region_id'])!=int(ir['tally_region_id']): region_fail+=1
        if val>0:
            positive+=1
            typ=row['position_status']
            if typ in ('segment','segment_to_boundary'):
                coords=[float(row[k]) for k in ('x','y','z','x_end','y_end','z_end')]
                refs=[float(ir[k]) for k in ('x_start','y_start','z_start','x_end','y_end','z_end')]
            elif typ=='secondary_point':
                coords=[float(row[k]) for k in ('x','y','z')]
                # For a dead secondary the scoring coordinate is the stopped
                # particle's location; in these records it equals parent end.
                refs=[float(ir[k]) for k in ('x_end','y_end','z_end')]
            else:
                position_fail+=1; continue
            pe=max(abs(x-y) for x,y in zip(coords,refs)); maxp=max(maxp,pe)
            if pe>a.position_tol:
                position_fail+=1
                if len(examples)<30: examples.append((row['addition_id'],row['event_id'],proc,'position',pe,line))
        elif row['position_status']!='none':
            # zero-value records have no score coordinate in the physical ledger
            if row['position_status']=='unresolved': position_fail+=1
    try: extra=next(it)
    except StopIteration: extra=None
    if extra: raise ValueError(f'internal trace has extra score calls after candidate ledger: line {extra[0]}')
    a.out.parent.mkdir(parents=True,exist_ok=True)
    with a.out.open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f); w.writerow(['metric','value']);
        for k,v in [('candidate_additions',cnt),('unresolved',unresolved),('energy_mismatch_over_tolerance',energy_fail),('region_mismatch',region_fail),('positive_deposits',positive),('positive_position_mismatch_over_tolerance',position_fail),('max_abs_energy_error_MeV',format(maxe,'.17g')),('max_abs_position_error_cm',format(maxp,'.17g')),('energy_tolerance_MeV',a.energy_tol),('position_tolerance_cm',a.position_tol)]: w.writerow([k,v])
        for k,v in sorted(byproc.items()):
            w.writerow([f'process_count_{k}',v[0]]); w.writerow([f'process_energy_match_count_{k}',v[1]]); w.writerow([f'process_candidate_energy_sum_{k}',format(v[2],'.17g')]); w.writerow([f'process_max_error_{k}',format(v[3],'.17g')])
    for x in examples: print('EXAMPLE',*x)
    print(f'additions={cnt} unresolved={unresolved} energy_fail={energy_fail} region_fail={region_fail} position_fail={position_fail} positive={positive} max_energy_error={maxe:.17g} max_position_error_cm={maxp:.17g}')
if __name__=='__main__': main()
