#!/usr/bin/env python3
"""Compare dead e+ child kinetic residuals against same-run PHITS tlv."""
import argparse,csv
from pathlib import Path

def main():
    p=argparse.ArgumentParser(); p.add_argument('--mapping',type=Path,required=True); p.add_argument('--raw',type=Path,required=True); p.add_argument('--out',type=Path,required=True); p.add_argument('--tol',type=float,default=1e-12); a=p.parse_args()
    internal={}
    with a.mapping.open(newline='',encoding='utf-8') as f:
        for r in csv.DictReader(f):
            if r['process']=='unsupported_secondary' and r['NCOL']=='13' and r['ITYP']=='13':
                internal[(int(r['NOCAS']),int(r['source_record_id']))]=float(r['PHITS_internal_tlv_MeV'])
    children={}
    with a.raw.open(newline='',encoding='utf-8') as f:
        for r in csv.DictReader(f):
            key=(int(r['NOCAS']),int(r['source_record_id']))
            if key in internal and r['child_ITYP']=='13' and r['child_state']=='-1':
                children[key]=r
    a.out.parent.mkdir(parents=True,exist_ok=True)
    n=ok=0; maxerr=totalerr=0.0
    with a.out.open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f); w.writerow(['NOCAS','source_record_id','QCLUSTS7_kinetic_MeV','QCLUSTS8_weight','rest_mass_GeV_not_scored','internal_tlv_MeV','abs_error_MeV','region','child_x_cm','child_y_cm','child_z_cm','pass_1e-12_MeV'])
        for key,internal_tlv in sorted(internal.items()):
            r=children.get(key)
            if r is None: continue
            q7=float(r['child_kinetic_energy_MeV']); err=abs(q7-internal_tlv); n+=1; ok+=err<=a.tol; totalerr+=err; maxerr=max(maxerr,err)
            w.writerow([*key,format(q7,'.17g'),r['child_weight'],r['child_rest_mass_GeV'],format(internal_tlv,'.17g'),format(err,'.17g'),r['region_before'],r['child_x_cm'],r['child_y_cm'],r['child_z_cm'],int(err<=a.tol)])
    print(f'dead_positron_secondaries={n} pass={ok} max_abs_error_MeV={maxerr:.17g} sum_abs_error_MeV={totalerr:.17g}')
if __name__=='__main__': main()
