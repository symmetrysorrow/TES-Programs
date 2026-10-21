#!/usr/bin/env python3
"""Summarize the nine histories that had unresolved positron rows before closure."""
import argparse,csv
from collections import defaultdict
from pathlib import Path

def main():
    p=argparse.ArgumentParser(); p.add_argument('--before',type=Path,required=True); p.add_argument('--after',type=Path,required=True); p.add_argument('--out',type=Path,required=True); a=p.parse_args()
    prior=defaultdict(dict)
    with a.before.open(newline='',encoding='utf-8') as f:
        for r in csv.DictReader(f):
            k=int(r['event_id_NOCAS']); prior[k]['rows']=int(r['unresolved_rows_in_history']); prior[k]['known']=float(r['known_deposit_lower_bound_MeV'])
    final={}
    with a.after.open(newline='',encoding='utf-8') as f:
        for r in csv.DictReader(f): final[int(r['NOCAS'])]=r
    a.out.parent.mkdir(parents=True,exist_ok=True)
    with a.out.open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f); w.writerow(['NOCAS','previous_unresolved_rows','previous_known_lower_bound_MeV','PHITS_internal_deposit_MeV','corrected_HDF5_deposit_MeV','absolute_difference_MeV','corrected_unresolved_rows','same_PHITS_bin'])
        for nocas in sorted(prior):
            r=final[nocas]
            w.writerow([nocas,prior[nocas]['rows'],format(prior[nocas]['known'],'.17g'),r['PHITS_internal_deposit_MeV'],r['HDF5_known_subtotal_MeV'],format(abs(float(r['internal_minus_known_MeV'])),'.17g'),r['unresolved_rows'],r['same_bin']])
    print(f'histories={len(prior)} prior_unresolved_rows={sum(x["rows"] for x in prior.values())}')
if __name__=='__main__': main()
