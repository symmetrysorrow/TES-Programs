#!/usr/bin/env python3
"""Evaluate PHITS 3.341's positron EGS5 stopping formula using raw inputs.

The internal trace is used only as a comparison oracle. Predictions use raw
Dumpall E, WT, endpoints, plus explicit PHITS thresholds and the local HATCH
coefficient table. The large trace is streamed.
"""
from __future__ import annotations
import argparse, csv, math
from collections import defaultdict
from pathlib import Path

FIELDS = ("label nocas nobch rcasc particle_id ityp jcoll ncol material_id "
          "tally_region_index tally_region_id region_before region_after oldwt E EC ein ecc edep tlv "
          "tally_unit species_slot lattice_slot time_slot x_start y_start z_start x_end y_end z_end tally_before").split()

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--raw-tracks',type=Path,required=True)
    ap.add_argument('--internal-map',type=Path,required=True)
    ap.add_argument('--mismatches',type=Path,required=True)
    ap.add_argument('--coefficients',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--emin-mev',type=float,default=0.001)
    ap.add_argument('--tolerance-mev',type=float,default=1e-12)
    a=ap.parse_args()
    # HATCH table for explicit medium 2 (Pb); PHITS uses integer truncation
    # of eke1*ln(E)+eke0 to select the table row.
    coeff={}
    with a.coefficients.open(newline='',encoding='utf-8-sig') as f:
        for r in csv.DictReader(f):
            if int(r['medium'])==2:
                coeff[int(r['lelke'])]=(float(r['pdedx0']),float(r['pdedx1']))
                eke0=float(r['eke0']); eke1=float(r['eke1'])
    targets={}
    with a.internal_map.open(newline='',encoding='utf-8') as f:
        for r in csv.DictReader(f):
            if r['NCOL'] in ('13','14') and r['ITYP']=='13' and r['process']=='unsupported_primary':
                targets[(int(r['NOCAS']),int(r['source_record_id']))]=('unresolved_primary',float(r['PHITS_internal_tlv_MeV']))
    with a.mismatches.open(newline='',encoding='utf-8') as f:
        for r in csv.DictReader(f):
            if r['NCOL'] in ('13','14') and r['ITYP']=='13' and r['JCOLL'] in ('2','14'):
                targets[(int(r['NOCAS']),int(r['source_record_id']))]=('exact_zero_mismatch',float(r['internal_tlv_MeV']))
    raw={}
    with a.raw_tracks.open(newline='',encoding='utf-8') as f:
        for r in csv.DictReader(f):
            key=(int(r['NOCAS']),int(r['source_record_id']))
            if key in targets and key not in raw:
                raw[key]=r
    a.out.parent.mkdir(parents=True,exist_ok=True)
    groups=defaultdict(lambda:[0,0,0.0,0.0,0.0])
    rows=0
    with a.out.open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f)
        w.writerow(['NOCAS','source_record_id','class','NCOL','ITYP','JCOLL','material_id','region','E_MeV','EC_MeV','WT','OLDWT','raw_chord_cm','PHITS_log_energy','PHITS_energy_index','pdedx_MeV_per_cm','predicted_tlv_MeV','internal_tlv_MeV','abs_error_MeV','relative_error','within_tolerance','raw_x_cm','raw_y_cm','raw_z_cm','raw_xc_cm','raw_yc_cm','raw_zc_cm'])
        for (nocas,srid),(kind,internal) in sorted(targets.items()):
            r=raw.get((nocas,srid))
            if r is None: continue
            E=float(r['E_MeV']); wt=float(r['OLDWT']); e=max(E,a.emin_mev); loge=math.log(e)
            idx=int(eke1*loge+eke0)
            if idx not in coeff: raise ValueError(f'coefficient index {idx} missing for {(nocas,srid)}')
            p0,p1=coeff[idx]; dedx=p1*loge+p0
            chord=math.sqrt(sum((float(r[k])-float(r[kk]))**2 for k,kk in [('x_cm','xc_cm'),('y_cm','yc_cm'),('z_cm','zc_cm')]))
            pred=wt*dedx*chord
            err=abs(pred-internal); rel=err/max(abs(internal),1e-300); ok=err<=a.tolerance_mev
            g=groups[kind]; g[0]+=1; g[1]+=int(ok); g[2]+=err; g[3]=max(g[3],err); g[4]+=pred
            w.writerow([nocas,srid,kind,r['NCOL'],r['ITYP'],r['JCOLL'],r['material_id'],r['region_before'],format(E,'.17g'),r['EC_MeV'],r['WT'],r['OLDWT'],format(chord,'.17g'),format(loge,'.17g'),idx,format(dedx,'.17g'),format(pred,'.17g'),format(internal,'.17g'),format(err,'.17g'),format(rel,'.17g'),int(ok),r['x_cm'],r['y_cm'],r['z_cm'],r['xc_cm'],r['yc_cm'],r['zc_cm']]); rows+=1
    print(f'targets={len(targets)} raw_matched={rows} missing_raw={len(targets)-rows}')
    for k,v in sorted(groups.items()): print(f'{k}: n={v[0]} within_tol={v[1]} sum_abs_error={v[2]:.17g} max_abs_error={v[3]:.17g} predicted_sum={v[4]:.17g}')

if __name__=='__main__': main()
