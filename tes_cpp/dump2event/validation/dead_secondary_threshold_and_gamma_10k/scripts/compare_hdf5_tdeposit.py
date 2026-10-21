#!/usr/bin/env python3
"""Compare HDF5 v2 per-history sums with PHITS unit=3 [t-deposit] bins."""
from __future__ import annotations
import argparse, csv, json, math
from collections import Counter, defaultdict
from pathlib import Path
import h5py

def phits_bins(path: Path):
    bins=[]; started=False
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.lstrip().startswith("h: n"):
            started=True; continue
        if not started or not line.strip() or line.lstrip().startswith("#"):
            continue
        fields=line.split()
        try: lo,hi,score=map(float,fields[:3])
        except (ValueError,IndexError):
            if bins: break
            continue
        if hi>lo: bins.append((lo,hi,score))
    if not bins: raise ValueError("PHITS tally has no histogram bins")
    return bins

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--hdf5",type=Path,required=True)
    ap.add_argument("--tally",type=Path,required=True)
    ap.add_argument("--histories",type=int,required=True)
    ap.add_argument("--output-csv",type=Path,required=True)
    ap.add_argument("--output-json",type=Path,required=True)
    a=ap.parse_args()
    bins=phits_bins(a.tally); total=defaultdict(float); invalid=0
    with h5py.File(a.hdf5,"r") as h:
        event_ids=[int(x) for x in h["events/event_id"][:]]
        d=h["deposits"]
        deposit_count=len(d["event_id"])
        for start in range(0,deposit_count,65536):
            stop=min(start+65536,deposit_count)
            ids=d["event_id"][start:stop]; energies=d["energy_MeV"][start:stop]
            valid=d["energy_valid"][start:stop]
            for eid,e,ok in zip(ids,energies,valid):
                if not bool(ok): invalid+=1; continue
                total[int(eid)]+=float(e)
    if len(event_ids)!=a.histories: raise ValueError(f"HDF5 events={len(event_ids)} expected={a.histories}")
    hcounts=Counter(); zero=under=over=0; values=[]; width=bins[0][1]-bins[0][0]
    for eid in event_ids:
        e=total.get(eid,0.0); values.append(e)
        if e==0: zero+=1
        elif e<bins[0][0]: under+=1
        elif e>=bins[-1][1]: over+=1
        else:
            i=math.floor((e-bins[0][0])/width)
            if i<0: under+=1
            elif i>=len(bins): over+=1
            else: hcounts[i]+=1
    diffs=[]; a.output_csv.parent.mkdir(parents=True,exist_ok=True)
    with a.output_csv.open("w",newline="",encoding="utf-8") as f:
        w=csv.writer(f); w.writerow(("bin","low_MeV","high_MeV","PHITS_counts","HDF5_counts","abs_diff","relative_diff"))
        for i,(lo,hi,score) in enumerate(bins):
            expected=score*a.histories; got=hcounts[i]; diff=got-expected; diffs.append(abs(diff))
            w.writerow((i,format(lo,".17g"),format(hi,".17g"),format(expected,".17g"),got,format(abs(diff),".17g"),format(diff/expected,".17g") if expected else ""))
    mean=sum(values)/a.histories
    variance=sum((v-mean)**2 for v in values)/(a.histories-1) if a.histories>1 else 0.0
    peak=max(range(len(bins)),key=lambda i:hcounts[i])
    result={"histories":a.histories,"hdf5_deposits":deposit_count,"invalid_deposit_rows":invalid,
      "zero_histories":zero,"underflow_histories":under,"overflow_histories":over,"bin_count":len(bins),
      "bin_width_MeV":width,"bin_count_L1_difference":sum(diffs),"max_abs_bin_count_difference":max(diffs,default=0.0),
      "mean_HDF5_deposit_MeV":mean,"sample_sd_MeV":math.sqrt(variance),"standard_error_MeV":math.sqrt(variance/a.histories),
      "peak_bin_low_MeV":bins[peak][0],"peak_bin_high_MeV":bins[peak][1],"max_event_deposit_MeV":max(values,default=0.0)}
    a.output_json.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2))
    return 0 if len(bins)==150 and invalid==0 else 1
if __name__=="__main__": raise SystemExit(main())
