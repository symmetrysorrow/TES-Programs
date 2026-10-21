#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,json,math
from pathlib import Path
from collections import Counter
import h5py

def main():
 p=argparse.ArgumentParser();p.add_argument('--tally',type=Path,required=True);p.add_argument('--hdf5',type=Path,required=True);p.add_argument('--output-csv',type=Path,required=True);p.add_argument('--output-json',type=Path,required=True);a=p.parse_args()
 rows=[];data=False;started=False
 for line in a.tally.read_text(encoding='utf-8',errors='replace').splitlines():
  if line.lstrip().startswith('h: n'): data=True;continue
  if not data or not line.strip() or line.lstrip().startswith('#'):continue
  s=line.split()
  try: rows.append((len(rows),float(s[0]),float(s[1]),float(s[2])));started=True
  except (ValueError,IndexError):
   if started:break
 bins=[(i,lo,hi,y) for i,lo,hi,y in rows if hi>lo]
 with h5py.File(a.hdf5,'r') as f:
  ids=f['events/event_id'][:];totals=f['events/deposit_total_MeV'][:];complete=f['events/complete'][:]
 if not bins:raise ValueError('no histogram data rows found')
 n=len(ids);width=bins[0][2]-bins[0][1];counts=Counter();under=over=zero=0
 for e in totals:
  if float(e)==0.0:
   zero+=1
   continue
  if e<bins[0][1]:under+=1
  elif e>=bins[-1][2]:over+=1
  else:
   i=math.floor((float(e)-bins[0][1])/width)
   if i<0:under+=1
   elif i>=len(bins):over+=1
   else:counts[i]+=1
 a.output_csv.parent.mkdir(parents=True,exist_ok=True);absdiff=[]
 with a.output_csv.open('w',newline='',encoding='utf-8') as f:
  w=csv.writer(f);w.writerow(['bin','low_MeV','high_MeV','PHITS_counts_from_unit3','HDF5_event_counts','difference_counts','relative_difference'])
  for i,(num,lo,hi,y) in enumerate(bins):
   phits=y*n;hv=counts[i];diff=hv-phits;absdiff.append(abs(diff));rel=diff/phits if phits else (0 if hv==0 else None)
   w.writerow([num,format(lo,'.9g'),format(hi,'.9g'),format(phits,'.9g'),hv,format(diff,'.9g'),'' if rel is None else format(rel,'.9g')])
 mean=sum(map(float,totals))/n if n else 0;variance=sum((float(x)-mean)**2 for x in totals)/(n-1) if n>1 else 0
 peak=max(range(len(bins)),key=lambda i:counts[i])
 summary={'histories':n,'all_hdf5_events_complete':bool(all(complete)),'bin_count':len(bins),'bin_width_MeV':width,'zero_deposit_histories_excluded':zero,'zero_deposit_histories_excluded':zero,'underflow_histories':under,'overflow_histories':over,'nonzero_hdf5_bins':sum(v>0 for v in counts.values()),'PHITS_count_total_reconstructed_from_unit3':sum(row[3]*n for row in bins),'HDF5_in_range_count':sum(counts.values()),'bin_count_L1_abs_difference':sum(abs(counts[i]-bins[i][3]*n) for i in range(len(bins))),'max_abs_count_difference_per_bin':max(absdiff,default=0),'mean_event_deposit_MeV':mean,'sample_sd_event_deposit_MeV':math.sqrt(variance),'standard_error_MeV':math.sqrt(variance/n) if n else None,'max_event_deposit_MeV':max(map(float,totals),default=None),'peak_bin_index':peak,'peak_bin_low_MeV':bins[peak][1],'peak_bin_high_MeV':bins[peak][2],'bins_at_phits_edges':[bins[0][1],bins[-1][2]]}
 a.output_json.write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
