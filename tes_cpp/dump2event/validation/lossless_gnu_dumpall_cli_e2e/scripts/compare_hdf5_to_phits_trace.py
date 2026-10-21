#!/usr/bin/env python3
"""Join Pb HDF5 v2 additions to PHITS's same-run raw sequence and internal trace."""
from __future__ import annotations
import argparse, csv, hashlib, json
from collections import Counter, defaultdict
from pathlib import Path
import h5py

def sha(path: Path) -> str:
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()
def rounded(v): return round(float(v), 13)
def raw_key(r):
    return (int(r['event_id']), int(r['particle_id']), int(r['ityp']), int(r['ncol']),
            rounded(r['energy_MeV']), rounded(r['x_cm']), rounded(r['y_cm']),
            rounded(r['z_cm']), int(r['region_id']))
def vec(row, names): return tuple(float(row[n]) for n in names)
def max_vec_delta(a,b): return max(abs(x-y) for x,y in zip(a,b))

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--run-dir',type=Path,required=True,help='same PHITS run directory')
    p.add_argument('--tracks-csv',type=Path,required=True,help='dumpall_transport_export tracks CSV')
    p.add_argument('--internal-csv',type=Path,required=True,help='paired PRE/POST trace CSV')
    p.add_argument('--hdf5',type=Path,required=True)
    p.add_argument('--output-csv',type=Path,required=True)
    p.add_argument('--output-json',type=Path,required=True)
    a=p.parse_args(); run=a.run_dir
    raw_tracks={}
    with a.tracks_csv.open(newline='',encoding='utf-8-sig') as f:
        for row in csv.DictReader(f): raw_tracks[int(row['source_record_id'])]=row
    # The PHITS DUMP sidecar gives a monotonically increasing record sequence.
    dump_by_key=defaultdict(list)
    with (run/'dump_record_sequence.log').open(encoding='ascii') as f:
        for line in f:
            s=line.split()
            if not s or s[0]!='DUMP': continue
            if len(s)<23: raise ValueError('short DUMP sidecar row')
            row={'event_id':s[2],'particle_id':s[5],'ityp':s[6],'ncol':s[7],
                 'energy_MeV':s[13],'x_cm':s[15],'y_cm':s[16],'z_cm':s[17],
                 'region_id':s[21]}
            dump_by_key[raw_key(row)].append(int(s[1]))
    line_to_seq={}; raw_ambiguous=0; raw_unmapped=0
    for line,row in raw_tracks.items():
        matches=dump_by_key.get(raw_key(row),[])
        if len(matches)==1: line_to_seq[line]=matches[0]
        elif matches: raw_ambiguous+=1
        else: raw_unmapped+=1
    seq=defaultdict(list); seq_order=[]
    with (run/'tally_dumpseq.log').open(encoding='ascii') as f:
        for line in f:
            s=line.split()
            if not s: continue
            if len(s)<17: raise ValueError('short TALLYSEQ row')
            serial=int(s[1]); seq_order.append(serial)
            seq[serial].append({'event_id':int(float(s[2])),'particle_id':int(s[5]),
                'ityp':int(s[6]),'jcoll':int(s[7]),'ncol':int(s[8]),
                'tally_region_id':int(s[10]),'weight':float(s[13]),
                'E_MeV':float(s[14]),'EC_MeV':float(s[15]),'tlv_MeV':float(s[16])})
    pre=list(csv.DictReader(a.internal_csv.open(newline='',encoding='utf-8-sig')))
    if len(pre)!=len(seq_order): raise ValueError('PRE/POST count differs from TALLYSEQ count')
    pre_by_seq=defaultdict(list); meta_mismatches=0; pre_increment_vs_tlv=[]
    for serial,row in zip(seq_order,pre):
        i=len(pre_by_seq[serial]); pre_by_seq[serial].append(row); q=seq[serial][i]
        for key,value in [('event_id',int(float(row['nocas'])),),('particle_id',int(row['particle_id'])),
                          ('ityp',int(row['ityp'])),('jcoll',int(row['jcoll'])),('ncol',int(row['ncol']))]:
            meta_mismatches += q[key]!=value
        pre_increment_vs_tlv.append(abs(float(row['tally_increment_MeV'])-q['tlv_MeV']))
    fields=['addition_id','source_line','dump_sequence','event_id','particle_id','process',
      'region_hdf5','region_internal','energy_hdf5_MeV','energy_internal_tlv_MeV',
      'energy_internal_increment_MeV','abs_diff_to_tlv_MeV','abs_diff_to_increment_MeV',
      'weight_internal','position_kind','position_max_abs_diff_cm','position_valid']
    a.output_csv.parent.mkdir(parents=True,exist_ok=True)
    counts=Counter(); occurrences=Counter(); energy_diffs=[]; increment_diffs=[]; position_diffs=[]
    event_h5=defaultdict(float); event_internal=defaultdict(float); event_counts_h5=Counter(); event_counts_internal=Counter()
    with a.output_csv.open('w',newline='',encoding='utf-8') as out, h5py.File(a.hdf5,'r') as h:
        w=csv.DictWriter(out,fieldnames=fields); w.writeheader(); d=h['deposits']; n=len(d['event_id'])
        chunk_fields=('addition_id','source_record_id','event_id','particle_id','process','region_id','energy_MeV','energy_valid','kind','start_position_valid','end_position_valid','score_position_valid','start_xyz_cm','end_xyz_cm','score_xyz_cm')
        for base in range(0,n,65536):
            arrays={name:d[name][base:min(base+65536,n)] for name in chunk_fields}
            for j in range(len(arrays['event_id'])):
                i=base+j
                sid=int(arrays['source_record_id'][j].decode() if isinstance(arrays['source_record_id'][j],bytes) else arrays['source_record_id'][j])
                serial=line_to_seq.get(sid)
                if serial is None: counts['unmapped_hdf5_source_line']+=1; continue
                occurrence=occurrences[serial]; occurrences[serial]+=1
                if occurrence>=len(seq[serial]) or occurrence>=len(pre_by_seq[serial]): counts['excess_hdf5_occurrence']+=1; continue
                q=seq[serial][occurrence]; r=pre_by_seq[serial][occurrence]
                event=int(arrays['event_id'][j]); pid=int(arrays['particle_id'][j]); region=int(arrays['region_id'][j])
                hv=float(arrays['energy_MeV'][j]) if bool(arrays['energy_valid'][j]) else None
                if hv is None: counts['invalid_hdf5_energy']+=1; continue
                diff=abs(hv-q['tlv_MeV']); dinc=abs(hv-float(r['tally_increment_MeV']))
                energy_diffs.append(diff); increment_diffs.append(dinc); event_h5[event]+=hv; event_internal[event]+=q['tlv_MeV']; event_counts_h5[event]+=1; event_counts_internal[event]+=1
                counts['region_mismatch']+=region!=q['tally_region_id']; counts['event_mismatch']+=event!=q['event_id']; counts['particle_mismatch']+=pid!=q['particle_id']
                kind=int(arrays['kind'][j]); posdiff=None; posvalid=False
                p0=vec(r,['x_start_cm','y_start_cm','z_start_cm']); p1=vec(r,['x_end_cm','y_end_cm','z_end_cm'])
                if kind==2 and bool(arrays['start_position_valid'][j]) and bool(arrays['end_position_valid'][j]):
                    h0=tuple(float(x) for x in arrays['start_xyz_cm'][j]); h1=tuple(float(x) for x in arrays['end_xyz_cm'][j])
                    posdiff=max(max_vec_delta(h0,p0),max_vec_delta(h1,p1)); posvalid=True
                elif bool(arrays['score_position_valid'][j]):
                    hp=tuple(float(x) for x in arrays['score_xyz_cm'][j]); posdiff=min(max_vec_delta(hp,p0),max_vec_delta(hp,p1)); posvalid=True
                if posvalid: position_diffs.append(posdiff)
                else: counts['position_unavailable']+=1
                w.writerow({'addition_id':arrays['addition_id'][j].decode() if isinstance(arrays['addition_id'][j],bytes) else arrays['addition_id'][j],
                  'source_line':sid,'dump_sequence':serial,'event_id':event,'particle_id':pid,
                  'process':arrays['process'][j].decode() if isinstance(arrays['process'][j],bytes) else arrays['process'][j],
                  'region_hdf5':region,'region_internal':q['tally_region_id'],'energy_hdf5_MeV':format(hv,'.17g'),
                  'energy_internal_tlv_MeV':format(q['tlv_MeV'],'.17g'),'energy_internal_increment_MeV':r['tally_increment_MeV'],
                  'abs_diff_to_tlv_MeV':format(diff,'.17g'),'abs_diff_to_increment_MeV':format(dinc,'.17g'),
                  'weight_internal':format(q['weight'],'.17g'),'position_kind':kind,
                  'position_max_abs_diff_cm':'' if posdiff is None else format(posdiff,'.17g'),'position_valid':posvalid})
    events=sorted(set(event_h5)|set(event_internal)); event_diff=[abs(event_h5[e]-event_internal[e]) for e in events]
    summary={'hdf5_rows':n,'matched_rows':len(energy_diffs),'raw_transport_rows':len(raw_tracks),
      'raw_rows_mapped_to_unique_dump_sequences':len(line_to_seq),'raw_mapping_ambiguous':raw_ambiguous,
      'raw_mapping_unmapped':raw_unmapped,'internal_tallyseq_rows':len(seq_order),
      'internal_unique_dump_sequences':len(seq),'internal_duplicate_sequence_rows':sum(len(v)-1 for v in seq.values()),
      'internal_PRE_POST_rows':len(pre),'pre_identity_mismatches':meta_mismatches,
      'pre_increment_vs_tlv_max_abs_MeV':max(pre_increment_vs_tlv,default=0),
      'hdf5_vs_internal_tlv_max_abs_MeV':max(energy_diffs,default=0),
      'hdf5_vs_internal_tlv_sum_abs_MeV':sum(energy_diffs),
      'hdf5_vs_internal_increment_max_abs_MeV':max(increment_diffs,default=0),
      'hdf5_vs_internal_increment_rows_over_1e-12_MeV':sum(x>1e-12 for x in increment_diffs),
      'hdf5_vs_internal_total_MeV':sum(event_h5.values()),'internal_tally_total_MeV':sum(event_internal.values()),
      'event_count':len(events),'event_count_mismatches':sum(event_counts_h5[e]!=event_counts_internal[e] for e in events),
      'event_sum_max_abs_diff_MeV':max(event_diff,default=0),'event_sum_abs_diff_sum_MeV':sum(event_diff),
      'event_total_abs_diff_MeV':abs(sum(event_h5.values())-sum(event_internal.values())),
      'position_rows_compared':len(position_diffs),'position_unavailable_rows':counts['position_unavailable'],
      'position_max_abs_diff_cm':max(position_diffs,default=0),
      'region_mismatches':counts['region_mismatch'],'event_mismatches':counts['event_mismatch'],
      'particle_mismatches':counts['particle_mismatch'],'unmapped_hdf5_source_line':counts['unmapped_hdf5_source_line'],
      'excess_hdf5_occurrence':counts['excess_hdf5_occurrence'],
      'dumpall_sha256':sha(run/'dumpall.dat'),'phits_input_sha256':sha(run/'trace.inp'),
      'internal_trace_sha256':sha(run/'tally_internal_trace.log'),'tally_dumpseq_sha256':sha(run/'tally_dumpseq.log'),
      'hdf5_sha256':sha(a.hdf5),'comparison_csv_sha256':sha(a.output_csv)}
    a.output_json.write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(summary,indent=2))
if __name__=='__main__': main()
