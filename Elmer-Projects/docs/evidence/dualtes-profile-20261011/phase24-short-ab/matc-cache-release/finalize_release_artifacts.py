from pathlib import Path
import json,hashlib,shutil
HERE=Path(__file__).resolve().parent; R=Path(r'C:\\Users\\Nao\\AppData\\Local\\Temp\\phase24-cache-runtime'); RESULT_DIR=R/'mesh_toml_95333801e1a3'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''): h.update(b)
 return h.hexdigest()
p=HERE/'release_runs.json'; d=json.loads(p.read_text())
for rec in d['runs']:
 label=rec['label']; sif=R/f'dualtes_release_{label}.sif'
 if sif.exists():
  rec['sif_sha256']=sha(sif); shutil.copy2(sif,HERE/f'dualtes_{label}.sif')
 outputs=[]
 for tes in (1,2):
  for kind in ('series','iterations'):
   src=(R/f'dualtes_release_{label}_{tes}_{kind}.csv')
   if not src.exists() and rec['mode']=='on':
    src=R/f'dualtes_cache_verify_{tes}_{kind}.csv'; rec['csv_telemetry_trial']='on_3'
   if src.exists():
    name=f'{label}_{tes}_{kind}.csv'; dst=HERE/name; shutil.copy2(src,dst); outputs.append({'file':name,'bytes':dst.stat().st_size,'sha256':sha(dst)})
 base=RESULT_DIR/f'dualtes_release_{label}.result'
 for rank in range(6):
  src=Path(str(base)+f'.{rank}')
  if src.exists(): outputs.append({'file':src.name,'rank':rank,'bytes':src.stat().st_size,'sha256':sha(src)})
 rec['outputs']=outputs
# include the corrected runner source digest and distinguish validation inputs
for label in d['order']:
 psrc=HERE/f'{label}.sif'
 if psrc.exists(): continue
 src=R/f'dualtes_release_{label}.sif'
 if src.exists(): shutil.copy2(src,psrc)
d['sif_source_copy_notes']='Each trial SIF is copied into this evidence directory. For ON telemetry CSV names from the legacy source SIF were shared; saved ON CSVs are from on_3; the corrected runner uses a unique output prefix per trial.'
d['result_location']='isolated Temp mesh alias; result shard hashes retained in runs.json and Release field comparison JSON'
p.write_text(json.dumps(d,indent=2)+'\n',encoding='utf-8')
