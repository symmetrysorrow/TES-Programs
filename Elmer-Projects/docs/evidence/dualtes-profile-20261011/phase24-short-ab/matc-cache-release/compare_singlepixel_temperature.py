from pathlib import Path
import hashlib,json
root=Path(r'C:\Users\Nao\AppData\Local\Temp\phase24-cache-runtime\mesh_hybrid_fullconf_h8')
out=Path(__file__).resolve().parent

def parse(p):
 lines=p.read_text(encoding='ascii',errors='replace').splitlines(); result={}; i=0
 while i<len(lines):
  if not lines[i].startswith('Time:'): i+=1; continue
  step=int(lines[i].split()[1]); i+=1; variable=lines[i].strip(); i+=1
  if not lines[i].startswith('Perm:'): raise ValueError(f'{p}: expected Perm')
  count=int(lines[i].split()[2]); i+=1; ids=[]
  for _ in range(count): ids.append(int(lines[i].split()[1])); i+=1
  vals=[]
  for _ in range(count): vals.append(float(lines[i].strip())); i+=1
  result[step]=(variable,dict(zip(ids,vals)))
 return result

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''): h.update(b)
 return h.hexdigest()
records=[]; hashes=[]
for rank in range(6):
 a=root/f'singlepixel_release_off.result.{rank}'; b=root/f'singlepixel_release_on.result.{rank}'; aa=parse(a); bb=parse(b)
 hashes.extend([{'file':a.name,'bytes':a.stat().st_size,'sha256':sha(a)},{'file':b.name,'bytes':b.stat().st_size,'sha256':sha(b)}])
 if aa.keys()!=bb.keys(): raise ValueError('step mismatch')
 for step in sorted(aa):
  av=aa[step][1]; bv=bb[step][1]
  if av.keys()!=bv.keys(): raise ValueError('node mapping mismatch')
  dif=[(abs(av[k]-bv[k]),av[k],bv[k]) for k in av]
  records.append({'rank':rank,'timestep':step,'node_count':len(av),'max_abs_K':max(x[0] for x in dif),'max_rel':max(x[0]/max(abs(x[1]),abs(x[2]),1e-300) for x in dif),'different_values':sum(x[0]!=0 for x in dif)})
summary={'ranks':6,'timesteps':4,'variable':'Temperature','max_abs_K':max(x['max_abs_K'] for x in records),'max_rel':max(x['max_rel'] for x in records),'total_different_values':sum(x['different_values'] for x in records),'per_rank_timestep':records,'result_shard_hashes':hashes}
(out/'singlepixel_temperature_comparison.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
print(json.dumps({k:v for k,v in summary.items() if k not in ('per_rank_timestep','result_shard_hashes')},indent=2))
