from pathlib import Path
import hashlib,json,sys
runtime=Path(sys.argv[1])
result=Path(sys.argv[2])

def parse(p):
 lines=p.read_text(encoding='ascii',errors='replace').splitlines(); out={}; i=0
 while i<len(lines):
  if not lines[i].startswith('Time:'): i+=1; continue
  timestep=int(lines[i].split()[1]); i+=1; variable=lines[i].strip(); i+=1
  if not lines[i].startswith('Perm:'): raise ValueError(f'{p}:{i}: expected Perm')
  count=int(lines[i].split()[2]); i+=1; ids=[]
  for _ in range(count): ids.append(int(lines[i].split()[1])); i+=1
  values=[]
  for _ in range(count): values.append(float(lines[i].strip())); i+=1
  out[timestep]=(variable,dict(zip(ids,values)))
 return out

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for block in iter(lambda:f.read(1<<20),b''): h.update(block)
 return h.hexdigest()
records=[]; hashes=[]
for rank in range(6):
 a_path=runtime/f'dualtes_release_off_1.result.{rank}'; b_path=runtime/f'dualtes_release_on_3.result.{rank}'
 a=parse(a_path); b=parse(b_path)
 hashes.extend([{'file':a_path.name,'sha256':sha(a_path),'bytes':a_path.stat().st_size},{'file':b_path.name,'sha256':sha(b_path),'bytes':b_path.stat().st_size}])
 for step in sorted(a):
  av=a[step][1]; bv=b[step][1]
  if av.keys()!=bv.keys(): raise ValueError('node index mismatch')
  diffs=[(abs(av[k]-bv[k]),k,av[k],bv[k]) for k in av]
  records.append({'rank':rank,'timestep':step,'nodes':len(av),'max_abs':max(x[0] for x in diffs),'max_rel':max(x[0]/max(abs(x[2]),abs(x[3]),1e-300) for x in diffs),'different_values':sum(x[0]!=0 for x in diffs)})
summary={'off':'dualtes_release_off_1','on':'dualtes_release_on_3','fields':'Temperature','ranks':6,'timesteps':4,'per_rank_timestep':records,'result_shard_hashes':hashes}
result.write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
