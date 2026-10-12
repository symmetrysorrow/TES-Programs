import hashlib,json
from pathlib import Path
here=Path(__file__).resolve().parent
meshroot=Path(r'C:\Users\Nao\AppData\Local\Temp\phase24-cache-runtime')
cases={'dualtes':'mesh_toml_95333801e1a3','singlepixel':'mesh_hybrid_fullconf_h8'}
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def parse(p):
 lines=p.read_text(encoding='ascii',errors='replace').splitlines();out={};i=0
 while i<len(lines):
  if not lines[i].startswith('Time:'): i+=1;continue
  step=int(lines[i].split()[1]);i+=1;var=lines[i].strip();i+=1
  if not lines[i].startswith('Perm:'): raise ValueError(f'{p}: missing Perm at {i}')
  n=int(lines[i].split()[2]);i+=1;ids=[int(lines[i+j].split()[1]) for j in range(n)];i+=n;vals=[float(lines[i+j]) for j in range(n)];i+=n
  out[step]=(var,dict(zip(ids,vals)))
 return out
comparisons={}
for case,meshname in cases.items():
 mesh=meshroot/meshname;per=[];shards=[]
 for rank in range(6):
  a=mesh/f'finalreplay_{case}_off.result.{rank}';b=mesh/f'finalreplay_{case}_on.result.{rank}'
  aa,bb=parse(a),parse(b);shards.extend([{'file':p.name,'bytes':p.stat().st_size,'sha256':sha(p)} for p in (a,b)])
  if aa.keys()!=bb.keys(): raise SystemExit(f'timestep mismatch {case}/{rank}')
  for step in aa:
   av,bv=aa[step][1],bb[step][1]
   if av.keys()!=bv.keys(): raise SystemExit(f'node mismatch {case}/{rank}/{step}')
   diffs=[abs(av[n]-bv[n]) for n in av]
   per.append({'rank':rank,'timestep':step,'nodes':len(av),'max_abs_K':max(diffs,default=0),'different_values':sum(d!=0 for d in diffs)})
 comparisons[case]={'ranks':6,'steps':4,'all_nodes':True,'per_rank_step':per,'max_abs_K':max(x['max_abs_K'] for x in per),'different_values':sum(x['different_values'] for x in per),'result_shard_hashes':shards}
out={'release_solver_sha256':'946ce238495957de1d95f6ca0e5f428ccaf966199df14cf7d83982ac22f40cae','release_heat_solve_sha256':'324c3168a16756f22ce1e80e85d0ac6939e2a8b751f651367f37636416ab8a12','comparisons':comparisons}
(here/'final_release_temperature_comparison.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print(json.dumps({k:{'max_abs_K':v['max_abs_K'],'different_values':v['different_values']} for k,v in comparisons.items()},indent=2))
