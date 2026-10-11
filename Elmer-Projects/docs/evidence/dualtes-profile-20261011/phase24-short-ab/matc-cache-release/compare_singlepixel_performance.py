from pathlib import Path
import csv,hashlib,json,re,sys
HERE=Path(__file__).resolve().parent; MESH=Path(r'C:\\Users\\Nao\\AppData\\Local\\Temp\\phase24-cache-runtime\\mesh_hybrid_fullconf_h8')
def field(p):
 l=p.read_text(encoding='ascii',errors='replace').splitlines();o={};i=0
 while i<len(l):
  if not l[i].startswith('Time:'):i+=1;continue
  t=int(l[i].split()[1]);i+=1;variable=l[i].strip();i+=1;n=int(l[i].split()[2]);i+=1;ids=[int(l[i+j].split()[1]) for j in range(n)];i+=n;values=[float(l[i+j].strip()) for j in range(n)];i+=n;o[t]=(variable,dict(zip(ids,values)))
 return o
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
rows=[];hashes=[];maxabs=maxrel=0.;neq=0
for rank in range(6):
 a=MESH/f'singlepixel_release_off.result.{rank}';b=MESH/f'singlepixel_release_on_perf.result.{rank}';x=field(a);y=field(b);hashes.extend([{'file':a.name,'sha256':sha(a),'bytes':a.stat().st_size},{'file':b.name,'sha256':sha(b),'bytes':b.stat().st_size}])
 for step in sorted(x):
  av=x[step][1];bv=y[step][1];d=[(abs(av[k]-bv[k]),av[k],bv[k]) for k in av];ma=max(z[0] for z in d);mr=max(z[0]/max(abs(z[1]),abs(z[2]),1e-300) for z in d);n=sum(z[0]!=0 for z in d);maxabs=max(maxabs,ma);maxrel=max(maxrel,mr);neq+=n;rows.append({'rank':rank,'timestep':step,'max_abs_K':ma,'max_rel':mr,'different_values':n})
checks={}
for kind in ('series','iterations'):
 a=HERE/f'singlepixel_{kind}_off.csv';b=HERE/f'singlepixel_{kind}_on_perf.csv';checks[kind+'_exact']=list(csv.reader(a.open()))==list(csv.reader(b.open()))
logs={}
for mode,name in [('off','singlepixel_off.solver.log'),('on_perf','singlepixel_on_perf.solver.log')]:
 s=(HERE/name).read_text(encoding='utf-8',errors='replace');logs[mode]={'result_norms':[float(z) for z in re.findall(r'HeatSolve:\s+Result Norm\s*:\s*([0-9.Ee+-]+)',s)],'hypre':[[int(a),float(b)] for a,b in re.findall(r'Phase24 PCG iterations=(\d+) final_rel_res=([0-9.Ee+-]+)',s)]}
checks['result_norms_exact']=logs['off']['result_norms']==logs['on_perf']['result_norms'];checks['hypre_pcg_exact']=logs['off']['hypre']==logs['on_perf']['hypre']
r={'ranks':6,'timesteps':4,'temperature_max_abs_K':maxabs,'temperature_max_rel':maxrel,'different_values':neq,'csv_and_solver_checks':checks,'solver_record_counts':{k:{'result_norms':len(v['result_norms']),'hypre_pcg':len(v['hypre'])} for k,v in logs.items()},'per_rank_timestep':rows,'result_shard_hashes':hashes}
(HERE/'singlepixel_performance_numerical_comparison.json').write_text(json.dumps(r,indent=2)+'\n',encoding='utf-8');print(json.dumps({k:v for k,v in r.items() if k not in ('per_rank_timestep','result_shard_hashes')},indent=2))
