"""Compare physical-body volumes and steady temperatures across Elmer meshes."""
from __future__ import annotations
import argparse, json, math, sys, hashlib
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.support.mesh_names import parse_mesh_names

def nodes(mesh: Path):
    return {int(a[0]): tuple(map(float,a[2:5])) for a in (line.split() for line in (mesh/'mesh.nodes').read_text().splitlines())}
def sha(path):
 h=hashlib.sha256();
 with path.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def volume(kind, p):
    def sub(a,b): return tuple(x-y for x,y in zip(a,b))
    def dot(a,b): return sum(x*y for x,y in zip(a,b))
    def cross(a,b): return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
    def tet(i,j,k,l): return abs(dot(sub(p[j],p[i]),cross(sub(p[k],p[i]),sub(p[l],p[i]))))/6
    if kind==504: return tet(0,1,2,3)
    if kind==706: return tet(0,1,2,3)+tet(1,2,3,4)+tet(2,3,4,5)
    return 0.0
def elements(mesh):
    xyz=nodes(mesh); out=[]
    for f in (mesh/'mesh.elements').read_text().splitlines():
        a=f.split(); kind=int(a[2]); ids=list(map(int,a[3:]));
        if kind in (504,706): out.append((int(a[1]),kind,ids,volume(kind,[xyz[x] for x in ids])))
    return out
def temps(result):
    l=result.read_text().splitlines(); i=next(i for i,x in enumerate(l) if x.startswith('Perm:')); n=int(l[i].split()[1]); perm=[int(x.split()[0]) for x in l[i+1:i+1+n]]
    vals=[float(x.replace('D','E')) for x in l[i+1+n:i+1+2*n]]
    # Elmer ASCII values are mesh-node order; permutation is metadata.
    return {node: vals[node-1] for node in perm}
def summary(mesh,result):
    names=parse_mesh_names(mesh/'mesh.names').bodies; rev={v:k for k,v in names.items()}; t=temps(result) if result else None; data={}
    for body,kind,ids,v in elements(mesh):
        d=data.setdefault(rev.get(body,f'body_{body}'),{'volume_m3':0.,'elements':{},'temperature_integral':0.,'min':math.inf,'max':-math.inf})
        d['volume_m3']+=v; d['elements'][str(kind)]=d['elements'].get(str(kind),0)+1
        if t:
            q=sum(t[x] for x in ids)/len(ids); d['temperature_integral']+=q*v; d['min']=min(d['min'],*(t[x] for x in ids)); d['max']=max(d['max'],*(t[x] for x in ids))
    for d in data.values():
        if t: d['temperature_mean_K']=d.pop('temperature_integral')/d['volume_m3']
        else: d.pop('temperature_integral'); d.pop('min'); d.pop('max')
    return data
def main():
 p=argparse.ArgumentParser();p.add_argument('tet_mesh',type=Path);p.add_argument('hybrid_mesh',type=Path);p.add_argument('--tet-result',type=Path);p.add_argument('--hybrid-result',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 tet=summary(a.tet_mesh,a.tet_result); hy=summary(a.hybrid_mesh,a.hybrid_result)
 out={'inputs':{'tet_mesh':str(a.tet_mesh.resolve()),'hybrid_mesh':str(a.hybrid_mesh.resolve()),'tet_result':str(a.tet_result.resolve()),'hybrid_result':str(a.hybrid_result.resolve()),'tet_result_sha256':sha(a.tet_result),'hybrid_result_sha256':sha(a.hybrid_result)},'all_tet':tet,'hybrid':hy,'body_volume_relative_difference':{k:(hy[k]['volume_m3']-v['volume_m3'])/v['volume_m3'] for k,v in tet.items() if k in hy},'interfaces':{}}
 for x,y in [('TES','Stycast'),('TES','Membrane_SiNx'),('Stycast','abs')]:
  if x in hy and y in hy: out['interfaces'][f'{x}-{y}']={'hybrid_mean_delta_K':hy[x].get('temperature_mean_K',0)-hy[y].get('temperature_mean_K',0),'tet_mean_delta_K':tet.get(x,{}).get('temperature_mean_K',0)-tet.get(y,{}).get('temperature_mean_K',0)}
 a.output.write_text(json.dumps(out,indent=2)+'\n');print(a.output)
if __name__=='__main__': main()
