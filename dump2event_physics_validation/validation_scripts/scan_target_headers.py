import time
from pathlib import Path
p=Path(r'H:\hata\1332_142_136\1332keV_6\dumpall.dat')
target_nocas={318209,119953}; target_rcasc={318209}
found={}; headers=[]; active=None; block=[]; n=0; started=time.time()
def parse(b):
    try:return [float(x.replace(b'D',b'E').replace(b'd',b'e')) for x in b.split()]
    except:return []
def save(k, rows):
    found[k]=rows
    out=Path('dump2event_physics_validation')/f'raw_history_NOCAS_{k}.txt'
    with out.open('w',encoding='ascii') as w:
        w.write(f'Source: {p}\nTarget NOCAS={k}\n')
        for ln,line in rows:w.write(f'{ln}: {line}\n')
    print(f'SAVED NOCAS={k} lines={len(rows)} file_lines={rows[0][0]}-{rows[-1][0]} elapsed={time.time()-started:.1f}s',flush=True)
with p.open('rb',buffering=8*1024*1024) as f:
    while True:
        raw=f.readline()
        if not raw: break
        n+=1; s=raw.strip()
        if s==b'4':
            nxt=f.readline()
            if not nxt: break
            n+=1; v=parse(nxt)
            is_header=(len(v)>=4 and int(v[0])>0 and int(v[1])>0 and v[1].is_integer())
            if is_header:
                nocas,nobch,rcasc,rsouin=v[:4]
                if int(nocas) in target_nocas or int(rcasc) in target_rcasc:
                    headers.append((n-1,int(nocas),int(nobch),rcasc,rsouin))
                    if len(headers)%1000==0: pass
                if active is not None:
                    block.extend([(n-1,'4'),(n,nxt.decode('ascii').rstrip())])
                    save(active,block); active=None; block=[]
                    if len(found)==len(target_nocas): break
                if int(nocas) in target_nocas and int(nobch)==1 and int(nocas) not in found:
                    active=int(nocas); block=[(n-1,'4'),(n,nxt.decode('ascii').rstrip())]
                continue
        if active is not None:block.append((n,raw.decode('ascii').rstrip()))
print(f'SCAN_DONE lines={n} elapsed={time.time()-started:.1f}s found={sorted(found)}',flush=True)
with Path('dump2event_physics_validation/header_identity_matches.csv').open('w',encoding='utf-8',newline='') as w:
    w.write('source_line,NOCAS,NOBCH,RCASC,RSOUIN\n')
    for row in headers:w.write(','.join(map(str,row))+'\n')
print('header matches',len(headers),headers[:20],flush=True)
