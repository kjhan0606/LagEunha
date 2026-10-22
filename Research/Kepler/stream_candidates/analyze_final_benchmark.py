from pathlib import Path
import json,re
import numpy as np
from audit_snapshot import read_snapshot

root=Path(__file__).resolve().parent
out={}
states=[]
for label in ('bench_final_cpu','bench_final_gpu'):
    d=root/label
    timing=json.loads((d/'timing.json').read_text())
    assert timing['exit_code']==0
    snap=next(p for p in d.glob('khout.*.dat') if not p.is_symlink())
    p,m=read_snapshot(snap)
    order=np.argsort(p['u4if']&np.uint64((1<<56)-1))
    states.append((p[order],m))
    rows=re.findall(r'\[CPU step (\d+)\] update=([0-9.]+)s force=([0-9.]+)s post=([0-9.]+)s final=([0-9.]+)s \| total=([0-9.]+)s',(d/'log').read_text())
    out[label]={'timing':timing,'time':m['t'],'native_steps':len(rows),
                'phase_median_seconds':dict(zip(('update','force','post','final','total'),np.median(np.array(rows[1:],dtype=float)[:,1:],axis=0).tolist()))}
p,m=states[0];q,n=states[1]
assert np.array_equal(p['u4if']&np.uint64((1<<56)-1),q['u4if']&np.uint64((1<<56)-1))
r=np.hypot(p['x']-12,p['y']-12)
out['wall_time_ratio_cpu_over_gpu']=out['bench_final_cpu']['timing']['elapsed_seconds']/out['bench_final_gpu']['timing']['elapsed_seconds']
out['caution']='15 sequential steps, shared H200 with old live baseline; not dedicated GPU scaling. CPU final areas have proven wrong cells, so differences from CPU are not all GPU errors.'
out['radial_differences']=[]
for lo,hi in ((0,.5),(.5,1),(1,2),(2,4),(4,8),(8,30)):
    sel=(r>=lo)&(r<hi)
    row={'Rlo':lo,'Rhi':hi,'n':int(sel.sum()),'fields':{}}
    for f in ('x','y','vx','vy','mass','volume','den','pressure','ie'):
        a,b=p[f][sel].astype(float),q[f][sel].astype(float)
        row['fields'][f]={'max_abs':float(abs(a-b).max()),'changed_cells':int(np.sum(a!=b)),
                          'relative_l2':float(np.linalg.norm(a-b)/max(np.linalg.norm(a),1e-300))}
    out['radial_differences'].append(row)
(root/'final_volume_gpu/benchmark.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k!='radial_differences'},indent=2))
