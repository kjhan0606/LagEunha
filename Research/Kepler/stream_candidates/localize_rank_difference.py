from pathlib import Path
import json,numpy as np
from audit_snapshot import read_snapshot
root=Path(__file__).resolve().parent
p,m=read_snapshot(root/'fixed_late_a/khout.009570.dat'); ids=(p['u4if']&np.uint64((1<<56)-1)).astype(int); order=np.argsort(ids);p=p[order];r=np.hypot(p['x']-12,p['y']-12)
fa=[]
for d in ('fixed_late_a','fixed_r4_late_a'):
 f=np.vstack([np.loadtxt(f,skiprows=1) for f in sorted((root/d).glob('face_dump_r*_call1.tsv'))]);f=f[np.lexsort((f[:,5],f[:,4]))];fa.append(f)
a,b=fa
assert np.array_equal(a[:,4:6],b[:,4:6])
fd=abs(a[:,7]-b[:,7]);fr=r[a[:,4].astype(int)]
st=[]
for d in ('fixed_late_a','fixed_r4_late_a'):
 q,_=read_snapshot(root/d/'khout.009571.dat'); qi=(q['u4if']&np.uint64((1<<56)-1)).astype(int);st.append(q[np.argsort(qi)])
dv=np.hypot(st[0]['vx']-st[1]['vx'],st[0]['vy']-st[1]['vy'])
out=[]
for lo,hi in [(0,.5),(.5,1),(1,2),(2,4),(4,8),(8,10),(10,12),(12,100)]:
 s=(r>=lo)&(r<hi);f=(fr>=lo)&(fr<hi)
 out.append({'Rlo':lo,'Rhi':hi,'max_dv':float(dv[s].max()),'cells_dv_gt_1e-12':int((dv[s]>1e-12).sum()),'max_dpstar':float(fd[f].max()),'face_ends_dpstar_gt_1e-12':int((fd[f]>1e-12).sum())})
Path('rank_difference_localization.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
