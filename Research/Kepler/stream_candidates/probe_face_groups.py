"""Inspect repeated ID pairs without assuming one edge record per side."""
from pathlib import Path
import json,sys
import numpy as np
case=Path(sys.argv[1]);f=np.loadtxt(case/'face_dump_r0_call1.tsv',skiprows=1)
a=f[:,4].astype(np.int64);b=f[:,5].astype(np.int64)
lo=np.minimum(a,b);hi=np.maximum(a,b);order=np.lexsort((hi,lo))
starts=np.r_[0,np.flatnonzero((np.diff(lo[order])!=0)|(np.diff(hi[order])!=0))+1]
counts=np.diff(np.r_[starts,len(f)])
S=np.add.reduceat(f[order,10:12],starts,axis=0)
force=np.add.reduceat(-f[order,7,None]*f[order,10:12],starts,axis=0)
work=np.add.reduceat(f[order,9],starts)
length=np.add.reduceat(np.linalg.norm(f[order,10:12],axis=1),starts)
worst=[]
for group in np.argsort(np.linalg.norm(force,axis=1))[-20:][::-1]:
    ix=order[starts[group]:starts[group]+counts[group]]
    worst.append({'ids':[int(lo[ix[0]]),int(hi[ix[0]])],'records':f[ix].tolist(),
                  'sum_S':S[group].tolist(),'sum_force':force[group].tolist(),'sum_work':float(work[group])})
out={'case':str(case),'geometry_group_residual_max':float(np.linalg.norm(S,axis=1).max()),
     'force_group_residual_abs_sum':float(np.linalg.norm(force,axis=1).sum()),
     'work_group_residual_abs_sum':float(abs(work).sum()),'groups_by_count':{},'worst_force_groups':worst}
for n in np.unique(counts):
    sel=counts==n
    out['groups_by_count'][str(n)]={'groups':int(sel.sum()),'length_sum':float(length[sel].sum()),
        'max_geometry_residual':float(np.linalg.norm(S[sel],axis=1).max()),
        'force_residual_abs_sum':float(np.linalg.norm(force[sel],axis=1).sum()),
        'work_residual_abs_sum':float(abs(work[sel]).sum())}
(case/'group_probe.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k!='worst_force_groups'},indent=2))
for z in worst[:3]:
    print('WORST',z['ids'],'force',z['sum_force'])
    for face in z['records']:print('END',int(face[4]),int(face[5]),'P',face[7],'S',face[10:12],'vertices',face[12:16])
