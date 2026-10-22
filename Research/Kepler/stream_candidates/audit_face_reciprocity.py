"""Audit both owned ends of every recorded internal/periodic face."""
from pathlib import Path
import json, sys
import numpy as np

root = Path(__file__).resolve().parent
out = {}
labels=sys.argv[1:]
output=root/'face_reciprocity_audit.json'
if labels and labels[0]=='--output':
    output=root/labels[1]
    labels=labels[2:]
for label in labels:
    d = root / label
    f = np.vstack([np.loadtxt(p,skiprows=1,ndmin=2) for p in sorted(d.glob('face_dump_r*_call1.tsv'))])
    if np.any(f[:,6] != 0):
        raise ValueError('Wall faces need a separate audit')
    a,b = f[:,4].astype(np.int64),f[:,5].astype(np.int64)
    keys = np.minimum(a,b)*65536 + np.maximum(a,b)
    order = np.argsort(keys,kind='stable')
    k = keys[order]
    starts = np.r_[0, np.flatnonzero(np.diff(k))+1]
    counts = np.diff(np.r_[starts,len(k)])
    pair = starts[counts==2]
    p,q = order[pair],order[pair+1]
    if not np.all(a[p]==b[q]) or not np.all(b[p]==a[q]):
        raise ValueError('Two records are not opposite ends')
    S1,S2 = f[p,10:12],f[q,10:12]
    area_scale = np.maximum(np.linalg.norm(S1,axis=1),np.linalg.norm(S2,axis=1))
    geom_abs = np.linalg.norm(S1+S2,axis=1)
    geom_rel = geom_abs/np.maximum(area_scale,1e-300)
    cm1=.5*(f[p,12:14]+f[p,14:16]); cm2=.5*(f[q,12:14]+f[q,14:16])
    dc=cm1-cm2;dc-=24*np.rint(dc/24)
    pressure_diff=abs(f[p,7]-f[q,7])
    pair_work=f[p,9]+f[q,9]
    pair_force=-(f[p,7,None]*S1+f[q,7,None]*S2)
    odd_groups=[np.arange(st,st+n) for st,n in zip(starts,counts) if n!=2]
    orphan=order[np.concatenate(odd_groups)] if odd_groups else np.empty(0,dtype=int)
    r=np.hypot(cm1[:,0]-12,cm1[:,1]-12)
    radial=[]
    for lo,hi in [(0,.5),(.5,1),(1,2),(2,4),(4,8),(8,30)]:
        sel=(r>=lo)&(r<hi)
        if not sel.any(): continue
        radial.append({'Rlo':lo,'Rhi':hi,'pairs':int(sel.sum()),
                       'geometry_rel_error_max':float(geom_rel[sel].max()),
                       'pairs_geometry_rel_error_gt_1e-8':int(np.sum(geom_rel[sel]>1e-8)),
                       'pressure_difference_max':float(pressure_diff[sel].max()),
                       'sum_abs_pair_work_residual':float(abs(pair_work[sel]).sum()),
                       'sum_abs_pair_force_residual':float(np.linalg.norm(pair_force[sel],axis=1).sum())})
    out[label]={'face_ends':len(f),'unique_pairs':len(starts),'counts_histogram':dict(zip(*[x.tolist() for x in np.unique(counts,return_counts=True)])),
                'orphan_face_ends':len(orphan),'orphan_ids':f[orphan,4:6].astype(int).tolist(),
                'geometry_abs_error_max':float(geom_abs.max()),'geometry_rel_error_max':float(geom_rel.max()),
                'centroid_periodic_distance_max':float(np.linalg.norm(dc,axis=1).max()),
                'pressure_difference_max':float(pressure_diff.max()),
                'sum_work_all_ends':float(f[:,9].sum()),
                'net_force_all_ends':(-np.sum(f[:,7,None]*f[:,10:12],axis=0)).tolist(),
                'radial':radial}
    print(label, 'face_ends',len(f),'unpaired',len(orphan),'geom_bad_pairs',int(np.sum(geom_rel>1e-8)), 'max_geom',geom_rel.max())
output.write_text(json.dumps(out,indent=2)+'\n')
