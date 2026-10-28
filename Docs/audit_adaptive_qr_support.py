"""Independent support-count candidate audit; no production solver changes."""
import json,sys
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read,sha
b=Path('/gpfs/kjhan/LagForce')
out={'target_kernel_sum':12.,'cases':{},'scope':'Independent full-cloud geometry audit. Not the unavailable failed smooth state, not production search coverage, not evolution or a rank guarantee.'}
for label,path in [('legacy_t3p4',b/'Codex_coupled_onset_to_t4_20261018/reference/khout.000743.dat'),('legacy_t4',b/'Codex_coupled_onset_to_t4_20261018/reference/khout.000984.dat')]:
 a,t=read(path);x=np.column_stack((a['x'],a['y']));tree=cKDTree(x)
 owners=np.flatnonzero(np.linalg.norm(x-3,axis=1)<1.9)
 distances,neighbors=tree.query(x[owners],k=129)
 ranks=[];conditions=[];radii=[];counts=[];residuals=[];oldcounts=[];worst=None
 for row,i in enumerate(owners):
  keep=neighbors[row]!=i;d=distances[row][keep];ids=neighbors[row][keep]
  assert len(d)==128 and np.all(d>0)
  # Root must lie inside the last queried distance so missing points have zero weight.
  lo=0.;hi=float(d[-1]);d2=d*d
  assert np.maximum(1-d2/(hi*hi),0).dot(np.maximum(1-d2/(hi*hi),0))>12
  for iteration in range(55):
   mid=(lo+hi)/2;z=np.maximum(1-d2/(mid*mid),0)
   if np.sum(z*z)<12:lo=mid
   else:hi=mid
  radius=(lo+hi)/2;assert radius<d[-1]
  select=d<radius;hs=np.sqrt(a['volume'][i]);xx,yy=((x[ids[select]]-x[i])/hs).T
  taper=np.maximum(1-d2[select]/radius**2,0)**2
  w=taper*hs/d[select]
  design=w[:,None]*np.column_stack((xx,yy,xx**2/2,xx*yy/np.sqrt(2),yy**2/2,xx**3/6,xx**2*yy/(2*np.sqrt(3)),xx*yy**2/(2*np.sqrt(3)),yy**3/6))
  singular=np.linalg.svd(design,compute_uv=False)
  rank=int(np.count_nonzero(singular>singular[0]*max(design.shape)*np.finfo(float).eps))
  condition=float(singular[0]/singular[-1]);ranks.append(rank);conditions.append(condition)
  radii.append(radius/hs);counts.append(int(select.sum()));residuals.append(abs(float(taper.sum())-12))
  oldcounts.append(len(tree.query_ball_point(x[i],2.5*hs))-1)
  if worst is None or condition>worst['condition']:worst={'owner_id':int(a['id'][i])&72057594037927935,'condition':condition,'support_count':int(select.sum()),'radius_over_sqrtV':radius/hs}
 def stats(v):return {'min':float(np.min(v)),'median':float(np.median(v)),'max':float(np.max(v))}
 out['cases'][label]={'time':t,'snapshot':str(path),'sha256':sha(path),'cells':len(owners),'minimum_rank':min(ranks),'rank_below9':sum(r<9 for r in ranks),'condition':stats(conditions),'radius_over_sqrtV':stats(radii),'support_count':stats(counts),'fixed_support_count':stats(oldcounts),'max_kernel_sum_residual':max(residuals),'worst_condition':worst}
dest=b/'adaptive_qr_support_audit_20261018.json'
with dest.open('x') as f:json.dump(out,f,indent=2)
print(json.dumps(out,indent=2))
