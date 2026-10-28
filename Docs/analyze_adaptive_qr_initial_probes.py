"""Validate actual exported smooth-QR operator and legacy force path."""
import json,sys
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read
b=Path('/gpfs/kjhan/LagForce');root=b/'Codex_adaptive_qr_periodic_probes_20261018';old=b/'Codex_coupled_response_forces_t3p4_20261018'
def rhs(p):
 a=np.loadtxt(p);return a[np.argsort(a[:,0])]
for name in ['off_plus','on_plus','on_minus']:assert (root/name/'exit.code').read_text().strip()=='0'
a=rhs(root/'off_plus/raw_faces.dat.rhs.rank0000');ref=rhs(old/'plus/raw_faces.dat.rhs.rank0000')
out={'legacy_off':{'all_RHS_columns_equal':bool(np.array_equal(a,ref)),'max_absolute_difference':float(abs(a-ref).max())},'smooth':{}}
for sign in ['plus','minus']:
 d=root/('on_'+sign);s=rhs(d/'raw_faces.dat.rhs.rank0000');a,t=read(d/'khout.000743.dat');n=len(a)
 assert np.array_equal(s[:,0],a['id'] & 72057594037927935)
 x=np.column_stack((a['x'],a['y']));v=np.column_stack((a['vx'],a['vy']));mask=np.ones(n,dtype=bool)
 for k,col in [('x',1),('y',2),('vx',3),('vy',4),('mass',7)]:assert np.array_equal(s[:,col],a[k])
 assert np.max(abs(s[mask,6]/a['volume'][mask]-1))<1e-10
 assert np.max(abs(s[mask,5]/a['pressure'][mask]-1))<1e-10
 op=np.fromfile(d/'qr_operator.bin',dtype=np.float64).reshape(-1,6)
 owner=op[:,0].astype(int);source=op[:,1].astype(int)
 assert np.all((owner>=0)&(owner<n)) and np.all((source>=0)&(source<n)) and np.isfinite(op).all()
 grad=np.column_stack([np.bincount(owner,weights=op[:,axis]*v[source,component],minlength=n) for component in [0,1] for axis in [2,3]])
 graderr=float(np.max(abs(grad[mask]-s[mask,21:25])))
 # Export is the unlimited QR map; RHS reports stress gradients AFTER
 # alpha_vx/alpha_vy limiting (exam.c). Compare these via that relationship.
 limiter_error=0;limited_components=0
 for c in [0,2]:
  g=grad[mask,c:c+2];h=s[mask,21+c:23+c];den=np.sum(g*g,axis=1)
  alpha=np.divide(np.sum(g*h,axis=1),den,out=np.zeros_like(den),where=den>1e-24)
  limiter_error=max(limiter_error,float(abs(h-alpha[:,None]*g).max()))
  assert np.all(alpha>=-1e-5) and np.all(alpha<=1+1e-5)
  limited_components+=int(np.count_nonzero(alpha<.99999))
 assert limiter_error<1e-5
 hs=np.sqrt(a['volume']);X=op[:,4]/hs[owner];Y=op[:,5]/hs[owner]
 maxerr=0
 for px,py in [(0,0),(1,0),(0,1),(2,0),(1,1),(0,2),(3,0),(2,1),(1,2),(0,3)]:
  monomial=X**px*Y**py
  for axis in [2,3]:
   val=np.bincount(owner,weights=op[:,axis]*hs[owner]*monomial,minlength=n)
   expected=float((axis==2 and (px,py)==(1,0)) or (axis==3 and (px,py)==(0,1)))
   maxerr=max(maxerr,float(abs(val[mask]-expected).max()))
 assert maxerr<1e-9
 order=np.argsort(owner,kind='stable');counts=np.bincount(owner,minlength=n);starts=np.r_[0,np.cumsum(counts)]
 tree=cKDTree(x,boxsize=6.0);mismatches=[];independent_error=0;minrank=9
 for i in np.flatnonzero(mask):
  actual=set(source[order[starts[i]:starts[i+1]]]);actual.discard(i)
  distances,near=tree.query(x[i],k=129)
  distances=distances[near!=i];near=near[near!=i]
  lo=0.;hi=float(distances[-1]);d2=distances**2
  assert np.sum(np.maximum(1-d2/hi**2,0)**2)>12
  for iteration in range(55):
   mid=(lo+hi)/2
   if np.sum(np.maximum(1-d2/mid**2,0)**2)<12:lo=mid
   else:hi=mid
  radius=(lo+hi)/2;assert radius<distances[-1]
  expected=set(near[distances<radius])
  if actual!=expected:mismatches.append({'owner':int(i),'missing':sorted(expected-actual),'extra':sorted(actual-expected)})
  ids=np.array(sorted(expected));delta=x[ids]-x[i];delta-=6*np.rint(delta/6);xx,yy=(delta/hs[i]).T;rr=np.hypot(xx,yy)
  w=np.maximum(1-(rr*hs[i]/radius)**2,0)**2/rr
  design=w[:,None]*np.column_stack((xx,yy,xx**2/2,xx*yy/np.sqrt(2),yy**2/2,xx**3/6,xx**2*yy/(2*np.sqrt(3)),xx*yy**2/(2*np.sqrt(3)),yy**3/6))
  fit,_,rank,_=np.linalg.lstsq(design,w[:,None]*(v[ids]-v[i]),rcond=None)
  minrank=min(minrank,int(rank));independent_error=max(independent_error,float(abs((fit[:2].T/hs[i]).reshape(-1)-grad[i]).max()))
 assert minrank==9 and independent_error<1e-7
 assert not mismatches, 'Production support differs from full periodic point cloud'
 out['smooth'][sign]={'time':t,'all_periodic_cells':int(mask.sum()),'max_exported_gradient_error':graderr,
  'max_limiter_collinearity_error':limiter_error,'limited_velocity_components':limited_components,
  'max_independent_smooth_fit_error':independent_error,'minimum_independent_rank':minrank,
  'max_normalized_cubic_moment_error':maxerr,'support_mismatch_count':len(mismatches),'support_mismatches_first10':mismatches[:10]}
out['scope']='Adaptive production terminal probes, actual QR export, independent kernel-sum radius and full periodic-box support for all cells (box length6). No time evolution or long-time stability claim.'
with (root/'analysis.json').open('x') as f:json.dump(out,f,indent=2)
print(json.dumps(out,indent=2))
