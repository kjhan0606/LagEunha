"""Certified lower-envelope segmentation and one-sided donor limits on moving weighted faces."""
import json, importlib.util, hashlib
from pathlib import Path
import numpy as np

def load(name,p):
    s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
base=Path('/gpfs/kjhan/LagForce');root=base/'Codex_checkpoint_weight_remap_20261018'
mod=load('mom',base/'Codex_weight_remap_moments_centered_20261017/audit.py')
geo=load('geo',base/'Codex_weighted_compression_geometry_20261017/audit.py')
reader=load('rd',base/'run_canonical_av_branches_20261018.py')
a,t=reader.read(base/'Codex_HG_constrained_20261017/omp8_t20/khout.000184.dat')
x=np.column_stack((a['x']-3,a['y']-3));vel=np.column_stack((a['vx'],a['vy']))
lift=np.load(root/'lift_coefficients.npz');ids=lift['ids'];rho=a['mass'][ids]/lift['volume']
lam=lift['coefficient'];coeff=lam/lift['feature_scale'];means=lift['cell_feature_mean']@lam
old=mod.polygons(x,np.zeros(len(x)),ids,geo.clip)
lo=np.array([p.min(0) for p in old]);hi=np.array([p.max(0) for p in old])
normals=[np.column_stack(((np.roll(p,-1,axis=0)-p)[:,1],-(np.roll(p,-1,axis=0)-p)[:,0])) for p in old]
target=np.load(base/'Codex_evolving_weight_transport_20261017/weights_0.npz')['weights']
from scipy.spatial import ConvexHull,cKDTree
tree=cKDTree(x)
results=[];all_faces=[]
for fraction in (.25,1.):
 weights=fraction*target;wdot=target.copy() # specified instantaneous rate, not an adopted evolution law
 hull=ConvexHull(np.column_stack((x,np.sum(x*x,axis=1)-weights)));tri=hull.simplices[hull.equations[:,2]<-1e-12];edges=np.unique(np.sort(np.concatenate((tri[:,[0,1]],tri[:,[0,2]],tri[:,[1,2]])),axis=1),axis=0);nb=[[] for _ in x]
 for i,j in edges:nb[i].append(j);nb[j].append(i)
 new=mod.polygons(x,weights,ids,geo.clip);local={int(i):k for k,i in enumerate(ids)};balance=np.zeros((len(ids),5));fluxes=[];excluded=0;kinerr=0.;algebra=0.;quaderr=0.;segments=0;eligible=0;max_certificate_error=0.;totalsegments=0;two_sided=0;maxsegments=0;outside=0;coverage_error=0.
 for k,p in enumerate(new):
  i=int(ids[k])
  if np.linalg.norm(x[i])>=.6:continue
  for aa,bb in zip(p,np.roll(p,-1,axis=0)):
   d=bb-aa;length=np.linalg.norm(d)
   if length<1e-10:continue
   mid=(aa+bb)/2;candidates=np.array(nb[i]);power=np.sum((x[candidates]-mid)**2,axis=1)-weights[candidates];j=int(candidates[np.argmin(abs(power-(np.sum((x[i]-mid)**2)-weights[i])))])
   if j not in local or i>=j:continue
   assert abs(np.sum((x[j]-mid)**2)-weights[j]-(np.sum((x[i]-mid)**2)-weights[i]))<1e-9
   eligible+=1;separation=x[j]-x[i];distance=np.linalg.norm(separation);n=separation/distance
   w0=(x[j]@vel[j]-x[i]@vel[i]-aa@(vel[j]-vel[i])+.5*(wdot[i]-wdot[j]))/distance;w1=-d@(vel[j]-vel[i])/distance
   # Differential weighted-bisector identity at both endpoints and midpoint.
   for u in (0.,.5,1.):
    pos=aa+u*d;normal_speed=w0+w1*u
    kinerr=max(kinerr,abs(2*distance*normal_speed+2*pos@(vel[j]-vel[i])-2*x[j]@vel[j]+2*x[i]@vel[i]-wdot[i]+wdot[j]))
   # Every nearest donor lies in this ball: a face owner supplies an upper
   # distance bound everywhere on the segment, via convexity of squared norm.
   bound=.5*length+max(np.linalg.norm(aa-x[i]),np.linalg.norm(bb-x[i]));cand=np.array(tree.query_ball_point(mid,bound+1e-12),dtype=int)
   offset=aa-x[cand];intercept=np.sum(offset*offset,axis=1);slope=2*(offset@d);tol=256*np.finfo(float).eps*max(1.,float(abs(intercept).max()),float(abs(slope).max()));cuts=[0.,1.]
   for l in range(len(cand)):
    for rr in range(l+1,len(cand)):
     den=slope[l]-slope[rr]
     if abs(den)>tol:
      u=(intercept[rr]-intercept[l])/den
      if 0<u<1:cuts.append(float(u))
   cuts=sorted(cuts);unique=[cuts[0]]
   for u in cuts[1:]:
    if u-unique[-1]>1e-12:unique.append(u)
   if unique[-1]<1:unique.append(1.)
   raw=[]
   for low,high in zip(unique[:-1],unique[1:]):
    if high<=low:continue
    u=.5*(low+high);values=intercept+slope*u;ties=np.flatnonzero(values<=values.min()+tol);pos=aa+u*d;normal_derivative=2*((pos-x[cand[ties]])@n)
    left=int(cand[ties[np.argmax(normal_derivative)]]);right=int(cand[ties[np.argmin(normal_derivative)]])
    for owner in (left,right):
     kk=int(np.flatnonzero(cand==owner)[0]);ends=np.array([low,high]);candidate_values=intercept[:,None]+slope[:,None]*ends;gap=float(np.max(candidate_values[kk]-candidate_values.min(axis=0)));max_certificate_error=max(max_certificate_error,gap);assert gap<=4*tol,(i,j,low,high,gap,tol)
    if raw and raw[-1][2:]==[left,right]:raw[-1][1]=high
    else:raw.append([low,high,left,right])
   covered=sum(row[1]-row[0] for row in raw);coverage_error=max(coverage_error,abs(covered-1));assert abs(covered-1)<1e-12
   totalsegments+=len(raw);two_sided+=sum(row[2]!=row[3] for row in raw);maxsegments=max(maxsegments,len(raw));outside+=sum(row[2] not in local or row[3] not in local for row in raw)
   all_faces.append({'fraction':fraction,'i':i,'j':j,'a':aa.tolist(),'b':bb.tolist(),'normal':n.tolist(),'normal_speed_coefficients':[float(w0),float(w1)],'segments':raw})
 row={'fraction':fraction,'unique_internal_faces':eligible,'segments':totalsegments,'two_sided_discontinuous_segments':two_sided,'max_segments_per_face':maxsegments,'donor_pairs_outside_patch':outside,'max_nearest_owner_certificate_error':max_certificate_error,'max_partition_coverage_error':coverage_error,'kinematic_identity_max_residual':kinerr};results.append(row);print(json.dumps(row),flush=True)
with (root/'certified_face_segments.json').open('x') as f:json.dump({'scope':__doc__,'rows':results,'faces':all_faces,'limits':['Finite floating tolerance, not exact arithmetic.','Geometry and one-sided states only; no Riemann flux or physical evolution.','Selected internal-face set only, not all domain boundary faces.']},f);f.write('\n')
