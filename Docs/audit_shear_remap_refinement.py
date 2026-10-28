"""Analytic strip-cell shear remap refinement; controlled spatial consistency only."""
from pathlib import Path
import json,numpy as np
root=Path('/gpfs/kjhan/LagForce/Codex_shear_remap_refinement_v2_20261018');root.mkdir(exist_ok=False)
# Unit-width 2D strips, v=(y,0), rho=1, y in[0,1]. q=[m,Y,YY,px,Ypx,Pixx].
def moments(a,b):return np.array([b-a,(b*b-a*a)/2,(b**3-a**3)/3])
def init(edges,kind):
 q=[]
 for a,b in zip(edges[:-1],edges[1:]):
  m,Y,YY=moments(a,b);c=Y/m
  q.append([m,Y,YY,m*c,Y*c,m*c*c] if kind=='point' else [m,Y,YY,Y,YY,YY])
 return np.array(q)
def remap(q,old,new,kind):
 out=np.zeros_like(q);minR=float('inf')
 for i,(a,b) in enumerate(zip(old[:-1],old[1:])):
  m,Y,YY,p,D,Pi=q[i];c=Y/m;v=p/m;S=YY-Y*c;C=D-Y*v;g=C/S if kind=='affine' else 0.;R=Pi-p*v-g*g*S;minR=min(minR,R)
  assert R > -1e-12
  # At constant projection, resetting D to Y*v discards cross moment.
  intercept=v-g*c;rho=m/(b-a);variance=R/m
  start=max(0,np.searchsorted(new,a,side='right')-1);stop=min(len(q),np.searchsorted(new,b,side='left')+1)
  for j in range(start,stop):
   low=max(a,new[j]);high=min(b,new[j+1])
   if high<=low:continue
   z0,z1,z2=rho*moments(low,high)
   out[j]+=np.array([z0,z1,z2,intercept*z0+g*z1,intercept*z1+g*z2,(intercept*intercept+variance)*z0+2*intercept*g*z1+g*g*z2])
 if kind=='constant':out[:,4]=out[:,1]*out[:,3]/out[:,0]
 return out,minR
rows=[]
for N in (16,32,64,128):
 h=1/N;original=np.linspace(0,1,N+1);X=(original[:-1]+original[1:])/2;delta=.15*h*(-1.)**np.arange(N+1);delta[[0,-1]]=0
 # Adjacent fixed generators obey face shift=(Wi-Wi+1)/(2h).
 weights=np.r_[0.,np.cumsum(-2*h*delta[1:-1])];face_residual=float(abs((weights[:-1]-weights[1:])/(2*h)-delta[1:-1]).max())
 for initial,reconstruction in [('point','constant'),('point','affine'),('exact_affine','affine')]:
  q=init(original,initial);q0=q.copy();old=original;minimum_R=float('inf')
  for fraction in np.r_[np.linspace(0,1,5)[1:],np.linspace(1,0,5)[1:]]:
   new=original+fraction*delta;q,mR=remap(q,old,new,reconstruction);minimum_R=min(minimum_R,mR);old=new
  m,Y,YY,p,D,Pi=q.T;c=Y/m;v=p/m;S=YY-Y*c;C=D-Y*v;g=C/S if reconstruction=='affine' else np.zeros(N);vg=v+g*(X-c);aff=.5*g*g*S;res=.5*(Pi-p*v-g*g*S)
  row={'N':N,'initial':initial,'reconstruction':reconstruction,'point_velocity_error_RMS':float(np.sqrt(np.average((vg-X)**2,weights=m))),'internal_K':float(.5*np.sum(Pi-p*v)),'initial_internal_K':float(.5*np.sum(q0[:,5]-q0[:,3]**2/q0[:,0])),'affine_K':float(aff.sum()),'residual_K':float(res.sum()),'mass_drift':float(m.sum()-q0[:,0].sum()),'momentum_drift':float(p.sum()-q0[:,3].sum()),'Pi_drift':float(Pi.sum()-q0[:,5].sum()),'represented_Lz_drift':float(-D.sum()+q0[:,4].sum()),'minimum_R':float(minimum_R),'weight_face_identity_error':face_residual}
  row['interior_point_velocity_error_RMS']=float(np.sqrt(np.mean((vg[(X>.25)&(X<.75)]-X[(X>.25)&(X<.75)])**2)))
  rows.append(row)
for kind in [('point','constant'),('point','affine')]:
 subset=[r for r in rows if (r['initial'],r['reconstruction'])==kind]
 for prev,r in zip(subset,subset[1:]):r['refinement_exponents']={k:float(np.log(prev[k]/r[k])/np.log(2)) for k in ('point_velocity_error_RMS','internal_K','residual_K')}
result={'scope':__doc__,'rows':rows,'limits':['Uniform-density rectangular strips, not irregular Kepler geometry or 2D vortex dynamics.','Weight amplitude scales as h^2; geometric shift as h. This holds deformation in cell units fixed, not absolute weights.','Point and exact-affine initializations differ in internal kinetic energy; exact-affine case is a separate consistency control, not an energy-free replacement of point data.','Pure remap, no physical time/pressure/gravity/AV, no CUDA/MPI validation.']}
(root/'result.json').write_text(json.dumps(result,indent=2)+'\n')
for r in rows:print(json.dumps(r))
