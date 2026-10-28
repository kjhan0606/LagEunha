"""Classify excluded donor-coincident faces and the equal-weight reduction boundary."""
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
from scipy.spatial import ConvexHull
results=[]
for fraction in (.25,1.):
 weights=fraction*target;wdot=target.copy() # specified instantaneous rate, not an adopted evolution law
 hull=ConvexHull(np.column_stack((x,np.sum(x*x,axis=1)-weights)));tri=hull.simplices[hull.equations[:,2]<-1e-12];edges=np.unique(np.sort(np.concatenate((tri[:,[0,1]],tri[:,[0,2]],tri[:,[1,2]])),axis=1),axis=0);nb=[[] for _ in x]
 for i,j in edges:nb[i].append(j);nb[j].append(i)
 new=mod.polygons(x,weights,ids,geo.clip);local={int(i):k for k,i in enumerate(ids)};balance=np.zeros((len(ids),5));fluxes=[];excluded=0;kinerr=0.;algebra=0.;quaderr=0.;segments=0;eligible=0;excluded_info=[]
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
   cand=np.flatnonzero(np.all(hi>=np.minimum(aa,bb)-1e-12,axis=1)&np.all(lo<=np.maximum(aa,bb)+1e-12,axis=1));intervals=[]
   for src in cand:
    low,high=0.,1.
    for normal,point in zip(normals[src],old[src]):
     den=normal@d;rhs=normal@(point-aa)
     if abs(den)<1e-18:
      if rhs < -1e-14:high=-1.;break
     elif den>0:high=min(high,rhs/den)
     else:low=max(low,rhs/den)
    if high-low>1e-10:intervals.append((low,high,src))
   intervals.sort();cursor=0.;ok=True
   for low,high,src in intervals:
    if abs(low-cursor)>1e-8:ok=False;break
    cursor=high
   if not ok or abs(cursor-1)>1e-8:
    excluded+=1
    jump=(vel[j]-vel[i])@n;psum=a['pressure'][i]+a['pressure'][j]
    excluded_info.append({'i':i,'j':j,'r':float(np.linalg.norm(mid)),'weight_difference_over_distance_squared':float((weights[i]-weights[j])/distance**2),'rate_difference_over_distance_squared':float((wdot[i]-wdot[j])/distance**2),'both_weights_exact_zero':bool(weights[i]==0 and weights[j]==0),'both_rates_exact_zero':bool(wdot[i]==0 and wdot[j]==0),'normal_velocity_jump':float(jump),'relative_pressure_jump':float(abs(a['pressure'][i]-a['pressure'][j])/psum),'intervals':len(intervals),'interval_length_sum':float(sum(b-a for a,b,_ in intervals))})
    continue
   segments+=len(intervals);orderflux=[]
   for order in (8,12):
    z,g=np.polynomial.legendre.leggauss(order);F=np.zeros(5);direct=np.zeros(5)
    for low,high,src in intervals:
     u=low+(z+1)*.5*(high-low);xy=aa+u[:,None]*d;dl=g*.5*(high-low)*length;v=vel[ids[src]];rel=xy-x[ids[src]];feature=rel@coeff[:2]+(rel[:,0]*v[1]-rel[:,1]*v[0])*coeff[2]+(mod.potential(xy)-mod.potential(x[ids[src]][None,:])[0])*coeff[3]-means[src];density=rho[src]*(1+rho[src]*feature);I0=dl@density;I1=dl@(density*u);fm=(v@n-w0)*I0-w1*I1;K=.5*(v@v);ud=a['ie'][ids[src]]/lift['volume'][src]
     FU=ud*length*((v@n-w0)*(high-low)-.5*w1*(high*high-low*low));F+=np.r_[fm,fm*v,fm*K,FU]
     relative=v@n-w0-w1*u;md=density*relative*dl;direct+=np.r_[md.sum(),md.sum()*v,md.sum()*K,ud*(relative@dl)]
    algebra=max(algebra,float(abs(F-direct).max()));orderflux.append(F)
   quaderr=max(quaderr,float(abs(orderflux[0]-orderflux[1]).max()));F=orderflux[1];balance[k]-=F;balance[local[j]]+=F;fluxes.append(F)
 row={'fraction':fraction,'specified_Wdot':'target weights per unit diagnostic time','eligible_unique_internal_faces':eligible,'retained_faces':len(fluxes),'excluded_ambiguous_faces':excluded,'segments':segments,'kinematic_identity_max_residual':float(kinerr),'two_moment_vs_direct_max_abs':algebra,'order8_vs12_max_abs':quaderr,'paired_accumulation_global_sum':balance.sum(0).tolist(),'flux_components':['mass','px','py','kinetic','thermal'],'sum_absolute_flux':np.sum(abs(np.array(fluxes)),axis=0).tolist()};row['excluded_classification']={'both_weights_and_rates_exact_zero':sum(r['both_weights_exact_zero'] and r['both_rates_exact_zero'] for r in excluded_info),'equal_weights_and_rates_scaled_1e12':sum(abs(r['weight_difference_over_distance_squared'])<1e-12 and abs(r['rate_difference_over_distance_squared'])<1e-12 for r in excluded_info),'overlapping_intervals':sum(r['interval_length_sum']>1+1e-8 for r in excluded_info),'undercovered_intervals':sum(r['interval_length_sum']<1-1e-8 for r in excluded_info),'r_min':min(r['r'] for r in excluded_info),'r_max':max(r['r'] for r in excluded_info),'normal_jump_abs_median':float(np.median([abs(r['normal_velocity_jump']) for r in excluded_info])),'normal_jump_abs_max':max(abs(r['normal_velocity_jump']) for r in excluded_info),'relative_pressure_jump_max':max(r['relative_pressure_jump'] for r in excluded_info),'largest_jump_examples':sorted(excluded_info,key=lambda r:abs(r['normal_velocity_jump']),reverse=True)[:5]};results.append(row);print(json.dumps(row['excluded_classification']),flush=True)
result={'scope':__doc__,'rows':results,'limits':['Specified Wdot is not a proposed weight evolution equation.','One owner per retained face, paired accumulation by construction; not independent MPI rank evaluation.','Ambiguous coincident donor boundaries excluded; not a complete closed mesh update.','Advective kinetic and thermal transport only; no pressure work, gravity, shock/Riemann solution, physical time integration or reduction proof.']}
with (root/'excluded_face_classification.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
