"""Affine velocity plus residual covariance under uniform spatial density projection."""
from pathlib import Path
import importlib.util,json,time
import numpy as np
from scipy.sparse import coo_matrix

def load(n,p):
 s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
base=Path('/gpfs/kjhan/LagForce');out=base/'Codex_remap_affine_path_20261018';out.mkdir(exist_ok=False)
mod=load('mom',base/'Codex_weight_remap_moments_centered_20261017/audit.py');geo=load('geo',base/'Codex_weighted_compression_geometry_20261017/audit.py');rd=load('rd',base/'run_canonical_av_branches_20261018.py')
a,t=rd.read(base/'Codex_HG_constrained_20261017/omp8_t20/khout.000184.dat');x=np.column_stack((a['x']-3,a['y']-3));ids=np.flatnonzero(np.linalg.norm(x,axis=1)<.7);vel=np.column_stack((a['vx'][ids],a['vy'][ids]));m=a['mass'][ids];U=a['ie'][ids]
# Uniform spatial lift per original cell, retaining exact m,p,Pi,U aggregates.
q0=np.column_stack((m,m[:,None]*vel,m*vel[:,0]**2,m*vel[:,0]*vel[:,1],m*vel[:,1]**2,U))
old=mod.polygons(x,np.zeros(len(x)),ids,geo.clip);vol=np.array([mod.area(p) for p in old]);lo=np.array([p.min(0) for p in old]);hi=np.array([p.max(0) for p in old]);target=np.load(base/'Codex_evolving_weight_transport_20261017/weights_0.npz')['weights'];rows=[]
def spatial(poly):
 xy,dv=mod.quadrature(poly,2);V=dv.sum();c=dv@xy/V;r=xy-c;S=(r.T*dv)@r/V
 return V,c,S

def initial(polys):
 q=[]
 for k,p in enumerate(polys):
  V,c,S=spatial(p);v=vel[k];D=m[k]*np.outer(c,v);Pi=m[k]*np.outer(v,v)
  q.append(np.r_[m[k],m[k]*v,D.ravel(),Pi[0,0],Pi[0,1],Pi[1,1],U[k]])
 return np.array(q)
qinit=initial(old)
def recover(polys,q):
 data=[];fail=[];min_eig=float('inf');max_error=0.
 for k,p in enumerate(polys):
  V,c,s=spatial(p);mass=q[k,0];v=q[k,1:3]/mass;D=q[k,3:7].reshape(2,2);Pi=np.array([[q[k,7],q[k,8]],[q[k,8],q[k,9]]]);C=D-mass*np.outer(c,v);S=mass*s;G=np.linalg.solve(S,C).T;R=Pi-mass*np.outer(v,v)-G@S@G.T;R=.5*(R+R.T);eig=np.linalg.eigvalsh(R)
  tol=1e-12*max(float(np.max(abs(Pi))),1e-20)+1e-20
  min_eig=min(min_eig,float(eig.min()))
  if eig.min() < -tol:fail.append({'id':int(ids[k]),'radius':float(np.linalg.norm(x[ids[k]])),'minimum_R_eigenvalue':float(eig.min()),'tolerance':tol,'thermal_U':float(q[k,10]),'K_affine':float(.5*np.trace(G@S@G.T)),'K_internal':float(.5*np.trace(Pi-mass*np.outer(v,v)))})
  rebuilt=mass*np.outer(c,v)+S@G.T;max_error=max(max_error,float(abs(rebuilt-D).max()))
  data.append((V,c,mass/V,v,G,R/mass,q[k,10]/V))
 return data,fail,min_eig,max_error

def advance(source,dest,q,data):
 lo=np.array([p.min(0) for p in source]);hi=np.array([p.max(0) for p in source]);post=np.zeros_like(q);coverage=np.zeros(len(ids));destcoverage=np.zeros(len(ids))
 for j,p in enumerate(dest):
  candidates=np.flatnonzero(np.all(hi>=p.min(0),axis=1)&np.all(lo<=p.max(0),axis=1))
  for i in candidates:
   overlap=source[i].copy()
   for b,c in zip(p,np.roll(p,-1,axis=0)):
    e=c-b;n=np.array([e[1],-e[0]]);overlap=geo.clip(overlap,n,n@b)
    if len(overlap)<3:break
   av=mod.area(overlap)
   if av<=1e-18:continue
   xy,dv=mod.quadrature(overlap,2);V,c,rho,v,G,res,ud=data[i];vv=v+(xy-c)@G.T;dm=rho*dv;D=np.einsum('n,ni,nj->ij',dm,xy,vv);Pi=np.einsum('n,ni,nj->ij',dm,vv,vv)+dm.sum()*res
   post[j]+=np.r_[dm.sum(),dm@vv,D.ravel(),Pi[0,0],Pi[0,1],Pi[1,1],ud*dv.sum()];coverage[i]+=dv.sum();destcoverage[j]+=dv.sum()
 sv=np.array([d[0] for d in data]);nv=np.array([spatial(p)[0] for p in dest]);err=max(float(abs(coverage/sv-1).max()),float(abs(destcoverage/nv-1).max()));assert err<1e-8
 return post,err

for steps in (1,2,4,8):
 source=old;q=qinit.copy();history=[];failed=False;started=time.monotonic()
 path=list(np.linspace(0,.25,steps+1)[1:])+list(np.linspace(.25,0,steps+1)[1:])
 for leg,fraction in enumerate(path):
  data,fail,mineig,derr=recover(source,q)
  if fail:
   history.append({'at_before_transition':leg,'failures':fail,'minimum_R_eigenvalue':mineig,'max_D_recovery_error':derr});failed=True;break
  dest=mod.polygons(x,fraction*target,ids,geo.clip);q,coverage=advance(source,dest,q,data);source=dest
  history.append({'transition':leg,'weight_fraction':float(fraction),'coverage_error':coverage,'max_aggregate_change':float(abs(q.sum(0)-qinit.sum(0)).max())})
 if not failed:
  data,fail,mineig,derr=recover(source,q)
  if fail:history.append({'at_final_recovery':True,'failures':fail,'minimum_R_eigenvalue':mineig,'max_D_recovery_error':derr});failed=True
 row={'subdivisions_per_leg':steps,'valid_through_full_path':not failed,'history':history,'seconds':time.monotonic()-started,'aggregate_change_m_p_D_Pi_U':(q.sum(0)-qinit.sum(0)).tolist()}
 if not failed:
  mean=q[:,1:3]/q[:,0,None];internal=.5*(q[:,7]+q[:,9])-.5*q[:,0]*np.sum(mean**2,axis=1)
  row.update(internal_K_over_original_U=float(internal.sum()/U.sum()),velocity_change_mass_RMS=float(np.sqrt(np.average(np.sum((mean-vel)**2,axis=1),weights=m))),represented_L_change=float(np.sum(q[:,4]-q[:,5])-np.sum(qinit[:,4]-qinit[:,5])),minimum_R_eigenvalue=mineig,max_D_recovery_error=derr)
 rows.append(row);np.savez(out/f'subdivisions{steps}.npz',ids=ids,state=q);print(json.dumps(row),flush=True)
result={'scope':__doc__,'rows':rows,'limits':['Uniform density is reimposed on every receiver; mass spatial moments are not preserved.','m,p,D,Pi,U retained; R is explicit nonthermal covariance. It is not a resolved spatial velocity realization.','Negative R beyond tolerance stops that case; no clipping or added energy.','No physical evolution, pressure/gravity, production checkpoint, MPI or backend parity test.']}
(out/'result.json').write_text(json.dumps(result,indent=2)+'\n')
