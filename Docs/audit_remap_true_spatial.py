"""First remap with actual transported spatial moments; uniform donor density."""
from pathlib import Path
import importlib.util,json,time
import numpy as np
from scipy.sparse import coo_matrix

def load(n,p):
 s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
base=Path('/gpfs/kjhan/LagForce');out=base/'Codex_remap_true_spatial_20261018';out.mkdir(exist_ok=False)
mod=load('mom',base/'Codex_weight_remap_moments_centered_20261017/audit.py');geo=load('geo',base/'Codex_weighted_compression_geometry_20261017/audit.py');rd=load('rd',base/'run_canonical_av_branches_20261018.py')
a,t=rd.read(base/'Codex_HG_constrained_20261017/omp8_t20/khout.000184.dat');x=np.column_stack((a['x']-3,a['y']-3));ids=np.flatnonzero(np.linalg.norm(x,axis=1)<.7);vel=np.column_stack((a['vx'][ids],a['vy'][ids]));m=a['mass'][ids];U=a['ie'][ids]
# Uniform spatial lift per original cell, retaining exact m,p,Pi,U aggregates.
q0=np.column_stack((m,m[:,None]*vel,m*vel[:,0]**2,m*vel[:,0]*vel[:,1],m*vel[:,1]**2,U))
old=mod.polygons(x,np.zeros(len(x)),ids,geo.clip);vol=np.array([mod.area(p) for p in old]);lo=np.array([p.min(0) for p in old]);hi=np.array([p.max(0) for p in old]);target=np.load(base/'Codex_evolving_weight_transport_20261017/weights_0.npz')['weights'];rows=[]
def recover(q,boost=None):
 mass=q[:,0];c=q[:,1:3]/mass[:,None];p=q[:,6:8].copy();D=q[:,8:12].reshape(-1,2,2).copy();Pi=np.array([[q[:,12],q[:,13]],[q[:,13],q[:,14]]]).transpose(2,0,1).copy()
 if boost is not None:
  b=np.asarray(boost);Pi+=p[:,:,None]*b[None,None,:]+b[None,:,None]*p[:,None,:]+mass[:,None,None]*np.outer(b,b);D+=q[:,1:3,None]*b[None,None,:];p+=mass[:,None]*b
 v=p/mass[:,None];Q=np.array([[q[:,3],q[:,4]],[q[:,4],q[:,5]]]).transpose(2,0,1);S=Q-mass[:,None,None]*c[:,:,None]*c[:,None,:];C=D-mass[:,None,None]*c[:,:,None]*v[:,None,:];T=Pi-mass[:,None,None]*v[:,:,None]*v[:,None,:]
 G=np.linalg.solve(S,C).transpose(0,2,1);R=T-G@S@G.transpose(0,2,1);R=.5*(R+R.transpose(0,2,1))
 return S,C,G,R
for steps in (1,2,4,8):
 fraction=.25/steps;new=mod.polygons(x,fraction*target,ids,geo.clip);q=np.zeros((len(ids),17));cover=np.zeros(len(ids));newcover=np.zeros(len(ids));start=time.monotonic()
 for j,p in enumerate(new):
  cand=np.flatnonzero(np.all(hi>=p.min(0),axis=1)&np.all(lo<=p.max(0),axis=1))
  for i in cand:
   poly=old[i].copy()
   for aa,bb in zip(p,np.roll(p,-1,axis=0)):
    e=bb-aa;n=np.array([e[1],-e[0]]);poly=geo.clip(poly,n,n@aa)
    if len(poly)<3:break
   if mod.area(poly)<=1e-18:continue
   xy,dv=mod.quadrature(poly,2);dm=m[i]/vol[i]*dv;v=vel[i];M=dm.sum();sp=dm@xy;Q=np.einsum('n,ni,nj->ij',dm,xy,xy);D=np.outer(sp,v);Pi=M*np.outer(v,v)
   q[j]+=np.r_[M,sp,Q[0,0],Q[0,1],Q[1,1],M*v,D.ravel(),Pi[0,0],Pi[0,1],Pi[1,1],U[i]/vol[i]*dv.sum()+.5*np.trace(Pi),0.]
   cover[i]+=dv.sum();newcover[j]+=dv.sum()
 nv=np.array([mod.area(p) for p in new]);err=max(float(abs(cover/vol-1).max()),float(abs(newcover/nv-1).max()));assert err<1e-8
 saved=np.load(base/'Codex_remap_affine_path_20261018'/f'subdivisions{steps}.npz')['state'];reduced=np.column_stack((q[:,0],q[:,6:12],q[:,12:15],q[:,15]-.5*(q[:,12]+q[:,14])))
 assert np.max(abs(reduced-saved))<1e-12
 S,C,G,R=recover(q);e=np.linalg.eigvalsh(R);tol=1e-12*np.max(abs(q[:,12:15]),axis=1)+1e-20;fail=e[:,0]<-tol
 Sb,Cb,Gb,Rb=recover(q,[7.,-3.])
 row={'subdivisions':steps,'fraction':fraction,'failed_R_cells':int(fail.sum()),'minimum_R_eigenvalue':float(e.min()),'minimum_S_eigenvalue':float(np.linalg.eigvalsh(S).min()),'saved_reduced_state_max_difference':float(abs(reduced-saved).max()),'boost_C_max_change':float(abs(Cb-C).max()),'boost_G_max_change':float(abs(Gb-G).max()),'boost_R_max_change':float(abs(Rb-R).max()),'coverage_max_relative_error':err,'minimum_U':float((q[:,15]-.5*(q[:,12]+q[:,14])).min()),'seconds':time.monotonic()-start}
 rows.append(row);np.savez(out/f'subdivisions{steps}.npz',ids=ids,state=q);print(json.dumps(row),flush=True)
result={'scope':__doc__,'rows':rows,'limits':['First remap only; spatial density reconstruction and repeated projection not provided.','True q,Q retained; agreement with failed reduced audit isolates spatial reset effect.','Raw-moment boost test has floating cancellation; no exact bitwise Galilean invariance claim.','No full Euler closure, production checkpoint conversion, MPI/backend parity or long-time Kepler claim.']}
(out/'result.json').write_text(json.dumps(result,indent=2)+'\n')
