"""One fixed weight path, subdivided; audit projection mixing and represented angular momentum."""
from pathlib import Path
import importlib.util,json,time
import numpy as np
from scipy.sparse import coo_matrix

def load(n,p):
 s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
base=Path('/gpfs/kjhan/LagForce');out=base/'Codex_remap_fixed_path_20261018';out.mkdir(exist_ok=False)
mod=load('mom',base/'Codex_weight_remap_moments_centered_20261017/audit.py');geo=load('geo',base/'Codex_weighted_compression_geometry_20261017/audit.py');rd=load('rd',base/'run_canonical_av_branches_20261018.py')
a,t=rd.read(base/'Codex_HG_constrained_20261017/omp8_t20/khout.000184.dat');x=np.column_stack((a['x']-3,a['y']-3));ids=np.flatnonzero(np.linalg.norm(x,axis=1)<.7);vel=np.column_stack((a['vx'][ids],a['vy'][ids]));m=a['mass'][ids];U=a['ie'][ids]
# Uniform spatial lift per original cell, retaining exact m,p,Pi,U aggregates.
q0=np.column_stack((m,m[:,None]*vel,m*vel[:,0]**2,m*vel[:,0]*vel[:,1],m*vel[:,1]**2,U))
old=mod.polygons(x,np.zeros(len(x)),ids,geo.clip);vol=np.array([mod.area(p) for p in old]);lo=np.array([p.min(0) for p in old]);hi=np.array([p.max(0) for p in old]);target=np.load(base/'Codex_evolving_weight_transport_20261017/weights_0.npz')['weights'];rows=[]
def centroid(poly):
 p=poly-poly[0];b=np.roll(p,-1,axis=0);cross=p[:,0]*b[:,1]-p[:,1]*b[:,0]
 return poly[0]+np.sum(cross[:,None]*(p+b),axis=0)/(3*cross.sum())
centers=np.array([centroid(p) for p in old])
L0=float(np.sum(centers[:,0]*q0[:,2]-centers[:,1]*q0[:,1]))
# A transported scalar L ledger is intentionally distinct from the spatially
# represented integral x cross rho*v; this comparison tests that distinction.
q0=np.column_stack((q0,centers[:,0]*q0[:,2]-centers[:,1]*q0[:,1]))
def transition(source,dest):
 sv=np.array([mod.area(p) for p in source]);dv=np.array([mod.area(p) for p in dest]);lo=np.array([p.min(0) for p in source]);hi=np.array([p.max(0) for p in source]);ii=[];jj=[];aa=[]
 for j,p in enumerate(dest):
  candidates=np.flatnonzero(np.all(hi>=p.min(0),axis=1)&np.all(lo<=p.max(0),axis=1))
  for i in candidates:
   overlap=source[i].copy()
   for b,c in zip(p,np.roll(p,-1,axis=0)):
    e=c-b;n=np.array([e[1],-e[0]]);overlap=geo.clip(overlap,n,n@b)
    if len(overlap)<3:break
   av=mod.area(overlap)
   if av>1e-18:ii.append(j);jj.append(i);aa.append(av)
 A=coo_matrix((aa,(ii,jj)),shape=(len(ids),len(ids))).tocsr()
 error=max(np.max(abs(np.asarray(A.sum(0)).ravel()/sv-1)),np.max(abs(np.asarray(A.sum(1)).ravel()/dv-1)));assert error<1e-8
 return A.multiply(1/sv[None,:]).tocsr(),A.T.multiply(1/dv[None,:]).tocsr(),error
for steps in (1,2,4,8):
 start=time.monotonic();source=old;q=q0.copy();reverse=[];error=0.
 for i in range(1,steps+1):
  dest=mod.polygons(x*(1.),(.25*i/steps)*target,ids,geo.clip)
  forward,back,err=transition(source,dest);q=forward@q;reverse.append(back);error=max(error,err);source=dest
 for back in reversed(reverse):q=back@q
 mean=q[:,1:3]/q[:,0,None];internal=.5*(q[:,3]+q[:,5])-.5*q[:,0]*np.sum(mean**2,axis=1)
 assert internal.min()>-1e-13 and q[:,6].min()>0
 represented=float(np.sum(centers[:,0]*q[:,2]-centers[:,1]*q[:,1]))
 row={'subdivisions_per_leg':steps,'total_fraction_path':.5,'internal_K_over_original_U':float(internal.sum()/U.sum()),'velocity_change_mass_RMS':float(np.sqrt(np.average(np.sum((mean-vel)**2,axis=1),weights=m))),'represented_L_change':represented-L0,'represented_L_relative_change':(represented-L0)/L0,'transported_L_ledger_change':float(q[:,7].sum()-L0),'aggregate_change_m_px_py_Pixx_Pixy_Piyy_U_Lledger':(q.sum(0)-q0.sum(0)).tolist(),'coverage_max_relative_error':float(error),'seconds':time.monotonic()-start}
 rows.append(row);np.savez(out/f'subdivisions{steps}.npz',ids=ids,state=q);print(json.dumps(row),flush=True)
result={'scope':__doc__,'source_time':t,'path':'0 -> .25 -> 0, one reversal; same path for all cases','initial_represented_L':L0,'rows':rows,'limits':['Uniform spatial lift, not fitted-density lift or exact original point Lz/Phi representation.','Piecewise constant projection only, not affine or production solver.','A scalar angular momentum ledger is not claimed to be a consistent spatial moment closure.','No pressure, gravity, physical evolution, MPI or CUDA test in this diagnostic.']}
(out/'result.json').write_text(json.dumps(result,indent=2)+'\n')
