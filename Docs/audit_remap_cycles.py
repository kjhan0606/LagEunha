"""Fixed-geometry weight cycles with piecewise constant projection, no physical evolution."""
from pathlib import Path
import importlib.util,json,time
import numpy as np
from scipy.sparse import coo_matrix

def load(n,p):
 s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
base=Path('/gpfs/kjhan/LagForce');out=base/'Codex_remap_cycles_20261018';out.mkdir(exist_ok=False)
mod=load('mom',base/'Codex_weight_remap_moments_centered_20261017/audit.py');geo=load('geo',base/'Codex_weighted_compression_geometry_20261017/audit.py');rd=load('rd',base/'run_canonical_av_branches_20261018.py')
a,t=rd.read(base/'Codex_HG_constrained_20261017/omp8_t20/khout.000184.dat');x=np.column_stack((a['x']-3,a['y']-3));ids=np.flatnonzero(np.linalg.norm(x,axis=1)<.7);vel=np.column_stack((a['vx'][ids],a['vy'][ids]));m=a['mass'][ids];U=a['ie'][ids]
# Uniform spatial lift per original cell, retaining exact m,p,Pi,U aggregates.
q0=np.column_stack((m,m[:,None]*vel,m*vel[:,0]**2,m*vel[:,0]*vel[:,1],m*vel[:,1]**2,U))
old=mod.polygons(x,np.zeros(len(x)),ids,geo.clip);vol=np.array([mod.area(p) for p in old]);lo=np.array([p.min(0) for p in old]);hi=np.array([p.max(0) for p in old]);target=np.load(base/'Codex_evolving_weight_transport_20261017/weights_0.npz')['weights'];rows=[]
for cycles in (1,4,16,64):
 start=time.monotonic();fraction=.25/cycles;new=mod.polygons(x,fraction*target,ids,geo.clip);nv=np.array([mod.area(p) for p in new]);ii=[];jj=[];aa=[]
 for dest,p in enumerate(new):
  candidates=np.flatnonzero(np.all(hi>=p.min(0),axis=1)&np.all(lo<=p.max(0),axis=1))
  for src in candidates:
   overlap=old[src].copy()
   for b,c in zip(p,np.roll(p,-1,axis=0)):
    e=c-b;n=np.array([e[1],-e[0]]);overlap=geo.clip(overlap,n,n@b)
    if len(overlap)<3:break
   area=mod.area(overlap)
   if area>1e-18:ii.append(dest);jj.append(src);aa.append(area)
 A=coo_matrix((aa,(ii,jj)),shape=(len(ids),len(ids))).tocsr()
 coverage=max(np.max(abs(np.asarray(A.sum(0)).ravel()/vol-1)),np.max(abs(np.asarray(A.sum(1)).ravel()/nv-1)));assert coverage<1e-8
 forward=A.multiply(1/vol[None,:]).tocsr();back=A.T.multiply(1/nv[None,:]).tocsr();q=q0.copy()
 for _ in range(cycles):q=back@(forward@q)
 mean=q[:,1:3]/q[:,0,None];k=.5*(q[:,3]+q[:,5]);internal=k-.5*q[:,0]*np.sum(mean**2,axis=1)
 assert internal.min()>-1e-13 and q[:,6].min()>0
 row={'cycles':cycles,'weight_fraction':fraction,'total_absolute_fraction_path':2*cycles*fraction,'coverage_max_relative_error':float(coverage),'aggregate_change_m_px_py_Pixx_Pixy_Piyy_U':(q.sum(0)-q0.sum(0)).tolist(),'internal_K_over_original_U':float(internal.sum()/U.sum()),'velocity_change_mass_RMS':float(np.sqrt(np.average(np.sum((mean-vel)**2,axis=1),weights=m))),'mass_change_L1_relative':float(abs(q[:,0]-m).sum()/m.sum()),'minimum_U':float(q[:,6].min()),'seconds':time.monotonic()-start}
 rows.append(row);np.savez(out/f'cycles{cycles}.npz',ids=ids,state=q);print(json.dumps(row),flush=True)
result={'scope':__doc__,'source_time':t,'rows':rows,'initial_thermal_sum':float(U.sum()),'reference':'Retaining original donor subcells without projection restores exactly the original field when geometry returns.','limitations':['Uniform spatial lift, unlike previous fitted-density lift; preserves m,p,Pi,U but not original point angular momentum or potential energy.','Projection is piecewise constant, not affine or production solver.','Total variation fixed but number of reversals changes; not time convergence along one fixed smooth weight path.','Internal kinetic energy is retained as covariance, not thermalized.','No pressure, gravity, physical evolution, backend parity or MPI claim.']}
(out/'result.json').write_text(json.dumps(result,indent=2)+'\n')
