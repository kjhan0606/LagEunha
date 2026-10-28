"""Cold constant-pressure density interface under one finite remap on perturbed grids."""
from pathlib import Path
import ast,importlib.util,json,time
import numpy as np
base=Path('/gpfs/kjhan/LagForce');out=base/'Codex_cold_interface_remap_20261018';out.mkdir(exist_ok=False)
def load(n,p):
 s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
for filename,name in [('audit_quadratic_affine_path.py','integrate'),('audit_remap_true_spatial.py','recover')]:
 tree=ast.parse((base/filename).read_text());nodes=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name==name];exec(compile(ast.Module(body=nodes,type_ignores=[]),filename,'exec'),globals())
mod=load('mom',base/'Codex_weight_remap_moments_centered_20261017/audit.py');geo=load('geo',base/'Codex_weighted_compression_geometry_20261017/audit.py');grav=load('grav',base/'audit_canonical_pair_history.py');rows=[]
for N in (16,32,64):
 start=time.monotonic();h=2/N;I,J=np.meshgrid(np.arange(N),np.arange(N));rng=np.random.default_rng(24680);x=np.column_stack((-1+(I.ravel()+.5)*h,-1+(J.ravel()+.5)*h))+rng.uniform(-.15,.15,(N*N,2))*h;ids=np.arange(len(x));r=np.linalg.norm(x,axis=1);vel=np.sqrt(grav.gravity_coefficient(r))[:,None]*np.column_stack((-x[:,1],x[:,0]));w=.03*h*h*np.sin(1.7*I.ravel()+.9*J.ravel())*np.maximum(1-(r/.8)**2,0)**2
 density=np.where(r<.5,.01,1.);old=mod.polygons(x,np.zeros(len(x)),ids,geo.clip);new=mod.polygons(x,w,ids,geo.clip);vol=np.array([mod.area(p) for p in old]);nv=np.array([mod.area(p) for p in new]);lo=np.array([p.min(0) for p in old]);hi=np.array([p.max(0) for p in old]);q0=[]
 for i,p in enumerate(old):
  xy,dv=mod.quadrature(p,2);q0.append(integrate(xy,dv,np.full(len(xy),density[i]),np.broadcast_to(vel[i],xy.shape),np.zeros((2,2)),1e-5))
 q0=np.array(q0);q=np.zeros_like(q0);cover=np.zeros(len(ids));newcover=np.zeros(len(ids));lowcover=np.zeros(len(ids));highcover=np.zeros(len(ids))
 for j,p in enumerate(new):
  cand=np.flatnonzero(np.all(hi>=p.min(0),axis=1)&np.all(lo<=p.max(0),axis=1))
  for i in cand:
   poly=old[i].copy()
   for a,b in zip(p,np.roll(p,-1,axis=0)):
    e=b-a;n=np.array([e[1],-e[0]]);poly=geo.clip(poly,n,n@a)
    if len(poly)<3:break
   if mod.area(poly)<=1e-18:continue
   xy,dv=mod.quadrature(poly,2);q[j]+=integrate(xy,dv,np.full(len(xy),density[i]),np.broadcast_to(vel[i],xy.shape),np.zeros((2,2)),1e-5);cover[i]+=dv.sum();newcover[j]+=dv.sum();lowcover[j]+=dv.sum() if density[i]<.5 else 0.;highcover[j]+=dv.sum() if density[i]>.5 else 0.
 err=max(float(abs(cover/vol-1).max()),float(abs(newcover/nv-1).max()));assert err<1e-8
 S,C,G,R=recover(q);mass=q[:,0];c=q[:,1:3]/mass[:,None];v=q[:,6:8]/mass[:,None];vg=v+np.einsum('nij,nj->ni',G,x-c);internal=.5*(q[:,12]+q[:,14])-.5*mass*np.sum(v*v,axis=1);res=.5*np.trace(R,axis1=1,axis2=2);regions={}
 for name,mask in [('inner_band',(r>.2)&(r<.5)),('outer_band',(r>=.5)&(r<.7))]:
  M=mass[mask].sum();regions[name]={'cells':int(mask.sum()),'internal_K_per_mass':float(internal[mask].sum()/M),'residual_K_per_mass':float(res[mask].sum()/M),'mean_velocity_change_RMS':float(np.sqrt(np.average(np.sum((v-vel)[mask]**2,axis=1),weights=mass[mask]))),'generator_velocity_change_RMS':float(np.sqrt(np.average(np.sum((vg-vel)[mask]**2,axis=1),weights=mass[mask])))}
 U=q[:,15]-.5*(q[:,12]+q[:,14]);expectedU=1e-5*nv;interface=(lowcover>1e-12*nv)&(highcover>1e-12*nv);pressure_metrics={}
 for label,mask in [('all_interior',r<.7),('mixed_density_cells',interface),('inner',(r>.2)&(r<.5))]:
  ratio=internal[mask]/expectedU[mask];pressure_metrics[label]={'cells':int(mask.sum()),'retained_covariance_pressure_relative_error_max':float(abs(U[mask]/expectedU[mask]-1).max()),'thermalized_pressure_relative_excess_max':float(ratio.max()),'thermalized_pressure_relative_excess_volume_RMS':float(np.sqrt(np.average(ratio**2,weights=nv[mask]))),'internal_K_over_thermal_sum':float(internal[mask].sum()/expectedU[mask].sum())}
 row={'pressure_metrics':pressure_metrics,'N':N,'h':h,'regions':regions,'coverage_error':err,'global_moment_drift_max':float(abs(q.sum(0)-q0.sum(0)).max()),'seconds':time.monotonic()-start};rows.append(row);print(json.dumps(row),flush=True)
for prev,row in zip(rows,rows[1:]):
 row['refinement_exponents']={region:{key:float(np.log(prev['regions'][region][key]/row['regions'][region][key])/np.log(2)) for key in ('internal_K_per_mass','residual_K_per_mass','generator_velocity_change_RMS')} for region in row['regions']}
(out/'result.json').write_text(json.dumps({'scope':__doc__,'rows':rows,'limits':['Perturbed grids with generator-based rho=.01 inside r=.5 and1 outside; interface follows donor cells, not exact circular cut cells. Constant thermal energy density1e-5.','Single remap, not fixed full path or physical orbit evolution.','W amplitude scales h^2; same dimensionless displacement family, not fixed absolute W field.','Fixed random seed, nonnested grids; three resolutions do not establish general convergence order.','Pointwise constant donor velocities, no affine energy added at initialization.']},indent=2)+'\n')
