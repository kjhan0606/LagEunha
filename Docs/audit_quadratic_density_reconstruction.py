"""Match actual transported spatial moments with quadratic density; continuous positivity gate."""
from pathlib import Path
import importlib.util,json
import numpy as np

def load(n,p):
 s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
base=Path('/gpfs/kjhan/LagForce');out=base/'Codex_quadratic_density_reconstruction_20261018';out.mkdir(exist_ok=False)
mod=load('mom',base/'Codex_weight_remap_moments_centered_20261017/audit.py');geo=load('geo',base/'Codex_weighted_compression_geometry_20261017/audit.py');rd=load('rd',base/'run_canonical_av_branches_20261018.py')
a,t=rd.read(base/'Codex_HG_constrained_20261017/omp8_t20/khout.000184.dat');x=np.column_stack((a['x']-3,a['y']-3));target=np.load(base/'Codex_evolving_weight_transport_20261017/weights_0.npz')['weights']
def feature(z):return np.column_stack((np.ones(len(z)),z[:,0],z[:,1],z[:,0]**2,z[:,0]*z[:,1],z[:,1]**2))
def minimum(poly,b):
 candidates=list(poly);H=np.array([[2*b[3],b[4]],[b[4],2*b[5]]]);g=b[1:3]
 for aa,bb in zip(poly,np.roll(poly,-1,axis=0)):
  d=bb-aa;den=d@H@d
  if den>0:
   u=-d@(g+H@aa)/den
   if 0<u<1:candidates.append(aa+u*d)
 if np.linalg.eigvalsh(H).min()>0:
  z=np.linalg.solve(H,-g);e=np.roll(poly,-1,axis=0)-poly;rel=z-poly
  if np.all(e[:,0]*rel[:,1]-e[:,1]*rel[:,0]>=-1e-12):candidates.append(z)
 return float(np.min(feature(np.array(candidates))@b))
rows=[]
for n in (1,2,4,8):
 saved=np.load(base/'Codex_remap_true_spatial_20261018'/f'subdivisions{n}.npz');ids=saved['ids'];q=saved['state'];polys=mod.polygons(x,.25/n*target,ids,geo.clip);fail=[];minratio=float('inf');maxerr=0.;maxcond=0.;coeff=[]
 for k,p in enumerate(polys):
  xy,dv=mod.quadrature(p,4);V=dv.sum();origin=x[ids[k]];h=np.sqrt(V);z=(xy-origin)/h;F=feature(z);gram=F.T@(dv[:,None]*F)/V;mass=q[k,0];c=q[k,1:3]/mass;Q=np.array([[q[k,3],q[k,4]],[q[k,4],q[k,5]]])/mass
  mean=(c-origin)/h;second=(Q-np.outer(origin,c)-np.outer(c,origin)+np.outer(origin,origin))/h**2
  rhs=np.array([1,*mean,second[0,0],second[0,1],second[1,1]]);b=np.linalg.solve(gram,rhs);res=float(abs(gram@b-rhs).max());assert res<1e-10
  low=minimum((p-origin)/h,b);minratio=min(minratio,low);maxerr=max(maxerr,res);maxcond=max(maxcond,float(np.linalg.cond(gram)));coeff.append(b)
  if low < -1e-10:fail.append({'id':int(ids[k]),'radius':float(np.linalg.norm(origin)),'minimum_density_over_mean':low})
 row={'subdivisions':n,'fraction':.25/n,'cells':len(ids),'negative_density_cells':len(fail),'minimum_density_over_mean':minratio,'max_normalized_moment_fit_residual':maxerr,'max_gram_condition':maxcond,'failures':fail};rows.append(row);np.savez(out/f'fit{n}.npz',ids=ids,coefficients=coeff)
 print(json.dumps({k:v for k,v in row.items() if k!='failures'}),flush=True)
result={'scope':__doc__,'rows':rows,'limits':['Signed polynomial fit with exact moment equations, not a positivity-constrained fit.','Quadratic minimum checked on vertices, edge stationary points and interior positive-definite stationary point.','No velocity/energy reconstruction, repeated remap, physical evolution or production integration.']};(out/'result.json').write_text(json.dumps(result,indent=2)+'\n')
