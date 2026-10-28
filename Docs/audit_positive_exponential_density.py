"""Positive exponential density fitting on saved rejected cells only."""
from pathlib import Path
import ast,importlib.util,json,time
import numpy as np
from scipy.optimize import least_squares
from scipy.special import logsumexp
base=Path('/gpfs/kjhan/LagForce');out=base/'Codex_positive_exponential_density_20261018';out.mkdir(exist_ok=False)
def load(n,p):
 s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
for filename,names in [('audit_quadratic_density_reconstruction.py',('feature','minimum'))]:
 tree=ast.parse((base/filename).read_text());functions=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in names];exec(compile(ast.Module(body=functions,type_ignores=[]),filename,'exec'),globals())
mod=load('mom',base/'Codex_weight_remap_moments_centered_20261017/audit.py');geo=load('geo',base/'Codex_weighted_compression_geometry_20261017/audit.py');rd=load('rd',base/'run_canonical_av_branches_20261018.py')
a,t=rd.read(base/'Codex_HG_constrained_20261017/omp8_t20/khout.000184.dat');x=np.column_stack((a['x']-3,a['y']-3));target=np.load(base/'Codex_evolving_weight_transport_20261017/weights_0.npz')['weights'];source=base/'Codex_quadratic_affine_path_20261018';failed=json.loads((source/'result.json').read_text());rows=[]
for run in failed['rows']:
 steps=run['subdivisions'];index=run['history'][-1]['state_index'];fraction=index*.25/steps;saved=np.load(source/f'steps{steps}.npz');ids=saved['ids'];state=saved['state'];polys=mod.polygons(x,fraction*target,ids,geo.clip)
 for f in run['history'][-1]['failures']:
  k=int(np.flatnonzero(ids==f['id'])[0]);p=polys[k];q=state[k];mass=q[0];V=mod.area(p);h=np.sqrt(V);origin=x[ids[k]];c=q[1:3]/mass;Q=np.array([[q[3],q[4]],[q[4],q[5]]])/mass;mean=(c-origin)/h;second=(Q-np.outer(origin,c)-np.outer(c,origin)+np.outer(origin,origin))/h**2;desired=np.r_[mean,second[0,0],second[0,1],second[1,1]]
  def grid(order):
   xy,dv=mod.quadrature(p,order);return feature((xy-origin)/h)[:,1:],np.log(dv/V)
  def stats(theta,F,logw):
   lw=logw+F@theta;Z=logsumexp(lw);prob=np.exp(lw-Z);mu=prob@F;delta=F-mu;cov=(delta.T*prob)@delta
   return mu,cov,Z
  theta=np.zeros(5);fits=[];started=time.monotonic()
  for order,check in ((8,16),(16,32),(32,48)):
   F,lw=grid(order)
   sol=least_squares(lambda b:stats(b,F,lw)[0]-desired,theta,jac=lambda b:stats(b,F,lw)[1],xtol=1e-13,ftol=1e-13,gtol=1e-13,max_nfev=200)
   theta=sol.x;mu,cov,Z=stats(theta,F,lw);Fc,lwc=grid(check);muc,_,Zc=stats(theta,Fc,lwc)
   # Check same normalization too, so the independent mass error is visible.
   ratio=np.exp(Zc-Z);rawcheck=np.r_[ratio,ratio*muc];ref=np.r_[1.,desired]
   lower=minimum((p-origin)/h,np.r_[0.,theta])-Z
   fits.append({'fit_order':order,'check_order':check,'solver_success':bool(sol.success),'nfev':sol.nfev,'fit_residual_max':float(abs(mu-desired).max()),'independent_mass_and_moment_residual_max':float(abs(rawcheck-ref).max()),'minimum_density_over_mean':float(np.exp(lower)),'log_normalizer':float(Z),'theta':theta.tolist()})
  row={'subdivisions':steps,'weight_fraction':fraction,'id':int(ids[k]),'fits':fits,'seconds':time.monotonic()-started};rows.append(row);print(json.dumps({'N':steps,'id':int(ids[k]),'checks':[v['independent_mass_and_moment_residual_max'] for v in fits],'minrho':fits[-1]['minimum_density_over_mean'],'seconds':row['seconds']}),flush=True)
(out/'result.json').write_text(json.dumps({'scope':__doc__,'rows':rows,'limits':['Five rejected cell states only, not the whole path or mesh.','Strict positivity from exponential form; normalization/moments verified numerically with independent higher-order quadrature.','No boundary transport, reconstruction uniqueness, pressure closure, physical evolution or production cost claim.']},indent=2)+'\n')
