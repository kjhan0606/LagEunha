"""Synthetic fixtures for the offline potential-constrained density API."""
import argparse,hashlib,json,sys
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce');modules=base/'Codex_contact_potential_modules_20261019'
sys.path[:0]=[str(modules),str(base/'Codex_compiled_midpoint_modules_20261018'),str(base/'Codex_coupled_relaxation_evolution_20261018')]
from contact_potential_density import ContactPotentialDensity,ReconstructionFailure
from contact_prior_density import ContactPriorDensity
from positive_radial_quadrature import rule
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();args.output.mkdir(exist_ok=False)
poly=np.array([[-.1,-.1],[.1,-.1],[.1,.1],[-.1,.1]])
points,w=rule(poly,48);ld=np.longdouble
known=np.exp(.4*(points[:,0]/.1)**4+.2*(points[:,1]/.1)**4);dm=w*known;dm/=dm.sum();m=1.
q=dm@points;Q=points.T@(dm[:,None]*points);prior=dict(flat=True,value=1.)
old=ContactPriorDensity(poly,m,q,Q,prior,12)
def polynomial(c,b,H):
 c=ld(c);b=np.asarray(b,dtype=ld);H=np.asarray(H,dtype=ld)
 def value(x):
  x=np.asarray(x,dtype=ld);return c+x@b+.5*np.einsum('ni,ij,nj->n',x,H,x)
 o=old.origin.astype(ld);s=ld(old.scale)
 reference=np.r_[value(o[None,:])[0],s*(b+H@o),.5*s*s*H[0,0],s*s*H[0,1],.5*s*s*H[1,1]]
 target=c*m+b@q.astype(ld)+.5*np.sum(H*Q.astype(ld))
 return value,reference,target
results={};passed=True
for label,c,b,H in (('constant',2.,[0.,0.],[[0.,0.],[0.,0.]]),('linear',.3,[.7,-.2],[[0.,0.],[0.,0.]]),('quadratic',.3,[.7,-.2],[[.4,.2],[.2,-.3]])):
 fun,ref,target=polynomial(c,b,H)
 try:
  fit=ContactPotentialDensity(poly,m,q,Q,prior,old,fun,target,ref,polynomial=True)
  assert fit.polynomial_dependent
  try:ContactPotentialDensity(poly,m,q,Q,prior,old,fun,target+ld('0.01'),ref,polynomial=True)
  except ReconstructionFailure as exc:assert exc.status=='inconsistent_dependent_potential_target'
  else:raise AssertionError('dependent target accepted')
  results[label]=dict(status='passed',residual=fit.residual)
 except Exception as exc:results[label]=dict(status='failed',error=repr(exc));passed=False
def nonpoly(x):
 x=np.asarray(x,dtype=ld);return x[:,0]**4+ld('.7')*x[:,1]**4+ld('.2')*x[:,0]**2*x[:,1]**2
target=np.sum(dm.astype(ld)*nonpoly(points));ref=np.zeros(6,dtype=ld)
try:
 fit=ContactPotentialDensity(poly,m,q,Q,prior,old,nonpoly,target,ref)
 pts,weights,rho0=fit.quadrature(48);rho=fit._value_with_prior(pts,rho0);wm=weights*rho
 target_error=float(np.sum(wm.astype(ld)*nonpoly(pts))-target)
 moment_error=max(abs(wm.sum()-m),abs(wm@pts-q).max(),abs(pts.T@(wm[:,None]*pts)-Q).max())
 assert abs(target_error)<=1e-14+1e-10*abs(float(target)) and moment_error<3e-13
 gauge=ld('.125');gref=ref.copy();gref[0]+=gauge
 shifted=ContactPotentialDensity(poly,m,q,Q,prior,old,lambda x:nonpoly(x)+gauge,target+gauge*m,gref)
 shiftedrho=shifted._value_with_prior(pts,rho0);difference=abs(shiftedrho-rho)
 assert np.all(difference<=1e-14+1e-10*abs(rho))
 results['nonpolynomial_and_gauge']=dict(status='passed',residual=fit.residual,independent_q48_Psi_error=target_error,independent_q48_moment_error=float(moment_error),sigma=float(fit.sigma),rank_error=float(fit.rank_error),condition=fit.condition,gauge_density_max_difference=float(difference.max()))
 bad_target=target+ld(100)*fit.sigma
 try:ContactPotentialDensity(poly,m,q,Q,prior,old,nonpoly,bad_target,ref)
 except ReconstructionFailure as exc:
  assert exc.status=='outside_quadrature_moment_hull';results['infeasible_target']=dict(status='passed',failure=exc.status,witness=exc.details)
 else:raise AssertionError('infeasible target accepted')
except Exception as exc:results['nonpolynomial_and_gauge']=dict(status='failed',error=repr(exc));passed=False
try:
 ContactPotentialDensity(poly,m,q,Q,prior,old,lambda x:np.ones(len(x),dtype=ld)+ld('1e-30')*nonpoly(x),ld(m),np.array([1.,0.,0.,0.,0.,0.],dtype=ld))
except ReconstructionFailure as exc:
 if exc.status!='unresolved_potential_rank':passed=False
 results['unresolved_rank']=dict(status='passed' if exc.status=='unresolved_potential_rank' else 'failed',failure=exc.status)
else:results['unresolved_rank']=dict(status='failed');passed=False
report=dict(status='passed' if passed else 'fixture_failed',results=results,source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),modules/'contact_potential_density.py',modules/'positive_radial_quadrature.py')},limits='Synthetic numerical realizability fixtures; no dynamics or proof of continuous moment hull membership.')
(args.output/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2));assert passed
