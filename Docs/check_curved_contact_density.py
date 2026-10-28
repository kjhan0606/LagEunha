"""Curved-domain m/q/Q/Psi and trace fixtures, no fluid provider."""
import argparse,hashlib,json,sys,time
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce');modules=base/'Codex_full_disk_density_modules_20261019';sys.path[:0]=[str(modules),str(base/'Codex_compiled_midpoint_modules_20261018'),str(base/'Codex_coupled_relaxation_evolution_20261018')]
from curved_disk_cell import DiskCell
from curved_contact_density import CurvedContactPriorDensity
from contact_potential_density import ContactPotentialDensity,features
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();args.output.mkdir(exist_ok=False)
cell=DiskCell([[-.1,-.1],[.1,-.1],[.1,.1],[-.1,.1]],[-1,-2,-3,-4],.1);ld=np.longdouble

def potential(x):
 x=np.asarray(x,dtype=ld);return x[:,0]**4+ld('.7')*x[:,1]**4+ld('.2')*x[:,0]**2*x[:,1]**2
rows={}
for label,meta in [('flat',dict(flat=True,value=1.)),('cut',dict(flat=False,normal=[1.,0.],offset=0.,low=.7,high=1.3))]:
 pieces=[(cell,1.)] if meta['flat'] else [(cell.clip_halfplane([1.,0.],0.),.7),(cell.clip_halfplane([-1.,0.],0.),1.3)]
 points=[];weights=[];prior=[]
 for p,rho in pieces:
  x,w=p.positive_quadrature(48);points.append(x);weights.append(w);prior.append(np.full(len(w),rho))
 x=np.concatenate(points);w=np.concatenate(weights);rho0=np.concatenate(prior);rho=rho0*np.exp(.2*(x[:,0]/.1)**4+.1*(x[:,1]/.1)**4)
 dm=w*rho;m=float(dm.sum());q=np.asarray(dm.astype(ld)@x.astype(ld),dtype=float);Q=np.asarray(np.einsum('n,ni,nj->ij',dm.astype(ld),x.astype(ld),x.astype(ld)),dtype=float);Psi=dm.astype(ld)@potential(x)
 old=CurvedContactPriorDensity(cell,m,q,Q,meta,12,())
 fit=ContactPotentialDensity(cell,m,q,Q,meta,old,potential,Psi,np.zeros(6,dtype=ld),24,())
 pts,ww,pp=fit.quadrature(48);mass=ww*fit._value_with_prior(pts,pp);A=features(pts,fit.origin,fit.scale)
 error=mass.astype(ld)@A-ld(m)*fit.target_A;psierror=mass.astype(ld)@potential(pts)-Psi
 assert abs(error).max()<3e-13 and abs(psierror)<1e-14 and fit.theta.shape==(6,)
 test=np.array([[0.,.03],[0.,-.02],[.04,.01],[-.04,.02]]);directions=np.array([[1.,0.],[-1.,0.],[0.,1.],[0.,1.]])
 values,branches=fit.value_trace_batch(test,directions,True);scalar=np.array([fit.value_trace(p,d)[0] for p,d in zip(test,directions)])
 assert np.all(abs(values-scalar)<=64*np.finfo(float).eps*abs(scalar))
 assert branches.tolist()==([0,0,0,0] if meta['flat'] else [1,-1,1,-1])
 rows[label]=dict(status='passed',independent_q48_moment_error=float(abs(error).max()),independent_q48_Psi_error=float(psierror),dual_residual=fit.residual,condition=fit.condition,sigma=float(fit.sigma),rank_error=float(fit.rank_error),density_min=float(mass.min()),scalar_batch_max_relative=float(abs(values/scalar-1).max()),prior_region_areas=[p.area for p in fit.parts],area=cell.area)
report=dict(status='passed',scope=__doc__,results=rows,source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),modules/'contact_potential_density.py',modules/'curved_contact_density.py',modules/'curved_disk_cell.py']},limits='Curved flat/cut density and trace fixtures; full-domain actual stage and long-time stability remain unverified.')
(args.output/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
