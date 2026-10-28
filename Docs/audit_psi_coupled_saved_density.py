"""q24 versus saved independent q48, and saved six-feature trace API audit.

Restore parameters only; zero Newton fits, providers, or evolution steps.
"""
import argparse,hashlib,importlib.util,json,sys,time
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce');modules=base/'Codex_psi_coupled_modules_20261019'
sys.path[:0]=[str(modules),str(base/'Codex_compiled_midpoint_modules_20261018'),str(base/'Codex_quadratic_thermal_audit_modules_20261018'),str(base/'Codex_quadratic_split_modules_20261018'),str(base/'Codex_coupled_relaxation_evolution_20261018'),str(base)]
from contact_potential_density import ContactPotentialDensity
from contact_prior_density import ContactPriorDensity,clip_halfplane,polygon_area
from positive_radial_quadrature import rule
from embedded_laguerre_geometry import EmbeddedLaguerreGeometry
from run_canonical_av_branches_20261018 import read
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--root',type=Path,required=True);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();args.output.mkdir(exist_ok=False);started=time.perf_counter()
def load(name,path):
 spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
potential=load('potential',base/'Codex_weight_remap_moments_centered_20261017/audit.py').potential
clip=load('clip',base/'Codex_weighted_compression_geometry_20261017/audit.py').clip
initial=np.load(args.root/'initial.npz');mass=initial['mass'];ids=initial['ids'];active=initial['psi_constrained_mask'];frozen=~active
ic,t=read(base/'Codex_canonical_initial_20261017/run/khout.000000.dat');builder=EmbeddedLaguerreGeometry(np.column_stack([ic['x']-3,ic['y']-3]),ids,clip)
report={};trace_count=0;eps=np.finfo(float).eps;ld=np.longdouble
for label,checkpoint in [('initial',args.root/'initial.npz'),('half',args.root/'half_update/candidate.npz'),('final',args.root/'accepted.npz')]:
 saved=np.load(checkpoint);a=saved['state'][:-4].reshape(-1,20);raw=np.column_stack([mass,a[:,3:]]);polygons,edges=builder(a[:,:2],a[:,2],0.)
 folder=args.root/(label+'_density');priors=json.loads((folder/'priors.json').read_text());params=[json.loads(line) for line in (folder/'parameters.jsonl').read_text().splitlines()];q48=np.load(args.root/(label+'_q48.npz'))['Psi']
 q24=[];abssums=[];maximum_trace_relative=0.;max_feature_omission_effect=0.
 for i,(p,r,meta,parameter) in enumerate(zip(polygons,raw,priors,params)):
  d=(ContactPotentialDensity if active[i] else ContactPriorDensity).__new__(ContactPotentialDensity if active[i] else ContactPriorDensity)
  d.mass=float(mass[i]);d.origin=np.array(parameter['origin']);d.scale=parameter['scale'];d.theta=np.array(parameter['theta']);d.logZ=parameter['logZ'];d.cut_metadata=meta
  if meta['flat']:parts=[(p,meta['value'])]
  else:
   n=np.array(meta['normal']);off=meta['offset'];parts=[(clip_halfplane(p,-n,-off),meta['low']),(clip_halfplane(p,n,off),meta['high'])]
  d.parts=[piece for piece,rho in parts if polygon_area(piece)>0];d.part_priors=[rho for piece,rho in parts if polygon_area(piece)>0]
  if active[i]:
   d.potential=potential;d.radii=(.25,.5);d.polynomial_dependent=False;d.reference=np.array(parameter['reference'],dtype=ld);d.remainder_projection=np.array(parameter['projection_remainder'],dtype=ld);d.sigma=ld(parameter['sigma'])
  val=ld(0);abssum=ld(0)
  for piece,rho in zip(d.parts,d.part_priors):
   pts,w=rule(piece,24,(.25,.5));dm=(w*d._value_with_prior(pts,np.full(len(w),rho))).astype(ld);Phi=potential(pts.astype(ld));val+=dm@Phi;abssum+=dm@abs(Phi)
  q24.append(val);abssums.append(abssum)
  if active[i]:
   points=(p+np.roll(p,-1,axis=0))/2;direction=np.mean(p,axis=0)-points
   batch,branches=d.value_trace_batch(points,direction,True)
   scalar=np.array([d.value_trace(point,inward)[0] for point,inward in zip(points,direction)])
   G=d._evaluation_features(points)
   bound=128*eps*(1+abs(G)@abs(d.theta)+abs(d.logZ))*abs(scalar)
   assert np.all(abs(batch-scalar)<=bound)
   maximum_trace_relative=max(maximum_trace_relative,float(np.max(abs(batch/scalar-1))))
   # The sixth coefficient is present in the actual evaluation exponent.
   without_six=batch*np.exp(-G[:,5]*d.theta[5])
   max_feature_omission_effect=max(max_feature_omission_effect,float(np.max(abs(without_six/batch-1))))
   trace_count+=len(points)
 q24=np.array(q24,dtype=ld);delta=q24-q48
 np.savez(args.output/(label+'.npz'),Psi_q24=q24,Psi_q48=q48,difference=delta,absolute_integrand_sum=np.array(abssums),active=active,ids=ids)
 report[label]=dict(active_q24_minus_q48_signed=float(delta[active].sum()),active_q24_minus_q48_absolute_sum=float(abs(delta[active]).sum()),active_q24_minus_q48_max=float(abs(delta[active]).max()),frozen_q24_minus_q48_signed=float(delta[frozen].sum()),maximum_scalar_batch_trace_relative_difference=maximum_trace_relative,maximum_omitting_sixth_feature_relative_effect=max_feature_omission_effect)
 print(json.dumps(dict(label=label,**report[label])),flush=True)
report=dict(status='passed',scope=__doc__,states=report,active_trace_points=trace_count,newton_fits=0,provider_calls=0,seconds=time.perf_counter()-started,source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),modules/'contact_potential_density.py',args.root/'manifest.json']},limits='q24/q48 differences are empirical quadrature estimates, not interval bounds; finite sample trace agreement does not prove continuous realization or physical accuracy.')
(args.output/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(dict(status='passed',seconds=report['seconds'])))
