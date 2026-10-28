"""Fixed four-state saved midpoint density derivative, not time evolution.

Only x/W/q/Q follow the saved actual stage direction. Mass and frozen rows
remain exact, factors/velocities/Pi remain untouched. No RHS or weight solve.
"""
import argparse,hashlib,importlib.util,json,sys,time
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce')
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
sys.path[:0]=[str(base/'Codex_compiled_midpoint_modules_20261018'),str(base/'Codex_moment_flux_batch_modules_20261018'),str(base/'Codex_coupled_relaxation_evolution_20261018'),str(base)]
from saved_potential_taylor import raw_p2,integrate_remainder
from contact_prior_density import ContactPriorDensity,local_contact_priors
from embedded_laguerre_geometry import EmbeddedLaguerreGeometry
from run_canonical_av_branches_20261018 import read
def load(name,path):
 spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
root=base/'Codex_compiled_midpoint_once_20261018';polynomial=base/'Codex_potential_polynomial_20261018'
data=np.load(root/'midpoint/inputs.npz');ref=np.load(root/'midpoint/reference.npz')
mass=data['mass'];frozen=data['frozen'].astype(bool);active=~frozen;ids=data['ids'];y0=data['state'];a0=y0[:-4].reshape(-1,20)
h=float(json.loads((root/'manifest.json').read_text())['dt']);epsilon=h/8
coeff=dict(np.load(polynomial/'coefficients.npz'));coeff.pop('ids')
base_integrals=np.load(polynomial/'midpoint_integrals.npz');basePsi=base_integrals['raw_P2']+base_integrals['integrals'][1,:,0]
clip_path=base/'Codex_weighted_compression_geometry_20261017/audit.py';potential_path=base/'Codex_weight_remap_moments_centered_20261017/audit.py'
clip=load('saved_clip',clip_path);potential=load('saved_potential',potential_path).potential
ic,t0=read(base/'Codex_canonical_initial_20261017/run/khout.000000.dat');assert t0==0.
builder=EmbeddedLaguerreGeometry(np.column_stack([ic['x']-3,ic['y']-3]),ids,clip.clip)
raw0=np.column_stack([mass,a0[:,3:]]);polys0,edges0=builder(a0[:,:2],a0[:,2],float(data['time']))
prior0=local_contact_priors(polys0,raw0,edges0)
def graph(edges):return {tuple(sorted((int(i),int(j)))) for i,j,*_ in edges if j>=0}
graph0=graph(edges0)
def branch(meta):return (meta['flat'],meta.get('reason'),meta.get('direction'),meta['stencil_cells'])
compiled_path=root/'midpoint/cuda.npz'
compiled_stage=np.load(compiled_path)['stage']
direction=np.column_stack([data['generator'],data['weight_rates'],compiled_stage[:,1:7]])
direction[frozen]=0.
args.output.mkdir(exist_ok=False)
manifest=dict(scope=__doc__,epsilon=epsilon,scales=[1.,-1.,.5,-.5],provider_calls=0,expected_newton_fits=4*len(ids),
    sources={str(p):sha(p) for p in (Path(__file__),clip_path,potential_path,polynomial/'coefficients.npz',polynomial/'midpoint_integrals.npz',root/'midpoint/inputs.npz',root/'midpoint/reference.npz',compiled_path)})
manifest['runtime_modules']={name:dict(path=str(m.__file__),sha256=sha(m.__file__)) for name,m in list(sys.modules.items()) if getattr(m,'__file__',None) and Path(m.__file__).suffix=='.py' and Path(m.__file__).is_relative_to(base)}
(args.output/'manifest.json').write_text(json.dumps(manifest,indent=2));results={};values={};errors={};fit_effect={};sum_error={};all_same=True
for label,scale in (('plus',1.),('minus',-1.),('plus_half',.5),('minus_half',-.5)):
 started=time.perf_counter();y=y0.copy();a=y[:-4].reshape(-1,20);a[:,:9]+=scale*epsilon*direction
 assert np.array_equal(a[frozen],a0[frozen]) and np.array_equal(a[:,9:],a0[:,9:])
 raw=np.column_stack([mass,a[:,3:]]);polys,edges=builder(a[:,:2],a[:,2],float(data['time']))
 priors=local_contact_priors(polys,raw,edges);densities=[ContactPriorDensity(p,r[0],r[1:3],r[3:7].reshape(2,2),m,12) for p,r,m in zip(polys,raw,priors)]
 offsets=np.r_[0,np.cumsum([len(p) for p in polys])]
 np.savez(args.output/(label+'_reconstruction.npz'),state=y,mass=mass,ids=ids,polygon_offsets=offsets,polygon_vertices=np.concatenate(polys))
 (args.output/(label+'_priors.json')).write_text(json.dumps(priors,indent=2))
 (args.output/(label+'_fits.json')).write_text(json.dumps([dict(theta=d.theta.tolist(),iterations=d.iterations,feature_residual=d.residual) for d in densities],indent=2))
 integral=integrate_remainder(densities,polys,coeff,potential);poly=raw_p2(raw,coeff);Psi=poly+integral[1,:,0]
 values[label]=Psi;errors[label]=abs(integral[1,:,0]-integral[0,:,0])
 fit_effect[label]=abs(integral[1,:,1]-poly)
 count=np.array([sum(len(part)*48**2 for part in d.parts) for d in densities])
 gamma=count*np.finfo(np.longdouble).eps/(1-count*np.finfo(np.longdouble).eps)
 sum_error[label]=gamma*integral[1,:,2]
 graphnew=graph(edges);changed=[k for k,(m,n) in enumerate(zip(prior0,priors)) if branch(m)!=branch(n)]
 same=graphnew==graph0 and not changed;all_same &= same
 np.savez(args.output/(label+'_integrals.npz'),raw_P2=poly,integrals=integral,Psi=Psi,
    P2_fit_constraint_error=integral[1,:,1]-poly)
 results[label]=dict(seconds=time.perf_counter()-started,newton_fits=len(densities),fit_residual_max=max(d.residual for d in densities),
    graph_added=len(graphnew-graph0),graph_removed=len(graph0-graphnew),prior_branch_changed=len(changed),prior_branch_changed_ids=ids[changed].tolist(),
    active_prior_branch_changed=int(active[changed].sum()),same_graph_and_prior_branch=same,
    quadrature_remainder_difference_absolute_sum=float(errors[label][active].sum()),
    P2_fit_constraint_error_absolute_sum=float(fit_effect[label][active].sum()),
    remainder_summation_roundoff_bound=float(sum_error[label][active].sum()),
    perturbation_roundoff_max=float(abs(a[:,:9]-a0[:,:9]-scale*epsilon*direction).max()))
 print(json.dumps(dict(label=label,report=results[label])),flush=True)

centered=(values['plus']-values['minus'])/(2*epsilon)
centered_half=(values['plus_half']-values['minus_half'])/epsilon
forward=(values['plus']-basePsi)/epsilon;backward=(basePsi-values['minus'])/epsilon
forward_half=(values['plus_half']-basePsi)/(epsilon/2);backward_half=(basePsi-values['minus_half'])/(epsilon/2)
transport_data=np.load(polynomial/'midpoint_transport.npz')
transport=transport_data['face_P2']+transport_data['volume_P2_raw']+transport_data['face_remainder']+transport_data['volume_remainder']
quad_derivative=(errors['plus_half']+errors['minus_half'])/epsilon
fit_derivative=(fit_effect['plus_half']+fit_effect['minus_half'])/epsilon
sum_derivative=(sum_error['plus_half']+sum_error['minus_half'])/epsilon
def summarize(v):
 v=v[active];return dict(signed_sum=float(v.sum()),absolute_sum=float(abs(v).sum()),maximum=float(abs(v).max(initial=0)))
derivatives={name:summarize(value) for name,value in (('centered',centered),('centered_half',centered_half),('forward',forward),('backward',backward),('forward_half',forward_half),('backward_half',backward_half),('scale_difference',centered_half-centered),('candidate_closure_defect',centered_half-transport),('quadrature_derivative_difference',quad_derivative),('P2_constraint_derivative_effect',fit_derivative),('remainder_summation_derivative_bound',sum_derivative))}
np.savez(args.output/'derivatives.npz',ids=ids,active=active,centered=centered,centered_half=centered_half,
    forward=forward,backward=backward,forward_half=forward_half,backward_half=backward_half,
    transport=transport,quad_derivative=quad_derivative,fit_derivative=fit_derivative,sum_derivative=sum_derivative)
report=dict(status='completed_derivative_requires_evaluation',manifest=manifest,states=results,
    all_same_graph_and_prior_branch=all_same,derivatives=derivatives,
    interpretation='Candidate closure derivative is unconfirmed if graph/prior branches change or two scales are unstable. No automatic epsilon sweep. P2 raw evaluation uses longdouble; remainder directly integrated with original binary64 density. Fit residual and P2 quadrature constraint defects are separately retained, not hidden by thermal correction.')
(args.output/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(derivatives),flush=True)
