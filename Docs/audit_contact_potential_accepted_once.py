"""One constrained fit per saved accepted active cell; no RHS/evolution."""
import argparse,collections,hashlib,importlib.util,json,sys,time
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce');modules=base/'Codex_contact_potential_modules_20261019'
sys.path[:0]=[str(modules),str(base/'Codex_compiled_midpoint_modules_20261018'),str(base/'Codex_face_state_batch_api_modules_20261018'),str(base/'Codex_moment_flux_batch_modules_20261018'),str(base/'Codex_coupled_relaxation_evolution_20261018'),str(base)]
from contact_potential_density import ContactPotentialDensity,ReconstructionFailure,features
from saved_contact_density_restore import restore_density
from saved_potential_taylor import taylor_coefficients
from embedded_laguerre_geometry import EmbeddedLaguerreGeometry
from audit_canonical_pair_history import gravity_coefficient
from run_canonical_av_branches_20261018 import read
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
def load(name,path):
 spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
args.output.mkdir(exist_ok=False);started=time.perf_counter()
root=base/'Codex_compiled_midpoint_once_20261018';archive=base/'Codex_numerical_Psi_endpoints_v2_20261018'
saved=np.load(archive/'accepted_reconstruction.npz');raw=np.column_stack([saved['mass'],saved['state'][:-4].reshape(-1,20)[:,3:]])
actual=np.load(root/'accepted.npz');assert np.array_equal(saved['state'],actual['state'])
mass=saved['mass'];ids=saved['ids'];N=len(ids);active=~actual['frozen'].astype(bool)
priors=json.loads((archive/'accepted_priors.json').read_text());fits=json.loads((archive/'accepted_fits.json').read_text())
offset=saved['polygon_offsets'];vertices=saved['polygon_vertices'];polys=[vertices[a:b] for a,b in zip(offset[:-1],offset[1:])]
initial_integrals=np.load(base/'Codex_potential_polynomial_20261018/initial_integrals.npz')
ledger=np.load(base/'Codex_saved_midpoint_Phi_flux_20261018/midpoint.npz')
h=float(json.loads((root/'manifest.json').read_text())['dt'])
ld=np.longdouble
Psi_initial=initial_integrals['integrals'][1,:,0]+initial_integrals['integrals'][1,:,1]
Psi_target=Psi_initial+ld(h)*ledger['Psi_evolution_masked'].astype(ld)
potential_path=base/'Codex_weight_remap_moments_centered_20261017/audit.py';potential=load('saved_potential',potential_path).potential
clip_path=base/'Codex_weighted_compression_geometry_20261017/audit.py';clip=load('saved_clip',clip_path)
ic,t0=read(base/'Codex_canonical_initial_20261017/run/khout.000000.dat');assert t0==0.
builder=EmbeddedLaguerreGeometry(np.column_stack([ic['x']-3,ic['y']-3]),ids,clip.clip)
a=actual['state'][:-4].reshape(-1,20);rebuilt,edges=builder(a[:,:2],a[:,2],float(actual['time']))
assert all(np.array_equal(p,q) for p,q in zip(polys,rebuilt))
gx,gw=np.polynomial.legendre.leggauss(8);gx=(gx+1)/2;gw/=2
trace_points=[];trace_owners=[];trace_weights=[];trace_normal=[];cell_nodes=[[] for _ in ids]
for edge_id,(i,j,p,q,n,length) in enumerate(edges):
 if j<0:continue
 delta=q-p;cuts=[0.,1.]
 for cell in (i,j):
  meta=priors[cell]
  if meta['flat']:continue
  normal=np.asarray(meta['normal']);den=float(delta@normal)
  if den!=0:
   t=(meta['offset']-float(p@normal))/den
   if 0<t<1:cuts.append(float(t))
 cuts=sorted(set(cuts))
 for lo,hi in zip(cuts[:-1],cuts[1:]):
  for g,w in zip(gx,gw):
   point=p+(lo+(hi-lo)*g)*delta;k=len(trace_points)
   trace_points.append(point);trace_owners.append((i,j));trace_weights.append(length*(hi-lo)*w);trace_normal.append(n)
   cell_nodes[i].append((k,0));cell_nodes[j].append((k,1))
trace_points=np.array(trace_points);trace_owners=np.array(trace_owners);trace_normal=np.array(trace_normal);trace_weights=np.array(trace_weights)
trace_old=np.full((len(trace_points),2),np.nan);trace_new=np.full_like(trace_old,np.nan)
manifest=dict(scope=__doc__,provider_calls=0,evolution_steps=0,active_cells=int(active.sum()),
 expected_constrained_fit_calls=int(active.sum()),fit_order=24,independent_order=48,potential_spline_radii=[.25,.5],
 target='Psi_initial numerical q48 + h * stored half-step Psi transport; raw18 untouched',
 input_sha256={str(p):sha(p) for p in (root/'accepted.npz',archive/'accepted_reconstruction.npz',archive/'accepted_priors.json',archive/'accepted_fits.json',base/'Codex_potential_polynomial_20261018/initial_integrals.npz',base/'Codex_saved_midpoint_Phi_flux_20261018/midpoint.npz')},
 source_sha256={str(p):sha(p) for p in (Path(__file__),modules/'contact_potential_density.py',modules/'positive_radial_quadrature.py',potential_path,clip_path)})
(args.output/'manifest.json').write_text(json.dumps(manifest,indent=2));np.savez(args.output/'targets.npz',ids=ids,active=active,Psi_initial=Psi_initial,Psi_target=Psi_target)
rows=[];parameters=[];errors=np.full((N,8),np.nan);attempts=0;counts=collections.Counter()
for k in range(N):
 r=raw[k];old=restore_density(polys[k],r,priors[k],fits[k]);items=cell_nodes[k]
 for node,side in items:trace_old[node,side]=old.value_trace(trace_points[node],(2*side-1)*trace_normal[node])[0]
 if not active[k]:continue
 attempts+=1
 try:
  taylor=taylor_coefficients(old.origin[None,:],potential,gravity_coefficient);s=ld(old.scale);H=taylor['H'][0].astype(ld)
  reference=np.r_[ld(taylor['Phi'][0]),s*taylor['gradient'][0].astype(ld),.5*s*s*H[0,0],s*s*H[0,1],.5*s*s*H[1,1]]
  fit=ContactPotentialDensity(polys[k],r[0],r[1:3],r[3:7].reshape(2,2),priors[k],old,potential,Psi_target[k],reference,24,(.25,.5))
  # Store sufficient parameters immediately, before independent integration.
  parameter=dict(id=int(ids[k]),theta=fit.theta.tolist(),logZ=fit.logZ,origin=fit.origin.tolist(),scale=fit.scale,
    reference=[str(x) for x in fit.reference],projection_remainder=[str(x) for x in fit.remainder_projection],sigma=str(fit.sigma),target=str(fit.Psi_target))
  parameters.append(parameter)
  with (args.output/'parameters.jsonl').open('a') as f:f.write(json.dumps(parameter)+'\n')
  points,weights,prior=fit.quadrature(48);density=fit._value_with_prior(points,prior);dm=weights*density
  A=features(points,fit.origin,fit.scale);Phi=potential(points.astype(ld));got=dm.astype(ld)@A
  expectation=np.asarray(fit.target_A)*ld(r[0]);Psi=dm.astype(ld)@Phi
  Psi_remainder=dm.astype(ld)@(Phi-A@fit.reference-A@fit.remainder_projection)
  expected_remainder=fit.Psi_target-ld(r[0])*(fit.target_A@fit.projection)
  target_equivalence=float((Psi-fit.Psi_target)-(Psi_remainder-expected_remainder))
  quadratic_error=got-expectation;Psi_error=Psi-fit.Psi_target
  # Existing component gate; normalized independent residuals are also reported.
  finite=bool(np.isfinite(dm).all() and np.all(density>0))
  gate=finite and bool(np.all(abs(quadratic_error)<=1e-14+1e-10*abs(expectation))) and abs(Psi_error)<=1e-14+1e-10*abs(fit.Psi_target)
  oldrho=old._value_with_prior(points,prior);old_dm=weights*oldrho
  oldPsi=old_dm.astype(ld)@Phi;oldA=old_dm.astype(ld)@A
  self_target_mean=oldPsi/np.sum(old_dm.astype(ld));self_quadratic_residual=float(abs(oldA/np.sum(old_dm.astype(ld))-fit.target_A).max())
  old_chi=Phi-A@fit.reference-A@fit.remainder_projection
  self_remainder_target=(self_target_mean-fit.target_A@fit.projection)/fit.sigma
  self_remainder_observed=np.sum(old_dm.astype(ld)*old_chi)/np.sum(old_dm.astype(ld))/fit.sigma
  old_entropy=float(np.sum(old_dm.astype(ld)*np.log((oldrho/prior).astype(ld))))
  trace_delta=[]
  for node,side in items:
   value=fit.value_trace(trace_points[node],(2*side-1)*trace_normal[node])[0];trace_new[node,side]=value
   trace_delta.append(value/trace_old[node,side]-1)
  normalized_error=float(abs(got/r[0]-fit.target_A).max())
  chi_residual=float(abs(Psi_remainder-expected_remainder)/(ld(r[0])*fit.sigma))
  errors[k]=[float(abs(quadratic_error).max()),float(Psi_error),normalized_error,chi_residual,fit.condition,float(fit.sigma),float(fit.rank_error),float(np.sum(weights*abs(density-oldrho))/r[0])]
  status='passed' if gate else 'independent_quadrature_gate_failed'
  row=dict(id=int(ids[k]),status=status,iterations=fit.iterations,dual_residual=fit.residual,condition=fit.condition,
    sigma=float(fit.sigma),rank_error=float(fit.rank_error),independent_q48_Psi_error=float(Psi_error),independent_q48_normalized_moment_error=normalized_error,
    independent_q48_scaled_remainder_residual=chi_residual,target_equivalence_residual=target_equivalence,
    self_control_reused_existing_fit=True,self_control_normalized_moment_residual=self_quadratic_residual,
    self_control_scaled_remainder_residual=float(abs(self_remainder_observed-self_remainder_target)),
    minimum_sampled_density=float(density.min()),entropy_change_q24_vs_old_q48=fit.entropy-old_entropy,
    density_relative_L1_mass=float(np.sum(weights*abs(density-oldrho))/r[0]),face_relative_change_max=float(max(abs(np.array(trace_delta)),default=0.)),
    old_Psi_minus_target=float(oldPsi-fit.Psi_target),initial_dual_residual=fit.initial_dual_residual)
 except ReconstructionFailure as exc:row=dict(id=int(ids[k]),status=exc.status,details=exc.details)
 except Exception as exc:row=dict(id=int(ids[k]),status='audit_exception',error=repr(exc))
 counts[row['status']]+=1;rows.append(row)
 with (args.output/'cells.jsonl').open('a') as f:f.write(json.dumps(row)+'\n')
 if attempts%100==0:print(json.dumps(dict(cells=attempts,counts=dict(counts),seconds=time.perf_counter()-started)),flush=True)
assert attempts==int(active.sum())
np.savez(args.output/'diagnostics.npz',ids=ids,active=active,errors=errors,trace_points=trace_points,trace_owners=trace_owners,trace_weights=trace_weights,trace_old=trace_old,trace_new=trace_new)
assert np.array_equal(actual['state'],saved['state']) and np.array_equal(mass,actual['mass'])
report=dict(status='completed_offline_reconstruction',manifest=manifest,constrained_fit_calls=attempts,counts=dict(counts),seconds=time.perf_counter()-started,
    raw_state_unchanged=True,maximum_independent_Psi_error=float(np.nanmax(abs(errors[:,1]))) if np.isfinite(errors[:,1]).any() else None,
    unresolved_or_failed_ids=[r['id'] for r in rows if r['status']!='passed'],
    limits='One constrained fit per active cell; old fit restored without Newton. No fallback. Independent q48 uses same prior and spline splits. Sampled positivity and numerical hull checks do not establish continuous realizability or evolved stability. Frozen cells are not assigned a replacement closure. Failed fit traces remain NaN, never filled from old density.')
(args.output/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2),flush=True)
