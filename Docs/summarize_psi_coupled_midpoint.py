"""Audit saved coupled connection, without a provider, evolution, or refit."""
import argparse,hashlib,json,sys
from pathlib import Path
import numpy as np
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--root',type=Path,required=True);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
root=args.root;summary=json.loads((root/'summary.json').read_text());assert summary['status']=='completed_connection_gate'
initial=np.load(root/'initial.npz');final=np.load(root/'accepted.npz');mass=initial['mass'];ids=initial['ids'];n=len(ids);active=initial['psi_constrained_mask'];frozen=~active;h=float(final['time']);ld=np.longdouble;eps=np.finfo(float).eps
assert summary['provider_calls']==2 and h==1e-4
assert np.array_equal(final['mass'].view(np.uint64),mass.view(np.uint64)) and np.array_equal(final['ids'],ids)
a0=initial['state'][:-4].reshape(-1,20);af=final['state'][:-4].reshape(-1,20)
assert np.array_equal(af[frozen].view(np.uint64),a0[frozen].view(np.uint64))
assert np.array_equal(final['factors'][frozen].view(np.uint64),initial['factors'][frozen].view(np.uint64))
assert np.array_equal(final['Psi'][frozen],initial['Psi'][frozen])
rates0=np.load(root/'stage0/Psi_rates.npz');ratesh=np.load(root/'midpoint/Psi_rates.npz');half=np.load(root/'half_Psi.npz')
assert np.array_equal(half['Psi'],initial['Psi']+ld(h/2)*rates0['rate'].astype(ld))
assert np.array_equal(final['Psi'],initial['Psi']+ld(h)*ratesh['rate'].astype(ld))

def totals(a):return np.r_[a[:,9:11].sum(0),a[:,19].sum(),(a[:,12]-a[:,13]).sum()]
def absolute(a):return np.r_[abs(a[:,9:11]).sum(0),abs(a[:,19]).sum(),(abs(a[:,12])+abs(a[:,13])).sum()]
ledger=final['state'][-4:]-initial['state'][-4:];residual=totals(af)-totals(a0)-ledger
cuda=np.load(root/'midpoint/cuda.npz');stage=cuda['stage'][:,:18];reported=json.loads((root/'midpoint/report.json').read_text())
absolute_rates=np.r_[abs(stage[:,7:9]).sum(0),abs(stage[:,17]).sum(),(abs(stage[:,10])+abs(stage[:,11])).sum()]
scale=absolute(af)+absolute(a0)+h*(absolute_rates+abs(cuda['ledgers'].reshape(3,4)).sum(0))
bound=(128+4*n)*eps*scale+h*np.array(reported['backend_comparison']['cuda']['global_ledger_bound'])
assert np.all(abs(residual)<bound)
reconstruction={}
for label in ('initial','half','final'):
 q=np.load(root/(label+'_q48.npz'));errors=q['errors'];params=[json.loads(line) for line in (root/(label+'_density')/'parameters.jsonl').read_text().splitlines()]
 selected=[p for p in params if p['active']];assert len(selected)==int(active.sum())
 assert all(len(p['theta'])==6 and float(p['sigma'])>float(p['rank_error']) and p['residual']<3e-13 for p in selected)
 reconstruction[label]=dict(maximum_q48_normalized_mqQ_error=float(abs(errors[:,1]).max()),maximum_q48_Psi_error=float(abs(errors[:,2]).max()),q48_Psi_signed_error=float(errors[:,2].sum()),q48_Psi_absolute_error=float(abs(errors[:,2]).sum()),maximum_q48_scaled_remainder_residual=float(abs(errors[:,3]).max()),maximum_dual_residual=max(p['residual'] for p in selected),maximum_condition=max(p['condition'] for p in selected),minimum_sigma=min(float(p['sigma']) for p in selected),maximum_rank_error=max(float(p['rank_error']) for p in selected))
q0=np.load(root/'initial_q48.npz')['Psi'];qf=np.load(root/'final_q48.npz')['Psi'];Psi0=initial['Psi'];Psif=final['Psi'];hydrodelta=np.sum((af[active,19]-a0[active,19]).astype(ld));boundary=ld(h)*(ld(ratesh['boundary_hydro'])+ld(ratesh['boundary_Phi']))
transported=hydrodelta+np.sum(Psif[active]-Psi0[active])-boundary
reconstructed=hydrodelta+np.sum(qf[active]-q0[active])-boundary
# Explicit envelope for stored finite sums and the measured q48 constraint
# discrepancy. Quadrature truncation remains a separate q24/48 diagnostic.
ops=128+4*n;gamma=ops*eps/(1-ops*eps)
energy_scale=float(np.sum(abs(af[active,19]))+np.sum(abs(a0[active,19]))+np.sum(abs(Psif[active]))+np.sum(abs(Psi0[active]))+abs(boundary))
roundoff_bound=gamma*energy_scale+h*float(ratesh['bound'])+h*float(np.array(reported['backend_comparison']['cuda']['global_ledger_bound'])[2])
constraint_uncertainty=np.sum(abs(q0[active]-Psi0[active]))+np.sum(abs(qf[active]-Psif[active]))
assert abs(transported)<=roundoff_bound
assert abs(reconstructed)<=ld(roundoff_bound)+constraint_uncertainty
physical=af[:,3:];Q=physical[:,2:6].reshape(-1,2,2);Pi=physical[:,12:16].reshape(-1,2,2)
assert np.array_equal(Q[:,0,1],Q[:,1,0]) and np.array_equal(Pi[:,0,1],Pi[:,1,0])
report=dict(status='passed_saved_connection_audit',scope=__doc__,accepted_time=h,provider_calls=summary['provider_calls'],mass_bitwise_unchanged=True,frozen_raw_factor_sidecar_unchanged=True,Q_Pi_symmetric=True,minimum_thermal=summary['minimum_thermal'],minimum_R=summary['minimum_R'],accepted_momentum_E_angular_ledger_residual=residual.tolist(),accepted_momentum_E_angular_ledger_bound=bound.tolist(),reconstruction=reconstruction,energy=dict(transported_residual=float(transported),reconstructed_q48_residual=float(reconstructed),stored_sum_roundoff_bound=float(roundoff_bound),measured_initial_and_final_q48_constraint_uncertainty=float(constraint_uncertainty),frozen_reconstructed_Psi_change=float(np.sum(qf[frozen]-q0[frozen])),frozen_transported_Psi_change=float(np.sum(Psif[frozen]-Psi0[frozen])),quadrature_truncation='pending q24 minus stored independent q48; not included in roundoff bound'),cost=summary['timers'],rhs_outer_cost=summary['rhs_outer'],physical_accuracy='pending each-current-polygon exact reference comparison',source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),root/'manifest.json',root/'summary.json',root/'initial.npz',root/'accepted.npz')},limits='Connection and numerical ledger audit only; no long-time stability, global circular domain, time order, or physical improvement verdict.')
args.output.write_text(json.dumps(report,indent=2));print(json.dumps({k:v for k,v in report.items() if k not in ('cost','rhs_outer_cost','source_sha256')},indent=2))
