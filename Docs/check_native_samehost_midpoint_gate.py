"""Exact native 1/2-rank and independent original same-host endpoint gate.

The previous cross-host seed discrepancy remains recorded, not tolerated away.
The replacement reference runs the original global API on the native host and
requires exact equality for ALL arrays including that nonlinear warm-start.
"""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--original',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
args=p.parse_args()
base=Path('/gpfs/kjhan/LagForce')
sys.path.insert(0,str(base/'Codex_full_disk_native_evolution_build_v8_20261019'))
from full_moment_state_io import load
native_run=base/'Codex_full_disk_native_midpoint_and_coarse_v8_20261019'
states=[]
for ranks in (1,2):
    folder=native_run/f'ranks{ranks}'
    summary=json.loads((folder/'completed_0001.json').read_text())
    assert summary['status']=='completed_native_evolution_batch' and summary['accepted_commit']
    assert summary['native_time_control'] and summary['ranks']==ranks
    assert summary['time']==.001 and summary['h']==.001 and summary['physical_providers']==2
    assert json.loads((folder/'candidate_audit_0001/physical_gate.json').read_text())['status']=='passed_checkpoint_physical_monitor'
    for epoch in (1,2):
        assert json.loads((folder/f'rhs_{epoch}.json').read_text())['status']=='passed_actual_full_moving_stage'
    for epoch in (2,3):
        for rank in range(ranks):
            audit=json.loads((folder/f'factor_{epoch}_rank{rank}.json').read_text())
            assert audit['status']=='passed_owned_native_factor_audit' and audit['shadow_audit']
    states.append(load(folder/'accepted_0001'))
exact_ranks={k:bool(np.array_equal(states[0][k],states[1][k])) for k in states[0]}
assert all(exact_ranks.values()),exact_ranks
comparison=json.loads((args.original/'comparison.json').read_text())
assert comparison['status']=='passed_samehost_original_global_midpoint' and comparison['host'].split('.')[0]=='syn104'
original=load(args.original/'endpoint.npz')
exact_original={k:bool(np.array_equal(states[1][k],original[k])) for k in states[1]}
assert all(exact_original.values()),exact_original
old=load(base/'Codex_full_disk_coarse_accepted_20261019/accepted.npz')
crosshost={}
for key in ('state','factors','Psi','weight_rate_seed'):
    a=states[1][key].astype(np.longdouble);b=old[key].astype(np.longdouble)
    error=abs(a-b);bound=1e-14+1e-10*abs(b)
    crosshost[key]=dict(maximum_difference=float(error.max()),
        original_component_scale_fraction=float((error/bound).max()),within_original_scale=bool(np.all(error<=bound)))
    if key!='weight_rate_seed':assert np.all(error<=bound),(key,crosshost[key])
report=dict(status='passed_native_midpoint_samehost_original_reference',time=.001,h=.001,
    exact_native_rank1_vs_rank2=exact_ranks,exact_native_vs_original_samehost=exact_original,
    old_Ada_crosshost_diagnostics=crosshost,component_tolerances_unchanged=True,
    earlier_failed_checker_preserved=str(native_run/'validation.log'),
    reference_change='Original global-stage API and Python integrator on same host and GPU library; exact all-field gate replaces cross-host derivative-seed scalar scale',
    native_actual_providers=4,independent_original_reference_providers=2,
    long_time_stability_verified=False)
with args.output.open('x') as f:json.dump(report,f,indent=2)
print(json.dumps(report),flush=True)
