"""Original native executable: 1/2-rank real midpoint and previous endpoint.

Exact rank reassembly, original component roundoff scale for the saved Ada
endpoint, and independent actual-cell physical checkpoint monitoring required.
This validates one step, never long-time Kepler stability.
"""
import argparse,json,sys,re
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,required=True);p.add_argument('--modules',type=Path,required=True);args=p.parse_args()
sys.path.insert(0,str(args.modules));from full_moment_state_io import load
base=Path('/gpfs/kjhan/LagForce');run=args.run;states=[];summaries=[];binary_sha={}
for ranks in (1,2):
    folder=run/f'ranks{ranks}';summary=json.loads((folder/'completed_0001.json').read_text())
    with (folder/'native.log').open() as log:
        first=next(line for line in log if '[FULL_MOMENT_NATIVE_RUNTIME]' in line)
    binary_sha[str(ranks)]=re.search(r'binary_sha=([0-9a-f]{64})',first).group(1)
    assert summary['status']=='completed_native_evolution_batch' and summary['native_time_control'] and summary['accepted_commit']
    assert summary['ranks']==ranks and summary['step']==1 and summary['time']==.001 and summary['h']==.001 and summary['physical_providers']==2
    state=load(folder/'accepted_0001');assert float(state['time'])==.001
    physical=json.loads((folder/'candidate_audit_0001/physical_gate.json').read_text());assert physical['status']=='passed_checkpoint_physical_monitor'
    for epoch in (2,3):
        for rank in range(ranks):
            audit=json.loads((folder/f'factor_{epoch}_rank{rank}.json').read_text())
            assert audit['status']=='passed_owned_native_factor_audit' and audit['shadow_audit'] is not None
    for epoch in (1,2):
        audit=json.loads((folder/f'rhs_{epoch}.json').read_text());assert audit['status']=='passed_actual_full_moving_stage'
    states.append(state);summaries.append(summary)
exact={key:bool(np.array_equal(states[0][key],states[1][key])) for key in states[0]};assert all(exact.values()),exact
previous=load(base/'Codex_full_disk_coarse_accepted_20261019/accepted.npz');old_comparison={}
for key in ('ids','mass','frozen','psi_constrained_mask'):
    assert np.array_equal(states[1][key],previous[key]),key
for key in ('state','factors','Psi','weight_rate_seed'):
    a=states[1][key].astype(np.longdouble);b=previous[key].astype(np.longdouble)
    error=abs(a-b);bound=1e-14+1e-10*abs(b)
    assert np.all(error<=bound),(key,float(error.max()),float((error/bound).max()))
    old_comparison[key]=dict(maximum_difference=float(error.max()),maximum_original_bound_fraction=float((error/bound).max()))
report=dict(status='passed_native_MPI_midpoint_one_step',native_time_control=True,cells=len(states[0]['ids']),
    time=.001,h=.001,exact_rank1_vs_rank2=exact,previous_Ada_endpoint=old_comparison,
    native_collective_commit=True,fixed_mass=True,longdouble_Psi_preserved=True,
    original_global_stage_factor_QR_geometry_density_and_physical_gates=True,
    validation_actual_providers=4,accepted_steps_each=1,native_MPI_evolution_verified=True,
    control_binary_sha256=binary_sha,same_control_binary=len(set(binary_sha.values()))==1,
    long_time_stability_verified=False,
    limits='One real native midpoint step; does not establish t1 or the historical t6–10 stability target. Saved cross-host endpoint uses existing 1e-14+1e-10*abs(component) audit scale; ranks require exact equality. If control binaries differ, this is not a controlled MPI speed benchmark; the preserved 1-rank result predates the two-phase resource-release correction.')
(run/'validation_report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
