"""One real restarted step vs saved uninterrupted native evolution."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np


def sha(file):return hashlib.sha256(Path(file).read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run',type=Path,required=True)
    p.add_argument('--source-run',type=Path,required=True)
    p.add_argument('--modules',type=Path,required=True)
    args=p.parse_args()
    sys.path.insert(0,str(args.modules))
    from full_moment_state_io import load,HEADER,FIELDS
    from full_moment_restart_provenance import verify
    base=Path('/gpfs/kjhan/LagForce');runtime=args.run/'runtime'
    old=args.source_run/'accepted_0002'
    ic=base/'Codex_full_curved_kepler_IC_20261019/initial.npz'
    header,record=verify(old,args.run/'input_manifest/manifest.json',ic,ic.parent/'report.json',.005)
    receipt=json.loads((runtime/'restart_import.json').read_text())
    assert receipt['status']=='verified_restart_provenance_no_state_update' and receipt['time']==.01
    assert receipt['old_run_token']==record['original_run_token'] and receipt['new_run_token']!=record['original_run_token']
    imported=Path(f"/dev/shm/eunha_run_{receipt['new_run_token']}/states/e1")
    new_header=HEADER.unpack((imported/'header.bin').read_bytes());expected=list(HEADER.unpack(header))
    expected[3]=1;expected[4]=receipt['new_run_token']
    assert list(new_header)==expected
    payload={key:sha(old/(key+'.bin'))==sha(imported/(key+'.bin')) for key in FIELDS}
    assert all(payload.values()),payload
    summary=json.loads((runtime/'completed_0003.json').read_text())
    assert summary['status']=='completed_native_evolution_batch' and summary['accepted_commit']
    assert summary['time']==.015 and summary['step']==3 and summary['h']==.005 and summary['physical_providers']==2
    assert summary['ranks']==2
    assert not (runtime/'accepted_0002').exists()
    restarted=load(runtime/'accepted_0003');continuous=load(args.source_run/'accepted_0003')
    exact={key:bool(np.array_equal(restarted[key],continuous[key])) for key in continuous}
    for label in ('restart_audit','candidate_audit_0003'):
        assert json.loads((runtime/label/'physical_gate.json').read_text())['status']=='passed_checkpoint_physical_monitor'
        assert json.loads((runtime/label/'density/report.json').read_text())['status']=='passed_full_stage_Psi_density_preparation'
    for epoch in (1,2):
        assert json.loads((runtime/f'rhs_{epoch}.json').read_text())['status']=='passed_actual_full_moving_stage'
    for epoch in (2,3):
        for rank in (0,1):
            row=json.loads((runtime/f'factor_{epoch}_rank{rank}.json').read_text())
            assert row['status']=='passed_owned_native_factor_audit' and row['shadow_audit']
    identity=json.loads((args.run/'runtime_identity.json').read_text())
    assert identity['status']=='passed_live_native_runtime_identity'
    report=dict(status='passed_native_restart_one_step_equivalence' if all(exact.values()) else 'failed_native_restart_one_step_equivalence',
        from_time=.01,to_time=.015,h=.005,new_physical_providers=2,new_accepted_steps=1,
        original_checkpoint_unchanged=True,import_payload_bytes_exact=payload,
        exact_value_restarted_vs_uninterrupted=exact,
        actual_loaded_runtime_identity_verified=True,initial_and_final_q48_physical_gates_passed=True,
        no_initial_interval_recalculation=True,long_time_stability_verified=False)
    with (args.run/'report.json').open('x') as f:json.dump(report,f,indent=2)
    assert all(exact.values()),exact
    print(json.dumps(report),flush=True)


if __name__=='__main__':main()
