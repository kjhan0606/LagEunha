"""Explicit change-of-integrator fork; no numerical update or acceptance.

Authenticate a saved fine .025 state and its completed original audits before
serializing all six numerical payloads. This is not a same-method restart.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import numpy as np
from full_moment_state_io import dump,load,HEADER,FIELDS
from full_moment_restart_provenance import sha,physics_sha,OPERATOR,OPERATOR_SHA

METHOD='self-consistent midpoint + Sylvester/exponential/positive QR'
SOURCE_METHOD='explicit midpoint + Sylvester/exponential/positive QR'
FORMAT='full_moment_explicit_fine_to_coupled_fork_v1'

def seal(checkpoint,initial,physics,output):
    checkpoint,initial,physics,output=map(lambda p:Path(p).resolve(),(checkpoint,initial,physics,output))
    saved=load(checkpoint/'accepted.npz');old=load(initial)
    assert checkpoint.name=='fine_0010' and float(saved['time'])==.025
    assert float(old['time'])==0 and len(saved['ids'])==19945
    assert np.array_equal(saved['ids'],old['ids']) and np.array_equal(saved['mass'],old['mass'])
    assert not saved['frozen'].any() and saved['psi_constrained_mask'].all()
    manifest=checkpoint.parent.parent/'manifest.json'
    source=json.loads(manifest.read_text())
    assert source['mode']=='fine' and source['steps']['fine_h']==.0025
    assert source['initial_sha256']==sha(initial) and source['fixed_mass'] and source['Re']=='infinity' and not source['steering']
    completed=json.loads((checkpoint.parent.parent/'result.json').read_text())
    assert completed['status']=='completed_comparison_branch'
    checks=json.loads((checkpoint/'state_checks.json').read_text())
    assert checks['mass_bitwise_initial'] and checks['minimum_thermal']>0
    assert all(abs(x)<=b for x,b in zip(checks['P_E_L_ledger_residual'],checks['P_E_L_ledger_bound']))
    assert abs(checks['combined_energy_change'])<=checks['combined_energy_roundoff_bound']
    gate=json.loads((checkpoint/'physical_gate.json').read_text())
    assert gate['status']=='passed_checkpoint_physical_monitor' and gate['time']==.025
    assert all(row['passed'] for row in gate['regions'].values())
    density=json.loads((checkpoint/'density/report.json').read_text())
    assert density['status']=='passed_full_stage_Psi_density_preparation'
    assert density['state_sha256']==sha(checkpoint/'accepted.npz')
    assert sha(OPERATOR)==OPERATOR_SHA
    evidence=[manifest,checkpoint.parent.parent/'result.json',checkpoint/'state_checks.json',
        checkpoint/'physical_gate.json',checkpoint/'density/report.json',checkpoint/'physics/report.json']
    record=dict(format=FORMAT,source=str(checkpoint/'accepted.npz'),source_sha256=sha(checkpoint/'accepted.npz'),
        source_integrator=SOURCE_METHOD,target_integrator=METHOD,source_h=.0025,time=.025,cells=len(saved['ids']),
        initial_state_sha256=sha(initial),physics_sha256=physics_sha(physics),operator_sha256=OPERATOR_SHA,
        fixed_mass=True,Re_infinity=True,steering=False,evidence={str(p):sha(p) for p in evidence},
        allowed_target_h=[.005,.0025],maximum_target_time=.05,first_gate_time=.03,
        physical_providers=0,accepted_time_steps=0)
    output.mkdir(exist_ok=False)
    with (output/'manifest.json').open('x') as f:
        json.dump(record,f,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
    return output/'manifest.json'

def verify_serialize(source,manifest,initial,physics,requested_h,output):
    source,manifest,initial,physics,output=map(lambda p:Path(p).resolve(),(source,manifest,initial,physics,output))
    record=json.loads(manifest.read_text())
    assert record['format']==FORMAT and record['source']==str(source) and record['source_sha256']==sha(source)
    assert record['source_integrator']==SOURCE_METHOD and record['target_integrator']==METHOD
    assert record['time']==.025 and record['source_h']==.0025 and record['cells']==19945
    assert record['fixed_mass'] and record['Re_infinity'] and record['steering'] is False
    assert math.isfinite(requested_h) and requested_h in record['allowed_target_h']==[.005,.0025]
    assert record['initial_state_sha256']==sha(initial) and record['physics_sha256']==physics_sha(physics)
    assert record['operator_sha256']==sha(OPERATOR)==OPERATOR_SHA
    assert record['evidence']
    for file,digest in record['evidence'].items():assert sha(file)==digest,f'Changed fork evidence: {file}'
    state=load(source);old=load(initial)
    assert float(state['time'])==record['time'] and len(state['ids'])==record['cells']
    assert np.array_equal(state['ids'],old['ids']) and np.array_equal(state['mass'],old['mass'])
    # Serialization identity only, explicitly not an old simulation token.
    token=int(record['source_sha256'][:16],16) or 1
    output.mkdir(exist_ok=False)
    target=dump(output/'serialized_source',state,1,token)
    header=(target/'header.bin').read_bytes()
    with (output/'expected_header.bin').open('xb') as f:f.write(header)
    observed=load(target)
    assert all(np.array_equal(observed[k],state[k]) for k in observed)
    receipt=dict(status='verified_explicit_to_coupled_fork_no_time_step',
        source=str(source),source_sha256=record['source_sha256'],manifest_sha256=sha(manifest),
        source_integrator=SOURCE_METHOD,target_integrator=METHOD,time=record['time'],
        cells=record['cells'],physics_sha256=record['physics_sha256'],requested_h=requested_h,
        serialized_source=str(target),serialization_token=token,
        payload_sha256={key:sha(target/(key+'.bin')) for key in FIELDS},
        all_six_payload_values_preserved=True,physical_providers=0,accepted_time_steps=0)
    with (output/'receipt.json').open('x') as f:json.dump(receipt,f,indent=2)
    return header,receipt

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('checkpoint','initial','physics','output'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();print(seal(a.checkpoint,a.initial,a.physics,a.output))
