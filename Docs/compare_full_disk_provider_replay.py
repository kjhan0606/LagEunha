"""Exact saved-stage comparison for the same-state optimized provider replay."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser(description=__doc__)
for k in ('baseline','candidate','output'):p.add_argument('--'+k,type=Path,required=True)
args=p.parse_args();args.output.mkdir(exist_ok=False)
a,b=[json.loads((root/'report.json').read_text()) for root in (args.baseline,args.candidate)]
assert a['status']==b['status']=='passed_actual_full_moving_stage'
assert a['time']==b['time'] and a['cells']==b['cells']
aa,bb=[np.load(root/'stage0.npz') for root in (args.baseline,args.candidate)]
assert set(aa.files)==set(bb.files)
checks={}
for key in aa.files:
    x,y=aa[key],bb[key]
    bitwise=x.dtype==y.dtype and x.shape==y.shape and x.tobytes()==y.tobytes()
    row=dict(bitwise_equal=bitwise,shape=list(x.shape),dtype=str(x.dtype))
    if not bitwise and np.issubdtype(x.dtype,np.number) and x.shape==y.shape:
        row['maximum_difference']=float(np.max(abs(x-y)))
    checks[key]=row
hash_checks={k:a[k]==b[k] for k in ('node_inputs_sha256','node_extra_sha256')}
passed=all(r['bitwise_equal'] for r in checks.values()) and all(hash_checks.values())
result=dict(status='passed_provider_replay_bitwise' if passed else 'failed_provider_replay_bitwise',
    time=a['time'],arrays=checks,compiled_input_hash_checks=hash_checks,
    baseline_wall_seconds=a['total_seconds'],candidate_wall_seconds=b['total_seconds'],
    observed_wall_ratio=b['total_seconds']/a['total_seconds'],baseline_outer=a['outer'],candidate_outer=b['outer'],
    cache=b.get('same_stage_face_cache'),
    source_sha256={str(root/'report.json'):hashlib.sha256((root/'report.json').read_bytes()).hexdigest() for root in (args.baseline,args.candidate)},
    new_unique_physical_provider_calls=0,same_state_provider_replays=1,
    limits='Same physical stage, source seed, existing prepared velocity/gravity and q24 density. Resource/timing differences do not alter physical state. No evolution or production exclusive-timing claim.')
(args.output/'report.json').write_text(json.dumps(result,indent=2));print(json.dumps({k:v for k,v in result.items() if k not in ('baseline_outer','candidate_outer')},indent=2))
assert passed,'Preserved mismatch report; no automatic replay or evolution'
