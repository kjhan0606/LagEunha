"""Bitwise-valued saved metadata/fit comparison for native prior integration."""
import argparse,hashlib,json,struct
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
for k in ('baseline','candidate','output'):p.add_argument('--'+k,type=Path,required=True)
args=p.parse_args();args.output.mkdir(exist_ok=False)
a,b=[json.loads((q/'report.json').read_text()) for q in (args.baseline,args.candidate)]
assert a['status']==b['status']=='passed_full_stage_Psi_density_preparation'
assert a['state_sha256']==b['state_sha256'] and a['cells']==b['cells']
def canonical(x):
    if isinstance(x,float):return ('binary64',struct.pack('>d',x).hex())
    if isinstance(x,list):return [canonical(v) for v in x]
    if isinstance(x,dict):return {k:canonical(v) for k,v in sorted(x.items())}
    return x
count=0;witnesses=0
with (args.baseline/'cells.jsonl').open() as af,(args.candidate/'cells.jsonl').open() as bf:
    for al,bl in zip(af,bf,strict=True):
        aa,bb=json.loads(al),json.loads(bl);assert aa['id']==bb['id'] and aa['status']==bb['status']=='passed'
        prior=dict(bb['prior']);witnesses+=int(prior.pop('exact_bisection_native_branch_witness',None) is not None)
        for k in aa:
            if k in ('timers','peak_RSS_KiB'):continue
            observed=prior if k=='prior' else bb[k]
            assert canonical(aa[k])==canonical(observed),(aa['id'],k)
        count+=1
assert count==a['cells']
result=dict(status='passed_full_disk_density_replay_bitwise_values',cells=count,native_witnesses=witnesses,
    state_sha256=a['state_sha256'],baseline_seconds=a['seconds'],candidate_seconds=b['seconds'],
    observed_wall_ratio=b['seconds']/a['seconds'],baseline_worker_wall_sums=a['worker_wall_timing_sums'],
    candidate_worker_wall_sums=b['worker_wall_timing_sums'],
    source_sha256={str(q/'report.json'):hashlib.sha256((q/'report.json').read_bytes()).hexdigest() for q in (args.baseline,args.candidate)},
    limits='Every existing per-cell physical/fit/quadrature field compared with float bit patterns, except resource/timing fields. Additional native witnesses intentionally differ; comparison performs no RHS or evolution.')
(args.output/'report.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
