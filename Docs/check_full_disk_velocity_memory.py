"""Compare anonymous shared gravity preparation against its saved-stage oracle."""
import argparse,hashlib,json,os,sys,time
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce');fast=base/'Codex_full_disk_fast_modules_v1_20261019'
sys.path.insert(0,str(fast))
from full_disk_velocity_memory import prepare
p=argparse.ArgumentParser(description=__doc__)
for k in ('state','geometry','density','baseline','output'):p.add_argument('--'+k,type=Path,required=True)
args=p.parse_args();args.output.mkdir(exist_ok=False);start=time.perf_counter()
data=np.load(args.state);native=json.loads((args.geometry/'full_exact_supports.json').read_text());mesh=np.load(args.geometry/'mesh.npz')
fits=[json.loads(l) for l in (args.density/'cells.jsonl').read_text().splitlines()]
sha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
state_sha=sha(args.state)
assert json.loads((args.density/'report.json').read_text())['state_sha256']==state_sha
assert json.loads((args.geometry/'report.json').read_text())['state_sha256']==state_sha
stage_key=state_sha+':'+sha(args.geometry/'full_exact_supports.json')
result=prepare(data,native,mesh,fits,stage_key,8)
baseline=np.load(args.baseline/'prepared.npz');checks={}
assert set(baseline.files)==set(result.prepared)
for key in baseline.files:
    a,b=baseline[key],result.prepared[key]
    checks[key]=dict(shape=list(a.shape),dtype=str(a.dtype),
        bitwise_equal=a.dtype==b.dtype and a.shape==b.shape and a.tobytes()==b.tobytes())
old=np.load(args.baseline/'gravity_nodes.npy',mmap_mode='r')
assert old.shape==result.gravity_nodes.shape
gravity_equal=True
for lo in range(0,len(old),250000):
    if not np.array_equal(old[lo:lo+250000],result.gravity_nodes[lo:lo+250000]):gravity_equal=False;break
passed=gravity_equal and all(row['bitwise_equal'] for row in checks.values())
report=dict(status='passed_in_memory_velocity_bitwise' if passed else 'failed_in_memory_velocity_bitwise',
    prepared_arrays=checks,gravity_nodes_bitwise_equal=gravity_equal,state_sha256=state_sha,
    preparation=result.report,baseline_seconds=json.loads((args.baseline/'report.json').read_text())['seconds'],
    seconds=time.perf_counter()-start,provider_calls=0,evolution_steps=0,
    source_sha256={str(q):sha(q) for q in (Path(__file__),fast/'full_disk_velocity_memory.py')},
    limits='Checks every prepared array and all gravity-node values. Native stage consumption reuses its existing contiguous-array ABI; no new physical RHS or evolution.')
(args.output/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
assert passed,'Preserved component mismatch; no production adoption'
