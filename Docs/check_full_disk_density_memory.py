"""Verify the production fit path against checkpoint-audited saved fits."""
import argparse,hashlib,json,sys,time,struct
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce');fast=base/'Codex_full_disk_fast_modules_v1_20261019'
sys.path.insert(0,str(fast))
from full_disk_density_memory import prepare
p=argparse.ArgumentParser(description=__doc__)
for k in ('state','geometry','baseline','output'):p.add_argument('--'+k,type=Path,required=True)
args=p.parse_args();args.output.mkdir(exist_ok=False)
data=np.load(args.state);native=json.loads((args.geometry/'full_exact_supports.json').read_text());mesh=np.load(args.geometry/'mesh.npz')
sha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
stage_key=sha(args.state)+':'+sha(args.geometry/'full_exact_supports.json')
result=prepare(data,native,mesh,stage_key,workers=8,independent_q48=False)
expected=[json.loads(l) for l in (args.baseline/'cells.jsonl').read_text().splitlines()]
def exact(x):
    if isinstance(x,float):return ('binary64',struct.pack('>d',x).hex())
    if isinstance(x,list):return [exact(v) for v in x]
    if isinstance(x,dict):return {k:exact(v) for k,v in sorted(x.items())}
    return x
checks=[]
for a,b in zip(expected,result.rows,strict=True):
    row=dict(id=a['id'],passed=all(exact(a[k])==exact(b[k]) for k in ('id','status','prior','parameter','nodes24')))
    checks.append(row)
    if not row['passed']:
        (args.output/'failure.json').write_text(json.dumps(dict(expected=a,actual=b),indent=2));raise RuntimeError(row)
report=dict(status='passed_production_density_fits_bitwise',cells=len(checks),preparation=result.report,
    matched_fields=['id','status','prior','parameter','nodes24'],checkpoint_q48_omitted=True,
    state_sha256=sha(args.state),provider_calls=0,evolution_steps=0,
    source_sha256={str(q):sha(q) for q in (Path(__file__),fast/'full_disk_density_memory.py')},
    limits='Every production fit matches the audited fit bitwise. Independent q48 remains required at specified physical checkpoints; omitting audit does not change q24 or Newton.')
(args.output/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
