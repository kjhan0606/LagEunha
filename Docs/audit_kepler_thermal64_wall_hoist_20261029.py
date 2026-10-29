"""Exact records and numerical ledgers, excluding interleaved timing output."""
import json
from pathlib import Path
import re
from audit_kepler_restart_bytes_20261029 import compare

root=Path('/gpfs/kjhan/LagForce/Kepler_Thermal64WallHoist_20261029')
cases=[('no_reject','Kepler_Thermal64WallNoReject_20261029','Kepler_Thermal64HoistNoReject_20261029'),
       ('forced_reject','Kepler_Thermal64WallSmoke_20261029','Kepler_Thermal64HoistReject_20261029')]
def ledgers(path):
    text=(path/'log').read_text()
    # stdout/stderr can be concatenated mid-line; timing is not a numerical ledger.
    pat=r'\[(?:RK4E|WALL_ACCEPT|WALL_ROLLBACK)\][^\n\[]*'
    return [x.rstrip() for x in re.findall(pat,text)]
results={}
for name,aa,bb in cases:
    a=Path('/gpfs/kjhan/LagForce')/aa;b=Path('/gpfs/kjhan/LagForce')/bb
    result=compare(a/'khout.000023.dat',b/'khout.000023.dat')
    reference=ledgers(a);candidate=ledgers(b)
    result['numerical_ledgers_exact']=reference==candidate and len(reference)==(3 if name=='forced_reject' else 2)
    result['exit_zero']=(b/'exit_status.txt').read_text().strip()=='0'
    result['passed']=result['passed'] and result['numerical_ledgers_exact'] and result['exit_zero']
    results[name]=result
out=root/'parity_normalized.json'
if out.exists():raise FileExistsError(out)
with out.open('x') as f:json.dump(results,f,indent=2);f.write('\n')
print(json.dumps(results,indent=2))
if not all(x['passed'] for x in results.values()):raise SystemExit(1)
