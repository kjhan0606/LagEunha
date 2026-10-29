"""Audit short H200 ABBA: full-record parity and observed wall timing."""
import hashlib
import json
from pathlib import Path
import re
import numpy as np
from audit_kepler_restart_bytes_20261029 import compare
from kepler_native_global_ids_20261029 import ordered_global, ID_MASK
from audit_kepler_q1_thermal64_promotion_20261029 import payload, TAG

root=Path('/gpfs/kjhan/LagForce/Kepler_Thermal64WallABBA_H200_20261029')
cfg=json.loads((root/'contract.json').read_text());out=root/'audit.json'
if out.exists():raise FileExistsError(out)
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
seed=root/'A1/khout.020790.dat'
if sha(seed)!=cfg['seed_sha256']:raise ValueError('Changed parent')
p0,_=ordered_global(seed);ids=p0['u4if']&ID_MASK;data={};logs={};endpoints={};checks={}
for e in cfg['order']:
    name=e['case'];d=root/name
    if sha(d/'eunha2.exe')!=e['binary_sha256'] or sha(d/'run.sh')!=e['run_sha256'] or sha(d/'inputs.sha256')!=e['inputs_sha256'] or sha(d/'khout.020790.dat')!=cfg['seed_sha256']:raise ValueError('Changed frozen run '+name)
    files=[x for x in d.glob('khout.*.dat') if x.name!='khout.020790.dat']
    if len(files)!=1:raise ValueError('Incomplete or ambiguous endpoint '+name)
    p,m=ordered_global(files[0],ids);q=payload(files[0],len(ids));text=(d/'log').read_text()
    checks[name+'_state']=(d/'exit_status.txt').read_text().strip()=='0' and m['t']==cfg['target_time'] and np.array_equal(p0['mass'],p['mass']) and np.all(p['w2']==0) and np.all(q['tag']==TAG) and np.all(np.isfinite(q['u'])) and np.all(q['u']>0) and np.all(np.isfinite(q['v'])) and np.all(q['v']>0)
    logs[name]=[s.rstrip() for s in re.findall(r'\[(?:RK4E|WALL_ACCEPT|WALL_ROLLBACK)\][^\n\[]*',text)]
    timers=re.findall(r'\[RK4 gap\][^\n\[]*',text)
    rows=[]
    for line in timers:
        step=int(re.search(r'\bstep=(\d+)',line)[1])
        if step>2:rows.append(float(re.search(r'\bwall=(\S+?)s(?:\s|$)',line)[1]))
    if len(rows)<10:raise ValueError('Insufficient warm step timings')
    checks[name+'_no_injection']=all(all(float(re.search(r'\b'+k+r'=(\S+)',line)[1])==0 for k in ('floor_cum','de_cum','sfl_cum')) for line in logs[name] if line.startswith('[RK4E]')) and any(line.startswith('[RK4E]') for line in logs[name])
    data[name]={'median_wall_seconds':float(np.median(rows)),'timed_steps':len(rows),'endpoint_sha256':sha(files[0])};endpoints[name]=files[0]
checks['numerical_ledgers_exact']=all(logs[n]==logs['A1'] for n in ('B1','B2','A2'))
parity={n:compare(endpoints['A1'],endpoints[n]) for n in ('B1','B2','A2')}
checks['all_nonpointer_records_and_footer_exact']=all(p['passed'] for p in parity.values())
baseline=float(np.median([data[k]['median_wall_seconds'] for k in ('A1','A2')]))
candidate=float(np.median([data[k]['median_wall_seconds'] for k in ('B1','B2')]))
result={'passed':all(bool(x) for x in checks.values()),'checks':{k:bool(v) for k,v in checks.items()},'cases':data,
        'baseline_abba_median_wall_seconds':baseline,'candidate_abba_median_wall_seconds':candidate,'observed_speedup':baseline/candidate,
        'full_restart_comparisons':parity,'scope':'Short same-GPU ABBA from t10; 2 warmup steps excluded. Timing is observational and does not prove long-time, collision or MPI correctness or scaling.'}
with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
if not result['passed']:raise SystemExit(1)
