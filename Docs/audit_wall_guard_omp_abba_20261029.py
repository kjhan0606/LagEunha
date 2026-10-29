"""Same-binary native thermal64 serial/OMP wall-guard ABBA qualification."""
import hashlib
import json
from pathlib import Path
import re
import statistics
import sys
import numpy as np
from audit_kepler_restart_bytes_20261029 import compare
from audit_kepler_wall_events_interleaved_20261029 import ledger_number
from kepler_native_global_ids_20261029 import ordered_global, ID_MASK

TAG=2026102603.0
def payload(path,n):
    dt=np.dtype({'names':['id','u','v','rhs','tag','anchor'],
                 'formats':['<u8']+['<f8']*5,
                 'offsets':[8,472,480,488,496,504],'itemsize':576})
    p=np.memmap(path,mode='r',dtype=dt,offset=16,shape=(n,))
    return p[np.argsort(p['id']&ID_MASK)]

root=Path(sys.argv[1]);stage=json.loads((root/'stage.json').read_text())
def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
initial,m0=ordered_global(Path(stage['initial']));ids=initial['u4if']&ID_MASK
checks={'initial_hash':sha(Path(stage['initial']))==stage['initial_sha256']}
cases={};endpoints={};ledgers={}
for name in stage['order']:
    d=root/name;text=(d/'log').read_text(errors='replace')
    files=[p for p in d.glob('khout.*.dat') if p.name!=Path(stage['initial']).name]
    if len(files)!=1:raise ValueError('expected one endpoint in '+name)
    endpoint=files[0];p,meta=ordered_global(endpoint,ids);q=payload(endpoint,len(p));endpoints[name]=endpoint
    records=[]
    for line in text.splitlines():
        for kind in ('[RK4E]','[WALL_ACCEPT]'):
            if kind in line:
                # Remove any unrelated log appended to the numeric record.
                fragment=line[line.index(kind):]
                fragment=re.split(r'(?=\[(?:CPU|GPU|UPDV)\b)',fragment,maxsplit=1)[0].strip()
                records.append(fragment)
    ledgers[name]=records
    rows=[]
    for line in text.splitlines():
        if '[RK4 gap]' in line:
            m=re.search(r'\bwall=([0-9.]+)s(?:\s|$)',line)
            if not m:raise ValueError('malformed RK4 timing')
            rows.append(float(m[1]))
    if len(rows)<10:raise ValueError('insufficient timing samples')
    injection=all(ledger_number(line,k)==0 for line in records if line.startswith('[RK4E]') for k in ('floor_cum','de_cum','sfl_cum'))
    expected_flag='1' if name.startswith('B') else '0'
    script=(d/'run.sh').read_text()
    checks[name+'_contract']=(sha(d/'eunha2.exe')==stage['binary_sha256'] and (d/'exit_status.txt').read_text().strip()=='0'
        and meta['t']==stage['target_time'] and np.array_equal(initial['mass'],p['mass']) and np.all(p['w2']==0)
        and all(np.isfinite(p[k]).all() and np.all(p[k]>0) for k in ('mass','volume','ie','den','pressure'))
        and all(np.isfinite(p[k]).all() for k in ('x','y','vx','vy'))
        and np.all(q['tag']==TAG) and np.all(np.isfinite(q['u'])) and np.all(q['u']>0)
        and np.all(np.isfinite(q['v'])) and np.all(q['v']>0)
        and np.array_equal(p['ie'],q['u'].astype(np.float32)) and np.array_equal(p['volume'],q['v'].astype(np.float32))
        and injection and len([r for r in records if r.startswith('[RK4E]')])==len(rows)
        and len([r for r in records if r.startswith('[WALL_ACCEPT]')])==len(rows)
        and '[WALL_EVENT_STOP]' not in text and 'GFS_WALL_OMP_VALIDATE='+expected_flag in script
        and 'GFS_THERMAL64=1' in script)
    cases[name]={'endpoint':str(endpoint),'endpoint_sha256':sha(endpoint),'steps':len(rows),
                 'median_step_seconds_after_two_warmups':statistics.median(rows[2:]),
                 'native_exec_wall_seconds':float((d/'elapsed_seconds.txt').read_text().strip())}
parity={name:compare(endpoints['A1'],endpoints[name]) for name in stage['order'][1:]}
checks['all_full_record_footer_parities']=all(v['passed'] for v in parity.values())
checks['all_numeric_ledgers_equal']=all(ledgers[name]==ledgers['A1'] for name in stage['order'][1:])
a=statistics.median(cases[n]['median_step_seconds_after_two_warmups'] for n in ('A1','A2'))
b=statistics.median(cases[n]['median_step_seconds_after_two_warmups'] for n in ('B1','B2'))
speed=a/b;checks['speedup_at_least_1p05']=speed>=stage['minimum_speedup']
result={'passed':all(checks.values()),'checks':{k:bool(v) for k,v in checks.items()},'speedup':speed,'cases':cases,'parity':parity,
        'scope':'H200 same-binary native Q1 thermal64 W0 short ABBA; one allocated GPU, node may be shared. This qualifies only the measured geometry and guard, not MPI, long-time stability or AREPO-comparable throughput.'}
with (root/'abba_audit.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
if not result['passed']:raise SystemExit(1)
