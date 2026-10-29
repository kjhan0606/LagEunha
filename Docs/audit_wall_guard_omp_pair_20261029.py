"""Check serial/OMP native guard candidates on the same checkpoint."""
import hashlib
import json
import re
import struct
import sys
from pathlib import Path
import numpy as np

root=Path(sys.argv[1]); retry=len(sys.argv)>2 and sys.argv[2]=='retry'
labels=('serial_retry','omp_retry') if retry else ('serial','omp')
step='000023' if retry else '000026'
paths=[root/label/('khout.'+step+'.dat') for label in labels]
def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def read(p):
    raw=p.read_bytes();header=raw[:16];magic,n,stride,stress=struct.unpack('<4i',header)
    if (magic,stride,stress)!=(-1,576,1) or len(raw)!=16+n*stride+52:raise ValueError('incomplete restart')
    return header,np.frombuffer(raw,dtype=np.uint8,count=n*stride,offset=16).reshape(n,stride),raw[-52:]
a,b=map(read,paths)
checks={'same_header':a[0]==b[0],'same_records_without_next_pointer':np.array_equal(a[1][:,8:],b[1][:,8:]),'same_footer':a[2]==b[2]}
for label in labels:
    checks[label+'_exit_zero']=(root/label/'exit_status.txt').read_text().strip()=='0'
    text=(root/label/'log').read_text(errors='replace')
    checks[label+'_thermal_and_cuda']='[FIXED_POLYGON] enabled' in text and '[GPU tess]' in text and '[WALL_ACCEPT]' in text
    checks[label+'_accepted_steps']=len(re.findall(r'^\[WALL_ACCEPT\]',text,re.M))==(1 if retry else 4)
    if retry:checks[label+'_exact_rollback']='payload_counts_dtold_ledgers_exact=1 device_cache_invalidated=1' in text
def ledgers(label):
    text=(root/label/'log').read_text(errors='replace')
    return [line for line in text.splitlines() if line.startswith(('[RK4E]','[WALL_ACCEPT]'))]
checks['numeric_ledgers_exact']=ledgers(labels[0])==ledgers(labels[1])
result={'passed':all(checks.values()),'checks':{k:bool(v) for k,v in checks.items()},
        'binary_sha256':digest(root/'project/eunha2.exe'),'input_sha256':digest(root/'serial/khout.000022.dat'),
        'output_sha256':[digest(x) for x in paths],'script_sha256':digest(Path(__file__)),
        'scope':('One-step forced rollback' if retry else 'Four-step no-impact')+' native Ada same-binary serial/OMP comparison on shared GPU. Excludes runtime next pointers, compares all other records, footer and printed numerical ledgers. Not a benchmark, physical wall-impact, multi-rank or long-time qualification.'}
with (root/('retry_audit.json' if retry else 'pair_audit_v2.json')).open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
if not result['passed']:raise SystemExit(1)
