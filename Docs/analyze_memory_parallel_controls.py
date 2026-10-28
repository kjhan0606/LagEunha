"""Actual saved RK4 state comparison; no whole-file pointer-field parity claim."""
import json,struct,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read
b=Path('/gpfs/kjhan/LagForce/Codex_memory_parallel_controls_20261018')
manifest=json.loads((b/'manifest.json').read_text());states={};aux={};cases={}
for name in manifest['cases']:
    d=b/name;assert (d/'exit.code').read_text().strip()=='0'
    p=d/'khout.000162.dat';a,t=read(p);states[name]=(a,t)
    data=Path(str(p)+'.rmem').read_bytes();assert struct.unpack('<d',data[48:56])[0]==t
    rec=np.frombuffer(data,dtype=[('id','<u8'),('x','<f8'),('y','<f8')],offset=56)
    assert np.array_equal(rec['id'],a['id']&np.uint64(72057594037927935))
    assert np.isfinite(rec['x']).all() and np.isfinite(rec['y']).all();aux[name]=rec
    lines=(d/'log').read_text().splitlines();accepted=[q for q in lines if q.startswith('[ACCEPTED_STATE]')]
    assert len(accepted)==1 and ' finite=1 ' in accepted[0] and ' retries=0 ' in accepted[0]
    cases[name]={'time':t,'wall_seconds':float((d/'wall_seconds').read_text()),'accepted':accepted[0]}
comparisons={};ok=True
for name in ('new_serial','parallel1','parallel4'):
    a,t=states['old'];c,tc=states[name];assert t==tc
    equal={key:bool(np.array_equal(a[key],c[key])) for key in a.dtype.names}
    zequal={key:bool(np.array_equal(aux['old'][key],aux[name][key])) for key in aux['old'].dtype.names}
    difference={key:float(np.max(abs(a[key].astype(float)-c[key].astype(float)))) for key in a.dtype.names if key!='id'}
    comparisons[name]={'physical_fields_bitwise_equal':equal,'auxiliary_fields_bitwise_equal':zequal,'physical_max_absolute_difference':difference}
    ok &= all(equal.values()) and all(zequal.values())
out={'cases':cases,'comparisons':comparisons,'all_bitwise_equal':ok,
    'scope':'Actual nonzero-memory restart and one full RK4 step with production neighbor search, row build and ordered force scatter. Short backend control, not long-time physical stability. Wall times include startup and concurrent load; not an isolated performance benchmark.'}
with (b/'analysis.json').open('x') as f:json.dump(out,f,indent=2)
print(json.dumps(out,indent=2));assert ok
