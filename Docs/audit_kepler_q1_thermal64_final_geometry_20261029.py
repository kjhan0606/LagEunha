"""Independently reconstruct final bounded Q1 areas for the thermal64 control."""
import hashlib
import json
from pathlib import Path
import numpy as np
from kepler_native_global_ids_20261029 import ordered_global, ID_MASK
from audit_kepler_q1_checkpoint_geometry_20261029 import independent
from audit_kepler_q1_thermal64_promotion_20261029 import payload, TAG

root=Path('/gpfs/kjhan/LagForce/Kepler_Q1_Thermal64Half_H200_20261029')
comparison=json.loads((root/'arithmetic_comparison.json').read_text())
seed=Path(json.loads((root/'contract.json').read_text())['source'])
ends=[p for p in root.glob('khout.*.dat') if p.name!=seed.name]
if len(ends)!=1:raise ValueError('Expected one endpoint')
endpoint=ends[0]
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
if sha(endpoint)!=comparison['thermal64_endpoint_sha256']:raise ValueError('Endpoint changed')
p0,_=ordered_global(seed);p,m=ordered_global(endpoint,p0['u4if']&ID_MASK);q=payload(endpoint,len(p))
if m['t']!=10.1 or not np.all(q['tag']==TAG):raise ValueError('Invalid tagged endpoint')
area,geometry=independent(np.column_stack((p['x'],p['y'])))
error=abs(q['v']/area-1.)
out={'passed':bool(np.max(error)<1e-6),'all_site_count':len(p),'max_relative_error':float(np.max(error)),
     'median_relative_error':float(np.median(error)),'geometry':geometry,'endpoint_sha256':sha(endpoint),
     'scope':'All-site independent fixed-polygon area verification at t=10.1; not a force, entropy or trajectory accuracy proof.'}
path=root/'independent_final_geometry.json'
if path.exists():raise FileExistsError(path)
with path.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
print(json.dumps(out,indent=2))
if not out['passed']:raise SystemExit(1)
