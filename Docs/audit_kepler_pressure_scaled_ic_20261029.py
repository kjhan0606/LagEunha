"""Validate a pressure-scaled diagnostic restart without altering its source."""
import hashlib
import json
import struct
import sys
from pathlib import Path
import numpy as np
from kepler_native_global_ids_20261029 import ordered_global,ID_MASK

source,x1,x4,out=map(Path,sys.argv[1:5])
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
expected='12ea9ff96ec70c4816f834618229c972d0a8e47814bcdc2a38a20a8e89aad50e'
if sha(source)!=expected or sha(x1)!=expected:raise ValueError('source or x1 mismatch')
a,ma=ordered_global(source);b,mb=ordered_global(x4,a['u4if']&ID_MASK)
if ma['t']!=30 or mb!=ma or len(a)!=55351:raise ValueError('header or time mismatch')
for name,factor in [('ie',4),('pressure',4),('mass',1),('x',1),('y',1),('vx',1),('vy',1),('den',1),('volume',1),('w2',1)]:
    if not np.array_equal(b[name],a[name]*factor):raise ValueError('unexpected '+name)
if np.any(b['ie']<=0) or np.any(b['pressure']<=0) or np.any(b['w2']!=0):raise ValueError('invalid physics')
offsets=[(108,4),(116,4),(124,4),(392,8),(400,8),(552,8)]
mask=np.ones(576,dtype=bool)
for offset,size in offsets:mask[offset:offset+size]=False
ra=np.memmap(source,mode='r',dtype=np.uint8,offset=16,shape=(len(a),576))
rb=np.memmap(x4,mode='r',dtype=np.uint8,offset=16,shape=(len(a),576))
if not np.array_equal(ra[:,mask],rb[:,mask]):raise ValueError('non-thermal record bytes changed')
gamma=float(np.float32(1.666667));P=(gamma-1)*b['ie'].astype(float)/b['volume'].astype(float)
error=abs(P/b['pressure'].astype(float)-1)
result={'source_sha256':expected,'x1_sha256':sha(x1),'x4_sha256':sha(x4),'time':ma['t'],
        'cells':len(a),'thermal_scaling':4,'eos_relative_error_max':float(error.max()),
        'eos_relative_error_p99':float(np.quantile(error,.99)),
        'changed_byte_ranges':offsets,'script_sha256':sha(Path(__file__)),
        'scope':'Counterfactual pressure-only diagnostic IC created from t30, not continuation of same physical solution; no evolved output yet.'}
with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
