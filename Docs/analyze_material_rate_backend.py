import json,hashlib
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce');root=base/'Codex_material_rate_backend_20261018'
ref=np.fromfile(base/'Codex_pair_hg_modes_t0_20261018/reynolds_source_expected.bin',dtype='<f8').reshape(-1,17)
out={'scope':'Force-free material source only, no complete time evolution.','cells':len(ref),'backends':{}}
for name in ['cpu1','cpu4','cuda']:
    a=np.fromfile(root/(name+'.bin'),dtype='<f8').reshape(-1,17)
    assert a.shape==ref.shape and np.isfinite(a).all()
    err=float(abs(a-ref).max());assert err==0
    out['backends'][name]={'max_abs_error_vs_reference':err}
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
out['header_sha256']=sha(base/'moment_material_drift.h');out['binary_sha256']=sha(root/'check.exe')
with (root/'analysis.json').open('x') as f:json.dump(out,f,indent=2);f.write('\n')
print(json.dumps(out,indent=2))
