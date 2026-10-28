"""Inspect completed C/OpenMP/CUDA face-flux outputs on lageunha."""
import hashlib
import json
from pathlib import Path
import numpy as np
root=Path('/gpfs/kjhan/LagForce/Codex_affine_flux_backend_20261018')
ref=np.fromfile(root/'expected.bin',dtype='<f8').reshape(-1,6)
out={'manifest':json.loads((root/'manifest.json').read_text()),'backends':{},
     'compiler_flags':'nvcc -O2 --fmad=false -Xcompiler=-fopenmp,-ffp-contract=off',
     'FMA_scope':'Default optimized/FMA-enabled build not checked.'}
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
out['inputs_sha256']=sha(root/'inputs.bin');out['binary_sha256']=sha(root/'check.exe')
out['header_sha256']=sha(root.parent/'affine_face_transport.h')
cpu=None
for name in ['cpu1','cpu4','cuda']:
    a=np.fromfile(root/(name+'.bin'),dtype='<f8').reshape(-1,6)
    assert a.shape==ref.shape and np.isfinite(a).all()
    if cpu is None:cpu=a
    err=abs(a-ref);scaled=err/(1+abs(ref))
    assert scaled.max()<1e-12,(name,scaled.max())
    equal=bool(np.array_equal(a.view(np.uint64),cpu.view(np.uint64)))
    assert equal,name
    out['backends'][name]={'bitwise_equal_to_cpu1':equal,
        'max_abs_error_vs_gauss_per_component':err.max(axis=0).tolist(),
        'max_scaled_error_vs_gauss':float(scaled.max())}
with (root/'analysis.json').open('x') as fp:json.dump(out,fp,indent=2);fp.write('\n')
print(json.dumps(out,indent=2))
