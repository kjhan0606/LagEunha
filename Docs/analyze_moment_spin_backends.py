"""CPU/OpenMP/CUDA spin ledger component parity, not full solver validation."""
import hashlib,json
from pathlib import Path
import numpy as np
root=Path('/gpfs/kjhan/LagForce/Codex_moment_spin_backends_20261018')
meta=json.loads((root/'cases.json').read_text());n=meta['valid'];ref=np.load(root/'reference.npy')
outputs={k:np.fromfile(root/(k+'.bin'),dtype='<f8').reshape(-1,9) for k in ('cpu1','cpu4','cuda')}
assert np.array_equal(outputs['cpu1'],outputs['cpu4'])
reports={}
for name,x in outputs.items():
    assert np.all(x[:n,0]==0) and np.isfinite(x).all()
    for j,g in enumerate(meta['guards']):assert x[n+j,0]==g['code'] and np.all(x[n+j,1:]==12345.)
    error=abs(x[:n,1:]-ref);scaled=error/(1e-11+1e-9*abs(ref))
    assert scaled.max()<1
    reports[name]=dict(maximum_absolute_error_by_field=error.max(0).tolist(),maximum_tolerance_fraction=float(scaled.max()))
out=dict(scope=__doc__,valid_cells=n,guards=meta['guards'],CPU1_OMP4_bitwise_equal=True,
    CUDA_CPU_max_absolute_difference=float(np.max(abs(outputs['cuda']-outputs['cpu1']))),backends=reports,
    fields=['I','L_internal','omega','Krot','Idot','Ldot','Krotdot','Udot'],
    sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [root/'input.bin',root/'moment_spin_ledger.h',root/'probe']},
    limits='Chord direction between two remap states is not a physical time derivative. No MPI, force scatter, or full RK4 validation.')
(root/'analysis.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
