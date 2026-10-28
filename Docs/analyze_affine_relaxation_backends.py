"""Read completed face-local CPU/OpenMP/CUDA relaxation outputs."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
parser=argparse.ArgumentParser();parser.add_argument('--root',type=Path,default=Path('/gpfs/kjhan/LagForce/Codex_affine_relaxation_backends_20261018'))
root=parser.parse_args().root
ref=np.fromfile(root/'expected.bin',dtype='<f8').reshape(-1,35)
manifest=json.loads((root/'manifest.json').read_text());valid=ref[:,0]==0
cpu=None;results={}
for name in ['cpu1','cpu4','cuda']:
    a=np.fromfile(root/(name+'.bin'),dtype='<f8').reshape(-1,35)
    assert a.shape==ref.shape and np.isfinite(a).all()
    assert np.array_equal(a[:,0],ref[:,0]) and np.array_equal(a[~valid],ref[~valid])
    if cpu is None:cpu=a
    maximum=float(np.max(abs(a[valid,1:]-ref[valid,1:])))
    assert maximum<1e-12
    if name=='cpu4':assert np.array_equal(a,cpu)
    da=a[valid,1:18];db=a[valid,18:35]
    energy=float(np.max(abs(da[:,15]+db[:,15])))
    angular=float(np.max(abs(da[:,9]-da[:,10]+db[:,9]-db[:,10])))
    heat=np.r_[da[:,15]-.5*(da[:,12]+da[:,14]),db[:,15]-.5*(db[:,12]+db[:,14])]
    assert energy<1e-12 and angular<1e-12 and heat.min()>-1e-14
    results[name]=dict(max_absolute_reference_error=maximum,
        max_absolute_cpu_difference=float(np.max(abs(a-cpu))),bitwise_equal_cpu=bool(np.array_equal(a.view(np.uint64),cpu.view(np.uint64))),
        max_pair_energy_rate_residual=energy,max_pair_angular_rate_residual=angular,minimum_thermal_rate=float(heat.min()))
sha=lambda name:hashlib.sha256((root/name).read_bytes()).hexdigest()
out=dict(scope=__doc__,valid_faces=int(valid.sum()),guard_cases=int((~valid).sum()),backends=results,
    hashes={name:sha(name) for name in ['inputs.bin','expected.bin','check.exe','affine_relaxation_rate.h']},
    compiler_flags='nvcc -O2 -std=c++14 -arch=sm_89 -Xcompiler=-fopenmp; default FMA',
    limits=['Face-local kernels only, no GPU scatter/assembly or whole RK4 evolution.',
            'Manufactured affine/curved traces on real patch faces; not long-time Kepler.',
            'MPI ownership and transport are untested for this operator.'])
(root/'analysis.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
