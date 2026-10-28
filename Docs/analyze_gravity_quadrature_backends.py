"""Compare actual CPU/OpenMP/CUDA with independent long-double reductions."""
import hashlib,json
from pathlib import Path
import numpy as np
root=Path('/gpfs/kjhan/LagForce/Codex_gravity_quadrature_backends_20261018')
samples=np.memmap(root/'samples.bin',dtype='<f8',mode='r').reshape(-1,7);offset=np.fromfile(root/'offsets.bin',dtype='<u8')
n=len(offset)-1;valid=n-2;refs=[]
for i in range(valid):
    r=np.asarray(samples[int(offset[i]):int(offset[i+1])],dtype=np.longdouble);x=r[:,:2];w=r[:,2];v=r[:,3:5];a=r[:,5:7]
    force=np.sum(w[:,None]*a,axis=0);D=x.T@(w[:,None]*a);VA=v.T@(w[:,None]*a);Pi=VA+VA.T
    refs.append(np.r_[force,D.ravel(),Pi[0,0],Pi[0,1],Pi[1,1],np.sum(w*np.sum(v*a,axis=1))])
reference=np.asarray(refs,dtype=float)
out={name:np.fromfile(root/f'{name}.bin',dtype='<f8').reshape(n,11) for name in ('cpu1','omp4','cuda')}
assert np.array_equal(out['cpu1'],out['omp4'])
rows=[]
for name,a in out.items():
    assert np.all(a[:valid,0]==0) and np.all(a[valid:,0]==1) and np.all(a[valid:,1:]==12345.)
    error=abs(a[:valid,1:]-reference);scaled=error/(1e-13+1e-10*abs(reference))
    assert np.isfinite(a).all() and np.max(scaled)<1
    rows.append(dict(backend=name,maximum_absolute_error=float(error.max()),maximum_tolerance_fraction=float(scaled.max()),
        maximum_thermal_source_error=float(np.max(abs(a[:valid,10]-.5*(a[:valid,7]+a[:valid,9])))),
        maximum_central_torque=float(np.max(abs(a[:valid,4]-a[:valid,5])))))
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return h.hexdigest()
report=dict(valid_cells=valid,invalid_cases=2,quadrature_samples=int(offset[valid]),cpu_openmp_bitwise=True,
    maximum_cpu_cuda_difference=float(np.max(abs(out['cpu1']-out['cuda']))),results=rows,
    sha256={p:sha(root/p) for p in ('samples.bin','offsets.bin','moment_gravity_quadrature.h','probe')},
    limits='Shared force/moment quadrature accumulation only. Acceleration, density, velocity reconstruction and nodes prepared on CPU. No full GPU hydro, MPI, RK4, or long-term stability claim.')
(root/'backend_report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
