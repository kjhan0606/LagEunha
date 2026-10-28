"""Compare captured raw18 node and ordered cell fluxes on a compute host."""
import hashlib
import json
from pathlib import Path
import numpy as np

base = Path('/gpfs/kjhan/LagForce')
capture = base/'Codex_kinetic_face_capture_20261018'
root = base/'Codex_kinetic_face_backends_20261018'
a = np.fromfile(capture/'inputs.bin', dtype='<f8').reshape(-1, 33)
expected = np.fromfile(capture/'expected.bin', dtype='<f8').reshape(-1, 18)
owners = np.load(capture/'owners.npy')
cell_expected = np.load(capture/'python_cell_flux.npy')
# Sum of absolute elementary terms: protects cancellation without using an
# arbitrary absolute floor. 64 eps covers this fixed small algebra expression.
x = abs(a[:, :2]); n = abs(a[:, 2:4]); scale = a[:, 4]
vl = a[:, 5:13].reshape(-1, 4, 2)
vr = a[:, 13:21].reshape(-1, 4, 2)
wl = abs(a[:, 21:25]); wr = abs(a[:, 25:29])
fm = wl.sum(axis=1)+wr.sum(axis=1)
fp = np.einsum('nk,nkj->nj', wl, abs(vl))+np.einsum('nk,nkj->nj', wr, abs(vr))+abs(a[:, 31, None])*n
pi = np.einsum('nk,nki,nkj->nij', wl, abs(vl), abs(vl))+np.einsum('nk,nki,nkj->nij', wr, abs(vr), abs(vr))
energy = (wl*(.5*(vl*vl).sum(axis=2)+abs(a[:, 29, None]))).sum(axis=1)+(wr*(.5*(vr*vr).sum(axis=2)+abs(a[:, 30, None]))).sum(axis=1)+abs(a[:, 31]*a[:, 32])
bound = np.column_stack([fm, x*fm[:, None], (x[:, :, None]*x[:, None, :]*fm[:, None, None]).reshape(-1, 4), fp, (x[:, :, None]*fp[:, None, :]).reshape(-1, 4), pi.reshape(-1, 4), energy])*scale[:, None]
eps = np.finfo(float).eps
node_tol = 64*eps*bound
cell_bound = np.zeros_like(cell_expected)
cell_node_tol = np.zeros_like(cell_expected)
counts = np.zeros(len(cell_expected), dtype=int)
for (i,j), b, t in zip(owners, bound, node_tol):
    cell_bound[i] += b; cell_node_tol[i] += t; counts[i] += 1
    if j >= 0:
        cell_bound[j] += b; cell_node_tol[j] += t; counts[j] += 1
gamma = counts*eps/(1-counts*eps)
cell_tol = cell_node_tol+2*gamma[:, None]*cell_bound
report = {'nodes':len(a), 'cells':len(cell_expected), 'tolerance':'node:64 eps sum(abs elementary terms); cell:sum(node tolerance)+2 gamma_count sum(node absolute bounds)', 'backends':{}}
arrays = {}
for name in ('cpu1','omp4','cuda'):
    path = root/(name+'.bin')
    output = np.fromfile(path, dtype='<f8').reshape(-1, 19)
    assert output.shape == (len(a), 19) and np.isfinite(output).all()
    assert not np.any(output[:, 18])
    flux = output[:, :18]; delta = abs(flux-expected)
    assert np.all(delta <= node_tol), (name, 'node', float(delta.max()))
    assert np.array_equal(flux[:,4], flux[:,5]) and np.array_equal(flux[:,14], flux[:,15])
    cells = np.zeros_like(cell_expected)
    for (i,j), value in zip(owners, flux):
        cells[i] -= value
        if j >= 0: cells[j] += value
    cell_delta = abs(cells-cell_expected)
    assert np.all(cell_delta <= cell_tol), (name, 'cell', float(cell_delta.max()))
    ratio = np.divide(delta, node_tol, out=np.zeros_like(delta), where=node_tol>0)
    cell_ratio = np.divide(cell_delta, cell_tol, out=np.zeros_like(cell_delta), where=cell_tol>0)
    report['backends'][name] = dict(passed=True, node_max_abs=float(delta.max()), node_max_tolerance_fraction=float(ratio.max()), cell_max_abs=float(cell_delta.max()), cell_max_tolerance_fraction=float(cell_ratio.max()), sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    arrays[name] = output
report['cpu_omp_bitwise_equal'] = np.array_equal(arrays['cpu1'].view(np.uint64), arrays['omp4'].view(np.uint64))
report['cpu_cuda_bitwise_equal'] = np.array_equal(arrays['cpu1'].view(np.uint64), arrays['cuda'].view(np.uint64))
report['scope'] = 'Node algebra and ordered cell face-flux assembly only; excludes geometry, reconstruction, massroot, pressure Pi/volume sources, gravity, MPI, evolution and Psi.'
report['provenance'] = {name:hashlib.sha256((root/name).read_bytes()).hexdigest() for name in ('moment_kinetic_face_harness.cu','moment_kinetic_face_flux.h','moment_state_layout.h','face_harness')}
(root/'comparison.json').write_text(json.dumps(report, indent=2))
print(json.dumps(report, indent=2))
