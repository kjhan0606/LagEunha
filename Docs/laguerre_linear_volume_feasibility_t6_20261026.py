"""One-time linearized Laguerre-volume feasibility probe at exact Kepler t=6.

This solves the frozen-topology power-cell volume Jacobian. It does not evolve
the fluid, construct the resulting nonlinear power diagram, or validate a
conservative weighted-face scheme.
"""

from pathlib import Path
import json
import sys
import time
import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import LinearOperator, cg
from scipy.spatial import Voronoi

ROOT = Path('/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020')
sys.path.insert(0, str(ROOT / 'rk4_diagnostics/varpress_w0_geom_cfl_long_20261026'))
from audit_varpress_checkpoint_20261026 import load_complete

p0, m0 = load_complete(ROOT / 'core_refinement/fine_integrated/khout.000001.dat')
if abs(float(m0['t'])) > 1e-12:
    raise ValueError('Initial snapshot is not t=0')
x0 = p0['x'].astype(float) - 12.
y0 = p0['y'].astype(float) - 12.
r0 = np.hypot(x0, y0)
omega = (r0*r0 + .047**2)**(-.75)
theta = 6.*omega
xy = np.column_stack((12. + x0*np.cos(theta) - y0*np.sin(theta),
                      12. + x0*np.sin(theta) + y0*np.cos(theta)))

start = time.monotonic()
vor = Voronoi(xy, qhull_options='Qbb Qc Qz')
active = np.flatnonzero((r0 >= 1.) & (r0 < 8.))
lookup = np.full(len(xy), -1, dtype=int)
lookup[active] = np.arange(len(active))
area = np.empty(len(active))
for k, i in enumerate(active):
    region = vor.regions[vor.point_region[i]]
    if not region or -1 in region:
        raise ValueError(f'Unbounded active cell: {i}')
    v = vor.vertices[region]
    area[k] = .5*abs(np.dot(v[:,0], np.roll(v[:,1],-1)) -
                       np.dot(v[:,1], np.roll(v[:,0],-1)))

target = p0['volume'][active].astype(float)
rhs = target - area
diag = np.zeros(len(active))
row = []
col = []
val = []
edges = []
for (i,j),vertices in zip(vor.ridge_points, vor.ridge_vertices):
    ii, jj = lookup[i], lookup[j]
    if ii < 0 and jj < 0:
        continue
    if -1 in vertices or len(vertices) != 2:
        raise ValueError(f'Non-finite face touching active region: {i}, {j}')
    length = np.linalg.norm(vor.vertices[vertices[0]] - vor.vertices[vertices[1]])
    distance = np.linalg.norm(xy[i] - xy[j])
    conductance = length/(2.*distance)
    if ii >= 0:
        diag[ii] += conductance
    if jj >= 0:
        diag[jj] += conductance
    if ii >= 0 and jj >= 0:
        row.extend((ii,jj))
        col.extend((jj,ii))
        val.extend((-conductance,-conductance))
    edges.append((i,j,distance*distance))
if np.any(diag <= 0):
    raise ValueError('Disconnected active cell')
jac = (coo_matrix((val,(row,col)),shape=(len(active),len(active))).tocsr() +
       diags(diag,format='csr'))
M = LinearOperator(jac.shape, matvec=lambda x:x/diag)
iterations = [0]
solve_start = time.monotonic()
weight, info = cg(jac,rhs,M=M,rtol=1e-10,atol=0.,maxiter=2000,
                  callback=lambda _:iterations.__setitem__(0,iterations[0]+1))
solve_seconds = time.monotonic()-solve_start
if info != 0:
    raise RuntimeError(f'PCG failed: info={info}, iterations={iterations[0]}')
relative_residual = np.linalg.norm(jac@weight-rhs)/np.linalg.norm(rhs)

jacobi = np.zeros(len(active))
jacobi_residual = {}
for step in range(1,129):
    jacobi += .8*(rhs-jac@jacobi)/diag
    if step in (1,2,4,8,16,32,64,128):
        jacobi_residual[step] = float(np.linalg.norm(rhs-jac@jacobi)/np.linalg.norm(rhs))

def weighted_summary(lo,hi):
    ids = np.flatnonzero((r0[active]>=lo)&(r0[active]<hi))
    volume = target[ids]
    return dict(rlo=lo,rhi=hi,cells=len(ids),
        volume_change_MW=float(np.sum(p0['mass'][active[ids]]*abs(rhs[ids])/volume)/np.sum(p0['mass'][active[ids]])),
        abs_weight_over_initial_area_p50_p95_p99=np.quantile(abs(weight[ids])/volume,[.5,.95,.99]).tolist())

face_ratios = np.array([abs((weight[lookup[i]] if lookup[i]>=0 else 0.)-
                             (weight[lookup[j]] if lookup[j]>=0 else 0.))/d2
                        for i,j,d2 in edges])
result = dict(n_all=len(xy),n_active=len(active),n_active_faces=len(edges),
    active_region='1 <= initial r < 8; outside weights fixed to zero',
    target='initial unweighted Voronoi area for every active cell',
    method='frozen-topology linearized power-volume Jacobian J_ii=sum Aij/(2dij), J_ij=-Aij/(2dij)',
    rhs_sum=float(rhs.sum()),rhs_L2=float(np.linalg.norm(rhs)),
    pcg_iterations=iterations[0],pcg_seconds=solve_seconds,
    pcg_relative_residual=float(relative_residual),
    jacobi_relative_residual=jacobi_residual,
    edge_abs_weight_difference_over_d2_quantiles=np.quantile(face_ratios,[.5,.95,.99,1.]).tolist(),
    edge_fraction_over_0p8=float(np.mean(face_ratios>.8)),
    bands=[weighted_summary(lo,hi) for lo,hi in [(1,2),(2,2.25),(2.25,2.5),(2.5,3),(3,4),(4,6),(6,8)]],
    total_wall_seconds=time.monotonic()-start,
    limitation='Linearization only; no nonlinear weighted-tessellation, conservation, or speed-of-evolution result')
output = Path(__file__).with_name('LAGUERRE_LINEAR_VOLUME_FEASIBILITY_T6_2026-10-26.json')
with output.open('x') as f:
    json.dump(result,f,indent=2)
    f.write('\n')
print(json.dumps(result,indent=2))
