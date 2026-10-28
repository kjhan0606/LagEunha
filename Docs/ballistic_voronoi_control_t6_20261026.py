"""One-off diagnostic: compare native t=6 density with exact-orbit Voronoi remeshing.

This does not evolve the fluid. It recomputes polygon areas once from saved
generators and once from their analytic circular-orbit positions.
"""

from pathlib import Path
import argparse
import json
import sys
import numpy as np
from scipy.spatial import Voronoi

ROOT = Path('/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020')
sys.path.insert(0, str(ROOT / 'rk4_diagnostics/varpress_w0_geom_cfl_long_20261026'))
from audit_varpress_checkpoint_20261026 import load_complete

INITIAL = ROOT / 'core_refinement/fine_integrated/khout.000001.dat'
parser = argparse.ArgumentParser()
parser.add_argument('--snapshot', type=Path, default=ROOT / 'rk4_diagnostics/varpress_w0_geom_cfl_long_20261026/C010/to6/khout.032479.dat')
parser.add_argument('--output', type=Path, default=Path(__file__).with_name('BALLISTIC_VORONOI_CONTROL_T6_2026-10-26.json'))
args = parser.parse_args()
ACTUAL = args.snapshot
OUTPUT = args.output

p0, m0 = load_complete(INITIAL)
p6, m6 = load_complete(ACTUAL)
if not np.array_equal(p0['u4if'], p6['u4if']):
    raise ValueError('Particle IDs differ')
if abs(float(m0['t'])) > 1e-12 or float(m6['t']) <= 0:
    raise ValueError('Unexpected times')
t = float(m6['t'])

x0 = p0['x'].astype(float) - 12.
y0 = p0['y'].astype(float) - 12.
r0 = np.hypot(x0, y0)
omega = (r0*r0 + .047**2)**(-.75)
theta = omega * t
ballistic = np.column_stack((12. + x0*np.cos(theta) - y0*np.sin(theta),
                             12. + x0*np.sin(theta) + y0*np.cos(theta)))
initial = np.column_stack((p0['x'], p0['y']))
actual = np.column_stack((p6['x'], p6['y']))
interior = np.flatnonzero((r0 >= 1.) & (r0 < 8.))

def polygon_areas(points):
    vor = Voronoi(points, qhull_options='Qbb Qc Qz')
    area = np.full(len(points), np.nan)
    for i in interior:
        region = vor.regions[vor.point_region[i]]
        if not region or -1 in region:
            raise ValueError(f'Unbounded interior cell {i}')
        xy = vor.vertices[region]
        area[i] = .5*abs(np.dot(xy[:, 0], np.roll(xy[:, 1], -1)) -
                            np.dot(xy[:, 1], np.roll(xy[:, 0], -1)))
    return area

a0 = polygon_areas(initial)
a6 = polygon_areas(actual)
ab = polygon_areas(ballistic)
mass = p0['mass'].astype(float)
rho0 = p0['den'].astype(float)
rho6 = p6['den'].astype(float)
rhob = mass / ab
gamma = float(np.float32(1.666667))
pb = p0['pressure'].astype(float) * (rhob/rho0)**gamma

def summary(lo, hi):
    s = (r0 >= lo) & (r0 < hi)
    w = mass[s] / mass[s].sum()
    mw = lambda a: float(np.sum(w*a[s]))
    d_actual = rho6/rho0-1.
    d_ball = rhob/rho0-1.
    d_pressure = p6['pressure'].astype(float)/p0['pressure'].astype(float)-1.
    d_pressure_ball = pb/p0['pressure'].astype(float)-1.
    corr = np.corrcoef(d_actual[s], d_ball[s])[0, 1]
    return dict(rlo=lo, rhi=hi, cells=int(s.sum()),
        native_volume_vs_scipy_t0_MW=mw(abs(p0['volume'].astype(float)/a0-1.)),
        native_volume_vs_scipy_t6_MW=mw(abs(p6['volume'].astype(float)/a6-1.)),
        actual_density_MW=mw(abs(d_actual)), ballistic_density_MW=mw(abs(d_ball)),
        density_residual_MW=mw(abs(d_actual-d_ball)),
        density_signed_correlation=float(corr),
        actual_pressure_MW=mw(abs(d_pressure)),
        ballistic_adiabatic_pressure_MW=mw(abs(d_pressure_ball)),
        pressure_residual_MW=mw(abs(d_pressure-d_pressure_ball)))

result = dict(initial=str(INITIAL), actual=str(ACTUAL), t=t,
    method='SciPy unweighted Voronoi polygon area, selected interior cells only; exact softened circular generator orbits',
    bins=[summary(lo,hi) for lo,hi in [(1,2),(2,2.25),(2.25,2.5),(2.5,3),(3,4),(4,6),(6,8)]])
with OUTPUT.open('x') as f:
    json.dump(result,f,indent=2)
    f.write('\n')
for b in result['bins']:
    print(f"{b['rlo']:g}-{b['rhi']:g} N={b['cells']} native/Scipy0={b['native_volume_vs_scipy_t0_MW']:.3g} "
          f"native/Scipy6={b['native_volume_vs_scipy_t6_MW']:.3g} "
          f"rho_actual={b['actual_density_MW']:.4g} rho_ball={b['ballistic_density_MW']:.4g} "
          f"rho_res={b['density_residual_MW']:.4g} corr={b['density_signed_correlation']:.4g} "
          f"P_res={b['pressure_residual_MW']:.4g}")
print(OUTPUT)
