"""Separate gravity and moving circular-reference terms of response power.

Together with saved pressure power these are part, not all, of d(half the
mass-weighted squared departure response)/dt. TCF/HG are not included.
"""
import argparse
import json
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, '/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read
from audit_canonical_pair_history import gravity_coefficient

parser = argparse.ArgumentParser()
parser.add_argument('--step', type=int, required=True)
args = parser.parse_args()
b = Path('/gpfs/kjhan/LagForce')
root = b/'Codex_coupled_onset_to_t4_20261018'
pressure = json.loads((root/f'pressure_response_{args.step:06d}.json').read_text())
ic, _ = read(b/'Codex_canonical_initial_20261017/run/khout.000000.dat')
states = [read(root/n/f'khout.{args.step:06d}.dat') for n in ['plus', 'minus']]
eps = pressure['epsilon']

def circular(x):
    r = np.linalg.norm(x, axis=1)
    return np.sqrt(gravity_coefficient(r))[:, None]*np.column_stack((-x[:, 1], x[:, 0]))

def terms(a):
    x = np.column_stack((a['x']-3, a['y']-3))
    v = np.column_stack((a['vx'], a['vy']))
    r = np.linalg.norm(x, axis=1)
    assert np.all(r > 0)
    k = gravity_coefficient(r)
    omega = np.sqrt(k)
    z = r/.5
    low = z <= .5
    mid = (z > .5) & (z < 1)
    dMdz = np.zeros_like(z)
    dMdz[low] = 32*z[low]**2-192*z[low]**4+192*z[low]**5
    dMdz[mid] = 64*z[mid]**2-192*z[mid]**3+192*z[mid]**4-64*z[mid]**5
    M = k*r**3
    slope = .5*(z*dMdz/M-3)  # r*Omega'/Omega
    jx = np.column_stack((-x[:, 1], x[:, 0]))
    jv = np.column_stack((-v[:, 1], v[:, 0]))
    reference_derivative = omega[:, None]*(jv+(slope*np.sum(x*v, axis=1)/r**2)[:, None]*jx)
    # Independent directional finite difference checks the analytic Jacobian.
    errors = []
    for h in [1e-5, 5e-6]:
        fd = (circular(x+h*v)-circular(x-h*v))/(2*h)
        errors.append(float(np.linalg.norm(fd-reference_derivative)/np.linalg.norm(reference_derivative)))
    assert errors[-1] < 1e-7
    return v-circular(x), -k[:, None]*x, reference_derivative, errors

for a,t in states:
    assert np.array_equal(a['id'], ic['id'])
    assert np.array_equal(a['mass'], ic['mass'])
    assert abs(t-pressure['time']) < 1e-12
p, n = [terms(a) for a,t in states]
du, dg, dr = [(p[i]-n[i])/(2*eps) for i in range(3)]
mass = ic['mass']
r0 = np.hypot(ic['x']-3, ic['y']-3)
masks = {'inner':r0<.5, 'shear_band':(r0>=.3)&(r0<.5),
         'central':r0<.2, 'edge':(r0>=.5)&(r0<.6), 'core':(r0>=.6)&(r0<1.9)}
rows = {}
for name, mask in masks.items():
    gravity = float(np.sum(mass[mask, None]*du[mask]*dg[mask]))
    reference = float(-np.sum(mass[mask, None]*du[mask]*dr[mask]))
    hydro = pressure['regions'][name]['actual']['departure_response_power']
    rows[name] = {'gravity':gravity, 'circular_reference_derivative':reference,
                  'gravity_plus_reference':gravity+reference,
                  'reversible_pressure':hydro,
                  'subtotal_without_TCF_HG':gravity+reference+hydro}
out = {'time':pressure['time'], 'epsilon':eps, 'regions':rows,
       'analytic_reference_Jacobian_relative_L2_checks':{'plus':p[3], 'minus':n[3]},
       'scope':'Finite central response on evolving hydro trajectories. Subtotal excludes TCF/HG; not net growth or a closed budget. Fixed initial material bins.'}
with (root/f'kinematic_response_{args.step:06d}.json').open('x') as f:
    json.dump(out, f, indent=2)
print(json.dumps(out, indent=2))
