"""Native fixed-mass Kepler gate at t=10 for the 54k-site profile IC."""

import argparse
import json
from pathlib import Path
import re
import sys

import numpy as np

sys.path.insert(0, '/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020')
from audit_snapshot import read_snapshot

N = 54312
EPS = .01175
CX = CY = 3.


def ordered(path):
    data, meta = read_snapshot(path)
    ids = (data['u4if'].astype(np.uint64) & np.uint64((1 << 56)-1)).astype(np.int64)
    if len(data) != N or not np.array_equal(np.sort(ids), np.arange(N)):
        raise ValueError(f'Invalid site IDs in {path}')
    return data[np.argsort(ids)], meta


def mean_abs(values, weights):
    return float(np.sum(abs(values)*weights, dtype=np.longdouble)/
                 np.sum(weights, dtype=np.longdouble))


def energy(p):
    x = p['x'].astype(float)-CX
    y = p['y'].astype(float)-CY
    r2 = x*x+y*y
    m = p['mass'].astype(float)
    kinetic = .5*np.sum(m*(p['vx'].astype(float)**2+p['vy'].astype(float)**2), dtype=np.longdouble)
    thermal = np.sum(p['ie'].astype(float), dtype=np.longdouble)
    potential = -np.sum(m/np.sqrt(r2+EPS*EPS), dtype=np.longdouble)
    return float(kinetic+thermal+potential)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('initial', type=Path)
    parser.add_argument('run', type=Path)
    args = parser.parse_args()
    p0, m0 = ordered(args.initial)
    if float(m0['t']) != 0:
        raise ValueError('Expected exact initial snapshot')
    matches = []
    for path in args.run.glob('khout.*.dat'):
        p, meta = ordered(path)
        if abs(float(meta['t'])-10.) < 1e-12:
            matches.append((path, p, meta))
    if len(matches) != 1:
        raise ValueError('Missing unique exact t=10 snapshot')
    path, p, meta = matches[0]
    if not (np.all(p0['ie']>0) and np.all(p['ie']>0) and
            np.all(p0['volume']>0) and np.all(p['volume']>0)):
        raise ValueError('Nonpositive thermodynamic state')
    x0 = p0['x'].astype(float)-CX
    y0 = p0['y'].astype(float)-CY
    r0 = np.hypot(x0, y0)
    disk = (r0 >= .5)&(r0 < 2.)
    mass = p0['mass'].astype(float)
    phase = 10.*(r0*r0+EPS*EPS)**(-.75)
    xb = CX+x0*np.cos(phase)-y0*np.sin(phase)
    yb = CY+x0*np.sin(phase)+y0*np.cos(phase)
    position_error = np.hypot(p['x']-xb, p['y']-yb)
    density_change = p['den'].astype(float)/p0['den'].astype(float)-1.
    gm1 = float(np.float32(1.666667))-1.
    entropy_log = (np.log(p['ie'].astype(float)/p0['ie'].astype(float))+
                   gm1*np.log(p['volume'].astype(float)/p0['volume'].astype(float)))
    x = p['x'].astype(float)-CX
    y = p['y'].astype(float)-CY
    r = np.hypot(x, y)
    annulus = (r>.8)&(r<1.2)
    if not np.any(annulus):
        raise ValueError('Empty GIZMO reference annulus')
    vphi = (x[annulus]*p['vy'][annulus]-y[annulus]*p['vx'][annulus])/r[annulus]
    vk = r[annulus]*(r[annulus]**2+EPS**2)**(-.75)
    l1_v = float(np.mean(abs(vphi-vk)/vk))
    l1_rho = float(np.mean(abs(p['den'][annulus].astype(float)-1.)))
    l1_rho_plateau = float(np.mean(abs(p['den'][annulus].astype(float)-1.01)/1.01))
    L0 = float(np.sum(mass[disk]*(x0[disk]*p0['vy'][disk]-y0[disk]*p0['vx'][disk]),
                      dtype=np.longdouble))
    Lt = float(np.sum(mass[disk]*(x[disk]*p['vy'][disk]-y[disk]*p['vx'][disk]),
                      dtype=np.longdouble))
    E0, Et = energy(p0), energy(p)
    log = (args.run/'log').read_text(errors='replace')
    lines = [line for line in log.splitlines() if '[RK4E]' in line]
    if not lines:
        raise ValueError('Missing native RK4 energy ledger')
    last = lines[-1]
    def number(key):
        match = re.search(rf'\b{re.escape(key)}=([0-9.eE+-]+)', last)
        if not match:
            raise ValueError(f'Missing {key} in final RK4 ledger')
        return float(match.group(1))
    floor_cum, de_cum, sfl_cum = number('floor_cum'), number('de_cum'), number('sfl_cum')
    disk_density = mean_abs(density_change[disk], mass[disk])
    disk_entropy = mean_abs(entropy_log[disk],mass[disk])
    result = {
        'initial': str(args.initial), 'endpoint': str(path), 'time': float(meta['t']),
        'n': N, 'disk_cells': int(np.sum(disk)),
        'disk_density_abs_mass_weighted': disk_density,
        'disk_position_error_abs_mass_weighted': mean_abs(position_error[disk],mass[disk]),
        'disk_entropy_log_signed_mass_weighted': float(np.sum(entropy_log[disk]*mass[disk],dtype=np.longdouble)/np.sum(mass[disk],dtype=np.longdouble)),
        'disk_entropy_log_abs_mass_weighted': disk_entropy,
        'disk_angular_momentum_relative': Lt/L0-1.,
        'total_energy_relative': (Et-E0)/abs(E0),
        'gizmo_annulus_cells': int(np.sum(annulus)),
        'gizmo_L1_velocity': l1_v,
        'gizmo_L1_density_literal': l1_rho,
        'gizmo_L1_density_plateau_corrected': l1_rho_plateau,
        'rk4_final_floor_cum': floor_cum,
        'rk4_final_de_cum': de_cum,
        'rk4_final_sfl_cum': sfl_cum,
        'minimum_internal_energy': float(np.min(p['ie'])),
        'minimum_volume': float(np.min(p['volume'])),
        'checks': {
            'finite': bool(all(np.isfinite(q[field]).all() for q in (p0,p) for field in
                               ('x','y','vx','vy','mass','volume','ie','den','pressure','w2'))),
            'mass_exact': bool(np.array_equal(p0['mass'],p['mass'])),
            'W0_exact': bool(np.all(p0['w2']==0) and np.all(p['w2']==0)),
            'positive_U_and_volume': bool(np.all(p['ie']>0) and np.all(p['volume']>0)),
            'no_floor_or_dual_energy_injection': bool(floor_cum==0 and de_cum==0 and sfl_cum==0),
            'gizmo_L1_velocity_below_0p01': bool(l1_v<.01),
            'gizmo_L1_density_below_0p25': bool(l1_rho<.25),
            'disk_density_change_below_0p05': bool(disk_density<.05),
            'disk_entropy_log_abs_below_0p01': bool(disk_entropy<.01),
            'disk_angular_momentum_drift_below_0p01': bool(abs(Lt/L0-1.)<.01),
            'total_energy_drift_below_0p001': bool(abs((Et-E0)/abs(E0))<.001),
        },
        'scope': 'Native CUDA+MPI+OpenMP RK4 fixed-mass W=0 t=10 gate; not the GIZMO t=600 endpoint',
    }
    with (args.run/'audit_t10.json').open('x') as stream:
        json.dump(result,stream,indent=2)
        stream.write('\n')
    print(json.dumps(result,indent=2))
    if not all(result['checks'].values()):
        raise SystemExit('t=10 native physics gate failed')


if __name__=='__main__':
    main()
