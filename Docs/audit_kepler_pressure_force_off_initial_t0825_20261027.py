"""Audit a causal pressure-force-off continuation from the accepted t0 state."""

import argparse
import json
from pathlib import Path
import re

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from audit_kepler_sparse_native_ballistic_20261027 import disk_geometry
from audit_kepler_sparse_state_consistency_20261027 import state_metrics
from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time


N=55351
T=0.825
INITIAL_SHA='8b561cdf05bfa06d2f42c5e7b782a2174719d5306967f4201b487d6113c4ab5b'
BINARY_SHA='0a65c373613a5b31f8449b456387ebcf80d72136a81e46eec10c9692a8bed5e7'
DOMAIN_AREA=.5*1024*2.75**2*np.sin(2*np.pi/1024)


def weighted_abs(values,mass):
    return float(np.sum(np.abs(values)*mass,dtype=np.longdouble)/
                 np.sum(mass,dtype=np.longdouble))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('run',type=Path)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    run=args.run.resolve()
    initial=run/'khout.000001.dat'
    standard_audit_path=run/'reference_native_t0825.json'
    reference=json.loads(standard_audit_path.read_text())
    if (reference['time']!=T or reference['initial_sha256']!=INITIAL_SHA or
            digest(initial)!=INITIAL_SHA or digest(Path(reference['initial']))!=INITIAL_SHA or
            digest(Path(reference['endpoint']))!=reference['endpoint_sha256'] or
            digest(run/'eunha2.exe')!=BINARY_SHA):
        raise ValueError('Changed baseline or diagnostic input')
    endpoints=[p for p in run.glob('khout.*.dat')
               if abs(snapshot_time(p,N)-T)<1e-10]
    if len(endpoints)!=1:
        raise ValueError(f'Expected exactly one diagnostic t={T} endpoint')
    endpoint=endpoints[0]
    p0,m0=ordered(initial,N)
    ps,ms=ordered(Path(reference['endpoint']),N)
    pd,md=ordered(endpoint,N)
    if abs(m0['t'])>1e-12 or abs(ms['t']-T)>1e-10 or abs(md['t']-T)>1e-10:
        raise ValueError('Wrong snapshot times')
    radius=np.hypot(p0['x'].astype(float)-3.,p0['y'].astype(float)-3.)
    disk=(radius>=.5)&(radius<2.)
    selected=np.flatnonzero(disk)
    exact_xy,exact_area,image_safety=disk_geometry(p0,T,selected)
    phase=T*(radius*radius+.01175**2)**(-.75)
    exact_vx=p0['vx'].astype(float)*np.cos(phase)-p0['vy'].astype(float)*np.sin(phase)
    exact_vy=p0['vx'].astype(float)*np.sin(phase)+p0['vy'].astype(float)*np.cos(phase)
    mass=p0['mass'][selected].astype(float)
    metrics={}
    for label,p in (('standard',ps),('pressure_off',pd)):
        pos=np.hypot(p['x'][selected].astype(float)-exact_xy[selected,0],
                     p['y'][selected].astype(float)-exact_xy[selected,1])
        vel=np.hypot(p['vx'][selected].astype(float)-exact_vx[selected],
                     p['vy'][selected].astype(float)-exact_vy[selected])
        density=p['den'][selected].astype(float)/p0['den'][selected].astype(float)-1.
        residual=exact_area/p['volume'][selected].astype(float)-1.
        metrics[label]={
            'disk_position_error_abs_mass_weighted':weighted_abs(pos,mass),
            'disk_velocity_error_abs_mass_weighted':weighted_abs(vel,mass),
            'disk_density_change_abs_mass_weighted':weighted_abs(density,mass),
            'disk_native_minus_ballistic_density_residual_abs_mass_weighted':weighted_abs(residual,mass),
        }
    delta_pos=np.hypot(pd['x'][selected].astype(float)-ps['x'][selected].astype(float),
                       pd['y'][selected].astype(float)-ps['y'][selected].astype(float))
    delta_vel=np.hypot(pd['vx'][selected].astype(float)-ps['vx'][selected].astype(float),
                       pd['vy'][selected].astype(float)-ps['vy'][selected].astype(float))
    log=(run/'log').read_text(errors='replace')
    ledgers=[line for line in log.splitlines() if '[RK4E]' in line]
    if not ledgers:
        raise ValueError('Missing RK4 energy ledger')
    no_injection=all(re.search(rf'\b{key}=0\.000e\+00\b',ledgers[-1])
                     for key in ('floor_cum','de_cum','sfl_cum'))
    s0,sd=state_metrics(p0),state_metrics(pd)
    checks={
        'exact_t0_and_t0825':True,
        'fixed_mass':bool(np.array_equal(p0['mass'],pd['mass'])),
        'W0':bool(np.all(p0['w2']==0) and np.all(pd['w2']==0)),
        'positive_states':bool(s0['positive'] and sd['positive']),
        'density_and_EOS':bool(max(sd['density_relative_max'],sd['eos_pressure_relative_max'])<1e-6),
        'polygon_area':bool(abs(sd['volume_sum']-DOMAIN_AREA)<1e-5),
        'zero_injection':bool(no_injection),
        'no_geometry_abort':bool('[GPU_TESS_INVALID]' not in log and '[GPU_CLIP_CAPACITY]' not in log),
        'no_wall_crossing':bool('[FIXED_POLYGON_CROSSING]' not in log),
    }
    result={
        'time':T,'initial':str(initial),'initial_sha256':INITIAL_SHA,
        'standard_audit':str(standard_audit_path),'standard_audit_sha256':digest(standard_audit_path),
        'standard_endpoint':reference['endpoint'],'standard_endpoint_sha256':reference['endpoint_sha256'],
        'pressure_off_endpoint':str(endpoint),'pressure_off_endpoint_sha256':digest(endpoint),
        'binary_sha256':BINARY_SHA,'disk_sites':len(selected),
        'minimum_ballistic_image_safety_ratio':image_safety,
        'metrics':metrics,
        'pressure_off_vs_standard_position_abs_mass_weighted':weighted_abs(delta_pos,mass),
        'pressure_off_vs_standard_velocity_abs_mass_weighted':weighted_abs(delta_vel,mass),
        'checks':checks,
        'scope':'Causal diagnostic only: pressure-force momentum term removed while geometric PdV and EOS remain. The force/work adjoint is broken; safety checks are not a physical hydro acceptance gate.'}
    with args.output.open('x') as stream:
        json.dump(result,stream,indent=2,allow_nan=False)
        stream.write('\n')
    print(json.dumps(result,indent=2))
    if not all(checks.values()):
        raise SystemExit('Diagnostic safety audit failed')


if __name__=='__main__':
    main()
