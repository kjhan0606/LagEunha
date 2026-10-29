"""Separate actual wall angular impulse from native angular-balance residual."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import numpy as np
from kepler_native_global_ids_20261029 import ordered_global, ID_MASK


def angular(p,mask=None):
    if mask is None:mask=np.ones(len(p),dtype=bool)
    m=p['mass'][mask].astype(float);x=p['x'][mask]-3.;y=p['y'][mask]-3.
    return float(np.sum(m*(x*p['vy'][mask]-y*p['vx'][mask]),dtype=np.longdouble))


def budget(initial,endpoint,logs):
    p0,m0=ordered_global(initial);ids=p0['u4if']&ID_MASK;p,meta=ordered_global(endpoint,ids)
    if not np.array_equal(p0['mass'],p['mass']):raise ValueError('Mass changed')
    pressure=np.zeros(3,dtype=np.longdouble);collision=pressure.copy();residual=np.longdouble(0)
    rows=[];impacts=0
    for log in logs:
        for line in log.read_text(errors='replace').splitlines():
            if '[WALL_ACCEPT]'not in line:continue
            start=float(re.search(r't_start=(\S+)',line)[1]);h=float(re.search(r'\bh=(\S+)',line)[1])
            rows.append((start,h));impacts+=int(re.search(r'impacts=(\d+)',line)[1])
            for name,total in [('wall_pressure_J',pressure),('wall_collision_J',collision)]:
                values=re.search(name+r'=\(([^,]+),([^,]+),([^\)]+)\)',line)
                total+=np.array([values[1],values[2],values[3]],dtype=np.longdouble)
            residual+=np.longdouble(re.search(r'angular_residual=(\S+)',line)[1])
    if not rows or rows[0][0]!=m0['t'] or abs(rows[-1][0]+rows[-1][1]-meta['t'])>1e-12:
        raise ValueError('Logs do not cover the complete snapshot interval')
    if any(abs(t+h-rows[i+1][0])>1e-12 for i,(t,h) in enumerate(rows[:-1])):
        raise ValueError('Logs contain a gap or duplicate physical interval')
    r0=np.hypot(p0['x']-3.,p0['y']-3.);disk=(r0>=.5)&(r0<2.)
    l0,lt=angular(p0),angular(p);d0,dt=angular(p0,disk),angular(p,disk)
    wall=float(pressure[2]+collision[2]);balance=(lt-l0)+wall
    mismatch=balance-float(residual)
    checks={'complete_accepted_interval':True,'fixed_mass_W0':bool(np.all(p0['w2']==0) and np.all(p['w2']==0)),
            'wall_budget_matches_logged_integrator_residual':abs(mismatch)<1e-10}
    return {'passed':all(checks.values()),'checks':checks,'initial_time':m0['t'],'time':meta['t'],
            'accepted_prefixes':len(rows),'impacts':impacts,'total_L_initial':l0,'total_L_final':lt,
            'total_fluid_L_change':lt-l0,'disk_L_initial':d0,'disk_L_final':dt,'disk_L_change':dt-d0,
            'disk_L_relative_change':dt/d0-1,'non_disk_L_change':(lt-l0)-(dt-d0),
            'pressure_wall_impulse_and_torque':pressure.astype(float).tolist(),
            'collision_wall_impulse_and_torque':collision.astype(float).tolist(),
            'wall_angular_impulse_over_initial_disk_L':wall/abs(d0),
            'fluid_plus_wall_angular_balance':balance,'sum_logged_integrator_angular_residual':float(residual),
            'budget_minus_sum_residual':mismatch,
            'source_hashes':{str(q):hashlib.sha256(q.read_bytes()).hexdigest() for q in [initial,endpoint,*logs]},
            'scope':'Wall torque itself and actual disk/non-disk angular changes are reported separately. Small corrected balance does not establish an isolated disk or GIZMO boundary equivalence. Non-disk is initial core plus exterior, not a measured pairwise inter-region flux.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--initial',type=Path,required=True);parser.add_argument('--endpoint',type=Path,required=True)
    parser.add_argument('--logs',type=Path,nargs='+',required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();result=budget(args.initial,args.endpoint,args.logs)
    with args.output.open('x') as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(result,indent=2))
    if not result['passed']:raise SystemExit(1)
