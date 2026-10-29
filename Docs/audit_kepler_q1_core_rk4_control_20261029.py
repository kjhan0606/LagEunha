"""Replay only core gravity ODEs with the actual accepted timestep sequence.

This is an independent numerical attribution control on 134 initial cells,
not a hydrodynamic rerun, steering, restart replacement or native promotion.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import time
import numpy as np
from kepler_native_global_ids_20261029 import ordered_global, ID_MASK
from audit_kepler_native_radial_20261026 import mean


def advance(q, h):
    def rhs(a):
        radius2=a[:,0]**2+a[:,1]**2+.01175**2
        f=-1/(radius2*np.sqrt(radius2))
        return np.column_stack((a[:,2],a[:,3],a[:,0]*f,a[:,1]*f))
    k1=h*rhs(q);k2=h*rhs(q+.5*k1);k3=h*rhs(q+.5*k2);k4=h*rhs(q+k3)
    return q+(k1+2*k2+2*k3+k4)/6


def analyze(initial, endpoint, log):
    started=time.monotonic()
    p0,m0=ordered_global(initial);ids=p0['u4if']&ID_MASK
    p,meta=ordered_global(endpoint,ids)
    x=p0['x']-3.;y=p0['y']-3.;r=np.hypot(x,y);core=r<.1
    if core.sum()!=134 or m0['t']!=0. or meta['t']!=.825:
        raise ValueError('Requires the canonical first Q1 core checkpoint')
    rows=[]
    for line in log.read_text(errors='replace').splitlines():
        if '[WALL_ACCEPT]' in line:
            rows.append((float(re.search(r't_start=(\S+)',line)[1]),float(re.search(r'\bh=(\S+)',line)[1])))
    if (not rows or rows[0][0]!=0 or abs(rows[-1][0]+rows[-1][1]-meta['t'])>1e-12
            or any(abs(t+h-rows[i+1][0])>1e-12 for i,(t,h) in enumerate(rows[:-1]))):
        raise ValueError('Accepted full-precision timestep sequence is incomplete')
    q0=np.column_stack((x[core],y[core],p0['vx'][core],p0['vy'][core]));mass=p0['mass'][core]
    omega=(r[core]**2+.01175**2)**(-.75);phase=meta['t']*omega
    exact=np.column_stack((x[core]*np.cos(phase)-y[core]*np.sin(phase),x[core]*np.sin(phase)+y[core]*np.cos(phase)))
    native=np.column_stack((p['x'][core]-3.,p['y'][core]-3.))
    native_error=mean(np.linalg.norm(native-exact,axis=1),mass)
    outputs=[];base=None
    for subdivision in [1,2,4]:
        q=q0.copy()
        for _,h in rows:
            for _ in range(subdivision):q=advance(q,h/subdivision)
        error=mean(np.linalg.norm(q[:,:2]-exact,axis=1),mass)
        outputs.append({'subdivision':subdivision,'mass_weighted_position_error':error,
                        'native_position_residual':mean(np.linalg.norm(q[:,:2]-native,axis=1),mass)})
        if subdivision==1:base=q
    residual=mean(np.linalg.norm(base[:,:2]-native,axis=1),mass)
    checks={'same_accepted_time_sequence':True,'all_control_states_finite':bool(np.isfinite(base).all()),
            'native_core_error_explained_to_5pct_by_gravity_RK4':residual<.05*native_error,
            'half_step_error_reduces_by_more_than_8':outputs[1]['mass_weighted_position_error']<outputs[0]['mass_weighted_position_error']/8,
            'quarter_step_error_reduces_again_by_more_than_8':outputs[2]['mass_weighted_position_error']<outputs[1]['mass_weighted_position_error']/8}
    return {'passed':all(checks.values()),'checks':checks,'cells':int(core.sum()),'time':meta['t'],
            'sources':{str(q):hashlib.sha256(q.read_bytes()).hexdigest() for q in [initial,endpoint,log]},
            'native_accepted_steps':len(rows),'native_core_position_error':native_error,
            'native_minus_gravity_control_position_residual':residual,'residual_over_native_orbit_error':residual/native_error,
            'gravity_only_RK4_controls':outputs,'analysis_seconds':time.monotonic()-started,
            'scope':'Only gravity ODE trajectories, with native accepted h and smaller h controls. An agreement attributes position error, not thermal entropy, density error, boundary effects, shocks or long-time gas stability. The native run is untouched.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['initial','endpoint','log','output']:parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args();result=analyze(args.initial,args.endpoint,args.log)
    with args.output.open('x') as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(result,indent=2))
    if not result['passed']:raise SystemExit(1)
