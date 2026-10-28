"""Propagate accepted t6.5 sites under only the exact softened point force."""

import argparse
import json
from pathlib import Path
import time

import numpy as np
from scipy.integrate import solve_ivp

from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest


ROOT = Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027')
PARENT_AUDIT = ROOT / 'jitter005_to65_a10_v1/audit_t65.json'
OFF_AUDIT = ROOT / 'pressure_force_off_jitter_to10_a10_v1/audit_t10_native.json'
STANDARD_AUDIT = ROOT / 'jitter005_to10_a10_v1/audit_t10.json'
EPS2 = .01175**2


def check_audit(path, t, allow_ballistic_failure=False):
    d = json.loads(path.read_text())
    failed = {key for key,value in d['checks'].items() if not value}
    allowed = {'disk_ballistic_residual'} if allow_ballistic_failure else set()
    if d['time'] != t or failed != allowed:
        raise ValueError(f'{path}: time or failed checks changed: {failed}')
    if digest(Path(d['endpoint'])) != d['endpoint_sha256']:
        raise ValueError(f'{path}: endpoint changed')
    return d


def weighted_abs(values, weights):
    return float(np.sum(np.abs(values)*weights,dtype=np.longdouble)
                 /np.sum(weights,dtype=np.longdouble))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    start=time.monotonic()
    parent=check_audit(PARENT_AUDIT,6.5)
    off=check_audit(OFF_AUDIT,10.,True)
    standard=check_audit(STANDARD_AUDIT,10.,True)
    p65,_=ordered(Path(parent['endpoint']),55351)
    p0,_=ordered(Path(parent['initial']),55351)
    poff,_=ordered(Path(off['endpoint']),55351)
    pstd,_=ordered(Path(standard['endpoint']),55351)
    ids=p65['u4if']
    if not (np.array_equal(ids,p0['u4if']) and np.array_equal(ids,poff['u4if'])
            and np.array_equal(ids,pstd['u4if'])):
        raise ValueError('Particle ID order changed')
    n=len(p65)
    y0=np.concatenate((p65['x'].astype(float)-3.,p65['y'].astype(float)-3.,
                       p65['vx'].astype(float),p65['vy'].astype(float)))

    def rhs(_, y):
        out=np.empty_like(y)
        x,ypos,vx,vy=np.split(y,4)
        r2=x*x+ypos*ypos+EPS2
        f=-1.0/(r2*np.sqrt(r2))
        out[:n]=vx
        out[n:2*n]=vy
        out[2*n:3*n]=x*f
        out[3*n:]=ypos*f
        return out

    solution=solve_ivp(rhs,(6.5,10.),y0,method='DOP853',rtol=1e-12,
                       atol=1e-14,max_step=.05,t_eval=[10.])
    if not solution.success:
        raise ValueError(solution.message)
    ode=solution.y[:,-1]
    sample=np.linspace(0,n-1,97,dtype=int)
    sample_y0=np.concatenate([y0[k*n+sample] for k in range(4)])
    count=len(sample)

    def sample_rhs(_,y):
        out=np.empty_like(y)
        x,ypos,vx,vy=np.split(y,4)
        r2=x*x+ypos*ypos+EPS2
        f=-1.0/(r2*np.sqrt(r2))
        out[:count]=vx;out[count:2*count]=vy
        out[2*count:3*count]=x*f;out[3*count:]=ypos*f
        return out

    tight=solve_ivp(sample_rhs,(6.5,10.),sample_y0,method='DOP853',
                    rtol=1e-13,atol=1e-15,max_step=.01,t_eval=[10.])
    if not tight.success:
        raise ValueError(tight.message)
    sample_max=float(max(np.max(np.abs(ode[k*n+sample]
                                      -tight.y[k*count:(k+1)*count,-1]))
                         for k in range(4)))
    if sample_max>1e-8:
        raise ValueError(f'ODE convergence check failed: {sample_max}')
    r0x=p0['x'].astype(float)-3.;r0y=p0['y'].astype(float)-3.
    r0=np.hypot(r0x,r0y)
    disk=(r0>=.5)&(r0<2.)
    mass=p0['mass'][disk].astype(float)
    phase=10.*(r0*r0+EPS2)**(-.75)
    x_exact=3.+r0x*np.cos(phase)-r0y*np.sin(phase)
    y_exact=3.+r0x*np.sin(phase)+r0y*np.cos(phase)
    ode_x=ode[:n]+3.;ode_y=ode[n:2*n]+3.

    def point_error(xa,ya,xb,yb):
        return weighted_abs(np.hypot(xa[disk]-xb[disk],ya[disk]-yb[disk]),mass)

    metrics={
        'ode_vs_exact_circular_position_mean':point_error(ode_x,ode_y,x_exact,y_exact),
        'pressure_off_vs_ode_position_mean':point_error(
            poff['x'].astype(float),poff['y'].astype(float),ode_x,ode_y),
        'pressure_off_vs_exact_circular_position_mean':point_error(
            poff['x'].astype(float),poff['y'].astype(float),x_exact,y_exact),
        'standard_vs_ode_position_mean':point_error(
            pstd['x'].astype(float),pstd['y'].astype(float),ode_x,ode_y),
        'pressure_off_vs_ode_speed_mean':weighted_abs(np.hypot(
            poff['vx'][disk].astype(float)-ode[2*n:3*n][disk],
            poff['vy'][disk].astype(float)-ode[3*n:][disk]),mass),
    }
    if abs(metrics['pressure_off_vs_exact_circular_position_mean']
           -off['disk_position_error_abs_mass_weighted'])>1e-10:
        raise ValueError('Pressure-off orbital metric differs from native auditor')
    result={'time_start':6.5,'time_end':10.,'n':n,'disk_sites':int(np.sum(disk)),
            'initial_parent_sha256':parent['endpoint_sha256'],
            'pressure_off_sha256':off['endpoint_sha256'],
            'standard_sha256':standard['endpoint_sha256'],
            'softening':.01175,'ode_method':'DOP853','ode_rtol':1e-12,
            'ode_atol':1e-14,'ode_max_step':.05,
            'ode_nfev':solution.nfev,'sample_tighter_max_abs_state_difference':sample_max,
            'metrics':metrics,'elapsed_seconds':time.monotonic()-start,
            'scope':'Independent point-force propagation of the already perturbed native t6.5 positions and velocities. Separates inherited orbital error from further native forces; does not reconstruct Voronoi volumes or hydro energy.'}
    with args.output.open('x') as out:
        json.dump(result,out,indent=2)
        out.write('\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
