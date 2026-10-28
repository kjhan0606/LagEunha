"""Evaluate a native polar IC at t=.825 against its exact t=0 snapshot."""

from pathlib import Path
import argparse
import json
import sys

import numpy as np

ROOT = Path('/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020')
sys.path.insert(0, str(ROOT))
from audit_snapshot import read_snapshot


def ordered(path):
    p, meta = read_snapshot(path)
    index = (p['u4if'].astype(np.uint64) & np.uint64((1<<56)-1)).astype(np.int64)
    if not np.array_equal(np.sort(index), np.arange(384*384)):
        raise ValueError('Snapshot ID set changed')
    return p[np.argsort(index)], meta


def weighted_mean(value, weight):
    return float(np.sum(value*weight, dtype=np.longdouble)/np.sum(weight, dtype=np.longdouble))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('initial', type=Path)
    parser.add_argument('run', type=Path)
    args = parser.parse_args()
    p0, m0 = ordered(args.initial)
    candidates = []
    for path in args.run.glob('khout.*.dat'):
        p, meta = ordered(path)
        if abs(float(meta['t'])-.825)<1e-12:
            candidates.append((path,p,meta))
    if len(candidates)!=1 or m0['t']!=0.:
        raise ValueError('Missing unique exact t=.825 or t=0 snapshot')
    path,p,m = candidates[0]
    dx0,dy0=p0['x']-3.,p0['y']-3.
    r0=np.hypot(dx0,dy0)
    omega=(r0*r0+.01175**2)**(-.75)
    angle=float(m['t'])*omega
    xb=3.+dx0*np.cos(angle)-dy0*np.sin(angle)
    yb=3.+dx0*np.sin(angle)+dy0*np.cos(angle)
    position_error=np.hypot(p['x']-xb,p['y']-yb)
    density_change=p['den'].astype(float)/p0['den'].astype(float)-1.
    pressure_change=p['pressure'].astype(float)/p0['pressure'].astype(float)-1.
    gm1=float(np.float32(1.666667))-1.
    entropy=p['ie'].astype(float)/(p0['ie'].astype(float)*(p0['volume'].astype(float)/p['volume'])**gm1)
    radius=np.hypot(p['x']-3.,p['y']-3.)
    vr=((p['x']-3.)*p['vx']+(p['y']-3.)*p['vy'])/np.maximum(radius,1e-300)
    bins=[]
    for lo,hi in ((0.,.5),(.5,.75),(.75,1.),(1.,2.),(.5,2.),(2.,2.35)):
        take=(r0>=lo)&(r0<hi)
        take_v=take&(r0>=.05)
        w=p0['mass'][take].astype(float)
        wv=p0['mass'][take_v].astype(float)
        bins.append(dict(rlo=lo,rhi=hi,cells=int(np.sum(take)),
                         density_error_mass_weighted=weighted_mean(abs(density_change[take]),w),
                         density_error_p95=float(np.quantile(abs(density_change[take]),.95)),
                         pressure_error_mass_weighted=weighted_mean(abs(pressure_change[take]),w),
                         entropy_error_mass_weighted=weighted_mean(abs(entropy[take]-1.),w),
                         position_error_mass_weighted=weighted_mean(position_error[take],w),
                         position_error_max=float(np.max(position_error[take])),
                         radial_speed_rms_over_circular=float(np.sqrt(weighted_mean((vr[take_v]/(r0[take_v]*omega[take_v]))**2,wv)))))
    checks=dict(finite=bool(all(np.isfinite(q[key]).all() for q in (p0,p) for key in
                                ('x','y','vx','vy','mass','volume','ie','den','pressure','w2'))),
                mass_exact=bool(np.array_equal(p0['mass'],p['mass'])),
                W0_exact=bool(np.all(p0['w2']==0) and np.all(p['w2']==0)),
                U_positive=bool(np.all(p['ie']>0)),
                volume_positive=bool(np.all(p['volume']>0)))
    result=dict(initial=str(args.initial), endpoint=str(path), time=float(m['t']),
                n=len(p), checks=checks, min_U=float(np.min(p['ie'])),
                min_volume=float(np.min(p['volume'])), bins=bins,
                scope='Native fixed-mass W0 RK4+CUDA evolution of Hopkins-profile polar disk; accuracy against exact circular orbits and initial fields. This is an early gate, not the paper final time.')
    with (args.run/'audit_t0825.json').open('x') as stream:
        json.dump(result,stream,indent=2)
        stream.write('\n')
    print(json.dumps(result,indent=2))
    if not all(checks.values()):
        raise SystemExit('Physical state gate failed')


if __name__=='__main__':
    main()
