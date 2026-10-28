"""Compare three saved 2D states at one time; never advances a solution.

Differences are reported per variable, without mixing their physical units.
Observed orders alone do not establish spatial or long-time physical accuracy.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np


def compare(states, timesteps):
    dt = np.asarray(timesteps, dtype=float)
    if dt.shape != (3,) or not np.isfinite(dt).all() or np.any(dt <= 0):
        raise ValueError('three finite positive timesteps required')
    ratios = dt[:-1] / dt[1:]
    if ratios[0] <= 1 or not np.allclose(ratios, ratios[0], rtol=1e-12, atol=0):
        raise ValueError('timesteps must decrease by a common refinement factor')
    if len(states) != 3:
        raise ValueError('three states required')
    potential_present=['potential_moments' in s for s in states]
    if any(potential_present) and not all(potential_present):
        raise ValueError('potential moments must be present in all states or none')
    times = np.array([s['t'] for s in states])
    if not np.isfinite(times).all() or np.ptp(times) > 16*np.finfo(float).eps*max(1., abs(times[0])):
        raise ValueError('states must have the same finite time')
    arrays = []
    for s in states:
        mass = np.asarray(s['mass'], dtype=float)
        if mass.ndim != 1 or not len(mass) or not np.isfinite(mass).all() or np.any(mass <= 0):
            raise ValueError('positive finite cell masses required')
        np.testing.assert_array_equal(mass, states[0]['mass'])
        y = np.asarray(s['state'], dtype=float)
        L = np.asarray(s['factors'], dtype=float)
        if y.shape != (20*len(mass)+4,) or L.shape != (len(mass),2,2):
            raise ValueError('2D state shape mismatch')
        if not all(np.isfinite(a).all() for a in (y,L)):
            raise ValueError('nonfinite state')
        a = y[:-4].reshape(-1,20)
        arrays.append(dict(generator=a[:,:2], weight=a[:,2], first_moment=a[:,3:5],
            second_moment=a[:,5:9], momentum=a[:,9:11], mixed_moment=a[:,11:15],
            velocity_second_moment=a[:,15:19], gas_energy=a[:,19],
            boundary_ledger=y[-4:],
            covariance=L@np.swapaxes(L,-1,-2)))
        if all(potential_present):
            U=np.asarray(s['potential_moments'],dtype=float)
            if U.shape!=mass.shape or not np.isfinite(U).all():
                raise ValueError('invalid potential moments')
            arrays[-1]['potential_moment']=U
    result = {}
    for key in arrays[0]:
        a,b,c = (s[key] for s in arrays)
        ab,bc = a-b,b-c
        norm1,norm2 = float(np.linalg.norm(ab)),float(np.linalg.norm(bc))
        scale = max(float(np.linalg.norm(v)) for v in (a,b,c))
        floor = 64*np.finfo(float).eps*max(scale,np.finfo(float).tiny)
        resolved = min(norm1,norm2) > floor
        result[key] = dict(coarse_fine_l2=norm1, fine_finest_l2=norm2,
            coarse_fine_max=float(np.max(abs(ab))), fine_finest_max=float(np.max(abs(bc))),
            roundoff_screen=floor, differences_resolved=bool(resolved),
            observed_order=float(np.log(norm1/norm2)/np.log(ratios[0])) if resolved else None,
            difference_alignment=float(np.vdot(ab,bc)/(norm1*norm2)) if resolved else None)
    return dict(t=float(times[0]),timesteps=dt.tolist(),variables=result,full_validation=False,
        scope='Same-time temporal self-convergence only; physical fields, spatial convergence and long duration remain separate.')


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('states',nargs=3)
    parser.add_argument('--dt',nargs=3,type=float,required=True)
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    data=[Path(p).read_bytes() for p in args.states]
    report=compare([json.loads(b) for b in data],args.dt)
    report['inputs']=[dict(path=str(Path(p).resolve()),sha256=hashlib.sha256(b).hexdigest()) for p,b in zip(args.states,data)]
    with Path(args.output).open('x') as f:json.dump(report,f,indent=2)
