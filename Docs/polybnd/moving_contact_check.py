"""Uniform translation of a cell-aligned contact and its enclosing boundary.

The analytical solution is used only for diagnostics, never for evolution.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from adaptive_coupled_mass import run,atomic
from adaptive_moment_batch import setup,diagnostic


def translated_state(y0,mass,velocity,t):
    out=y0.copy();a=out[:-4].reshape(-1,20);b=y0[:-4].reshape(-1,20)
    shift=t*np.asarray(velocity)
    a[:,:2]+=shift
    a[:,3:5]+=mass[:,None]*shift
    a[:,5:9]=(b[:,5:9].reshape(-1,2,2)
        +np.einsum('ni,j->nij',b[:,3:5],shift)
        +np.einsum('i,nj->nij',shift,b[:,3:5])
        +mass[:,None,None]*np.outer(shift,shift)).reshape(-1,4)
    a[:,11:15]=(b[:,11:15].reshape(-1,2,2)
        +np.einsum('i,nj->nij',shift,b[:,9:11])).reshape(-1,4)
    return out


def main():
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,required=True)
    p.add_argument('--speed',type=float,required=True);args=p.parse_args()
    if not json.loads(Path('predecessor_completion.json').read_text()).get('passed'):
        raise RuntimeError('predecessor failed')
    if Path('failure.json').exists():raise RuntimeError('previous failure needs investigation')
    try:
        path=Path('moving_contact/result.json')
        if path.exists():
            result=json.loads(path.read_text())
            for key,value in dict(kind='moving_contact',n=args.n,amplitude=args.speed,duration=.05,mass_budget=5e-11).items():
                if result.get(key)!=value:raise RuntimeError('configuration mismatch')
        else:result=run('moving_contact',n=args.n,amplitude=args.speed)
        mass,y0=setup('moving_contact',args.n,args.speed)
        c=json.loads(Path('moving_contact/checkpoint.json').read_text())
        np.testing.assert_allclose(c['mass'],mass,rtol=0,atol=1e-15)
        if c['t']!=.05:raise RuntimeError('endpoint incomplete')
        y=np.asarray(c['state']);u=np.array([args.speed,-.5*args.speed])
        expected=translated_state(y0,mass,u,.05)
        a=y[:-4].reshape(-1,20);e=expected[:-4].reshape(-1,20)
        velocity=float(np.max(abs(a[:,9:11]/mass[:,None]-u)))
        state=float(np.max(abs(a-e)/np.maximum(1.,abs(e))))
        d=diagnostic(.05,y,mass,y0)
        exact_gate=bool(velocity<1e-10 and state<1e-10 and d['pressure_departure']<1e-10)
        summary=dict(result=result,translation_gate=exact_gate,
            max_velocity_error=velocity,max_scaled_state_error=state,
            max_site_error=float(np.max(abs(a[:,:2]-e[:,:2]))),
            pressure_departure=d['pressure_departure'],velocity=u.tolist(),
            passed=exact_gate and all(result.get(k) is True for k in
                ('state_gate','mass_path_gate','generator_path_gate')),
            note='Translating material mesh and domain; does not test transport through a stationary mesh.')
        atomic('case_summary.json',summary)
    except Exception as exc:
        atomic('failure.json',dict(error=repr(exc)));raise
    raise SystemExit(0 if summary['passed'] else 2)


if __name__=='__main__':main()
