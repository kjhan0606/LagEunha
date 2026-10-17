"""Fixed-step temporal sensitivity of density error, with checkpoint resume."""
import argparse
import json
from pathlib import Path
import numpy as np
from adaptive_coupled_mass import atomic,trial
from adaptive_moment_batch import setup,diagnostic
from factor_moment_evolve import rebuild
from shear_interface_audit import audit


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--n',type=int,required=True)
    parser.add_argument('--steps',type=int,choices=(32,64,128),required=True);args=parser.parse_args()
    old=json.loads(Path('reference.checkpoint.json').read_text())
    if (old['config']['kind'],old['config']['n'],old['config']['amplitude'],old['t'])!=('shear_contact',args.n,100.,.01):
        raise RuntimeError('adaptive reference mismatch')
    if not json.loads(Path('predecessor_completion.json').read_text())['reference_valid']:
        raise RuntimeError('interface reference invalid')
    cfg=dict(kind='shear_contact',n=args.n,amplitude=100.,end=.01,fixed_steps=args.steps)
    mass,y0=setup('shear_contact',args.n,100.)
    np.testing.assert_allclose(old['mass'],mass,rtol=0,atol=1e-15)
    L=[np.zeros((2,2)) for _ in mass];y0=rebuild(y0,mass,L)
    cp=Path('inputs/current/checkpoint.json');cp.parent.mkdir(parents=True,exist_ok=True)
    if Path('failure.json').exists():raise RuntimeError('previous failure requires diagnosis')
    if cp.exists():
        saved=json.loads(cp.read_text())
        if saved['config']!=cfg:raise RuntimeError('checkpoint mismatch')
        y=np.asarray(saved['state']);L=[np.asarray(f) for f in saved['factors']]
        start=saved['steps'];peak=saved['max_ledger_error'];minimum=saved['min_pressure']
    else:y=y0.copy();start=0;peak=0.;minimum=1e300
    try:
        dt=.01/args.steps
        for k in range(start,args.steps):
            record,new,newL=trial(k*dt,y,L,mass,dt)
            d=diagnostic((k+1)*dt,new,mass,y0)
            if not np.all(np.isfinite(new)) or d['min_pressure']<=0 or d['ledger_error']>=1e-10:
                raise RuntimeError('finite/pressure/ledger gate failed')
            y,L=new,newL;peak=max(peak,d['ledger_error']);minimum=min(minimum,d['min_pressure'])
            atomic(cp,dict(config=cfg,t=(k+1)*dt,steps=k+1,state=y,mass=mass,factors=L,
                max_ledger_error=peak,min_pressure=minimum))
            if (k+1)%16==0:print(json.dumps(dict(step=k+1,total=args.steps,t=(k+1)*dt)),flush=True)
        measured=audit('current')
        reference=np.asarray(old['state'])
        difference=float(np.max(abs(y-reference)/np.maximum(1.,abs(reference))))
        atomic('result.json',dict(n=args.n,steps=args.steps,interface=measured,
            max_scaled_state_difference_from_adaptive=difference,
            max_ledger_error=peak,min_pressure=minimum,
            note='Fixed-step sensitivity only; dense-path gates not evaluated for these coarse runs.'))
    except Exception as exc:
        atomic('failure.json',dict(error=repr(exc)));raise


if __name__=='__main__':main()
