"""Contact equilibrium and combined wave extension, with resumable evolution."""
import argparse
import json
from pathlib import Path
import numpy as np
from adaptive_coupled_mass import run, atomic
from adaptive_moment_batch import setup, diagnostic


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--kind',choices=('contact','wave'),required=True)
    p.add_argument('--n',type=int,required=True)
    args=p.parse_args()
    predecessor=json.loads(Path('predecessor_completion.json').read_text())
    if not predecessor.get('passed'):raise RuntimeError('temporal predecessor failed')
    amplitude=0. if args.kind=='contact' else .1
    if Path('failure.json').exists():raise RuntimeError('investigate previous failure before retry')
    try:
        output=Path(args.kind)/'result.json'
        if output.exists():
            result=json.loads(output.read_text())
            for key,value in dict(kind=args.kind,n=args.n,amplitude=amplitude,duration=.05,mass_budget=5e-11).items():
                if result.get(key)!=value:raise RuntimeError('completed configuration mismatch')
        else:result=run(args.kind,n=args.n,amplitude=amplitude)
        passed=all(result.get(k) is True for k in ('state_gate','mass_path_gate','generator_path_gate'))
        summary=dict(passed=passed,result=result)
        if args.kind=='contact':
            mass,y0=setup('contact',args.n,0.)
            checkpoint=json.loads((Path('contact')/'checkpoint.json').read_text())
            np.testing.assert_allclose(checkpoint['mass'],mass,rtol=0,atol=1e-15)
            y=np.asarray(checkpoint['state'])
            a=y[:-4].reshape(-1,20);a0=y0[:-4].reshape(-1,20)
            d=diagnostic(.05,y,mass,y0)
            velocity=float(np.max(abs(a[:,9:11]/mass[:,None])))
            drift=float(np.max(abs(y-y0)))
            equilibrium=bool(velocity<1e-10 and drift<1e-10 and d['pressure_departure']<1e-10)
            summary.update(equilibrium_gate=equilibrium,max_velocity=velocity,
                max_state_drift=drift,pressure_departure=d['pressure_departure'],
                max_site_drift=float(np.max(abs(a[:,:2]-a0[:,:2]))),
                note='Cell-aligned stationary density jump 100:1, fixed wall, P=.6, v=0; not an advected contact or shock test.')
            summary['passed']=passed and equilibrium
        atomic('case_summary.json',summary)
    except Exception as exc:
        atomic('failure.json',dict(error=repr(exc)));raise
    raise SystemExit(0 if summary['passed'] else 2)


if __name__=='__main__':main()
