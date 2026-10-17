"""One independent, resumable extension case, with predecessor science gate."""
import argparse
import json
from pathlib import Path
from adaptive_coupled_mass import run,atomic,configuration


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--kind',choices=('wave','shear'),default='wave')
    parser.add_argument('--n',type=int,required=True)
    parser.add_argument('--amplitude',type=float,required=True)
    args=parser.parse_args()
    parameters=dict(n=args.n,amplitude=args.amplitude,end=.05,budget=5e-11)
    cfg=configuration(args.kind,**parameters)
    predecessor=json.loads(Path('predecessor_completion.json').read_text())
    if not predecessor.get('passed') or not all(predecessor['cases'][c]['passed'] for c in ('shear','wave')):
        atomic('gate.json',dict(status='blocked_by_predecessor'))
        raise SystemExit(2)
    atomic('gate.json',dict(status='passed',config=cfg))
    output=Path(args.kind)/'result.json'
    if output.exists():
        result=json.loads(output.read_text())
        if any(result.get(k)!=v for k,v in dict(kind=args.kind,n=args.n,amplitude=args.amplitude,duration=.05,mass_budget=5e-11).items()):
            raise RuntimeError('completed case configuration differs')
    else:
        # Time-limit interruption leaves no failure file and resumes normally.
        # A numerical exception is preserved for diagnosis, not blindly retried.
        if Path('failure.json').exists():
            raise RuntimeError('previous numerical failure requires investigation')
        try:result=run(args.kind,**parameters)
        except Exception as exc:
            atomic('failure.json',dict(error=repr(exc),config=cfg))
            raise
    passed=result.get('status')=='completed' and all(result.get(k) is True for k in
        ('state_gate','mass_path_gate','generator_path_gate'))
    atomic('case_summary.json',dict(passed=passed,result=result))
    raise SystemExit(0 if passed else 2)


if __name__=='__main__':main()
