"""Deforming density contact: affine velocity/pressure diagnostics.

No claim that sheared material polygons remain a Euclidean power diagram.
Velocity-moment identities are measured on the actual numerical cells.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from adaptive_coupled_mass import run,atomic
from adaptive_moment_batch import setup,diagnostic


def affine_errors(y,mass):
    a=y[:-4].reshape(-1,20);q=a[:,3:5];Q=a[:,5:9].reshape(-1,2,2)
    p=a[:,9:11];D=a[:,11:15].reshape(-1,2,2);Pi=a[:,15:19].reshape(-1,2,2)
    u=np.array([.2,0.]);B=np.array([[0.,1.1],[0.,0.]])
    bq=q@B.T
    ep=mass[:,None]*u+bq
    ed=np.einsum('ni,j->nij',q,u)+Q@B.T
    epi=mass[:,None,None]*np.outer(u,u)+np.einsum('ni,j->nij',bq,u)+np.einsum('i,nj->nij',u,bq)+B@Q@B.T
    return dict(p=float(np.max(abs(p-ep)/mass[:,None])),
        D=float(np.max(abs(D-ed)/mass[:,None,None])),
        Pi=float(np.max(abs(Pi-epi)/mass[:,None,None])))


def main():
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,required=True)
    p.add_argument('--ratio',type=float,required=True);args=p.parse_args()
    if not json.loads(Path('predecessor_completion.json').read_text()).get('passed'):
        raise RuntimeError('predecessor failed')
    if Path('failure.json').exists():raise RuntimeError('previous failure needs investigation')
    try:
        output=Path('shear_contact/result.json')
        if output.exists():
            result=json.loads(output.read_text())
            for key,value in dict(kind='shear_contact',n=args.n,amplitude=args.ratio,duration=.01,mass_budget=5e-11).items():
                if result.get(key)!=value:raise RuntimeError('completed configuration mismatch')
        else:result=run('shear_contact',n=args.n,amplitude=args.ratio,end=.01)
        mass,y0=setup('shear_contact',args.n,args.ratio)
        c=json.loads(Path('shear_contact/checkpoint.json').read_text())
        np.testing.assert_allclose(c['mass'],mass,rtol=0,atol=1e-15)
        if c['t']!=.01:raise RuntimeError('incomplete endpoint')
        y=np.asarray(c['state']);errors=affine_errors(y,mass);d=diagnostic(.01,y,mass,y0)
        physics_gate=bool(max(errors.values())<1e-8 and d['pressure_departure']<1e-8)
        summary=dict(result=result,affine_moment_errors=errors,pressure_departure=d['pressure_departure'],
            affine_pressure_gate=physics_gate,
            passed=physics_gate and all(result.get(k) is True for k in ('state_gate','mass_path_gate','generator_path_gate')),
            note='Short deforming-contact diagnostic. Does not establish exact material-cell geometry or convergence of density interface.')
        atomic('case_summary.json',summary)
    except Exception as exc:
        atomic('failure.json',dict(error=repr(exc)));raise
    raise SystemExit(0 if summary['passed'] else 2)


if __name__=='__main__':main()
