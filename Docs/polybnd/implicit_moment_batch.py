"""Implicit midpoint diagnostic for compatible pressure/transport moments.
Nonlinear solver iterates are not physical stages; converged midpoints are.
No eigenvalue clipping. Failure is recorded, never replaced with a fallback.
"""
import json,time
from pathlib import Path
import numpy as np
from adaptive_moment_batch import setup,diagnostic
from compatible_moment_rhs import rhs


def advance(t,y,m,dt):
    trial=y.copy();calls=0
    for iteration in range(80):
        mid=(y+trial)/2;dy,d=rhs(t+dt/2,mid,m);calls+=1
        target=y+dt*dy;residual=float(np.max(abs(target-trial)))
        if residual<2e-14:
            # Return trial itself: its midpoint generated the checked residual.
            return trial,mid,d,calls,residual
        trial=.2*trial+.8*target
    raise RuntimeError('midpoint nonlinear solve did not converge')


def run(kind,amplitude,dt,end=.02):
    m,y=setup(kind,4,amplitude);y0=y.copy();start=time.monotonic();history=[];calls=0;maxres=0.
    tag=f'{kind}_a{amplitude:g}_dt{dt:g}'
    try:
        for step in range(round(end/dt)):
            t=step*dt;new,mid,d,n,res=advance(t,y,m,dt);calls+=n;maxres=max(maxres,res)
            dm=diagnostic(t+dt/2,mid,m,y0);df=diagnostic(t+dt,new,m,y0)
            history.append(dict(t=t+dt,midpoint_min_R=dm['min_R'],endpoint_min_R=df['min_R'],pressure_error=df['pressure_departure'],ledger_error=df['ledger_error']))
            y=new
        result=dict(kind=kind,amplitude=amplitude,dt=dt,end=end,status='completed',calls=calls,seconds=time.monotonic()-start,
                    max_nonlinear_residual=maxres,min_midpoint_R=min(r['midpoint_min_R'] for r in history),
                    min_endpoint_R=min(r['endpoint_min_R'] for r in history),final=diagnostic(end,y,m,y0))
        Path(tag+'.history.json').write_text(json.dumps(history,indent=2))
        Path(tag+'.checkpoint.json').write_text(json.dumps(dict(t=end,mass=m.tolist(),state=y.tolist()),indent=2))
    except Exception as e:result=dict(kind=kind,amplitude=amplitude,dt=dt,status='failed',error=repr(e),calls=calls,seconds=time.monotonic()-start)
    return result

if __name__=='__main__':
    for kind,a in [('wave',.001),('wave',.1),('shear',0.)]:
        for dt in (.001,.0005):print(json.dumps(run(kind,a,dt)),flush=True)
