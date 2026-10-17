"""2D convergence tests with a primary PSD covariance factor at every stage.
Affine donor transport model; residual-stress option deliberately excluded.
"""
import json,time
from pathlib import Path
import numpy as np
from compatible_2d_suite import initial,compare,gresho
from adaptive_moment_batch import diagnostic
from compatible_moment_rhs import rhs
from factor_moment_evolve import rebuild,advance_factors
from periodic2d_hllc import quadrature
import moment_ale_evolve as model


def run(kind,n,end,dt=.001):
    mass,y=initial(kind,n);a=y[:-4].reshape(n*n,20);polys,_=model.geometry(a[:,:2],a[:,2],0.)
    factors=[]
    for i,poly in enumerate(polys):
        if kind=='wave':factors.append(np.zeros((2,2)));continue
        r,w=quadrature(poly);v,_=gresho(r);s=model.decode(poly,np.r_[mass[i],a[i,3:]])
        residual=v-s['v']-(r-s['c'])@s['B'].T
        _,L=np.linalg.qr(np.sqrt(w)[:,None]*residual,mode='reduced');factors.append(L.T)
    y=rebuild(y,mass,factors);y0=y.copy();initial_error=compare(kind,0.,y,mass);start=time.monotonic();history=[];calls=0;max_mass_flux=0.
    def f(t,z,L):
        nonlocal calls,max_mass_flux;calls+=1
        value,diag=rhs(t,z,mass,linear_pressure=True,covariance_factors=L,return_covariance_terms=True)
        max_mass_flux=max(max_mass_flux,float(diag['mass_flux']))
        return value,diag
    try:
        for k in range(round(end/dt)):
            t=k*dt;dy,d=f(t,y,factors)
            hmL=advance_factors(factors,d,dt/2);hm=rebuild(y+dt*dy/2,mass,hmL)
            dh,dd=f(t+dt/2,hm,hmL)
            nextL=advance_factors(factors,dd,dt);ny=rebuild(y+dt*dh,mass,nextL)
            mid_diag=diagnostic(t+dt/2,hm,mass,y0);last_diag=diagnostic(t+dt,ny,mass,y0)
            minR=min(float(np.linalg.eigvalsh(L@L.T).min()/m) for L,m in zip(hmL+nextL,list(mass)*2))
            history.append(dict(t=t+dt,factor_min_R=minR,min_pressure=min(mid_diag['min_pressure'],last_diag['min_pressure']),ledger_error=last_diag['ledger_error']))
            y,factors=ny,nextL
        err=compare(kind,end,y,mass)
        result=dict(kind=kind,n=n,end=end,dt=dt,status='completed',max_mass_flux=max_mass_flux,calls=calls,seconds=time.monotonic()-start,
                    errors=err,min_factor_R=min(h['factor_min_R'] for h in history),min_pressure=min(h['min_pressure'] for h in history),
                    kinetic_fraction=err['kinetic']/initial_error['kinetic'] if initial_error['kinetic'] else None,
                    final=diagnostic(end,y,mass,y0))
        tag=f'{kind}_n{n}_dt{dt:g}'
        Path(tag+'.history.json').write_text(json.dumps(history,indent=2))
        Path(tag+'.checkpoint.json').write_text(json.dumps(dict(t=end,mass=mass.tolist(),state=y.tolist(),factors=[L.tolist() for L in factors]),indent=2))
    except Exception as e:result=dict(kind=kind,n=n,end=end,status='failed',error=repr(e),calls=calls,seconds=time.monotonic()-start)
    return result

if __name__=='__main__':
    import sys
    if '--kinetic' in sys.argv:
        from kinetic_factor_rhs import rhs
    for kind,end in [('wave',.05),('gresho',.02)]:
        for n in (4,6,8):print(json.dumps(run(kind,n,end)),flush=True)
