"""PSD residual covariance integration via factors and exponential transport.
R is a primary variable, Pi is derived. No negative-eigenvalue clipping.
Default: donor/compatible-pressure closure WITHOUT residual stress.
--kinetic: experimental sample-wise transport with residual stress. Its time
order must be measured; the right-factor loss extension is not proved order two.
"""
import json,time
from pathlib import Path
import numpy as np
from scipy.linalg import expm
from adaptive_moment_batch import setup,diagnostic
from compatible_moment_rhs import rhs


def rebuild(y,mass,factors):
    y=y.copy();a=y[:-4].reshape(len(mass),20)
    for i,m in enumerate(mass):
        q=a[i,3:5];Q=a[i,5:9].reshape(2,2);p=a[i,9:11];D=a[i,11:15].reshape(2,2)
        S=Q-np.outer(q,q)/m;C=D-np.outer(q,p)/m
        if np.linalg.eigvalsh((S+S.T)/2).min()<=0:raise RuntimeError('nonpositive spatial covariance')
        B=np.linalg.solve(S,C).T
        Pi=np.outer(p,p)/m+B@S@B.T+factors[i]@factors[i].T
        a[i,15:19]=Pi.ravel()
    return y


def advance_factors_exponential(factors,diag,dt):
    result=[]
    for i,(L,A,Q) in enumerate(zip(factors,diag['cov_A'],diag['cov_Qfactor'])):
        loss=expm(-.5*dt*diag['cov_K'][i]) if 'cov_K' in diag else np.eye(L.shape[1])
        U=np.concatenate([expm(dt*A)@L@loss,np.sqrt(dt)*(expm(.5*dt*A)@Q)],axis=1)
        _,r=np.linalg.qr(U.T,mode='reduced');result.append(r.T)
    return result


advance_factors = advance_factors_exponential

def run(kind,amplitude,dt,end=.02):
    mass,y=setup(kind,4,amplitude);n=len(mass)
    # All supported initial fields here have exact R=0 (rest or global affine shear).
    factors=[np.zeros((2,2)) for _ in mass];initial_adjustment=0.
    new=rebuild(y,mass,factors);initial_adjustment=float(np.max(abs(new-y)));y=new;y0=y.copy()
    history=[];start=time.monotonic();calls=0;max_mass_flux=0.
    def f(t,state,L):
        nonlocal calls,max_mass_flux;calls+=1
        value,diag=rhs(t,state,mass,linear_pressure=True,covariance_factors=L,return_covariance_terms=True)
        max_mass_flux=max(max_mass_flux,float(diag['mass_flux']))
        return value,diag
    try:
        for k in range(round(end/dt)):
            t=k*dt;dy,d=f(t,y,factors)
            halfL=advance_factors(factors,d,dt/2);half=rebuild(y+dt*dy/2,mass,halfL)
            dh,diag=f(t+dt/2,half,halfL)
            nextL=advance_factors(factors,diag,dt);nexty=rebuild(y+dt*dh,mass,nextL)
            hm=diagnostic(t+dt/2,half,mass,y0);hn=diagnostic(t+dt,nexty,mass,y0)
            min_gram=min(float(np.linalg.eigvalsh(L@L.T).min()/m) for L,m in zip(halfL+nextL,list(mass)*2))
            history.append(dict(t=t+dt,min_factor_R=min_gram,raw_subtraction_R=hn['min_R'],min_pressure=min(hm['min_pressure'],hn['min_pressure']),ledger_error=hn['ledger_error']))
            y,factors=nexty,nextL
        result=dict(kind=kind,amplitude=amplitude,dt=dt,end=end,status='completed',max_mass_flux=max_mass_flux,calls=calls,seconds=time.monotonic()-start,
                    initial_moment_roundoff_adjustment=initial_adjustment,
                    min_factor_R=min(h['min_factor_R'] for h in history),min_pressure=min(h['min_pressure'] for h in history),
                    final=diagnostic(end,y,mass,y0))
        tag=f'{kind}_a{amplitude:g}_dt{dt:g}'
        Path(tag+'.history.json').write_text(json.dumps(history,indent=2))
        Path(tag+'.checkpoint.json').write_text(json.dumps(dict(t=end,mass=mass.tolist(),state=y.tolist(),factors=[a.tolist() for a in factors]),indent=2))
    except Exception as e:result=dict(kind=kind,amplitude=amplitude,dt=dt,status='failed',error=repr(e),calls=calls,seconds=time.monotonic()-start)
    return result

if __name__=='__main__':
    import sys
    if '--kinetic' in sys.argv:
        from kinetic_factor_rhs import rhs
    for kind,a in [('wave',.001),('wave',.1),('shear',0.),('contact',0.)]:
        for dt in (.001,.0005):print(json.dumps(run(kind,a,dt)),flush=True)


def advance_factors_sylvester(factors,diag,dt):
    """Put the kinetic loss into a bounded symmetric Lyapunov generator.
    A_loss R+R A_loss^T=-L K L^T, including rank-deficient L.
    SVD avoids forming/subtracting a raw covariance or clipping eigenvalues.
    With midpoint coefficients, exp(dt*A/2) also damps incoming covariance.
    Order at rank changes still needs an evolution test.
    """
    corrected=dict(diag);matrices=[]
    for i,(L,A) in enumerate(zip(diag.get('cov_linearization_factors',factors),diag['cov_A'])):
        if 'cov_K' not in diag:matrices.append(A);continue
        U,s,Vh=np.linalg.svd(L)
        if s.max()==0:loss=np.zeros_like(A)
        else:
            scaled=s/s.max();den=scaled[:,None]**2+scaled[None,:]**2
            ratio=np.divide(scaled[:,None]*scaled[None,:],den,out=np.zeros_like(den),where=den>0)
            loss=U@((Vh@diag['cov_K'][i]@Vh.T)*ratio)@U.T
        matrices.append(A-loss)
    corrected['cov_A']=matrices;corrected.pop('cov_K',None)
    return advance_factors_exponential(factors,corrected,dt)
