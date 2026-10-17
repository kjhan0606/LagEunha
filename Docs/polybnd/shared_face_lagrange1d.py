"""1D shared-face geometry/work gate. Acoustic Riemann + explicit midpoint.
Fixed cell masses, physically advected generators, reconstructible power weights.
Piecewise constant FV state: no claim of 2D subcell-moment closure.
"""
import json
from pathlib import Path
import numpy as np
GAMMA=5/3


def decode(y,m):
    f,x,p,E=y;V=np.diff(f);u=p/m;U=E-p*p/(2*m)
    if np.any(V<=0) or np.any(U<=0) or not all(np.isfinite(a).all() for a in y):
        raise ValueError('nonpositive geometry/thermal state')
    rho=m/V;P=(GAMMA-1)*U/V
    return V,u,rho,P


def power_weights(f,x):
    if np.any(np.diff(x)<=0):raise ValueError('generator crossing')
    w=np.r_[0.,np.cumsum(np.diff(x*x)-2*np.diff(x)*f[1:-1])]
    reconstructed=(np.diff(x*x)-np.diff(w))/(2*np.diff(x))
    return w,float(np.max(abs(reconstructed-f[1:-1])))


def rhs(y,m,boundary):
    f,x,p,E=y;V,u,rho,P=decode(y,m);Z=np.sqrt(GAMMA*rho*P)
    uf=np.empty(len(f));pf=np.empty(len(f));s=Z[:-1]+Z[1:]
    uf[1:-1]=(Z[:-1]*u[:-1]+Z[1:]*u[1:]+P[:-1]-P[1:])/s
    pf[1:-1]=(Z[1:]*P[:-1]+Z[:-1]*P[1:]+Z[:-1]*Z[1:]*(u[:-1]-u[1:]))/s
    uf[0],uf[-1]=boundary
    pf[0]=P[0]+Z[0]*(boundary[0]-u[0])
    pf[-1]=P[-1]+Z[-1]*(u[-1]-boundary[1])
    if min(pf)<=0:raise ValueError('nonpositive acoustic star pressure')
    work=pf*uf
    return (uf,u,-np.diff(pf),-np.diff(work)),np.array([pf[0]-pf[-1],work[0]-work[-1]])


def add(y,k,a):return tuple(v+a*d for v,d in zip(y,k))


def step(y,m,dt,boundary):
    k,_=rhs(y,m,boundary);half=add(y,k,dt/2);decode(half,m)
    kh,ledger=rhs(half,m,boundary);out=add(y,kh,dt);decode(out,m)
    w,power_error=power_weights(out[0],out[1])
    if np.any(out[1]<=out[0][:-1]) or np.any(out[1]>=out[0][1:]):
        raise ValueError('generator left its represented cell')
    # This identity is local, not just a global volume check.
    gcl=float(np.max(abs(np.diff(out[0])-np.diff(y[0])-dt*np.diff(kh[0]))))
    return out,dt*ledger,dict(gcl=gcl,power_error=power_error,weights=w)


def initial(n,kind,boost=0.):
    f=np.linspace(0,1,n+1);x=(f[:-1]+f[1:])/2;V=np.diff(f);P0=.6;c=np.sqrt(GAMMA*P0)
    if kind=='acoustic':
        eps=1e-6;avg=np.diff(np.sin(2*np.pi*f))/(2*np.pi*V)
        rho=1+eps/c**2*avg;P=P0+eps*avg;u=np.full(n,boost)
    elif kind=='compression':rho=np.ones(n);P=np.full(n,P0);u=-.2*(x-.5)+boost
    elif kind=='contact':rho=np.where(x<.5,1.,100.);P=np.full(n,P0);u=np.full(n,boost)
    else:raise ValueError(kind)
    m=rho*V;E=P*V/(GAMMA-1)+.5*m*u*u
    return (f,x,m*u,E),m


def run(n,kind,cfl=.4,end=.1,boost=0.):
    y,m=initial(n,kind,boost);initial_y=tuple(a.copy() for a in y)
    boundary=(boost+.1,boost-.1) if kind=='compression' else (boost,boost)
    ledger=np.zeros(2);initial_tot=np.array([sum(y[2]),sum(y[3])]);t=0.;steps=rejects=0
    max_gcl=max_power=max_ledger=0.;minP=minV=np.inf;history=[]
    while t<end-1e-15:
        V,u,rho,P=decode(y,m);dt=min(cfl*np.min(V/np.sqrt(GAMMA*P/rho)),end-t)
        for attempt in range(21):
            try:trial,dl,d=step(y,m,dt,boundary);break
            except ValueError:
                if attempt==20:raise
                dt*=.5;rejects+=1
        y=trial;ledger+=dl;t+=dt;steps+=1
        V,u,rho,P=decode(y,m)
        err=float(np.max(abs(np.array([sum(y[2]),sum(y[3])])-initial_tot-ledger)))
        max_gcl=max(max_gcl,d['gcl']);max_power=max(max_power,d['power_error']);max_ledger=max(max_ledger,err)
        minP=min(minP,float(min(P)));minV=min(minV,float(min(V)))
        history.append(dict(t=t,ledger_error=err,minP=float(min(P)),minV=float(min(V))))
        if steps>100000:raise RuntimeError('step bound exceeded')
    f,x,p,E=y;V,u,rho,P=decode(y,m)
    if kind=='acoustic':
        z=f-boost*t;avgcos=np.diff(np.sin(2*np.pi*z))/(2*np.pi*V)
        avgsin=-np.diff(np.cos(2*np.pi*z))/(2*np.pi*V);c=1.;eps=1e-6
        Pref=.6+eps*avgcos*np.cos(2*np.pi*c*t)
        uref=boost+eps/c*avgsin*np.sin(2*np.pi*c*t)
        physical_error=float(np.sum(V*(abs(P-Pref)+abs(u-uref)))/eps)
    elif kind=='compression':
        Pref=.6*(1-.2*t)**(-GAMMA)
        physical_error=float(np.max(abs(P/Pref-1)))
    else:
        physical_error=float(max(np.max(abs(P-.6)),np.max(abs(u-boost)),np.max(abs(f-initial_y[0]-boost*t))))
    result=dict(n=n,kind=kind,cfl=cfl,end=end,boost=boost,steps=steps,rejects=rejects,
        gcl_error=max_gcl,power_error=max_power,ledger_error=max_ledger,min_pressure=minP,min_volume=minV,
        physical_error=physical_error,max_weight=float(np.max(abs(d['weights']))),
        mass_error=0.,pass_algebra=bool(max(max_gcl,max_power,max_ledger)<2e-11 and minP>0 and minV>0))
    tag=f'{kind}_n{n}_cfl{cfl:g}_boost{boost:g}'
    Path(tag+'.history.json').write_text(json.dumps(history,indent=2))
    Path(tag+'.checkpoint.json').write_text(json.dumps(dict(t=t,faces=f.tolist(),generators=x.tolist(),mass=m.tolist(),momentum=p.tolist(),energy=E.tolist(),weights=d['weights'].tolist(),ledger=ledger.tolist()),indent=2))
    return result

if __name__=='__main__':
    results=[]
    cases=[(n,'acoustic',cfl,0.) for n in (16,32,64) for cfl in (.4,.2)]
    cases += [(n,'compression',.4,0.) for n in (16,32,64)] + [(32,'contact',.4,.7),(32,'acoustic',.4,3.)]
    for n,kind,cfl,boost in cases:
        r=run(n,kind,cfl=cfl,boost=boost);results.append(r);print(json.dumps(r),flush=True)
        if not r['pass_algebra']:raise RuntimeError('geometry/conservation failed')
    for kind in ('acoustic','compression'):
        rr=[r for r in results if r['kind']==kind and r['cfl']==.4 and r['boost']==0.]
        if not all(b['physical_error']<a['physical_error'] for a,b in zip(rr,rr[1:])):
            raise RuntimeError('physical error does not decrease under refinement')
    if next(r for r in results if r['kind']=='contact')['physical_error']>2e-11:raise RuntimeError('contact failed')
