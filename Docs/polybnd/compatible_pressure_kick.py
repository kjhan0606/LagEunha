"""Fixed-geometry pressure substep algebra gate, not a coupled hydro solver.
Positive quadrature states; affine Galerkin force lift; shared midpoint work.
Pressure is frozen over this substep. No entropy/positivity theorem is claimed.
"""
from dataclasses import dataclass,replace
import numpy as np
from shapely.geometry.polygon import orient
from positive_transport_patch import rule,GAMMA

@dataclass
class Cloud:
    polygon: object
    x: np.ndarray
    m: np.ndarray
    u: np.ndarray  # thermal energy per quadrature point, not density
    v: np.ndarray


def cloud(poly,rho,pressure,velocity):
    x,w=rule(poly)
    return Cloud(poly,x,w*np.broadcast_to(rho,(len(w),)),w*np.asarray(pressure(x))/(GAMMA-1),np.asarray(velocity(x)))


def fit(s):
    mass=s.m.sum();c=s.m@s.x/mass;mean=s.m@s.v/mass;dx=s.x-c
    S=dx.T@(s.m[:,None]*dx);C=dx.T@(s.m[:,None]*(s.v-mean))
    B=np.linalg.solve(S,C).T;res=s.v-mean-dx@B.T
    R=res.T@(s.m[:,None]*res)
    return mass,c,mean,S,B,R


def face_velocity(f,x): return f[2]+f[4]@(x-f[1])


def faces(sites,weights,cells):
    out=[]
    for i,cell in enumerate(cells):
        p=np.asarray(orient(cell,1).exterior.coords[:-1]);delta=sites-sites[i]
        bound=np.sum(delta*delta,axis=1)+weights[i]-weights
        for a,b in zip(p,np.roll(p,-1,axis=0)):
            length=np.linalg.norm(b-a)
            if length<1e-14: continue
            normal=np.array([b[1]-a[1],a[0]-b[0]])/length
            r=abs(2*delta@((a+b)/2-sites[i])-bound);r[i]=np.inf;j=int(np.argmin(r))
            if r[j]>1e-10:j=-1
            if j<0 or i<j: out.append((i,j,a,b,normal,length))
    return out


def totals(states):
    p=np.zeros(2);L=K=U=0.
    for s in states:
        p+=s.m@s.v;L+=np.sum(s.m*(s.x[:,0]*s.v[:,1]-s.x[:,1]*s.v[:,0]))
        K+=.5*np.sum(s.m*np.sum(s.v*s.v,axis=1));U+=s.u.sum()
    return np.r_[p,L,K+U]


def kick(states,edges,dt,boundary='wall',pressure_face=None):
    """Return new states or reject without mutating inputs.
    pressure_face(x) is a manufactured reference override for static tests only.
    boundary='trace' uses reconstructed midpoint fluid velocity at the boundary;
    boundary='wall' uses zero normal boundary work. Neither moves geometry.
    """
    if dt<=0 or boundary not in ('wall','trace'): raise ValueError('invalid kick options')
    if any(np.any(s.m<=0) or np.any(s.u<=0) for s in states): raise ValueError('nonpositive state')
    fs=[fit(s) for s in states];N=len(states)
    PV=np.array([(GAMMA-1)*s.u.sum() for s in states])
    P=np.array([pv/s.polygon.area for pv,s in zip(PV,states)])
    force=np.zeros((N,2));moment=PV[:,None,None]*np.eye(2)[None,:,:]
    records=[];momentum_boundary=np.zeros(2);angular_boundary=0.
    for i,j,a,b,n,length in edges:
        for t in (.5-.5/np.sqrt(3),.5+.5/np.sqrt(3)):
            x=a+t*(b-a);weight=length/2
            pf=float(pressure_face(x)) if pressure_face is not None else (P[i]+P[j])/2 if j>=0 else P[i]
            traction=pf*weight*n
            force[i]-=traction;moment[i]-=np.outer(x-fs[i][1],traction)
            if j>=0:
                force[j]+=traction;moment[j]+=np.outer(x-fs[j][1],traction)
            else:
                momentum_boundary-=dt*traction
                angular_boundary-=dt*(x[0]*traction[1]-x[1]*traction[0])
            records.append((i,j,x,n,pf*weight))
    velocity_states=[];dK=[]
    for i,s in enumerate(states):
        mass,c,_,S,_,_=fs[i]
        A=np.linalg.solve(S,moment[i]).T
        dv=dt*(force[i]/mass+(s.x-c)@A.T)
        velocity_states.append(replace(s,v=s.v+dv))
        # Midpoint identity gives the exact discrete kinetic energy change.
        dK.append(float(np.sum(s.m*np.sum((s.v+.5*dv)*dv,axis=1))))
    mid=[fit(replace(s,v=(s.v+v.v)/2)) for s,v in zip(states,velocity_states)]
    dE=np.zeros(N);energy_boundary=0.
    for i,j,x,n,pw in records:
        vi=face_velocity(mid[i],x)
        vf=(vi+face_velocity(mid[j],x))/2 if j>=0 else vi if boundary=='trace' else np.zeros(2)
        work=dt*pw*np.dot(vf,n);dE[i]-=work
        if j>=0: dE[j]+=work
        else: energy_boundary-=work
    result=[]
    for i,(s,v) in enumerate(zip(states,velocity_states)):
        U=s.u.sum();newU=U+dE[i]-dK[i]
        if newU<=0 or not np.isfinite(newU): raise RuntimeError('pressure step rejected: nonpositive internal energy')
        # Keep the positive thermal profile; no residual kinetic energy is removed.
        result.append(replace(v,u=s.u*(newU/U)))
    ledger=np.r_[momentum_boundary,angular_boundary,energy_boundary]
    residual=totals(result)-totals(states)-ledger
    R_error=max(np.max(abs(fit(a)[5]-fit(b)[5])) for a,b in zip(states,result))
    return result,dict(ledger=ledger.tolist(),conservation_error=float(np.max(abs(residual))),
        residual_covariance_change=float(R_error),min_thermal_energy=min(s.u.min() for s in result),
        kinetic_change=sum(dK),thermal_change=sum(s.u.sum()-a.u.sum() for s,a in zip(result,states)))
