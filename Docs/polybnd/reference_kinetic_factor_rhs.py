"""Positive discrete-velocity face transport plus compatible pressure force.
Four symmetric velocity samples represent each residual covariance. Each sample
is upwinded by its own relative normal speed. HLLC supplies pressure/work only.
This is an experimental moment closure, not a proved general Euler solver.
"""
import numpy as np
import moment_ale_evolve as base

ALPHA=np.array([[np.sqrt(2),0.],[-np.sqrt(2),0.],[0.,np.sqrt(2)],[0.,-np.sqrt(2)]])


def transport(velL,velR,rhoL,rhoR,n,wn):
    sl=velL@n-wn;sr=velR@n-wn
    left=rhoL*.25*np.maximum(sl,0.);right=rhoR*.25*np.minimum(sr,0.)
    derivative=rhoL*.25*np.where(sl>0,1.,np.where(sl<0,0.,.5)).sum()+rhoR*.25*np.where(sr<0,1.,np.where(sr>0,0.,.5)).sum()
    return left,right,float(derivative)


def evaluate_records(records,N,z):
    r=np.zeros(N);L=np.zeros((N,N));values=[]
    for rec in records:
        i,j,scale,wn0,d,vl,vr,rl,rr,n,*_=rec
        wn=wn0+(z[i]-z[j])/(2*d) if j>=0 else wn0
        a,b,rho=transport(vl,vr,rl,rr,n,wn);fm=sum(a)+sum(b)
        values.append((a,b,wn));r[i]+=scale*fm
        if j>=0:
            r[j]-=scale*fm;c=scale*rho/(2*d)
            L[i,i]+=c;L[j,j]+=c;L[i,j]-=c;L[j,i]-=c
    return r,L,values


def solve(records,N):
    def evaluate(z):return evaluate_records(records,N,z)
    z=np.zeros(N)
    for iteration in range(50):
        r,L,values=evaluate(z);norm=float(max(abs(r)))
        if norm<2e-13:return z,values,norm
        # Vacuum intervals may leave several independent weight gauges.
        step=np.r_[np.linalg.lstsq(L[:-1,:-1],r[:-1],rcond=1e-13)[0],0.]
        for k in range(30):
            trial=z+step*2.**(-k)
            if max(abs(evaluate(trial)[0]))<norm:z=trial;break
        else:raise RuntimeError('kinetic mass constraint line search failed')
    raise RuntimeError('kinetic mass constraint did not converge')


def rhs(t,flat,mass,linear_pressure=True,covariance_factors=None,return_covariance_terms=True,
        audit_geometry_rates=None):
    # Optional rates ONLY evaluate the residual of an existing dense path.
    # Evolution callers leave this unset, retaining unsteered material sites
    # and the instantaneous mass solve. Do not use the audit mode to advance.
    if covariance_factors is None:raise ValueError('positive covariance factors required')
    N=len(mass);y=flat[:-4].reshape(N,20);x=y[:,:2];w=y[:,2]
    polys,edges=base.geometry(x,w,t);raw=np.column_stack([mass,y[:,3:]])
    states=[base.decode(p,r) for p,r in zip(polys,raw)]
    for s,L in zip(states,covariance_factors):s['R']=L@L.T
    material_gen=np.array([base.velocity(s,p) for s,p in zip(states,x)])
    gen=material_gen
    if audit_geometry_rates is not None:
        gen=np.asarray(audit_geometry_rates[0]);prescribed_z=np.asarray(audit_geometry_rates[1])
        if gen.shape!=(N,2) or prescribed_z.shape!=(N,) or not np.all(np.isfinite(gen)) or not np.all(np.isfinite(prescribed_z)):
            raise ValueError('invalid diagnostic geometry rates')
    centers=np.array([s['density'].origin for s in states]);grad=np.zeros((N,2))
    if linear_pressure:
        neighbors=[set() for _ in range(N)]
        for i,j,*_ in edges:
            if j>=0:neighbors[i].add(j);neighbors[j].add(i)
        for i,nn in enumerate(neighbors):
            ids=list(nn);dr=centers[ids]-centers[i];d=np.linalg.norm(dr,axis=1)
            dp=np.array([states[j]['P']-states[i]['P'] for j in ids])
            g,_,rank,_=np.linalg.lstsq(dr/d[:,None],dp/d,rcond=None)
            if rank<2:raise RuntimeError('pressure stencil deficient')
            grad[i]=g
    def pressure(i,r):
        p=states[i]['P']+grad[i]@(r-centers[i])
        if p<=0:raise RuntimeError('nonpositive reconstructed pressure')
        return p
    records=[];gx,gw=np.polynomial.legendre.leggauss(4);gx=(gx+1)/2;gw/=2
    for i,j,a,b,n,length in edges:
        d=np.linalg.norm(x[j]-x[i]) if j>=0 else 0.
        for g,q in zip(gx,gw):
            point=a+g*(b-a);vl=base.velocity(states[i],point);rl=float(states[i]['density'].value(point)[0]);pl=pressure(i,point)
            vls=vl+ALPHA@covariance_factors[i].T/np.sqrt(mass[i])
            if j>=0:
                vr=base.velocity(states[j],point);rr=float(states[j]['density'].value(point)[0]);pr=pressure(j,point)
                vrs=vr+ALPHA@covariance_factors[j].T/np.sqrt(mass[j])
                wn=((x[j]-x[i])@gen[j]-(gen[j]-gen[i])@(point-x[i]))/d
            else:
                wn=base.boundary_velocity(point)@n;vr=vl-2*(vl@n-wn)*n;rr=rl;pr=pl
                vrs=vls-2*((vls@n-wn)[:,None]*n)
            fan=base.fan(base.state(rl,vl,pl),base.state(rr,vr,pr),n)
            sm=fan[0][1];ps=float((fan[2][1][1:3]-fan[1][1][1:3]*sm)@n)
            records.append((i,j,length*q,wn,d,vls,vrs,rl,rr,n,point,pl,pr,ps,sm))
    if audit_geometry_rates is None:
        z,values,mass_residual=solve(records,N)
    else:
        z=prescribed_z
        residual,_,values=evaluate_records(records,N,z)
        mass_residual=float(np.max(abs(residual)))
    dr=np.zeros((N,18));ledger=np.zeros(4);F=np.zeros((N,2));H=np.zeros((N,2,2));K=np.zeros((N,2,2));incoming=[[] for _ in range(N)]
    for rec,(wl,wr,wn) in zip(records,values):
        i,j,scale,_,_,vl,vr,rl,rr,n,point,pl,pr,ps,sm=rec
        fm=sum(wl)+sum(wr);fp=wl@vl+wr@vr+ps*n
        Pi=vl.T@(wl[:,None]*vl)+vr.T@(wr[:,None]*vr)
        fe=wl@(.5*np.sum(vl*vl,axis=1)+pl/((base.GAMMA-1)*rl))+wr@(.5*np.sum(vr*vr,axis=1)+pr/((base.GAMMA-1)*rr))+ps*sm
        ff=scale*np.r_[fm,point*fm,(np.outer(point,point)*fm).ravel(),fp,np.outer(point,fp).ravel(),Pi.ravel(),fe]
        dr[i]-=ff;traction=scale*ps*n;F[i]-=traction;H[i]-=np.outer(point-states[i]['c'],traction)
        if j>=0:dr[j]+=ff;F[j]+=traction;H[j]+=np.outer(point-states[j]['c'],traction)
        else:ledger-=scale*np.r_[fp,fe,point[0]*fp[1]-point[1]*fp[0]]
        K[i]+=scale/mass[i]*(ALPHA.T@(wl[:,None]*ALPHA))
        if j>=0:K[j]+=scale/mass[j]*(ALPHA.T@((-wr)[:,None]*ALPHA))
        # Incoming samples are positive contributions to the receiver residual.
        for receiver,vel,rates in [(i,vr,-scale*wr),(j,vl,scale*wl)]:
            if receiver<0:continue
            delta=vel-base.velocity(states[receiver],point)
            incoming[receiver].append((np.sqrt(rates)[:,None]*delta).T)
    for i,s in enumerate(states):
        r=raw[i];p=r[7:9];D=r[9:13].reshape(2,2);Pi=r[13:17].reshape(2,2);PV=s['P']*s['density'].area
        dr[i,1:3]+=p;dr[i,3:7]+=(D+D.T).ravel();dr[i,9:13]+=(Pi+PV*np.eye(2)).ravel()
        S=r[3:7].reshape(2,2)-np.outer(r[1:3],r[1:3])/mass[i];H[i]+=PV*np.eye(2)
        A=np.linalg.solve(S,H[i]).T;a0=F[i]/mass[i]-A@s['c']
        dr[i,13:17]+=(np.outer(a0,p)+np.outer(p,a0)+A@D+D.T@A.T).ravel()
    diag=dict(cov_linearization_factors=covariance_factors,mass_flux=mass_residual,cov_A=[-s['B'] for s in states],cov_K=list(K),
        cov_Qfactor=[np.concatenate(a,axis=1) if a else np.zeros((2,0)) for a in incoming],
        min_pressure=min(s['P'] for s in states),mass_rates=dr[:,0].copy(),
        generator_velocity_defect=float(np.max(abs(gen-material_gen))))
    return np.r_[np.column_stack([gen,z,dr[:,1:]]).ravel(),ledger],diag
