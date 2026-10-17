"""Experimental compatible pressure lift and donor moment transport.
Shared HLLC mass/pressure samples; donor affine velocity and residual covariance.
Finite-R flux closure remains a research model, not proven Euler consistency.
"""
import numpy as np
import moment_ale_evolve as base


def kinetic_flux(v,R,fm):
    """Third central velocity moment set to zero; R is covariance per mass."""
    second=np.outer(v,v)+R
    # n enters separately below; this helper returns the noncentral moments.
    return v*fm,second*fm,.5*(v@v+np.trace(R))*fm


def gaussian_flux(v,R,rho,n,fm):
    """Flux moments of a symmetric positive velocity distribution at a face."""
    p,Pi,E=kinetic_flux(v,R,fm);sn=rho*(R@n)
    return p+sn,Pi+np.outer(v,sn)+np.outer(sn,v),E+v@sn


def rhs(t,flat,mass,residual_stress=False,linear_pressure=False,covariance_factors=None,return_covariance_terms=False):
    count=len(mass);y=flat[:-4].reshape(count,20);x=y[:,:2];w=y[:,2]
    polys,edges=base.geometry(x,w,t);raw=np.column_stack([mass,y[:,3:]])
    states=[base.decode(poly,r) for poly,r in zip(polys,raw)]
    if covariance_factors is not None:
        for s,L in zip(states,covariance_factors):s['R']=L@L.T
    gen=np.array([base.velocity(s,point) for s,point in zip(states,x)])
    centers=np.array([s['density'].origin for s in states]);grad=np.zeros((count,2))
    if linear_pressure:
        neighbors=[set() for _ in range(count)]
        for i,j,*_ in edges:
            if j>=0:neighbors[i].add(j);neighbors[j].add(i)
        for i,nn in enumerate(neighbors):
            ids=list(nn);dr=centers[ids]-centers[i];length=np.linalg.norm(dr,axis=1)
            delta=np.array([states[j]['P']-states[i]['P'] for j in ids])
            fit,_,rank,_=np.linalg.lstsq(dr/length[:,None],delta/length,rcond=None)
            if rank<2:raise RuntimeError('pressure gradient stencil rank deficient')
            grad[i]=fit
    def pressure_at(i,point):
        pressure=states[i]['P']+grad[i]@(point-centers[i])
        if pressure<=0:raise RuntimeError('linear pressure reconstruction became nonpositive')
        return pressure
    records=[];data=[];gx,gw=np.polynomial.legendre.leggauss(4);gx=(gx+1)/2;gw/=2
    for i,j,a,b,n,length in edges:
        distance=np.linalg.norm(x[j]-x[i]) if j>=0 else 0.
        for u,weight in zip(gx,gw):
            point=a+u*(b-a);vi=base.velocity(states[i],point);rhoi=float(states[i]['density'].value(point)[0]);pi=pressure_at(i,point)
            if j>=0:
                vj=base.velocity(states[j],point);rhoj=float(states[j]['density'].value(point)[0]);pj=pressure_at(j,point)
                wn=((x[j]-x[i])@gen[j]-(gen[j]-gen[i])@(point-x[i]))/distance
            else:
                wn=base.boundary_velocity(point)@n;vj=vi-2*(vi@n-wn)*n;rhoj=rhoi;pj=pi
            fan=base.fan(base.state(rhoi,vi,pi),base.state(rhoj,vj,pj),n)
            records.append((i,j,length*weight,wn,distance,fan));data.append((point,n))
    z,values,mass_residual=base.solve(records,count)
    derivative=np.zeros((count,18));ledger=np.zeros(4);F=np.zeros((count,2));H=np.zeros((count,2,2))
    outgoing=np.zeros(count);incoming=[[] for _ in range(count)]
    for rec,value,info in zip(records,values,data):
        i,j,length,_,_,fan=rec;flux,branch,wn=value;point,n=info
        uf=fan[1][branch];vf=uf[1:3]/uf[0];fm=flux[0]
        pressure=float((flux[1:3]-vf*fm)@n)
        donor=i if fm>=0 or j<0 else j;s=states[donor]
        v=base.velocity(s,point);rho=float(s['density'].value(point)[0]);R=s['R']/mass[donor]
        if return_covariance_terms:
            if residual_stress or covariance_factors is None:
                raise ValueError('factor covariance transport requires no residual stress and explicit factors')
            outcell=i if fm>=0 else j;receiver=j if fm>=0 else i;rate=length*abs(fm)
            if outcell>=0:outgoing[outcell]+=rate
            if receiver>=0 and rate>0:
                delta=v-base.velocity(states[receiver],point)
                incoming[receiver].append(np.sqrt(rate)*np.column_stack([covariance_factors[donor]/np.sqrt(mass[donor]),delta]))
        if residual_stress:
            advp,advPi,advE=gaussian_flux(v,R,rho,n,fm)
        else:
            advp,advPi,advE=kinetic_flux(v,R,fm)
        fp=advp+pressure*n
        fe=advE+pressure_at(donor,point)/((base.GAMMA-1)*rho)*fm+pressure*(vf@n)
        ff=length*np.r_[fm,point*fm,(np.outer(point,point)*fm).ravel(),fp,
                         np.outer(point,fp).ravel(),advPi.ravel(),fe]
        derivative[i]-=ff
        traction=length*pressure*n;F[i]-=traction;H[i]-=np.outer(point-states[i]['c'],traction)
        if j>=0:
            derivative[j]+=ff;F[j]+=traction;H[j]+=np.outer(point-states[j]['c'],traction)
        else:ledger-=length*np.r_[fp,fe,point[0]*fp[1]-point[1]*fp[0]]
    for i,s in enumerate(states):
        r=raw[i];p=r[7:9];D=r[9:13].reshape(2,2);Pi=r[13:17].reshape(2,2)
        area=s['density'].area;PV=s['P']*area
        derivative[i,1:3]+=p;derivative[i,3:7]+=(D+D.T).ravel()
        derivative[i,9:13]+=(Pi+PV*np.eye(2)).ravel()
        S=r[3:7].reshape(2,2)-np.outer(r[1:3],r[1:3])/mass[i]
        H[i]+=PV*np.eye(2)
        A=np.linalg.solve(S,H[i]).T;a0=F[i]/mass[i]-A@s['c']
        # Same affine acceleration acts on all positive velocity samples.
        pressure_Pi=np.outer(a0,p)+np.outer(p,a0)+A@D+D.T@A.T
        derivative[i,13:17]+=pressure_Pi.ravel()
    diag=dict(mass_flux=mass_residual,min_pressure=min(s['P'] for s in states),
        density_fit_residual=max(s['density'].residual for s in states),
        density_fit_iterations=max(s['density'].iterations for s in states),
        residual_covariance_min=min(float(np.linalg.eigvalsh(s['R']).min()/mass[i]) for i,s in enumerate(states)),
        residual_covariance_max=max(float(np.max(abs(s['R']))/mass[i]) for i,s in enumerate(states)))
    if return_covariance_terms:
        diag['cov_A']=[-s['B']-.5*outgoing[i]/mass[i]*np.eye(2) for i,s in enumerate(states)]
        diag['cov_Qfactor']=[np.concatenate(a,axis=1) if a else np.zeros((2,0)) for a in incoming]
    return np.r_[np.column_stack([gen,z,derivative[:,1:]]).ravel(),ledger],diag
