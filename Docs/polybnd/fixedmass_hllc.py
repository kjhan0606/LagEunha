"""Contact-resolving four-cell ALE research model; lab-frame HLLC fan.

For a fixed Riemann fan each ALE flux branch is F_k - w_n U_k.
The mass residual Jacobian is minus a positive-density graph Laplacian.
Solve its piecewise-linear root with damped Newton, gauge wdot[-1]=0.
No steering, mass projection, pressure floors, or Rusanov fallback.
"""
import json
import numpy as np
from fixedmass_ale_evolve import GAMMA, DOMAIN, initial, unpack, faces, evolve
from weighted_mass_patch import power_cells


def state(rho, velocity, pressure):
    return np.r_[rho, rho*np.asarray(velocity), pressure/(GAMMA-1)+.5*rho*np.dot(velocity, velocity)]


def fan(ul, ur, normal):
    """Lab-frame Davis bounds and HLLC Rankine-Hugoniot star states."""
    if len(ul)!=len(ur) or len(ul)!=len(normal)+2:
        raise ValueError('Riemann state dimension mismatch')
    def primitive(u):
        rho=u[0]; v=u[1:-1]/rho; vn=np.dot(v,normal)
        p=(GAMMA-1)*(u[-1]-.5*rho*np.dot(v,v))
        if rho<=0 or p<=0: raise RuntimeError('nonpositive HLLC input')
        f=np.r_[rho*vn, u[1:-1]*vn+p*normal, (u[-1]+p)*vn]
        return rho,v,vn,p,np.sqrt(GAMMA*p/rho),f
    rl,vl,ql,pl,cl,fl=primitive(ul); rr,vr,qr,pr,cr,fr=primitive(ur)
    sl=min(ql-cl,qr-cr); sr=max(ql+cl,qr+cr)
    al=rl*(sl-ql); ar=rr*(sr-qr)
    sm=(pr-pl+al*ql-ar*qr)/(al-ar)
    if not sl<sm<sr: raise RuntimeError('unordered HLLC waves')
    ps=pl+al*(sm-ql)
    if ps<=0: raise RuntimeError('nonpositive HLLC star pressure')
    def star(u,r,v,q,p,s,f):
        rs=r*(s-q)/(s-sm)
        us=np.r_[rs, rs*(v+(sm-q)*normal), ((s-q)*u[-1]-p*q+ps*sm)/(s-sm)]
        if rs<=0 or us[-1]-.5*np.dot(us[1:-1],us[1:-1])/rs<=0:
            raise RuntimeError('nonpositive HLLC star state')
        return us, f+s*(us-u)
    uls,fls=star(ul,rl,vl,ql,pl,sl,fl); urs,frs=star(ur,rr,vr,qr,pr,sr,fr)
    return np.array([sl,sm,sr]), np.array([ul,uls,urs,ur]), np.array([fl,fls,frs,fr])


def sample(riemann, wn):
    speeds,states,fluxes=riemann
    branch=int(np.searchsorted(speeds,wn,side='right'))
    return fluxes[branch]-wn*states[branch], states[branch,0], branch


def solve_weights(records, tolerance=2e-13):
    count=1+max(max(r[0],r[1]) for r in records)
    def evaluate(wdot):
        residual=np.zeros(count); lap=np.zeros((count,count)); samples=[]
        for i,j,length,base,distance,f in records:
            wn=base+(wdot[i]-wdot[j])/(2*distance) if j>=0 else 0.
            flux,rho,branch=sample(f,wn); samples.append((flux,branch,wn))
            residual[i]+=length*flux[0]
            if j>=0:
                residual[j]-=length*flux[0]; c=length*rho/(2*distance)
                lap[i,i]+=c; lap[j,j]+=c; lap[i,j]-=c; lap[j,i]-=c
        return residual,lap,samples
    wdot=np.zeros(count)
    for iteration in range(40):
        residual,lap,samples=evaluate(wdot); norm=np.max(abs(residual))
        if norm<tolerance: return wdot,samples,iteration+1,norm
        step=np.r_[np.linalg.solve(lap[:-1,:-1],residual[:-1]),0.]
        for backtrack in range(25):
            trial=wdot+step*2.**(-backtrack)
            rnew,_,_=evaluate(trial)
            if np.max(abs(rnew))<norm or np.max(abs(rnew))<tolerance:
                wdot=trial; break
        else: raise RuntimeError('HLLC weight line search failed')
    raise RuntimeError('HLLC weight Newton failed')


def rhs(y,mass):
    sites,weights,momentum,energy=unpack(y)
    cells=power_cells(sites,weights,DOMAIN); area=np.array([c.area for c in cells])
    if np.any(area<=0): raise RuntimeError('empty power cell')
    velocity=momentum/mass[:,None]; rho=mass/area
    pressure=(GAMMA-1)*(energy-.5*np.sum(momentum**2,axis=1)/mass)/area
    if np.any(pressure<=0) or not np.isfinite(pressure).all(): raise RuntimeError('nonpositive pressure')
    states=np.column_stack((rho,momentum/area[:,None],energy/area)); records=[]; pairs=[]
    for i,j,mid,n,length in faces(sites,weights,cells):
        if j>=0:
            distance=np.linalg.norm(sites[j]-sites[i])
            base=(np.dot(sites[j],velocity[j])-np.dot(sites[i],velocity[i])-np.dot(velocity[j]-velocity[i],mid))/distance
            other=states[j]
        else:
            distance=0.; base=0.
            other=states[i].copy(); other[1:3]-=2*np.dot(other[1:3],n)*n
        records.append((i,j,length,base,distance,fan(states[i],other,n)))
        pairs.append((i,j))
    wdot,samples,iterations,residual=solve_weights(records)
    change=np.zeros((4,4)); impulse=np.zeros(2); max_edge_mass=0.
    for record,(flux,branch,wn) in zip(records,samples):
        i,j,length,*_=record; f=length*flux; change[i]-=f
        if j>=0: change[j]+=f; max_edge_mass=max(max_edge_mass,abs(f[0]))
        else: impulse-=f[1:3]
    return np.r_[velocity.ravel(),wdot,change[:,1:3].ravel(),change[:,3],impulse], dict(
        mass_flux=float(np.max(abs(change[:,0]))),pressure_min=float(pressure.min()),
        pressure_max=float(pressure.max()),area_min=float(area.min()),wave_iterations=iterations,
        edge_mass_flux=float(max_edge_mass),signature=tuple((i,j,s[1]) for (i,j),s in zip(pairs,samples)),
        topology=tuple(pairs))


def initial_case(kind):
    if kind not in ('contact100','mixed'): return initial(kind)
    mass,y=initial('contact' if kind=='contact100' else 'acoustic')
    sites,w,p,e=unpack(y); area=np.array([c.area for c in power_cells(sites,w,DOMAIN)])
    mass=area*np.array([1.,100.,100.,1.] if kind=='contact100' else [1.,.7,.7,1.])
    return mass,y


def run(kind,dt):
    counts=dict(branch_changes=0,topology_changes=0,max_edge_mass_flux=0.)
    previous=[None,None]
    def observed(y,mass):
        dy,d=rhs(y,mass)
        for index,(key,field) in enumerate((('branch_changes','signature'),('topology_changes','topology'))):
            if previous[index] is not None and previous[index]!=d[field]: counts[key]+=1
            previous[index]=d[field]
        counts['max_edge_mass_flux']=max(counts['max_edge_mass_flux'],d['edge_mass_flux'])
        return dy,d
    result=evolve(kind,dt,rhs_fn=observed,initial_fn=initial_case)
    mass,y0=initial_case(kind); sites,w,_,_=unpack(np.array(result['final_state']))
    a=np.array([c.area for c in power_cells(sites,w,DOMAIN)])
    x0,w0,_,_=unpack(y0); a0=np.array([c.area for c in power_cells(x0,w0,DOMAIN)])
    result.update(counts,final_density=(mass/a).tolist(),density_change=float(np.max(abs(mass/a-mass/a0))))
    return result


if __name__=='__main__':
    for kind in ('uniform','contact','contact100','acoustic','mixed'):
        for dt in (.01,.005,.0025,.00125):
            print(json.dumps(run(kind,dt)),flush=True)
